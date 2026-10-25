"""make_stations.py -- station files for your site from the FCC's public database.

The radar needs, per illuminator: centre frequency, tower latitude /
longitude / radiation-centre height, and the signal type. The FCC's FM and
TV query services (transition.fcc.gov/fcc-bin/fmq and tvq) return all of
it for every licensed US station within a radius of a point, as
pipe-delimited rows. This turns those rows into tools/stations/*.json in
the shape every tool's --station= expects.

Privacy: the query centre is sent to the FCC. By default it is the site's
observer position ROUNDED to 0.1 degree (about 10 km) with the radius
padded by 15 km, so the request says which city, not which house. --exact
sends the real position.

Heights: the FCC's RCAMSL (radiation centre above mean sea level) goes in
as tx_alt_m when present; otherwise HAAT + 100 m, flagged as an estimate.
Coordinates are as filed (NAD27 for older records: tens of metres off
WGS84, under one 8-VSB range cell).

What the FCC does not say: which TV channels run ATSC 3.0. Every DTV
channel is written as 8-VSB; --atsc3 25,36 marks channels you know are
ATSC 3.0 (rabbitears.info lists them). Co-channel licensees are ranked by
ERP over distance squared from the centre; the others are recorded under
"ambiguous_with" (another mast) or "shares_channel_with" (the same mast:
channel sharing after the repack). Mast heights are the median of the
FCC records on that mast when three or more agree.

Usage:
  python tools/make_stations.py --list                     # table only, from the site's observer
  python tools/make_stations.py --radius-km 80 --out tools/stations
  python tools/make_stations.py --centre 39.7,-105.0 --radius-km 60 --out lab/stations --atsc3 24
  options: --fm-only / --tv-only, --min-erp-kw FM TV (default 1 10), --include-lptv, --include-cp, --exact
"""
import argparse
import json
import math
import os
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import radar_site                                              # noqa: E402

FCC = "https://transition.fcc.gov/fcc-bin/"
FS_VSB = 7174825.174825176          # vsb_ref FS3 (10/8 of the 8-VSB symbol rate... 2/3 of 10.76 MS/s)
FS_ATSC3 = 6.912e6
FS_FM = 250e3


def tv_fc(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    if rf <= 6:
        lo = 54 + (rf - 2) * 6 if rf <= 4 else 76 + (rf - 5) * 6
    return (lo + 3.0) * 1e6


def dms(cols, i):
    """Columns i..i+3 are N/S, deg, min, sec (or E/W, deg, min, sec)."""
    hemi, d, m, s = cols[i].strip(), float(cols[i + 1]), float(cols[i + 2]), float(cols[i + 3])
    v = d + m / 60 + s / 3600
    return -v if hemi in ("S", "W") else v


def num(s, default=None):
    m = re.search(r"[-+]?\d*\.?\d+", s or "")
    return float(m.group()) if m else default


def query(kind, lat, lon, radius_km, lo, hi):
    """kind 'fm' (lo/hi MHz) or 'tv' (lo/hi RF channel). Returns a list of column lists."""
    ns = "N" if lat >= 0 else "S"; ew = "W" if lon < 0 else "E"
    la, lo_ = abs(lat), abs(lon)
    pos = dict(dlat2=int(la), mlat2=int((la * 60) % 60), slat2=round((la * 3600) % 60, 1), NS=ns,
               dlon2=int(lo_), mlon2=int((lo_ * 60) % 60), slon2=round((lo_ * 3600) % 60, 1), EW=ew)
    if kind == "fm":
        q = dict(state="", call="", arn="", city="", freq=lo, fre2=hi, serv="FM", vac=3, facid="", asrn="",
                 **{"class": ""}, list=4, dist=int(radius_km), **pos, size=9)
        url = FCC + "fmq?"
    else:
        q = dict(state="", call="", arn="", city="", chan=lo, cha2=hi, serv="", type=3, facid="", asrn="",
                 list=4, dist=int(radius_km), **pos, size=9)
        url = FCC + "tvq?"
    url += "&".join(f"{k}={v}" for k, v in q.items())
    text = fetch(url)
    return [ln.split("|") for ln in text.splitlines() if ln.startswith("|")]


def fetch(url):
    """GET as text. The FCC front end refuses some client/agent pairings (a
    python client with a browser-like agent got 403 here); urllib with its
    own agent passed, and curl always does, so curl is the fallback."""
    import shutil
    import subprocess
    try:
        with urllib.request.urlopen(urllib.request.Request(url), timeout=90) as r:
            return r.read().decode("latin-1")
    except urllib.error.HTTPError as e:
        if not shutil.which("curl"):
            raise
        print(f"[fcc] urllib got HTTP {e.code}; retrying with curl")
        r = subprocess.run(["curl", "-s", "-m", "120", "-A", "curl/8.0", url], capture_output=True)
        if r.returncode != 0 or not r.stdout:
            raise
        return r.stdout.decode("latin-1")


def parse_rows(rows, kind):
    out = []
    for c in rows:
        if len(c) < 32:
            continue
        call = c[1].strip()
        service = c[3].strip()
        status = c[9].strip()
        rec = dict(call=call, service=service, status=status, city=c[10].strip(), state=c[11].strip(),
                   facility_id=c[18].strip(), erp_kw=num(c[14], 0.0), haat_m=num(c[16]),
                   lat=dms(c, 19), lon=dms(c, 23), licensee=re.sub(r"\s+", " ", c[27].strip()),
                   rcamsl_m=num(c[31]))
        if kind == "fm":
            rec["freq_mhz"] = num(c[2])
            rec["fm_class"] = c[7].strip()
        else:
            rec["rf"] = int(num(c[4], 0))
        out.append(rec)
    return out


def dedupe(recs, include_cp):
    """One record per facility: a licence beats a construction permit beats an application."""
    rank = {"LIC": 0, "LIC-M": 0, "CP": 1, "CP MOD": 1, "APP": 2}
    best = {}
    for r in recs:
        if r["status"] not in rank and not include_cp:
            continue
        if r["status"] != "LIC" and not include_cp:
            continue
        k = r["facility_id"] or (r["call"], r.get("freq_mhz"), r.get("rf"))
        if k not in best or rank.get(r["status"], 3) < rank.get(best[k]["status"], 3):
            best[k] = r
    return list(best.values())


def strength(r, lat, lon):
    """Expected received power, ERP over distance squared, from the (rounded)
    query centre. Ranks co-channel licensees: the one you receive, not the
    strongest one in the radius (89.7 at 135 km: a 27 kW station 120 km away
    is not the 2.5 kW one 15 km away)."""
    d = math.hypot((r["lat"] - lat) * 111.0, (r["lon"] - lon) * 111.0 * math.cos(math.radians(lat)))
    return r["erp_kw"] / max(d, 2.0) ** 2


def slug(call):
    return re.sub(r"[^a-z0-9]+", "_", call.lower()).strip("_")


def site_key(r):
    """Records within ~500 m are one mast (channel sharers, candelabras, antenna farms)."""
    return (round(r["lat"] / 0.005), round(r["lon"] / 0.005))


def fix_heights(recs):
    """tx height per record. The FCC's RCAMSL is right most of the time and
    wrong by a kilometre now and then (WDCW on RF15: 1532 m in Washington,
    every other record on that mast 310 m). Where three or more records
    share a mast, the mast's height is their median and a record more than
    300 m off it is replaced and flagged."""
    masts = {}
    for r in recs:
        if r["rcamsl_m"] and r["rcamsl_m"] > 0:
            masts.setdefault(site_key(r), []).append(float(r["rcamsl_m"]))
    for r in recs:
        hs = sorted(masts.get(site_key(r), []))
        med = hs[len(hs) // 2] if len(hs) >= 3 else None
        if r["rcamsl_m"] and r["rcamsl_m"] > 0:
            if med is not None and abs(r["rcamsl_m"] - med) > 300:
                r["alt_m"], r["alt_note"] = med, f"mast median of {len(hs)} FCC records; this record said {r['rcamsl_m']:.0f} m"
            else:
                r["alt_m"], r["alt_note"] = float(r["rcamsl_m"]), "FCC RCAMSL"
        elif med is not None:
            r["alt_m"], r["alt_note"] = med, f"mast median of {len(hs)} FCC records (none on this one)"
        elif r["haat_m"] is not None:
            r["alt_m"], r["alt_note"] = float(r["haat_m"]) + 100.0, "ESTIMATE: HAAT + 100 m ground"
        else:
            r["alt_m"], r["alt_note"] = 300.0, "ESTIMATE: no height on file"


def alt_m(r):
    return r["alt_m"], r["alt_note"]


def fm_station(r, others):
    alt, note = alt_m(r)
    st = {
        "name": f"{r['call']} {r['freq_mhz']:.1f} {r['city']}, {r['state']} (FCC public record)",
        "short_name": f"{r['freq_mhz']:.1f} FM / {r['call']}",
        "fc_hz": round(r["freq_mhz"] * 1e6, 1),
        "fs_hz": FS_FM,
        "tx_lat": round(r["lat"], 6), "tx_lon": round(r["lon"], 6), "tx_alt_m": alt, "tx_alt_note": note,
        "erp_kw": r["erp_kw"], "fm_class": r["fm_class"], "facility_id": r["facility_id"],
        "max_range_m": 120000.0, "r_gate_m": 2200.0, "f_gate_hz": 8.0, "notch_hz": 5.0,
        "identity": "FCC licensee on this frequency with the highest ERP/distance^2 from the query centre; RDS not checked",
        "illuminator": "FM broadcast, analog, 200 kHz; reference = reconstructed",
        "source": "FCC fmq, tools/make_stations.py",
    }
    same = [o for o in others if site_key(o) == site_key(r)]
    diff = [o for o in others if site_key(o) != site_key(r)]
    if same:
        st["shares_mast_with"] = [f"{o['call']} ({o['erp_kw']:g} kW)" for o in same]
    if diff:
        st["ambiguous_with"] = [f"{o['call']} {o['city']}, {o['state']} ({o['erp_kw']:g} kW)" for o in diff]
    return st


def tv_station(r, atsc3, others=()):
    alt, note = alt_m(r)
    is3 = r["rf"] in atsc3
    same = [o for o in others if site_key(o) == site_key(r)]
    diff = [o for o in others if site_key(o) != site_key(r)]
    extra = {}
    if same:
        extra["shares_channel_with"] = [f"{o['call']} ({o['erp_kw']:g} kW)" for o in same]
    if diff:
        extra["ambiguous_with"] = [f"{o['call']} {o['city']}, {o['state']} ({o['erp_kw']:g} kW)" for o in diff]
    return {**extra,
        "name": f"{r['call']} RF{r['rf']} {r['city']}, {r['state']} (FCC public record)",
        "short_name": f"RF{r['rf']} / {r['call'].replace('-TV', '').replace('-DT', '')} ({'ATSC 3.0' if is3 else '8-VSB'})",
        "fc_hz": tv_fc(r["rf"]),
        "fs_hz": FS_ATSC3 if is3 else FS_VSB,
        "tx_lat": round(r["lat"], 6), "tx_lon": round(r["lon"], 6), "tx_alt_m": alt, "tx_alt_note": note,
        "erp_kw": r["erp_kw"], "service": r["service"], "facility_id": r["facility_id"],
        "max_range_m": 60000.0, "notch_hz": 15.0,
        "illuminator": ("ATSC 3.0 OFDM; pilot-synthesized reference (track_hunt.py)" if is3
                        else "ATSC 1.0 8-VSB; reference rebuilt from symbol decisions (vsb_ref.py)"),
        "source": "FCC tvq, tools/make_stations.py" + ("; ATSC 3.0 per --atsc3" if is3 else ""),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--centre", default=None, help="lat,lon of the query centre (default: the site's observer)")
    ap.add_argument("--exact", action="store_true", help="send the exact centre (default rounds to 0.1 deg)")
    ap.add_argument("--radius-km", type=float, default=80.0)
    ap.add_argument("--fm-only", action="store_true")
    ap.add_argument("--tv-only", action="store_true")
    ap.add_argument("--min-erp-kw", nargs=2, type=float, default=(1.0, 10.0), metavar=("FM", "TV"))
    ap.add_argument("--include-lptv", action="store_true", help="keep low-power TV (LD/LPD) records")
    ap.add_argument("--include-cp", action="store_true", help="keep construction permits and applications")
    ap.add_argument("--atsc3", default="", help="comma list of RF channels that run ATSC 3.0")
    ap.add_argument("--list", action="store_true", help="print the table, write nothing")
    ap.add_argument("--out", default=str(HERE / "stations"), help="TV files go here, FM files in <out>/fm")
    a = ap.parse_args()

    if a.centre:
        lat, lon = (float(v) for v in a.centre.split(","))
    else:
        lat, lon, _ = radar_site.observer()
    radius = a.radius_km
    if not a.exact:
        lat, lon = round(lat, 1), round(lon, 1)
        radius += 15.0
    atsc3 = {int(v) for v in a.atsc3.split(",") if v.strip()}
    print(f"[fcc] query centre {'exact' if a.exact else 'rounded to 0.1 deg'} (withheld), radius {radius:.0f} km")

    written = []
    fm_recs, tv_recs, fm_rows, tv_rows = [], [], [], []
    if not a.tv_only:
        fm_rows = query("fm", lat, lon, radius, 87.5, 108.1)
        fm_recs = [r for r in dedupe(parse_rows(fm_rows, "fm"), a.include_cp) if r["erp_kw"] >= a.min_erp_kw[0]]
    if not a.fm_only:
        tv_rows = query("tv", lat, lon, radius, 2, 36)
        tv_recs = parse_rows(tv_rows, "tv")
        if not a.include_lptv:
            tv_recs = [r for r in tv_recs if r["service"] in ("DTV", "DT", "DC", "DCA") or
                       (r["service"].startswith("D") and "L" not in r["service"])]
        tv_recs = [r for r in dedupe(tv_recs, a.include_cp) if r["erp_kw"] >= a.min_erp_kw[1]]
    fix_heights(fm_recs + tv_recs)
    if not a.tv_only:
        rows, recs = fm_rows, fm_recs
        by_f = {}
        for r in recs:
            by_f.setdefault(round(r["freq_mhz"], 1), []).append(r)
        print(f"[fcc] FM: {len(rows)} rows -> {len(recs)} licensed stations >= {a.min_erp_kw[0]:g} kW on {len(by_f)} frequencies")
        out_fm = Path(a.out) / "fm"
        for f, rs in sorted(by_f.items()):
            rs.sort(key=lambda r: -strength(r, lat, lon))
            r, others = rs[0], rs[1:]
            st = fm_station(r, others)
            amb = f"  AMBIG with {', '.join(o['call'] for o in others)}" if st.get("ambiguous_with") else ""
            print(f"  {f:6.1f}  {r['call']:8s} {r['erp_kw']:6.1f} kW  {r['city']}, {r['state']}  alt {st['tx_alt_m']:.0f} m{amb}")
            if not a.list:
                out_fm.mkdir(parents=True, exist_ok=True)
                p = out_fm / f"{slug(r['call'])}_{f:.1f}.json"
                json.dump(st, open(p, "w"), indent=1)
                written.append(p)
    if not a.fm_only:
        rows, recs = tv_rows, tv_recs
        by_rf = {}
        for r in recs:
            by_rf.setdefault(r["rf"], []).append(r)
        print(f"[fcc] TV: {len(rows)} rows -> {len(recs)} licensed full-power digital stations >= {a.min_erp_kw[1]:g} kW on {len(by_rf)} channels")
        for rf, rs in sorted(by_rf.items()):
            # the record to name the channel after: full power before class A, then ERP
            rs.sort(key=lambda r: (r["service"] not in ("DTV", "DT"), -strength(r, lat, lon)))
            r = rs[0]
            st = tv_station(r, atsc3, rs[1:])
            amb = ""
            if st.get("shares_channel_with"):
                amb += f"  shares with {', '.join(o.split(' (')[0] for o in st['shares_channel_with'])}"
            if st.get("ambiguous_with"):
                amb += f"  AMBIG with {', '.join(o.split(' ')[0] for o in st['ambiguous_with'])}"
            print(f"  RF{rf:<3d} {r['call']:10s} {r['erp_kw']:7.1f} kW  {r['city']}, {r['state']}  alt {st['tx_alt_m']:.0f} m  "
                  f"{'ATSC 3.0' if rf in atsc3 else '8-VSB'}{amb}")
            if not a.list:
                Path(a.out).mkdir(parents=True, exist_ok=True)
                p = Path(a.out) / f"{slug(r['call'].replace('-TV', '').replace('-DT', ''))}_{rf}.json"
                json.dump(st, open(p, "w"), indent=1)
                written.append(p)
    if written:
        print(f"[fcc] wrote {len(written)} station files under {a.out}")
        print("      these name your city's towers, not your receiver; the observer position stays in lab/site.json")


if __name__ == "__main__":
    main()
