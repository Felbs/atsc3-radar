"""adsb_history.py -- ADS-B truth after the fact, from adsb.lol's daily history.

OpenSky's anonymous API has a daily credit cap (we hit it 2026-10-25
00:28Z, mid-session) and 10 s resolution. adsb.lol publishes every
aircraft's full trace for each UTC day as a GitHub release
(adsblol/globe_history_<year>), ~3.5 h after the day ends: a 'prod'
release (~4 GB, split tar) and an 'mlatonly' release (~0.3 GB) of
aircraft positioned only by multilateration, i.e. NOT broadcasting
ADS-B positions -- the population the unmatched-track goal is about.

This turns the traces inside a lat/lon box and a time window into the
same .adsb.jsonl rows truth_capture writes (OpenSky state-vector order),
so truth_compare / truth_trajectory / unmatched_tracks run unchanged.
position_source is 0 for prod and 2 for mlatonly.

Downloads are untrusted data: each release goes into its own empty
directory and only traces/*.json(.gz) members are read, never extracted
as files.

Usage:
  python tools/adsb_history.py fetch 2026-10-24 [--kind prod|mlatonly]      # download the split tar
  python tools/adsb_history.py truth <meta.json> [<meta.json> ...] [--pad 300] [--kinds prod,mlatonly]
        -> writes <meta stem>.hist.adsb.jsonl and points a copy of the meta at it
"""
import gzip
import io
import json
import math
import os
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent / "lab"
HIST = LAB / "adsb_history"
FT = 0.3048
KT = 0.514444
FPM = 0.00508

sys.path.insert(0, str(HERE))


def release_tag(day, kind):
    return f"v{day.replace('-', '.')}-planes-readsb-{kind}-0"


def fetch(day, kind="prod"):
    year = day[:4]
    tag = release_tag(day, kind)
    d = HIST / tag
    d.mkdir(parents=True, exist_ok=True)
    print(f"[hist] downloading {tag} -> {d}", flush=True)
    r = subprocess.run(["gh", "release", "download", tag, "-R", f"adsblol/globe_history_{year}",
                        "-D", str(d), "--skip-existing"], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-800:])
        return None
    parts = sorted(p for p in d.iterdir() if p.name.startswith(tag) and ".tar" in p.name)
    print(f"[hist] {len(parts)} part(s): " + ", ".join(f"{p.name} {p.stat().st_size/1e6:.0f} MB" for p in parts))
    return parts


class Concat(io.RawIOBase):
    """Stream several split-tar parts as one file."""
    def __init__(self, paths):
        self.files = [open(p, "rb") for p in paths]
        self.i = 0

    def readable(self):
        return True

    def readinto(self, b):
        while self.i < len(self.files):
            n = self.files[self.i].readinto(b)
            if n:
                return n
            self.i += 1
        return 0


def iter_traces(parts):
    """Yield (icao, trace_json) for every traces member of the release."""
    raw = io.BufferedReader(Concat(parts), buffer_size=1 << 22)
    with tarfile.open(fileobj=raw, mode="r|*") as tf:
        for m in tf:
            if not m.isfile() or "trace_full_" not in m.name:
                continue
            f = tf.extractfile(m)
            if f is None:
                continue
            data = f.read()
            try:
                if data[:2] == b"\x1f\x8b":
                    data = gzip.decompress(data)
                yield m.name, json.loads(data)
            except Exception:
                continue


def rows_from_trace(t, box, t0, t1, source):
    """OpenSky-ordered state vectors for the trace points inside box and window."""
    lamin, lamax, lomin, lomax = box
    base = float(t.get("timestamp", 0))
    icao = t.get("icao", "")
    cs = ""
    out = []
    for p in t.get("trace", []):
        ts = base + float(p[0])
        if ts < t0 or ts > t1:
            continue
        lat, lon = p[1], p[2]
        if lat is None or lon is None or not (lamin <= lat <= lamax and lomin <= lon <= lomax):
            continue
        alt = p[3]
        on_ground = alt == "ground"
        if len(p) > 8 and isinstance(p[8], dict) and p[8].get("flight"):
            cs = p[8]["flight"].strip()
        geo = p[10] if len(p) > 10 and isinstance(p[10], (int, float)) else None
        baro = None if on_ground or alt is None else float(alt) * FT
        gs = p[4] * KT if isinstance(p[4], (int, float)) else None
        vr = p[7] * FPM if len(p) > 7 and isinstance(p[7], (int, float)) else None
        out.append([icao, cs, "", ts, ts, lon, lat, baro, on_ground, gs,
                    p[5] if isinstance(p[5], (int, float)) else None, vr, None,
                    geo * FT if geo is not None else None, None, False, source])
    return out


def truth(meta_paths, pad=300.0, kinds=("prod", "mlatonly")):
    """Truth for one or many captures; each day's release is read once for all of them."""
    import truth_compare as tc
    import datetime as dt
    if isinstance(meta_paths, (str, Path)):
        meta_paths = [meta_paths]
    lat, lon, alt = tc.observer()
    blat, blon = 0.6, 0.75
    import radar_site
    box = radar_site.truth_box(lat, lon, blat, blon)          # receiver and every known tower
    jobs = []
    for mp in meta_paths:
        meta = json.load(open(mp))
        t0 = float(meta["t_launch_unix"]) - pad
        t1 = float(meta.get("t_done_unix", meta["t_launch_unix"] + meta["secs"])) + pad
        day = dt.datetime.fromtimestamp(t0, dt.timezone.utc).strftime("%Y-%m-%d")
        jobs.append(dict(path=mp, meta=meta, t0=t0, t1=t1, day=day, rows=[], n_ac={}))
    for day in sorted({j["day"] for j in jobs}):
        J = [j for j in jobs if j["day"] == day]
        for kind in kinds:
            d = HIST / release_tag(day, kind)
            parts = sorted(p for p in d.glob(f"{release_tag(day, kind)}*")) if d.exists() else []
            if not parts:
                print(f"[hist] {kind} {day}: not downloaded (python tools/adsb_history.py fetch {day} --kind {kind})")
                continue
            src = 2 if kind == "mlatonly" else 0
            for j in J:
                j["n_ac"][kind] = 0
            for name, t in iter_traces(parts):
                for j in J:
                    r = rows_from_trace(t, box, j["t0"], j["t1"], src)
                    if r:
                        j["rows"] += r
                        j["n_ac"][kind] += 1
            print(f"[hist] {kind} {day}: read for {len(J)} capture(s)", flush=True)
    return [write_truth(j) for j in jobs]


def write_truth(j):
    meta_path, meta, rows, n_ac, day = j["path"], j["meta"], j["rows"], j["n_ac"], j["day"]
    rows.sort(key=lambda r: r[3])
    stem = Path(meta_path).name.replace(".meta.json", "")
    out = Path(meta_path).with_name(stem + ".hist.adsb.jsonl")
    with open(out, "w") as f:
        # one 'poll' per second so the file has the same shape as truth_capture's
        cur, buf = None, []
        for r in rows:
            sec = int(r[3])
            if cur is not None and sec != cur:
                f.write(json.dumps(dict(poll_utc=cur, opensky_time=cur, n=len(buf), states=buf,
                                        source="adsb.lol history")) + "\n")
                buf = []
            cur = sec
            buf.append(r)
        if buf:
            f.write(json.dumps(dict(poll_utc=cur, opensky_time=cur, n=len(buf), states=buf,
                                    source="adsb.lol history")) + "\n")
    m2 = dict(meta)
    m2["truth"] = str(out)
    m2["truth_source"] = f"adsb.lol globe_history {day} ({', '.join(f'{k}: {v}' for k, v in n_ac.items())})"
    m2_path = Path(meta_path).with_name(stem + ".hist.meta.json")
    json.dump(m2, open(m2_path, "w"), indent=1)
    print(f"[hist] {Path(meta_path).name}: {sum(n_ac.values())} aircraft, {len(rows)} trace points -> {out.name}")
    return m2_path


if __name__ == "__main__":
    import truth_compare as tc
    tc.station_arg()
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opt = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    if not args:
        print(__doc__)
    elif args[0] == "fetch":
        fetch(args[1], opt.get("--kind", "prod"))
    elif args[0] == "truth":
        truth(args[1:], float(opt.get("--pad", 300)), tuple(opt.get("--kinds", "prod,mlatonly").split(",")))
