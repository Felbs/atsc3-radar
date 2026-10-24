"""unmatched_tracks.py -- the matcher inverted: tracks that NO ADS-B aircraft
explains.

A positive match is strong because the truth is specific. An unmatched
track rests on an absence, so it is held to a higher bar, all of it
written here before the first run (2026-10-24):

  candidate = a track_hunt track with
    - every point >= 5 km excess range and outside the +-15 Hz notch
    - >= MIN_PTS points
    - NO ADS-B prediction within a DOUBLED gate (2 x 490 m, 60 Hz) at
      any of its points (truth includes MLAT-positioned aircraft when
      OpenSky has them; position_source counts are printed)
    - walk rate consistent with its own Doppler (|obs - pred| <= 25 % of
      pred, or both under 10 m/s)
    - not a DOPPLER TWIN of a matched aircraft: a second-path echo of a
      real aircraft sits at a longer range with the same Doppler. A
      candidate whose Doppler is within TWIN_F of any matched prediction
      at the same dwells, 0-4 km beyond it, is a ghost suspect.
    - not at a fixed cell listed by gain_pair.py (if a *_gainpair.json
      is given)
  and it is reported with the capture's own miss rate beside it, because
  a radar that misses 80 % of known aircraft produces unmatched tracks
  for free.

Cross-window: candidates in different windows are linked when one
extrapolates to the other (range within 1 km after walking at its own
rate, Doppler within 30 Hz). A candidate seen in >= 2 windows >= 60 s
apart is the only grade worth a sentence.

Usage: python tools/unmatched_tracks.py <meta.json> <tracks_t5.json> <tracks_t30.json> <tracks_t90.json> [--gainpair=<json>]
"""
import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402

MIN_PTS = 8
R_MIN = 5e3
NOTCH = 15.0
R_GATE2 = 2 * (3 * tc.BIN_M + 360)
F_GATE2 = 60.0
TWIN_F = 10.0
TWIN_R = (0.0, 4e3)


def load_truth_with_source(path):
    """icao -> samples, plus a Counter of position_source."""
    truth = tc.load_truth(path)
    src = Counter()
    for line in open(path):
        row = json.loads(line)
        for s in row.get("states", []):
            src[s[16] if len(s) > 16 else None] += 1
    return truth, src


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    gp = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--gainpair=")]
    meta = json.load(open(args[0]))
    truth_path = meta["truth"]
    if not os.path.exists(truth_path):
        truth_path = os.path.join(os.path.dirname(os.path.abspath(args[0])),
                                  os.path.basename(truth_path))
    truth, src = load_truth_with_source(truth_path)
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    T = tc.ecef(tc.TX_LAT, tc.TX_LON, tc.TX_ANT_AMSL_M)
    dRT = float(np.linalg.norm(T - R))
    t_launch = float(meta["t_launch_unix"])
    fixed = []
    if gp:
        fixed = json.load(open(gp[0]))["common_cells"]

    print("=" * 76)
    print(f"unmatched_tracks -- {Path(args[0]).stem}  rfgain {meta.get('rfgain')}")
    print(f"  truth: {len(truth)} aircraft; position_source "
          f"{ {str(k): v for k, v in src.items()} }  (0 ADS-B, 1 ASTERIX, 2 MLAT, 3 FLARM)")
    print(f"  rules: R_b >= {R_MIN / 1e3:.0f} km, |f| > {NOTCH:.0f} Hz, >= {MIN_PTS} pts, "
          f"no truth within {R_GATE2:.0f} m / {F_GATE2:.0f} Hz, walk within 25 %, "
          f"no Doppler twin ({TWIN_F:.0f} Hz, 0-{TWIN_R[1] / 1e3:.0f} km beyond)")
    print("=" * 76)

    all_cands = []
    for path in args[1:]:
        trk = json.load(open(path))
        w = Path(path).stem
        det = trk["detections"]
        times = sorted({d["t_s"] for d in det})
        # predictions per dwell, every aircraft, wide map (<= 70 km)
        pred = {ts: [] for ts in times}
        for icao, samples in truth.items():
            cs = samples[0][7] or icao
            for ts in times:
                pr = tc.predict(samples, t_launch + ts, T, R, dRT)
                if pr and 0 < pr[0] < 70e3 and abs(pr[1]) < 1200:
                    pred[ts].append((icao, cs, pr[0], pr[1]))
        # matched detections (normal gate) for the twin test and the miss rate
        R_GATE = 3 * tc.BIN_M + 360
        matched = []
        for d in det:
            for icao, cs, rb, fd in pred[d["t_s"]]:
                if abs(d["range_m"] - rb) < R_GATE and abs(d["doppler_hz"] - fd) < 30:
                    matched.append((d["t_s"], rb, fd, cs))
                    break
        in_map = {icao for ts in times for icao, _, rb, fd in pred[ts] if rb < 60e3 and abs(fd) < 1000}
        seen = {cs for _, _, _, cs in matched}
        print(f"\n[{w}] {len(det)} detections, {len(matched)} matched; "
              f"{len(in_map)} aircraft in map, {len(seen)} seen "
              f"(miss rate {1 - len(seen) / max(len(in_map), 1):.0%}); {len(trk['tracks'])} tracks")

        for ti, t in enumerate(trk["tracks"]):
            pts = t["points"]
            n = len(pts)
            why = []
            if n < MIN_PTS:
                why.append(f"{n} pts")
            if min(p["range_m"] for p in pts) < R_MIN:
                why.append("inside 5 km")
            if min(abs(p["doppler_hz"]) for p in pts) <= NOTCH:
                why.append("in notch")
            # truth within the doubled gate at any point
            expl = Counter()
            for p in pts:
                for icao, cs, rb, fd in pred.get(p["t_s"], []):
                    if abs(p["range_m"] - rb) < R_GATE2 and abs(p["doppler_hz"] - fd) < F_GATE2:
                        expl[cs] += 1
            if expl:
                top = expl.most_common(1)[0]
                why.append(f"ADS-B {top[0]} explains {top[1]}/{n} pts")
            # walk consistency
            obs, prd = t["observed_rate_ms"], t["predicted_rate_ms"]
            if not (abs(obs - prd) <= 0.25 * abs(prd) or (abs(obs) < 10 and abs(prd) < 10)):
                why.append(f"walk {obs:+.0f} vs {prd:+.0f} m/s")
            # Doppler twin of a matched aircraft
            twins = Counter()
            for p in pts:
                for ts, rb, fd, cs in matched:
                    if abs(ts - p["t_s"]) < 1e-3 and abs(p["doppler_hz"] - fd) < TWIN_F \
                            and TWIN_R[0] <= p["range_m"] - rb <= TWIN_R[1]:
                        twins[cs] += 1
            if twins and twins.most_common(1)[0][1] >= 0.5 * n:
                why.append(f"Doppler twin of {twins.most_common(1)[0][0]}")
            # fixed cell
            for c in fixed:
                if abs(c["r_m"] - t["range0_m"]) < 2 * tc.BIN_M and abs(c["f_hz"] - t["doppler_mean_hz"]) < 6:
                    why.append("fixed cell (gain pair)")
            cand = dict(window=w, idx=ti, n=n, t0=t["t0"], t1=t["t1"],
                        r0=t["range0_m"], r1=t["range1_m"], fd=t["doppler_mean_hz"],
                        obs=obs, prd=prd, snr=max(t["snr_db"]), why=why)
            tag = "CANDIDATE" if not why else "no"
            if not why or (n >= MIN_PTS and min(p["range_m"] for p in pts) >= R_MIN):
                print(f"  track {ti:2d}: {n:3d} pts t={t['t0']:6.2f}-{t['t1']:6.2f}s  "
                      f"{t['range0_m'] / 1e3:6.2f}->{t['range1_m'] / 1e3:6.2f} km {t['doppler_mean_hz']:+6.0f} Hz  "
                      f"walk {obs:+5.0f}/{prd:+5.0f} m/s  {max(t['snr_db']):4.1f} dB  "
                      f"{tag:9s} {'; '.join(why)}")
            if not why:
                all_cands.append(cand)

    # cross-window linking
    print(f"\nCANDIDATES (all rules): {len(all_cands)}")
    links = []
    for i, a in enumerate(all_cands):
        for b in all_cands[i + 1:]:
            if a["window"] == b["window"]:
                continue
            dt = b["t0"] - a["t1"]
            if dt <= 0:
                a, b = b, a
                dt = b["t0"] - a["t1"]
            r_ext = a["r1"] + a["obs"] * dt
            if abs(r_ext - b["r0"]) < 1e3 and abs(a["fd"] - b["fd"]) < 30:
                links.append((a, b, dt))
    for a, b, dt in links:
        print(f"  LINKED across windows ({dt:.0f} s apart): {a['window']} #{a['idx']} -> {b['window']} #{b['idx']}  "
              f"{a['r0'] / 1e3:.2f}->{b['r1'] / 1e3:.2f} km, {a['fd']:+.0f}/{b['fd']:+.0f} Hz")
    if not links:
        print("  none persist across windows")
    out = Path(args[1]).with_name(f"unmatched_{Path(args[0]).stem.replace('.meta', '')}.json")
    json.dump(dict(candidates=all_cands, links=[(a["window"], a["idx"], b["window"], b["idx"]) for a, b, _ in links],
                   position_source={str(k): v for k, v in src.items()}), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
