"""truth_trajectory.py -- the mandatory trajectory test, across windows.

The 8/21 confirming protocol: a matched aircraft's detection must MOVE
with its own ADS-B prediction, window to window, and the cell it leaves
must empty. Per-dwell match counts cannot show that; this lists, per
aircraft and per window, every matched detection against its prediction,
the within-window walk of both, and which aircraft are tracked across
>= 2 windows (the only way an R_b < 5 km match earns a vote).

Usage: python tools/truth_trajectory.py <meta.json> <tracks_A.json> [tracks_B.json ...]
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402


def main():
    meta = json.load(open(sys.argv[1]))
    truth_path = meta["truth"]
    if not os.path.exists(truth_path):
        truth_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])),
                                  os.path.basename(truth_path))
    truth = tc.load_truth(truth_path)
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    T = tc.ecef(tc.TX_LAT, tc.TX_LON, tc.TX_ANT_AMSL_M)
    dRT = float(np.linalg.norm(T - R))
    t_launch = float(meta["t_launch_unix"])
    R_GATE = 3 * tc.BIN_M + 180 * 2.0
    F_GATE = 30.0

    per_ac = {}          # icao -> window -> list of rows
    cs_of = {}
    windows = []
    for path in sys.argv[2:]:
        trk = json.load(open(path))
        det = trk["detections"]
        wname = Path(path).stem
        windows.append(wname)
        times = sorted({d["t_s"] for d in det})
        for d in det:
            best = None
            for icao, samples in truth.items():
                pr = tc.predict(samples, t_launch + d["t_s"], T, R, dRT)
                if pr is None:
                    continue
                rb, fd = pr
                if abs(d["range_m"] - rb) < R_GATE and abs(d["doppler_hz"] - fd) < F_GATE:
                    err = abs(d["range_m"] - rb) / R_GATE + abs(d["doppler_hz"] - fd) / F_GATE
                    if best is None or err < best[0]:
                        best = (err, icao, rb, fd, samples[0][7])
            if best:
                _, icao, rb, fd, cs = best
                cs_of[icao] = cs
                per_ac.setdefault(icao, {}).setdefault(wname, []).append(
                    (d["t_s"], d["range_m"], d["doppler_hz"], d.get("snr_db", float("nan")), rb, fd))

    print("=" * 78)
    print("truth_trajectory -- matched detections per aircraft, per window")
    print(f"  windows: {', '.join(windows)}   gate |dR| < {R_GATE:.0f} m, |df| < {F_GATE:.0f} Hz")
    print("=" * 78)
    multi = []
    order = sorted(per_ac, key=lambda k: -sum(len(v) for v in per_ac[k].values()))
    for icao in order:
        wins = per_ac[icao]
        n = sum(len(v) for v in wins.values())
        tag = "TRACKED ACROSS %d WINDOWS" % len(wins) if len(wins) >= 2 else "one window"
        if len(wins) >= 2:
            multi.append(icao)
        print(f"\n{icao} {cs_of[icao]:8s} {n} matched dwells  [{tag}]")
        for w in windows:
            rows = sorted(wins.get(w, []))
            if not rows:
                # where was it predicted during this window, if in map?
                continue
            t0, t1 = rows[0][0], rows[-1][0]
            dr_det = rows[-1][1] - rows[0][1]
            dr_pre = rows[-1][4] - rows[0][4]
            snr = [r[3] for r in rows]
            print(f"  {w:16s} {len(rows):3d} dwells t={t0:6.2f}..{t1:6.2f}s  "
                  f"det {rows[0][1]/1e3:6.2f}->{rows[-1][1]/1e3:6.2f} km {np.mean([r[2] for r in rows]):+6.0f} Hz  "
                  f"pred {rows[0][4]/1e3:6.2f}->{rows[-1][4]/1e3:6.2f} km {np.mean([r[5] for r in rows]):+6.0f} Hz  "
                  f"snr {np.nanmin(snr):4.1f}-{np.nanmax(snr):4.1f} dB")
            if len(rows) >= 3 and t1 > t0:
                print(f"  {'':16s}     walk det {dr_det/(t1-t0):+6.1f} m/s  pred {dr_pre/(t1-t0):+6.1f} m/s  "
                      f"mean |dR| {np.mean([abs(r[1]-r[4]) for r in rows]):5.0f} m  "
                      f"mean |df| {np.mean([abs(r[2]-r[5]) for r in rows]):4.1f} Hz")
        # the vacated-cell test: in windows where it was NOT matched, was it predicted in map?
        for w in windows:
            if w in wins:
                continue
            path = [p for p in sys.argv[2:] if Path(p).stem == w][0]
            det = json.load(open(path))["detections"]
            ts = sorted({d["t_s"] for d in det})
            prs = [tc.predict(truth[icao], t_launch + t, T, R, dRT) for t in ts]
            prs = [p for p in prs if p and 0 < p[0] < 60e3 and abs(p[1]) < 1000]
            if prs:
                print(f"  {w:16s}   not matched; predicted in map at "
                      f"{min(p[0] for p in prs)/1e3:.2f}-{max(p[0] for p in prs)/1e3:.2f} km "
                      f"{min(p[1] for p in prs):+.0f}..{max(p[1] for p in prs):+.0f} Hz")
            else:
                print(f"  {w:16s}   not matched; not predicted in map (outside 60 km / 1 kHz or no truth)")

    print(f"\nAircraft tracked across >= 2 windows: {len(multi)}: "
          + ", ".join(f"{i} {cs_of[i]}" for i in multi))
    out = {i: {w: len(v) for w, v in per_ac[i].items()} for i in per_ac}
    outp = Path(sys.argv[2]).with_name("truth_trajectory.json")
    json.dump(dict(windows=windows, per_aircraft=out, tracked_multi=multi,
                   callsigns=cs_of), open(outp, "w"), indent=1)
    print(f"-> {outp}")


if __name__ == "__main__":
    main()
