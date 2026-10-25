"""fix_position.py -- a POSITION from one receiver, by using several
illuminators at once.

One illuminator gives an ellipsoid (constant bistatic range) and a
Doppler. Two or three illuminators in the SAME capture give two or three
ellipsoids through the same aircraft at the same instant, and ellipsoids
intersect at a point. An 8 MHz FM capture holds 40 channels, so three
stations on three different towers (e.g. 101.5 / 104.1 / 106.7 in the
102-108 MHz span) are three illuminators for free -- the multistatic fix
that otherwise needs a second SDR.

Input: N tracks files from fm_radar.py run on the same capture with N
different --station files (so their dwell times coincide). For every
dwell, every combination of one detection per station is solved by
Gauss-Newton for the point P that reproduces all N bistatic ranges; a
triplet is kept when every residual is inside one range cell and the
N Doppler values are consistent with one velocity (least-squares
velocity fit, residual under F_TOL). Fixes are linked over time into
tracks and written as lat/lon/alt so they can be compared DIRECTLY with
ADS-B positions -- no bistatic prediction step in between.

Controls, same as everywhere else: a fix must persist (>= MIN_PTS dwells),
move at a plausible speed, and the whole thing is run with one station's
detections shifted by one dwell (a null that must produce ~nothing).

Usage: python tools/fix_position.py lab/sweep_X/tracks_wbqb_101.5.json lab/sweep_X/tracks_wlno_104.1.json \
           lab/sweep_X/tracks_wjfk_fm_106.7.json [--truth meta.json] [--out fixes.json]
"""
import itertools
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402

C = 2.998e8
R_TOL_CELLS = 1.0        # per-station range residual allowed, in that station's cells
F_TOL = 6.0              # Hz, Doppler residual after a one-velocity fit
MIN_PTS = 5
ALT_MIN, ALT_MAX = 200.0, 14000.0
GRID_KM = 15.0           # initial guesses on a grid around the receiver
GRID_R = 90e3


def ecef_to_lla(p):
    x, y, z = p
    a, e2 = 6378137.0, 6.69437999014e-3
    lon = math.atan2(y, x)
    r = math.hypot(x, y)
    lat = math.atan2(z, r * (1 - e2))
    for _ in range(6):
        n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = r / math.cos(lat) - n
        lat = math.atan2(z, r * (1 - e2 * n / (n + h)))
    n = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    h = r / math.cos(lat) - n
    return math.degrees(lat), math.degrees(lon), h


def load_station_tracks(path):
    d = json.load(open(path))
    st = d["station"]
    T = tc.ecef(st["tx_lat"], st["tx_lon"], st.get("tx_alt_m", 200.0))
    fs = d["fs_hz"]
    lam = C / d["fc_hz"]
    by_t = {}
    for x in d["detections"]:
        by_t.setdefault(round(x["t_s"], 3), []).append(x)
    return dict(name=st.get("short_name", Path(path).stem), T=T, cell=C / fs, lam=lam, by_t=by_t, path=path)


def solve_point(stations, ranges, R, P0, iters=12):
    """Gauss-Newton for P given bistatic excess ranges to each station."""
    P = P0.astype(float).copy()
    for _ in range(iters):
        f = np.zeros(len(stations))
        J = np.zeros((len(stations), 3))
        for i, s in enumerate(stations):
            dT = P - s["T"]; dR = P - R
            nT, nR = np.linalg.norm(dT), np.linalg.norm(dR)
            f[i] = nT + nR - np.linalg.norm(s["T"] - R) - ranges[i]
            J[i] = dT / nT + dR / nR
        # regularize toward a plausible altitude when under-determined (2 stations)
        if len(stations) < 3:
            up = R / np.linalg.norm(R)
            h = float(np.dot(P - R, up))
            target = min(max(h, ALT_MIN + 1500), ALT_MAX - 1500)
            f = np.append(f, 0.05 * (h - target)); J = np.vstack([J, 0.05 * up])
        step, *_ = np.linalg.lstsq(J, -f, rcond=None)
        P += step
        if np.linalg.norm(step) < 0.5:
            break
    res = []
    for i, s in enumerate(stations):
        res.append(np.linalg.norm(P - s["T"]) + np.linalg.norm(P - R) - np.linalg.norm(s["T"] - R) - ranges[i])
    return P, np.array(res)


def doppler_velocity(stations, P, R, dops):
    """Least-squares velocity from N bistatic Dopplers at P; returns (v, residual_hz)."""
    A, b = [], []
    for s, fd in zip(stations, dops):
        uT = (P - s["T"]) / np.linalg.norm(P - s["T"])
        uR = (P - R) / np.linalg.norm(P - R)
        A.append(-(uT + uR) / s["lam"]); b.append(fd)      # fd = -(d/dt)(|P-T|+|P-R|)/lam
    A, b = np.array(A), np.array(b)
    if len(A) < 3:
        # constrain vertical velocity ~ 0
        up = R / np.linalg.norm(R)
        A = np.vstack([A, 0.02 * up]); b = np.append(b, 0.0)
    v, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = A[:len(dops)] @ v - b[:len(dops)]
    return v, float(np.sqrt(np.mean(res ** 2)))


def main():
    tc.station_arg()
    argv = sys.argv[1:]
    args, opt, i = [], {}, 0
    while i < len(argv):                      # accepts --key=value and --key value
        a = argv[i]
        if a.startswith("--") and "=" in a:
            k, v = a.split("=", 1); opt[k] = v
        elif a.startswith("--") and i + 1 < len(argv) and not argv[i + 1].startswith("--"):
            opt[a] = argv[i + 1]; i += 1
        elif not a.startswith("--"):
            args.append(a)
        i += 1
    shift = int(opt.get("--null-shift", "0"))          # dwells to shift station 0 by (the null)
    stations = [load_station_tracks(p) for p in args]
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    up = R / np.linalg.norm(R)
    east = np.cross([0, 0, 1], up); east /= np.linalg.norm(east)
    north = np.cross(up, east)
    grid = [R + up * 6000 + east * gx + north * gy
            for gx in np.arange(-GRID_R, GRID_R + 1, GRID_KM * 1e3)
            for gy in np.arange(-GRID_R, GRID_R + 1, GRID_KM * 1e3)]
    times = sorted(set.intersection(*[set(s["by_t"]) for s in stations]))
    if shift:
        keys = sorted(stations[0]["by_t"])
        stations[0]["by_t"] = {keys[(i + shift) % len(keys)]: stations[0]["by_t"][k] for i, k in enumerate(keys)}
        times = sorted(set.intersection(*[set(s["by_t"]) for s in stations]))
    print(f"fix_position -- {len(stations)} illuminators: " + ", ".join(s["name"] for s in stations))
    print(f"  common dwells: {len(times)}; null shift {shift} dwells")
    fixes = []
    n_tried = 0
    for t in times:
        lists = [s["by_t"][t] for s in stations]
        for combo in itertools.product(*lists):
            n_tried += 1
            ranges = [d["range_m"] for d in combo]
            dops = [d["doppler_hz"] for d in combo]
            best = None
            for P0 in grid:
                P, res = solve_point(stations, ranges, R, P0)
                if np.all(np.abs(res) < R_TOL_CELLS * np.array([s["cell"] for s in stations])):
                    h = float(np.dot(P - R, up))
                    if not (ALT_MIN < h < ALT_MAX):
                        continue
                    v, fres = doppler_velocity(stations, P, R, dops)
                    if fres > F_TOL or np.linalg.norm(v) > 350:
                        continue
                    score = float(np.sum(res ** 2)) + 1e4 * fres ** 2
                    if best is None or score < best[0]:
                        best = (score, P, res, v, fres)
            if best:
                _, P, res, v, fres = best
                la, lo, h = ecef_to_lla(P)
                fixes.append(dict(t_s=t, lat=la, lon=lo, alt_m=h, speed_ms=float(np.linalg.norm(v)),
                                  range_res_m=[float(r) for r in res], doppler_res_hz=fres,
                                  snr_db=[d["snr_db"] for d in combo]))
    # link fixes into tracks (same aircraft: within 2 km per second of gap)
    fixes.sort(key=lambda f: f["t_s"])
    tracks = []
    for f in fixes:
        P = tc.ecef(f["lat"], f["lon"], f["alt_m"])
        for tr in tracks:
            g = tr[-1]
            dt = f["t_s"] - g["t_s"]
            if 0 < dt <= 6 and np.linalg.norm(P - tc.ecef(g["lat"], g["lon"], g["alt_m"])) < 350 * dt + 1500:
                tr.append(f); break
        else:
            tracks.append([f])
    tracks = [tr for tr in tracks if len(tr) >= MIN_PTS]
    print(f"  combinations tried {n_tried}, fixes {len(fixes)}, fix tracks (>= {MIN_PTS} dwells) {len(tracks)}")
    for i, tr in enumerate(tracks[:12]):
        a, b = tr[0], tr[-1]
        dist = np.linalg.norm(tc.ecef(a["lat"], a["lon"], a["alt_m"]) - tc.ecef(b["lat"], b["lon"], b["alt_m"]))
        print(f"  fix track {i}: {len(tr)} dwells t={a['t_s']:.1f}-{b['t_s']:.1f}s  "
              f"alt {np.median([f['alt_m'] for f in tr]):6.0f} m  speed {np.median([f['speed_ms'] for f in tr]):5.0f} m/s  "
              f"moved {dist/1e3:5.1f} km in {b['t_s']-a['t_s']:.0f} s  (lat/lon withheld here; in the json)")
    out = opt.get("--out") or Path(args[0]).with_name("fixes.json")
    json.dump(dict(stations=[s["name"] for s in stations], null_shift=shift, fixes=fixes,
                   tracks=[[dict(f) for f in tr] for tr in tracks]), open(out, "w"), indent=1)
    print(f"-> {out}")
    # truth comparison, position against position
    if "--truth" in opt:
        meta = json.load(open(opt["--truth"]))
        truth = tc.load_truth(meta["truth"])
        tl = float(meta["t_launch_unix"])
        hits = 0; errs = []
        for tr in tracks:
            best_ac, best_e = None, 1e12
            for icao, samples in truth.items():
                e = []
                for f in tr:
                    p = tc.interp(samples, tl + f["t_s"])
                    if p is None:
                        e = None; break
                    e.append(np.linalg.norm(p - tc.ecef(f["lat"], f["lon"], f["alt_m"])))
                if e and np.mean(e) < best_e:
                    best_ac, best_e = (icao, samples[0][7]), float(np.mean(e))
            if best_ac and best_e < 5000:
                hits += 1; errs.append(best_e)
                print(f"  fix track matches {best_ac[1] or best_ac[0]}: mean position error {best_e/1e3:.2f} km over {len(tr)} dwells")
        print(f"  TRUTH: {hits} of {len(tracks)} fix tracks within 5 km of an ADS-B aircraft; "
              f"median error {np.median(errs)/1e3 if errs else float('nan'):.2f} km")


if __name__ == "__main__":
    main()
