"""fuse_fix.py -- multistatic position fix by fusing range-Doppler maps BEFORE detection.

fix_position.py intersects detections: every station must pass its own CFAR
threshold on the same dwell, and a weak station (one that cancels 22 dB
instead of 31) throws most aircraft away. Here nothing is thresholded per
station. For every point P of a 3-D grid around the receiver the bistatic
range to each illuminator is fixed geometry, so each station contributes the
strongest SNR it has at that range (over all Doppler outside the notch);
the per-station profiles are robustly normalised, summed over stations and
over W consecutive dwells (an aircraft moves < 1 km in 3 s), and the grid is
searched for peaks.

Verification, in order:
  1. Null: the same search with ONE station's maps shifted by 60 s (120
     dwells, circular), for each station in turn. Shifting by one dwell is
     not a null -- an aircraft moves 100 m in 0.5 s, inside one 1.2 km cell.
     The detection threshold is the highest peak score any null produces,
     so the nulls give zero detections by construction and the real count
     above it is the result.
  2. Doppler: at a peak, the three Doppler profiles at the three ranges
     must be explained by ONE horizontal velocity (residual < F_TOL Hz).
     Range only enters the score, so this is an independent check.
  3. Persistence: peaks linked over windows into tracks that move at an
     aircraft speed.

Positions go to the json (lab/, gitignored); stdout prints counts and
errors only.

Usage: python tools/fuse_fix.py <tracks_A.json> <tracks_B.json> <tracks_C.json> [--window 6]
       (tracks files from fm_radar.py --save-maps on the same capture)
"""
import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402
from fix_position import doppler_velocity, ecef_to_lla         # noqa: E402

C = 2.998e8
NULL_SHIFT_S = 60.0
F_TOL = 3.0          # Hz; Doppler cells are 0.59 Hz
NOTCH = 5.0
MIN_BIN = 3       # first 3.6 km of bistatic range excluded, as the single-station detector does


def load_station(path):
    d = json.load(open(path))
    st = d["station"]
    maps_dir = Path(path).parent / "maps" / Path(path).stem
    rax, _, dax = np.load(maps_dir / "axes.npy", allow_pickle=True)
    files = sorted(glob.glob(str(maps_dir / "rd_*.npy")))
    t = np.array([float(Path(f).stem[3:]) for f in files])
    maps = np.stack([np.load(f).astype(np.float32) for f in files])   # (n_t, n_r, n_d) dB over CFAR noise
    keep = np.abs(np.asarray(dax, float)) > NOTCH
    prof = maps[:, :, keep].max(axis=2)                                # (n_t, n_r) best SNR at each range
    # robust per-dwell normalisation so a strong station does not dominate
    med = np.median(prof, axis=1, keepdims=True)
    mad = np.median(np.abs(prof - med), axis=1, keepdims=True) * 1.4826 + 1e-6
    z = (prof - med) / mad
    # remove what never changes: residual direct path / static clutter sits in the
    # first range bins in every dwell and, summed over three stations, outscored
    # every aircraft (and every null) at grid points near the receiver
    z = z - np.median(z, axis=0, keepdims=True)
    z[:, :MIN_BIN] = -10.0
    return dict(name=st.get("short_name", Path(path).stem), T=tc.ecef(st["tx_lat"], st["tx_lon"], st.get("tx_alt_m", 200.0)),
                lam=C / d["fc_hz"], rax=np.asarray(rax, float), dax=np.asarray(dax, float), t=t, maps=maps, z=z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tracks", nargs="+")
    ap.add_argument("--window", type=int, default=6, help="dwells summed per score (hop 0.5 s)")
    ap.add_argument("--extent-km", type=float, default=100.0)
    ap.add_argument("--step-km", type=float, default=1.0)
    ap.add_argument("--alts", default="1000,2500,4000,6000,8500,11000")
    ap.add_argument("--peaks", type=int, default=8, help="peaks examined per window")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    S = [load_station(p) for p in a.tracks]
    n_t = min(len(s["t"]) for s in S)
    for s in S:
        assert np.allclose(s["t"][:n_t], S[0]["t"][:n_t]), "dwell times differ: not the same capture/settings"
    t = S[0]["t"][:n_t]
    hop = float(np.median(np.diff(t)))
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    up = R / np.linalg.norm(R)
    east = np.cross([0, 0, 1], up); east /= np.linalg.norm(east)
    north = np.cross(up, east)
    xs = np.arange(-a.extent_km, a.extent_km + 1e-9, a.step_km) * 1e3
    alts = np.array([float(v) for v in a.alts.split(",")])
    X, Y, H = np.meshgrid(xs, xs, alts, indexing="ij")
    G = R + X[..., None] * east + Y[..., None] * north + H[..., None] * up     # (nx, ny, na, 3)
    # fixed geometry: fractional range bin of every grid point for every station
    idx = []
    for s in S:
        rb = np.linalg.norm(G - s["T"], axis=-1) + np.linalg.norm(G - R, axis=-1) - np.linalg.norm(s["T"] - R)
        cell = s["rax"][1] - s["rax"][0]
        f = (rb - s["rax"][0]) / cell
        lo = np.clip(np.floor(f).astype(int), 0, len(s["rax"]) - 2)
        valid = (f >= 0) & (f <= len(s["rax"]) - 1)
        idx.append((lo, valid, rb))
    print(f"fuse_fix -- {len(S)} illuminators: " + ", ".join(s["name"] for s in S))
    print(f"  {n_t} dwells, hop {hop:.2f} s, window {a.window} dwells; grid {len(xs)}x{len(xs)}x{len(alts)} "
          f"= {X.size} points at {a.step_km:.1f} km")

    def scores(shifts):
        """per-window peak list for given per-station dwell shifts"""
        zs = [np.roll(s["z"][:n_t], sh, axis=0) for s, sh in zip(S, shifts)]
        out = []
        for w0 in range(0, n_t - a.window + 1, a.window // 2):
            tot = np.zeros(X.shape, np.float32)
            for zi, (lo, valid, _) in zip(zs, idx):
                prof = zi[w0:w0 + a.window].sum(axis=0)                      # integrate over the window
                v = np.maximum(prof[lo], prof[lo + 1])
                tot += np.where(valid, v, -10.0)
            # local maxima over x, y (any altitude), strongest first
            flat = tot.max(axis=2)
            from scipy.ndimage import maximum_filter
            pk = (flat == maximum_filter(flat, size=5)) & np.isfinite(flat)
            ii = np.argwhere(pk)
            vals = flat[pk]
            order = np.argsort(-vals)[:a.peaks]
            for o in order:
                ix, iy = ii[o]
                ia = int(np.argmax(tot[ix, iy]))
                out.append(dict(w0=int(w0), score=float(vals[o]), ix=int(ix), iy=int(iy), ia=ia))
        return out

    real = scores([0] * len(S))
    sh = int(round(NULL_SHIFT_S / hop))
    nulls = {}
    for k in range(len(S)):
        shifts = [0] * len(S); shifts[k] = sh
        nulls[S[k]["name"]] = scores(shifts)
    thr = max(p["score"] for v in nulls.values() for p in v)
    print(f"  null threshold (highest peak in {len(S)} nulls, one station shifted {NULL_SHIFT_S:.0f} s each): {thr:.1f}")
    for k, v in nulls.items():
        sc = sorted((p["score"] for p in v), reverse=True)
        print(f"    null {k:22s} top peaks {', '.join(f'{x:.1f}' for x in sc[:3])}")
    rs = sorted((p["score"] for p in real), reverse=True)
    print(f"    real                        top peaks {', '.join(f'{x:.1f}' for x in rs[:5])}")
    det = [p for p in real if p["score"] > thr]
    print(f"  REAL peaks above the null threshold: {len(det)} (nulls: 0 by construction) in "
          f"{len(set(p['w0'] for p in det))} of {len(set(p['w0'] for p in real))} windows")

    # Doppler check: one velocity must explain all three Doppler profiles
    fixes = []
    for p in det:
        P = G[p["ix"], p["iy"], p["ia"]]
        w = slice(p["w0"], p["w0"] + a.window)
        cand = []
        for s, (lo, _, rb) in zip(S, idx):
            b = int(round((rb[p["ix"], p["iy"], p["ia"]] - s["rax"][0]) / (s["rax"][1] - s["rax"][0])))
            b = int(np.clip(b, 1, len(s["rax"]) - 2))
            prof = s["maps"][w, b - 1:b + 2, :].max(axis=1).mean(axis=0)
            prof[np.abs(s["dax"]) <= NOTCH] = -99
            top = np.argsort(-prof)[:4]
            cand.append([(float(s["dax"][j]), float(prof[j])) for j in top])
        best = None
        import itertools
        for combo in itertools.product(*cand):
            v, fres = doppler_velocity(S, P, R, [c[0] for c in combo])
            sp = float(np.linalg.norm(v))
            if fres < F_TOL and 40 < sp < 320 and (best is None or fres < best[0]):
                best = (fres, v, combo)
        la, lo_, h = ecef_to_lla(P)
        rec = dict(t_s=float(t[p["w0"]] + hop * a.window / 2), score=p["score"], lat=la, lon=lo_, alt_m=h,
                   doppler_ok=best is not None)
        if best:
            rec.update(doppler_res_hz=best[0], speed_ms=float(np.linalg.norm(best[1])),
                       dopplers_hz=[c[0] for c in best[2]])
        fixes.append(rec)
    ok = [f for f in fixes if f["doppler_ok"]]
    print(f"  Doppler-consistent (one velocity, residual < {F_TOL} Hz, 40-320 m/s): {len(ok)} of {len(fixes)}")

    # persistence: link Doppler-consistent fixes
    tracks = []
    for f in sorted(ok, key=lambda f: f["t_s"]):
        P = tc.ecef(f["lat"], f["lon"], f["alt_m"])
        for tr in tracks:
            g = tr[-1]; dt = f["t_s"] - g["t_s"]
            if 0 < dt <= 10 and np.linalg.norm((P - tc.ecef(g["lat"], g["lon"], g["alt_m"]))[:2]) < 320 * dt + 2500:
                tr.append(f); break
        else:
            tracks.append([f])
    tracks = [tr for tr in tracks if len(tr) >= 3]
    print(f"  fix tracks (>= 3 windows): {len(tracks)}")
    for i, tr in enumerate(sorted(tracks, key=len, reverse=True)[:10]):
        A, B = tr[0], tr[-1]
        d = np.linalg.norm(tc.ecef(A["lat"], A["lon"], A["alt_m"]) - tc.ecef(B["lat"], B["lon"], B["alt_m"]))
        span = B["t_s"] - A["t_s"]
        print(f"    track {i}: {len(tr)} windows over {span:.0f} s, moved {d/1e3:.1f} km "
              f"({d/max(span,1e-9):.0f} m/s by position) vs Doppler speed {np.median([f['speed_ms'] for f in tr]):.0f} m/s, "
              f"alt {np.median([f['alt_m'] for f in tr]):.0f} m")
    out = a.out or str(Path(a.tracks[0]).parent / "fuse_fix.json")
    json.dump(dict(stations=[s["name"] for s in S], window=a.window, threshold=thr,
                   null_top={k: sorted((p["score"] for p in v), reverse=True)[:10] for k, v in nulls.items()},
                   real_top=rs[:20], fixes=fixes, tracks=tracks), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
