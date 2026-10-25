"""fuse_fix.py -- multistatic position fix by fusing range-Doppler maps BEFORE detection,
with Doppler inside the score.

Version 1 (commit c709c56) scored a grid point by each station's best SNR at
its bistatic range over ANY Doppler. That searches ~1e6 range-cell triples,
so its noise maximum was as high as an aircraft (real 165.7 vs null 162.5).

Version 2. For a grid point P (east, north, altitude) and a horizontal
velocity v, station i sees bistatic range r_i(P) and Doppler
    f_i = a_i(P) . v,     a_i = -(u_T,i + u_R)_horizontal / lambda_i
(u = unit vectors from tower and receiver to P; +Hz = closing). Two stations'
Doppler peaks at their ranges fix v (2 equations, 2 unknowns); v then
PREDICTS the third station's Doppler to within one 0.59 Hz cell. Score =
the WEAKEST station's value at its cell (min over stations, declared
before running -- a fix needs all three): noise rarely puts peaks on two
maps whose implied velocity lands a third peak exactly where predicted.
Every pairing (A,B->C), (A,C->B), (B,C->A) over each station's top K Doppler
peaks at its range is tried; the best is kept per point.

Maps are summed over W dwells (3 s at 0.5 s hop: an aircraft moves < 1 km,
under one 1.2 km cell), in linear SNR, with the first MIN_BIN range bins and
each bin's capture median removed (static clutter).

Nulls: the same search with ONE station's maps shifted by 60 s (circular),
for each station. Threshold = highest null peak; real peaks above it are
the result. Then persistence: peaks linked across windows when the second
position is where the first one's velocity says it should be.

Positions go to the json (lab/, gitignored); stdout prints counts only.

HELD-OUT CHECK (--heldout <tracks...>): stations not used in the search. A
peak's position and velocity predict its (range, Doppler) on each of them
along its own path (+-15 s); the map there is ranked against 100 random
paths (same altitude and speed, random place and heading). A held-out
station "sees" the peak when it beats all 100. Null peaks get the same
check, which calibrates it: in a busy sky the 60 s-shift null pairs real
aircraft from different stations, so it is not empty, but a chance pairing
puts no echo on a fourth tower. The 10-25 fix (102-108 MHz): 4 of 4 held-out
towers see it at both fix times.

Usage: python tools/fuse_fix.py <tracks_A.json> <tracks_B.json> <tracks_C.json> [--window 6] [--k 4]
                                [--heldout <tracks_D.json> ...]
       (tracks files from fm_radar.py --save-maps on the same capture)
"""
import argparse
import glob
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import numba

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402
from fix_position import ecef_to_lla                           # noqa: E402

C = 2.998e8
NULL_SHIFT_S = 60.0
NOTCH = 5.0
MIN_BIN = 3
V_MIN, V_MAX = 40.0, 320.0


def load_station(path):
    d = json.load(open(path))
    st = d["station"]
    maps_dir = Path(path).parent / "maps" / Path(path).stem
    rax, _, dax = np.load(maps_dir / "axes.npy", allow_pickle=True)
    files = sorted(glob.glob(str(maps_dir / "rd_*.npy")))
    t = np.array([float(Path(f).stem[3:]) for f in files])
    db = np.stack([np.load(f).astype(np.float32) for f in files])      # (n_t, n_r, n_d) dB over CFAR noise
    lin = 10 ** (np.clip(db, -20, 40) / 10)
    lin -= np.median(lin, axis=0, keepdims=True)                        # what never changes
    lin[:, :MIN_BIN, :] = 0
    dax = np.asarray(dax, float)
    lin[:, :, np.abs(dax) <= NOTCH] = 0
    return dict(name=st.get("short_name", Path(path).stem),
                T=tc.ecef(st["tx_lat"], st["tx_lon"], st.get("tx_alt_m", 200.0)),
                lam=C / d["fc_hz"], rax=np.asarray(rax, float), dax=dax, t=t, lin=lin)


def path_value(s, R, P0, v, tc_, span=15.0):
    """Mean map value (linear SNR, median-removed) at a path's predicted cells on one station."""
    dr = s["rax"][1] - s["rax"][0]; df = s["dax"][1] - s["dax"][0]
    vals = []
    for j, tj in enumerate(s["t"]):
        if abs(tj - tc_) > span:
            continue
        P = P0 + v * (tj - tc_)
        uT = (P - s["T"]) / np.linalg.norm(P - s["T"]); uR = (P - R) / np.linalg.norm(P - R)
        rr = np.linalg.norm(P - s["T"]) + np.linalg.norm(P - R) - np.linalg.norm(s["T"] - R)
        b = int(round((rr - s["rax"][0]) / dr))
        k = int(round((float(-(uT + uR) @ v / s["lam"]) - s["dax"][0]) / df))
        if MIN_BIN <= b < s["lin"].shape[1] and 1 <= k < s["lin"].shape[2] - 1:
            vals.append(s["lin"][j, b, k - 1:k + 2].max())
    return float(np.mean(vals)) if vals else float("nan")


def heldout_seen(H, R, east, north, up, P0, v, tc_, alt, rng, n_rand=100, extent=80e3):
    """How many held-out stations see the path: its value beats all n_rand random paths."""
    sp = float(np.linalg.norm(v))
    seen = 0
    for s in H:
        real = path_value(s, R, P0, v, tc_)
        rnd = []
        for _ in range(n_rand):
            Pn = R + rng.uniform(-extent, extent) * east + rng.uniform(-extent, extent) * north + alt * up
            h = rng.uniform(0, 2 * np.pi)
            rnd.append(path_value(s, R, Pn, sp * (np.sin(h) * east + np.cos(h) * north), tc_))
        rnd = np.array(rnd); rnd = rnd[np.isfinite(rnd)]
        seen += bool(np.isfinite(real) and len(rnd) and real > rnd.max())
    return seen


@numba.njit(parallel=True, cache=False)
def score_grid(M, pk, rb, A, f0, df, K):
    """M: (3, n_r, n_d) window maps; pk: (3, n_r, K) top Doppler bins per range;
    rb: (nP, 3) range bin (-1 invalid); A: (nP, 3, 2) Doppler-per-velocity.
    Returns best score, its velocity, per point."""
    nP = rb.shape[0]
    n_d = M.shape[2]
    best = np.zeros(nP, np.float32)
    bv = np.zeros((nP, 2), np.float32)
    for p in numba.prange(nP):
        ok = True
        for i in range(3):
            if rb[p, i] < 0:
                ok = False
        if not ok:
            continue
        for c in range(3):               # c = predicted station
            a_ = (c + 1) % 3
            b_ = (c + 2) % 3
            ra = rb[p, a_]; rbb = rb[p, b_]; rc = rb[p, c]
            det = A[p, a_, 0] * A[p, b_, 1] - A[p, a_, 1] * A[p, b_, 0]
            if abs(det) < 1e-6:
                continue
            for ka in range(K):
                ia = pk[a_, ra, ka]
                fa = f0 + ia * df
                va_ = M[a_, ra, ia]
                if va_ <= best[p]:
                    continue
                for kb in range(K):
                    ib = pk[b_, rbb, kb]
                    fb = f0 + ib * df
                    vb_ = M[b_, rbb, ib]
                    if vb_ <= best[p]:
                        continue
                    vx = (fa * A[p, b_, 1] - fb * A[p, a_, 1]) / det
                    vy = (A[p, a_, 0] * fb - A[p, b_, 0] * fa) / det
                    sp = np.sqrt(vx * vx + vy * vy)
                    if sp < 40.0 or sp > 320.0:
                        continue
                    fc = A[p, c, 0] * vx + A[p, c, 1] * vy
                    ic = int(round((fc - f0) / df))
                    if ic < 1 or ic > n_d - 2:
                        continue
                    vc = max(M[c, rc, ic - 1], max(M[c, rc, ic], M[c, rc, ic + 1]))
                    s = min(va_, min(vb_, vc))
                    if s > best[p]:
                        best[p] = s
                        bv[p, 0] = vx
                        bv[p, 1] = vy
    return best, bv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tracks", nargs=3)
    ap.add_argument("--window", type=int, default=6)
    ap.add_argument("--k", type=int, default=4, help="Doppler peaks per station per range bin")
    ap.add_argument("--extent-km", type=float, default=100.0)
    ap.add_argument("--step-km", type=float, default=1.0)
    ap.add_argument("--alts", default="1000,2500,4000,6000,8500,11000")
    ap.add_argument("--peaks", type=int, default=6, help="peaks examined per window")
    ap.add_argument("--heldout", nargs="*", default=[], help="tracks files of stations NOT used in the search")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    S = [load_station(p) for p in a.tracks]
    n_t = min(len(s["t"]) for s in S)
    for s in S:
        assert np.allclose(s["t"][:n_t], S[0]["t"][:n_t]), "dwell times differ"
        assert np.allclose(s["dax"], S[0]["dax"]), "Doppler axes differ"
    t = S[0]["t"][:n_t]
    hop = float(np.median(np.diff(t)))
    dax = S[0]["dax"]; f0 = float(dax[0]); df = float(dax[1] - dax[0])
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    up = R / np.linalg.norm(R)
    east = np.cross([0, 0, 1], up); east /= np.linalg.norm(east)
    north = np.cross(up, east)
    xs = np.arange(-a.extent_km, a.extent_km + 1e-9, a.step_km) * 1e3
    alts = np.array([float(v) for v in a.alts.split(",")])
    X, Y, H = np.meshgrid(xs, xs, alts, indexing="ij")
    G = (R + X[..., None] * east + Y[..., None] * north + H[..., None] * up).reshape(-1, 3)
    nP = len(G)
    rb = np.full((nP, 3), -1, np.int64)
    A = np.zeros((nP, 3, 2), np.float32)
    uR = (G - R) / np.linalg.norm(G - R, axis=1, keepdims=True)
    for i, s in enumerate(S):
        uT = (G - s["T"]) / np.linalg.norm(G - s["T"], axis=1, keepdims=True)
        rr = np.linalg.norm(G - s["T"], axis=1) + np.linalg.norm(G - R, axis=1) - np.linalg.norm(s["T"] - R)
        b = np.round((rr - s["rax"][0]) / (s["rax"][1] - s["rax"][0])).astype(np.int64)
        valid = (b >= MIN_BIN) & (b < len(s["rax"]))
        rb[valid, i] = b[valid]
        g = -(uT + uR) / s["lam"]
        A[:, i, 0] = g @ east
        A[:, i, 1] = g @ north
    print(f"fuse_fix v2 -- " + ", ".join(s["name"] for s in S))
    print(f"  {n_t} dwells, hop {hop:.2f} s, window {a.window}, top {a.k} Doppler peaks/station/range; "
          f"grid {X.shape[0]}x{X.shape[1]}x{len(alts)} = {nP} points; score = weakest station")

    windows = list(range(0, n_t - a.window + 1, a.window // 2))

    def search(shifts):
        out = []
        for w0 in windows:
            M = np.stack([np.roll(s["lin"][:n_t], sh, axis=0)[w0:w0 + a.window].sum(axis=0) / a.window
                          for s, sh in zip(S, shifts)]).astype(np.float32)
            pk = np.argsort(-M, axis=2)[:, :, :a.k].astype(np.int64)
            best, bv = score_grid(M, pk, rb, A, f0, df, a.k)
            sc = best.reshape(X.shape).max(axis=2)
            from scipy.ndimage import maximum_filter
            loc = (sc == maximum_filter(sc, size=5)) & (sc > 0)
            ii = np.argwhere(loc)
            order = np.argsort(-sc[loc])[:a.peaks]
            for o in order:
                ix, iy = ii[o]
                ia = int(np.argmax(best.reshape(X.shape)[ix, iy]))
                flat = np.ravel_multi_index((ix, iy, ia), X.shape)
                out.append(dict(w0=int(w0), score=float(sc[ix, iy]), p=int(flat), v=bv[flat].tolist()))
        return out

    tw = time.time()
    real = search([0, 0, 0])
    print(f"  real search {time.time()-tw:.0f} s")
    sh = int(round(NULL_SHIFT_S / hop))
    nulls = {}
    for k in range(3):
        shifts = [0, 0, 0]; shifts[k] = sh
        nulls[S[k]["name"]] = search(shifts)
    thr = max(p["score"] for v in nulls.values() for p in v)
    print(f"  null threshold (highest peak, one station shifted {NULL_SHIFT_S:.0f} s, each in turn): {thr:.2f}")
    for k, v in nulls.items():
        sc = sorted((p["score"] for p in v), reverse=True)
        print(f"    null {k:22s} top {', '.join(f'{x:.2f}' for x in sc[:4])}")
    rs = sorted((p["score"] for p in real), reverse=True)
    print(f"    real                        top {', '.join(f'{x:.2f}' for x in rs[:8])}")
    det = [p for p in real if p["score"] > thr]
    print(f"  REAL peaks above the null threshold: {len(det)} in {len(set(p['w0'] for p in det))} of {len(windows)} windows")

    # FOLLOW: project a fix along its own velocity and count the dwells (+-30 s)
    # where >= 2 stations have a CFAR detection at the predicted range and
    # Doppler. Calibrated by doing the same for every null's top peaks in the
    # null's own (shifted) data. The 10-25 fix: WLNO and WJFK detections sat
    # on its projected path for 25 s.
    dets_by = [json.load(open(p_))["detections"] for p_ in a.tracks]

    def follow(P0, v, tc_, shifts):
        n = 0
        for dt in np.arange(-30.0, 30.01, hop):
            P = P0 + v * dt
            hits = 0
            for s, dd, shf in zip(S, dets_by, shifts):
                tt = tc_ + dt - shf * hop               # a shifted station's data comes from shf dwells earlier
                uT = (P - s["T"]) / np.linalg.norm(P - s["T"]); uRv = (P - R) / np.linalg.norm(P - R)
                rr = np.linalg.norm(P - s["T"]) + np.linalg.norm(P - R) - np.linalg.norm(s["T"] - R)
                fd = float(-(uT + uRv) @ v / s["lam"])
                key = round(float(tt) * 2) / 2
                if any(abs(x["range_m"] - rr) < 2500 and abs(x["doppler_hz"] - fd) < 4.0 for x in s.setdefault("_by", {}).get(key, [])):
                    hits += 1
            n += hits >= 2
        return int(n)

    for s, dd in zip(S, dets_by):
        s["_by"] = {}
        for x in dd:
            s["_by"].setdefault(round(x["t_s"] * 2) / 2, []).append(x)
    null_support = []
    for k, (name, v_) in enumerate(nulls.items()):
        for p_ in sorted(v_, key=lambda q: -q["score"])[:10]:
            P = G[p_["p"]]; vv = p_["v"][0] * east + p_["v"][1] * north
            # in the null, station k's maps were rolled by +sh: its dwell j holds original dwell j - sh
            null_support.append(follow(P, vv, float(t[p_["w0"]] + hop * a.window / 2), [sh if i == k else 0 for i in range(3)]))
    print(f"  follow support, null top peaks (dwells with >= 2 stations on the path, of {int(60/hop)+1}): "
          f"max {max(null_support)}, median {int(np.median(null_support))}")

    fixes = []
    for p in det:
        P = G[p["p"]]
        la, lo_, h = ecef_to_lla(P)
        vv = p["v"][0] * east + p["v"][1] * north
        tcen = float(t[p["w0"]] + hop * a.window / 2)
        fixes.append(dict(t_s=tcen, score=p["score"], lat=la, lon=lo_, alt_m=h,
                          ve_ms=p["v"][0], vn_ms=p["v"][1], speed_ms=float(np.hypot(*p["v"])),
                          heading_deg=float(np.degrees(np.arctan2(p["v"][0], p["v"][1])) % 360),
                          follow_support=follow(P, vv, tcen, [0, 0, 0]), _P=P))
    for f in sorted(fixes, key=lambda f: f["t_s"]):
        print(f"    fix t={f['t_s']:6.1f} s score {f['score']:.2f}  {f['speed_ms']:.0f} m/s hdg {f['heading_deg']:.0f}  "
              f"alt {f['alt_m']:.0f} m  follow support {f['follow_support']}")
    # persistence: the next window's fix must sit where this one's velocity puts it
    tracks = []
    for f in sorted(fixes, key=lambda f: f["t_s"]):
        placed = False
        for tr in tracks:
            g = tr[-1]; dt = f["t_s"] - g["t_s"]
            dh = abs((f["heading_deg"] - g["heading_deg"] + 180) % 360 - 180)
            if 0 < dt <= 120:
                pred = g["_P"] + dt * (g["ve_ms"] * east + g["vn_ms"] * north)
                miss = float(np.linalg.norm(f["_P"] - pred))
                if miss < 2000 and abs(f["speed_ms"] - g["speed_ms"]) < 30 and dh < 15:
                    f["link_miss_m"] = miss
                    tr.append(f); placed = True; break
        if not placed:
            tracks.append([f])
    tracks = [tr for tr in tracks if len(tr) >= 2]
    print(f"  fix tracks (>= 2 independent fixes, the later where the earlier's velocity predicts, gaps <= 120 s): {len(tracks)}")
    for tr in tracks:
        print(f"    {len(tr)} fixes over {tr[-1]['t_s']-tr[0]['t_s']:.1f} s, projection misses "
              f"{', '.join(f'{f_['link_miss_m']/1e3:.2f}' for f_ in tr[1:])} km")
    heldout = None
    if a.heldout:
        H = [load_station(p) for p in a.heldout]
        rng = np.random.default_rng(0)
        def check(p, shift_k=None):
            P = G[p["p"]]; vv = p["v"][0] * east + p["v"][1] * north
            tcen = float(t[p["w0"]] + hop * a.window / 2)
            return heldout_seen(H, R, east, north, up, P, vv, tcen, float(ecef_to_lla(P)[2]), rng)
        tw = time.time()
        null_seen = [check(p_) for v_ in nulls.values() for p_ in sorted(v_, key=lambda q: -q["score"])[:10]]
        real_top = sorted(real, key=lambda q: -q["score"])[:20]
        real_seen = [check(p_) for p_ in real_top]
        print(f"  held-out check ({len(H)} stations: {', '.join(h['name'] for h in H)}; {time.time()-tw:.0f} s)")
        print(f"    null top peaks seen by k held-out stations: " + " ".join(str(x) for x in sorted(null_seen, reverse=True)))
        print(f"    real top peaks:                             " + " ".join(str(x) for x in real_seen))
        for p_, k_ in zip(real_top, real_seen):
            for f in fixes:
                if f["t_s"] == float(t[p_["w0"]] + hop * a.window / 2) and abs(f["score"] - p_["score"]) < 1e-6:
                    f["heldout_seen"] = k_
        heldout = dict(stations=[h["name"] for h in H], null_seen=null_seen,
                       real_top=[dict(t_s=float(t[p_["w0"]] + hop * a.window / 2), score=p_["score"], seen=k_,
                                      **dict(zip(("lat", "lon", "alt_m"), map(float, ecef_to_lla(G[p_["p"]])))),
                                      ve_ms=p_["v"][0], vn_ms=p_["v"][1])
                                 for p_, k_ in zip(real_top, real_seen)])
    for f in fixes:
        f.pop("_P")
    out = a.out or str(Path(a.tracks[0]).parent / "fuse_fix_v2.json")
    json.dump(dict(stations=[s["name"] for s in S], window=a.window, k=a.k, threshold=thr,
                   null_top={k: sorted((p["score"] for p in v), reverse=True)[:10] for k, v in nulls.items()},
                   real_top=rs[:20], null_follow_support=null_support, fixes=fixes, tracks=tracks,
                   heldout=heldout), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
