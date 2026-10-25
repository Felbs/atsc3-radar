"""two_channel_radar.py -- passive radar on a RECORDED reference: any two coherent channels.

The rest of this repo rebuilds the reference from the broadcast, because it
has one receiver. A coherent two-channel SDR (USRP B210, KrakenSDR, RSPduo,
a clock-shared RTL pair) records the reference directly on an antenna
pointed at the tower. This runs the same canceller, CAF, CFAR and track
linker on such a recording, whatever the illuminator: DVB-T, DAB, FM,
ATSC -- the waveform never has to be understood.

Input: a .npy of shape (2, N) complex64 (channel 0 = reference, 1 =
surveillance; --swap if not), loaded memory-mapped with allow_pickle=False,
or two cs16 files (--surv). The sample rate is not in a .npy; give --fs.

Without --fc the carrier is unknown, so a track's Doppler cannot predict its
range walk. The fit runs the other way: walk = -lambda x Doppler, so every
track gives a wavelength, and real tracks agree on one. --fc then checks it.

Usage:
  python tools/two_channel_radar.py lab/ext/x.npy --fs 8e6 [--fc 594e6] [--max-range-km 3]
         [--dwell 0.5] [--hop 0.5] [--doppler 128] [--taps auto] [--notch-hz 1.5] [--jobs 16] [--out ...]
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pr_core as P                                            # noqa: E402

C = 2.998e8
_SRC = {}


def _open(a):
    if a.capture.endswith(".npy"):
        X = np.load(a.capture, mmap_mode="r", allow_pickle=False)
        if X.ndim != 2 or X.shape[0] != 2:
            raise SystemExit(f"expected shape (2, N), got {X.shape}")
        r, s = (1, 0) if a.swap else (0, 1)
        return (lambda o, n: np.asarray(X[r, o:o + n], np.complex64),
                lambda o, n: np.asarray(X[s, o:o + n], np.complex64), X.shape[1])
    n = P.cs16_len(a.capture)
    return (lambda o, k: P.load_cs16(a.capture, k, o), lambda o, k: P.load_cs16(a.surv, k, o), n)


def _dwell(job):
    (cap, surv, swap, fs, t, dwell, taps, max_r, dop, n_dop, pfa, notch, min_bins, save_dir) = job
    key = (cap, surv, swap)
    if key not in _SRC:
        _SRC[key] = _open(argparse.Namespace(capture=cap, surv=surv, swap=swap))
    ref_at, surv_at, _ = _SRC[key]
    o, n = int(t * fs), int(dwell * fs)
    ref, sv = ref_at(o, n), surv_at(o, n)
    ref = ref - ref.mean(); sv = sv - sv.mean()
    clean = P.cancel_direct(sv, ref, taps)
    canc = P.cancellation_db(sv, clean)
    rd, rax, dax = P.caf(ref, clean, fs, max_r, dop, n_dop)
    d = P.detect(rd, rax, dax, pfa=pfa)
    bin_m = C / fs
    keep = []
    for q in d:
        if abs(q["doppler_hz"]) > notch and q["range_m"] >= min_bins * bin_m:
            keep.append(dict(t_s=round(t + dwell / 2, 4), range_m=float(q["range_m"]), doppler_hz=float(q["doppler_hz"]),
                             snr_db=float(q["snr_db"])))
    if save_dir:
        noise, _ = P.cfar_threshold_map(rd, pfa=pfa)
        np.save(os.path.join(save_dir, f"rd_{t + dwell/2:09.4f}.npy"),
                (10 * np.log10(rd + 1e-30) - 10 * np.log10(noise + 1e-30)).astype(np.float16))
    return t, canc, keep, (np.asarray(rax), np.asarray(dax))


def link(dets, bin_m, max_dr_bins=2.0, max_df=3.0, max_gap=2, min_n=4):
    """Nearest-neighbour tracks in (range, Doppler), tolerating max_gap missed dwells.
    No walk-vs-Doppler gate here: with the carrier unknown, that is the test, not the filter."""
    ts = sorted({d["t_s"] for d in dets})
    by_t = {t: sorted((d for d in dets if d["t_s"] == t), key=lambda d: -d["snr_db"]) for t in ts}
    idx = {t: i for i, t in enumerate(ts)}
    live, done = [], []
    for t in ts:
        used = set()
        for tr in sorted(live, key=lambda tr: -tr["pts"][-1]["snr_db"]):
            last = tr["pts"][-1]
            best = None
            for j, d in enumerate(by_t[t]):
                if j in used:
                    continue
                if abs(d["range_m"] - last["range_m"]) <= max_dr_bins * bin_m * (idx[t] - idx[last["t_s"]]) and \
                        abs(d["doppler_hz"] - last["doppler_hz"]) <= max_df:
                    best = j; break
            if best is not None:
                used.add(best); tr["pts"].append(by_t[t][best])
        for j, d in enumerate(by_t[t]):
            if j not in used:
                live.append(dict(pts=[d]))
        still = []
        for tr in live:
            (still if idx[t] - idx[tr["pts"][-1]["t_s"]] <= max_gap else done).append(tr)
        live = still
    done += live
    out = []
    for tr in done:
        p = tr["pts"]
        if len(p) < min_n:
            continue
        tt = np.array([q["t_s"] for q in p]); rr = np.array([q["range_m"] for q in p]); ff = np.array([q["doppler_hz"] for q in p])
        rate = float(np.polyfit(tt, rr, 1)[0]) if tt[-1] > tt[0] else 0.0
        fm = float(np.mean(ff))
        out.append(dict(n=len(p), t0=float(tt[0]), t1=float(tt[-1]), span_s=float(tt[-1] - tt[0]),
                        range0_m=float(rr[0]), range1_m=float(rr[-1]), walk_bins=float(abs(rr[-1] - rr[0]) / bin_m),
                        doppler_mean_hz=fm, observed_rate_ms=rate,
                        lambda_fit_m=(-rate / fm) if abs(fm) > 1e-6 else None,
                        snr_db=[q["snr_db"] for q in p], points=p))
    return sorted(out, key=lambda tr: -tr["n"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--surv", default=None, help="surveillance cs16 when capture is the reference cs16")
    ap.add_argument("--swap", action="store_true", help=".npy channel 1 is the reference")
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--fc", type=float, default=None, help="carrier, Hz (optional: tracks fit it)")
    ap.add_argument("--max-range-km", type=float, default=3.0)
    ap.add_argument("--dwell", type=float, default=0.5)
    ap.add_argument("--hop", type=float, default=0.5)
    ap.add_argument("--doppler", type=float, default=128.0)
    ap.add_argument("--n-doppler", type=int, default=None, help="default: 2 x doppler x dwell (one bin = 1/dwell)")
    ap.add_argument("--taps", default="auto", help="ECA taps (auto = max range)")
    ap.add_argument("--pfa", type=float, default=1e-6)
    ap.add_argument("--notch-hz", type=float, default=1.5)
    ap.add_argument("--min-range-bins", type=int, default=2)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--secs", type=float, default=None)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--save-maps", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    _, _, n_total = _open(a)
    fs = a.fs
    bin_m = C / fs
    secs = a.secs or (n_total / fs - a.t0)
    taps = int(np.ceil(a.max_range_km * 1e3 / bin_m)) + 4 if a.taps == "auto" else int(a.taps)
    n_dop = a.n_doppler or int(round(2 * a.doppler * a.dwell))
    out = a.out or os.path.splitext(a.capture)[0] + "_tracks.json"
    save_dir = None
    if a.save_maps:
        save_dir = os.path.join(os.path.dirname(out), "maps", os.path.splitext(os.path.basename(out))[0])
        os.makedirs(save_dir, exist_ok=True)
    print(f"[2ch] {os.path.basename(a.capture)}: {n_total/fs:.1f} s at {fs/1e6:.3f} MS/s; range cell {bin_m:.1f} m to "
          f"{a.max_range_km} km, ECA {taps} taps; dwell {a.dwell} s ({1/a.dwell:.2f} Hz cells) +-{a.doppler:.0f} Hz; "
          f"notch +-{a.notch_hz} Hz", flush=True)
    jobs = []
    t = a.t0
    while t + a.dwell <= a.t0 + secs + 1e-9:
        jobs.append((a.capture, a.surv, a.swap, fs, round(t, 6), a.dwell, taps, a.max_range_km * 1e3, a.doppler,
                     n_dop, a.pfa, a.notch_hz, a.min_range_bins, save_dir))
        t += a.hop
    tw = time.time()
    res = []
    with ProcessPoolExecutor(a.jobs) as ex:
        for r in ex.map(_dwell, jobs):
            res.append(r)
            if len(res) % 20 == 0:
                print(f"   {len(res)}/{len(jobs)} dwells, {time.time()-tw:.0f} s", flush=True)
    res.sort(key=lambda r: r[0])
    canc = [r[1] for r in res]
    dets = [d for r in res for d in r[2]]
    if save_dir:
        rax, dax = res[0][3]
        np.save(os.path.join(save_dir, "axes.npy"), np.array([rax, None, dax], dtype=object), allow_pickle=True)
    tracks = link(dets, bin_m)
    print(f"[2ch] {len(res)} dwells in {time.time()-tw:.0f} s; cancellation median {np.median(canc):.1f} dB "
          f"(range {min(canc):.1f}-{max(canc):.1f}); {len(dets)} detections -> {len(tracks)} tracks of >= 4 dwells")
    lam_fc = C / a.fc if a.fc else None
    for i, tr in enumerate(tracks[:15]):
        lf = tr["lambda_fit_m"]
        extra = ""
        if lf is not None and tr["walk_bins"] >= 3:
            extra = f"  lambda fit {lf:+.3f} m" + (f" -> {C/lf/1e6:.0f} MHz" if lf > 0 else "")
        if lam_fc and tr["walk_bins"] >= 3:
            pred = -lam_fc * tr["doppler_mean_hz"]
            extra += f" | walk {tr['observed_rate_ms']:+.1f} vs predicted {pred:+.1f} m/s"
        print(f"  track {i:2d}: {tr['n']:3d} dwells {tr['t0']:5.1f}-{tr['t1']:5.1f} s  {tr['range0_m']:7.0f} -> {tr['range1_m']:7.0f} m "
              f"({tr['walk_bins']:4.1f} cells)  fd {tr['doppler_mean_hz']:+6.1f} Hz  max {max(tr['snr_db']):4.1f} dB{extra}")
    walkers = [tr for tr in tracks if tr["walk_bins"] >= 3 and tr["lambda_fit_m"] and tr["lambda_fit_m"] > 0]
    if walkers:
        lf = np.array([tr["lambda_fit_m"] for tr in walkers]); w = np.array([tr["n"] for tr in walkers], float)
        med = float(np.median(lf))
        print(f"[2ch] wavelength from {len(walkers)} tracks that walk >= 3 cells with the right sign: median {med:.3f} m "
              f"({C/med/1e6:.0f} MHz), dwell-weighted {np.average(lf, weights=w):.3f} m; "
              f"{int(np.sum(np.abs(lf/med - 1) < 0.15))} within 15 % of the median")
    neg = [tr for tr in tracks if tr["walk_bins"] >= 3 and tr["lambda_fit_m"] is not None and tr["lambda_fit_m"] <= 0]
    print(f"[2ch] walking tracks whose range moves the WRONG way for their Doppler: {len(neg)} (noise or clutter links)")
    json.dump(dict(capture=os.path.abspath(a.capture), fs_hz=fs, fc_hz=a.fc, dwell_s=a.dwell, hop_s=a.hop, taps=taps,
                   cancellation_db=canc, n_detections=len(dets), tracks=tracks, detections=dets), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
