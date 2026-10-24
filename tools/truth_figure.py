"""truth_figure.py -- one picture of a truth run: a range-Doppler dwell with
the ADS-B predictions on it, and the tracked aircraft's range and Doppler
against their own telemetry across the whole capture.

Usage: python tools/truth_figure.py <meta.json> <t_dwell_s> <out.png> <tracks_A.json> [tracks_B.json ...]
       [--focus ICAO,ICAO,...]
"""
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_compare as tc                                     # noqa: E402

# Okabe-Ito subset; validated (dataviz six checks) with direct labels as the
# secondary encoding; text stays in ink colours.
COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
MARKERS = ["o", "s", "D", "^"]
INK, INK2, GRID = "#1f2328", "#59636e", "#d0d7de"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--focus")]
    focus = [a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--focus=")]
    focus = focus[0].split(",") if focus else []
    meta = json.load(open(args[0]))
    t_dwell = float(args[1])
    out = args[2]
    track_files = args[3:]

    truth_path = meta["truth"]
    if not os.path.exists(truth_path):
        truth_path = os.path.join(os.path.dirname(os.path.abspath(args[0])),
                                  os.path.basename(truth_path))
    truth = tc.load_truth(truth_path)
    lat, lon, alt = tc.observer()
    R = tc.ecef(lat, lon, alt)
    T = tc.ecef(tc.TX_LAT, tc.TX_LON, tc.TX_ANT_AMSL_M)
    dRT = float(np.linalg.norm(T - R))
    t_launch = float(meta["t_launch_unix"])
    cs_of = {k: v[0][7] or k for k, v in truth.items()}
    R_GATE = 3 * tc.BIN_M + 180 * 2.0
    F_GATE = 30.0

    # all detections, all windows
    dets, windows = [], []
    for p in track_files:
        trk = json.load(open(p))
        ts = [d["t_s"] for d in trk["detections"]]
        windows.append((Path(p).stem, min(ts), max(ts)))
        dets += trk["detections"]

    # the dwell map: nearest saved dwell to t_dwell among the windows
    best = None
    for p in track_files:
        for f in glob.glob(os.path.join(os.path.dirname(p), "maps", Path(p).stem, "rd_*.npy")):
            t = float(Path(f).stem[3:])
            if best is None or abs(t - t_dwell) < abs(best[0] - t_dwell):
                best = (t, f, os.path.join(os.path.dirname(f), "axes.npy"))
    t_map, map_f, ax_f = best
    m = np.load(map_f).astype(np.float32)
    axes = np.load(ax_f, allow_pickle=True)
    rng_m, dop_hz = np.asarray(axes[0], float), np.asarray(axes[2], float)

    fig = plt.figure(figsize=(16, 9.5), dpi=100, facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1],
                          hspace=0.32, wspace=0.22, left=0.055, right=0.99,
                          top=0.90, bottom=0.07)
    axm = fig.add_subplot(gs[:, 0])
    axr = fig.add_subplot(gs[0, 1])
    axd = fig.add_subplot(gs[1, 1], sharex=axr)

    # --- panel A: the dwell
    RMAX, DMAX = 20e3, 320.0
    ri = rng_m <= RMAX
    di = np.abs(dop_hz) <= DMAX
    sub = m[np.ix_(ri, di)]
    axm.imshow(sub, origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=20,
               extent=[dop_hz[di][0], dop_hz[di][-1], rng_m[ri][0] / 1e3, rng_m[ri][-1] / 1e3],
               interpolation="nearest")
    dd = [d for d in dets if abs(d["t_s"] - t_map) < 1e-3]
    axm.scatter([d["doppler_hz"] for d in dd], [d["range_m"] / 1e3 for d in dd],
                s=90, facecolors="none", edgecolors=INK, linewidths=1.2, label="detection (CFAR)", zorder=4)
    # predictions at this dwell
    px, py, pl = [], [], []
    for icao, samples in truth.items():
        pr = tc.predict(samples, t_launch + t_map, T, R, dRT)
        if pr and 0 < pr[0] < RMAX and abs(pr[1]) < DMAX:
            px.append(pr[1]); py.append(pr[0] / 1e3); pl.append(icao)
    axm.scatter(px, py, marker="x", s=36, c=INK2, linewidths=1.0, label="ADS-B prediction (any aircraft)", zorder=5)
    for k, icao in enumerate(focus):
        if icao not in truth:
            continue
        pr = tc.predict(truth[icao], t_launch + t_map, T, R, dRT)
        if not pr or not (0 < pr[0] < RMAX and abs(pr[1]) < DMAX):
            continue
        c = COLORS[k % 4]
        axm.scatter([pr[1]], [pr[0] / 1e3], marker=MARKERS[k % 4], s=150, facecolors="none",
                    edgecolors=c, linewidths=2.2, zorder=6)
        axm.annotate(cs_of[icao], (pr[1], pr[0] / 1e3), xytext=(10, 6), textcoords="offset points",
                     color=INK, fontsize=10, fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=c, lw=1.2))
    axm.axvspan(-15, 15, color="#ffffff", alpha=0.55, lw=0)
    axm.text(0, RMAX / 1e3 * 0.985, "clutter notch ±15 Hz", ha="center", va="top", fontsize=8, color=INK2)
    axm.set_xlabel("bistatic Doppler (Hz)   + = closing", color=INK)
    axm.set_ylabel("excess bistatic range (km)", color=INK)
    axm.set_title(f"one dwell, t = {t_map:.2f} s   (shade = dB over local noise, 0 to 20)",
                  fontsize=11, color=INK, loc="left")
    axm.legend(loc="upper right", fontsize=9, frameon=True)
    for s in axm.spines.values():
        s.set_color(GRID)

    # --- panels B, C: tracked aircraft across the capture
    tgrid = np.arange(0.0, meta["secs"] + 0.5, 1.0)
    for ax in (axr, axd):
        for (_, a, b) in windows:
            ax.axvspan(a, b, color="#f0f3f6", lw=0, zorder=0)
        ax.grid(True, color=GRID, lw=0.6)
        for s in ax.spines.values():
            s.set_color(GRID)
    for k, icao in enumerate(focus):
        if icao not in truth:
            continue
        c = COLORS[k % 4]
        prs = [tc.predict(truth[icao], t_launch + t, T, R, dRT) for t in tgrid]
        tt = [t for t, p in zip(tgrid, prs) if p]
        rr = [p[0] / 1e3 for p in prs if p]
        ff = [p[1] for p in prs if p]
        axr.plot(tt, rr, color=c, lw=2, zorder=2)
        axd.plot(tt, ff, color=c, lw=2, zorder=2)
        # matched detections: nearest prediction within gate
        mt, mr, mf = [], [], []
        for d in dets:
            pr = tc.predict(truth[icao], t_launch + d["t_s"], T, R, dRT)
            if pr and abs(d["range_m"] - pr[0]) < R_GATE and abs(d["doppler_hz"] - pr[1]) < F_GATE:
                mt.append(d["t_s"]); mr.append(d["range_m"] / 1e3); mf.append(d["doppler_hz"])
        axr.scatter(mt, mr, marker=MARKERS[k % 4], s=34, facecolors="white", edgecolors=c, linewidths=1.4, zorder=3)
        axd.scatter(mt, mf, marker=MARKERS[k % 4], s=34, facecolors="white", edgecolors=c, linewidths=1.4, zorder=3)
        if tt:
            axr.annotate(f"{cs_of[icao]}  ({len(mt)} dwells)", (tt[-1], rr[-1]), xytext=(4, 0),
                         textcoords="offset points", color=INK, fontsize=9, va="center")
    axr.set_ylabel("excess range (km)", color=INK)
    axr.set_title("tracked aircraft: line = own ADS-B prediction, markers = our detections, shade = processed windows", fontsize=10, color=INK, loc="left")
    axd.set_ylabel("bistatic Doppler (Hz)", color=INK)
    axd.set_xlabel("time into capture (s)", color=INK)
    axr.set_xlim(0, meta["secs"] + 30)
    plt.setp(axr.get_xticklabels(), visible=False)

    fig.suptitle(f"radarTuna truth run {Path(args[0]).stem.split('_')[-1].replace('.meta', '')}  "
                 f"RF25 / WNUV (97 km)  {meta['ant']}  rfgain {meta['rfgain']}  bias-T {meta['biast']}",
                 fontsize=13, color=INK, x=0.055, ha="left")
    fig.savefig(out, facecolor="white")
    print("->", out, "dwell", t_map, "detections at dwell", len(dd), "predictions in view", len(px))


if __name__ == "__main__":
    main()
