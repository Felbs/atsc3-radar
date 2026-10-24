"""fm_radar.py -- passive radar off a broadcast FM station, same receiver,
same CAF/ECA/CFAR, the opposite trade to ATSC 3.0.

FM is 150-200 kHz wide, so a range cell is 1.2-2 km instead of 43 m -- but
a target sits in one cell for seconds, so the dwell is a full second and
Doppler resolves to 1 Hz at a 3 m wavelength. ATSC 3.0 sees WHERE; FM
sees HOW FAST. Hams have done FM for decades; this is that, on the rig
that just did ATSC 3.0, so the two can be compared on the same sky with
the same ADS-B truth pipeline.

Reference: no decoding. The direct path dominates the capture, so the
capture itself is a usable reference (--ref raw). Better is the
RECONSTRUCTED reference: FM-discriminate, low-pass the multiplex, and
re-modulate at constant amplitude. That strips every echo and all
amplitude fading out of the reference, which otherwise put a ghost of
every target into the reference and hence into the map (--ref recon,
the default).

Range-Doppler per dwell: ECA cancels the direct path and static clutter
(n_taps delays), the batches-CAF makes the map, CA-CFAR detects, and a
dwell-to-dwell linker builds tracks in track_hunt's JSON shape so
truth_compare / truth_trajectory / truth_figure work unchanged with a
--station=<json> describing the illuminator.

Usage: python tools/fm_radar.py <capture.cs16> --station=lab/station_wfls.json
       [--dwell 1.0] [--hop 0.5] [--max-range-km 120] [--doppler 300]
       [--ref recon|raw] [--taps 32] [--save-maps] [--out lab/tracks_FM.json]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.signal import resample_poly, firwin, filtfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pr_core as P                                            # noqa: E402

C = 2.998e8


def load_block(path, fs_in, t0, secs):
    n0 = int(t0 * fs_in)
    n = int(secs * fs_in)
    raw = np.fromfile(path, dtype=np.int16, count=2 * n, offset=4 * n0)
    if len(raw) < 2 * n:
        return None
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def decimate(x, decim, cutoff_frac=0.42):
    """Sharp FIR decimation. HD Radio (IBOC) digital sidebands live at
    +-130-200 kHz, 14-15 dB below the analog carrier, and the FM
    reconstruction cannot model them; resample_poly's default filter let
    enough of them through to cap cancellation at 23 dB. A 385-tap FIR cut
    at 0.42 fs_out (105 kHz at 250 kS/s) gives 29.6 dB on the same dwell;
    narrower bands lose cancellation again because the constant-envelope
    model of the analog signal breaks as it is band-limited harder."""
    fs_out = 1.0  # cutoff is expressed relative to the OUTPUT rate
    h = firwin(decim * 48 + 1, cutoff_frac / decim * 2)   # fraction of input Nyquist
    return resample_poly(x, 1, decim, window=h * decim).astype(np.complex64)


def reconstruct(x, fs, mpx_hz=60e3):
    """FM discriminate -> low-pass the multiplex -> re-modulate, unit amplitude."""
    ph = np.unwrap(np.angle(x))
    finst = np.diff(ph, prepend=ph[0]) * fs / (2 * np.pi)
    taps = firwin(129, mpx_hz / (fs / 2))
    finst = filtfilt(taps, 1.0, finst)      # zero-phase: no group delay between ref and surv
    # the remodulated phase, re-centred so the carrier offset is kept
    ph2 = np.cumsum(finst) * 2 * np.pi / fs
    return np.exp(1j * ph2).astype(np.complex64)


def link_tracks(dets, bin_m, lam, max_dr_bins=1.5, max_df=8.0, min_n=3):
    """Dwell-to-dwell linker: nearest detection within gates, else new track."""
    tracks = []
    by_t = {}
    for d in dets:
        by_t.setdefault(d["t_s"], []).append(d)
    times = sorted(by_t)
    open_tr = []
    for ts in times:
        used = set()
        for tr in open_tr:
            last = next(p for p in reversed(tr) if p)
            best = None
            for i, d in enumerate(by_t[ts]):
                if i in used:
                    continue
                dr = abs(d["range_m"] - last["range_m"]) / bin_m
                df = abs(d["doppler_hz"] - last["doppler_hz"])
                if dr <= max_dr_bins and df <= max_df:
                    sc = dr + df / max_df
                    if best is None or sc < best[0]:
                        best = (sc, i)
            if best:
                tr.append(by_t[ts][best[1]])
                used.add(best[1])
            else:
                tr.append(None)
        # close tracks that missed twice
        still = []
        for tr in open_tr:
            if tr[-1] is None and (len(tr) < 2 or tr[-2] is None):
                pts = [p for p in tr if p]
                if len(pts) >= min_n:
                    tracks.append(pts)
            else:
                still.append(tr)
        open_tr = still
        for i, d in enumerate(by_t[ts]):
            if i not in used:
                open_tr.append([d])
    for tr in open_tr:
        pts = [p for p in tr if p]
        if len(pts) >= min_n:
            tracks.append(pts)
    out = []
    for pts in tracks:
        t0, t1 = pts[0]["t_s"], pts[-1]["t_s"]
        r0, r1 = pts[0]["range_m"], pts[-1]["range_m"]
        fd = float(np.mean([p["doppler_hz"] for p in pts]))
        obs = (r1 - r0) / (t1 - t0) if t1 > t0 else 0.0
        prd = -fd * lam                 # +Hz = closing = range shrinking
        out.append(dict(n=len(pts), t0=t0, t1=t1, span_s=t1 - t0, range0_m=r0, range1_m=r1,
                        doppler_mean_hz=round(fd, 1), observed_rate_ms=round(obs, 1),
                        predicted_rate_ms=round(prd, 1),
                        agreement_pct=round(100 * abs(obs - prd) / max(abs(prd), 1e-9), 1),
                        walk_bins=round(abs(r1 - r0) / bin_m, 2),
                        snr_db=[p["snr_db"] for p in pts], points=pts))
    out.sort(key=lambda t: -t["n"])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--station", required=True, help="json: fc_hz, tx_lat, tx_lon, tx_alt_m, ...")
    ap.add_argument("--fs-in", type=float, default=None, help="capture rate (default: sidecar)")
    ap.add_argument("--decim", type=int, default=8)
    ap.add_argument("--dwell", type=float, default=1.0)
    ap.add_argument("--hop", type=float, default=0.5)
    ap.add_argument("--max-range-km", type=float, default=120.0)
    ap.add_argument("--doppler", type=float, default=300.0)
    ap.add_argument("--n-doppler", type=int, default=1024)
    ap.add_argument("--ref", default="recon", choices=("recon", "raw"))
    ap.add_argument("--taps", type=int, default=32)
    ap.add_argument("--pfa", type=float, default=1e-6)
    ap.add_argument("--notch-hz", type=float, default=5.0)
    ap.add_argument("--min-range-bins", type=int, default=2)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--secs", type=float, default=None)
    ap.add_argument("--save-maps", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    st = json.load(open(a.station))
    side_p = os.path.splitext(a.capture)[0] + ".json"
    side = json.load(open(side_p)) if os.path.exists(side_p) else {}
    fs_in = a.fs_in or float(side.get("rate_hz", 2e6))
    fc = float(st.get("fc_hz") or side.get("centre_hz"))
    fs = fs_in / a.decim
    lam = C / fc
    bin_m = C / fs
    n_total = os.path.getsize(a.capture) // 4
    secs = a.secs or (n_total / fs_in - a.t0)
    out = a.out or os.path.join(os.path.dirname(a.capture), "tracks_FM.json")
    maps_dir = None
    if a.save_maps:
        maps_dir = os.path.join(os.path.dirname(out), "maps", os.path.splitext(os.path.basename(out))[0])
        os.makedirs(maps_dir, exist_ok=True)

    print(f"[fm  ] {os.path.basename(a.capture)}  fc {fc/1e6:.3f} MHz  fs {fs_in/1e6:.3f} -> {fs/1e3:.0f} kS/s  "
          f"bin {bin_m:.0f} m  lambda {lam:.2f} m  dwell {a.dwell} s hop {a.hop} s  ref={a.ref} taps={a.taps}")
    print(f"[fm  ] map: {a.max_range_km:.0f} km x +-{a.doppler:.0f} Hz, {a.n_doppler} Doppler bins "
          f"({2*a.doppler/a.n_doppler:.2f} Hz each); notch +-{a.notch_hz} Hz; "
          f"expected CFAR false alarms/dwell ~ {a.pfa * int(a.max_range_km*1e3/bin_m) * a.n_doppler:.2f}")
    dets = []
    canc = []
    t = a.t0
    nd = 0
    tw = time.time()
    while t + a.dwell <= a.t0 + secs + 1e-6:
        x = load_block(a.capture, fs_in, t, a.dwell)
        if x is None:
            break
        x = decimate(x, a.decim)
        x -= x.mean()
        ref = reconstruct(x, fs) if a.ref == "recon" else x.copy()
        surv = P.cancel_direct(x, ref, a.taps)
        canc.append(P.cancellation_db(x, surv))
        rd, rax, dax = P.caf(ref, surv, fs, a.max_range_km * 1e3, a.doppler, a.n_doppler)
        rd_db = 10 * np.log10(rd + 1e-30)
        d = P.detect(rd, rax, dax, pfa=a.pfa)
        keep = [q for q in d if abs(q["doppler_hz"]) > a.notch_hz and q["range_m"] >= a.min_range_bins * bin_m]
        for q in keep:
            q["t_s"] = round(t + a.dwell / 2, 3)
            q["range_m"] = float(q["range_m"]); q["doppler_hz"] = float(q["doppler_hz"])
            q["snr_db"] = float(q["snr_db"]); q["power"] = float(q["power"]); q["threshold"] = float(q["threshold"])
            q.pop("bin_idx", None)
        dets += keep
        if maps_dir:
            noise, alpha = P.cfar_threshold_map(rd, pfa=a.pfa)
            np.save(os.path.join(maps_dir, f"rd_{t + a.dwell/2:09.4f}.npy"),
                    (10 * np.log10(rd + 1e-30) - 10 * np.log10(noise + 1e-30)).astype(np.float16))
            if nd == 0:
                np.save(os.path.join(maps_dir, "axes.npy"), np.array([rax, None, dax], dtype=object), allow_pickle=True)
        best = max(keep, key=lambda q: q["snr_db"]) if keep else None
        print(f"   d{nd:3d} t={t:7.2f}s canc {canc[-1]:5.1f} dB  det {len(d):2d} -> {len(keep):2d}"
              + (f"  best {best['snr_db']:5.1f} dB @ {best['range_m']/1e3:6.1f} km {best['doppler_hz']:+7.1f} Hz" if best else ""))
        nd += 1
        t += a.hop
    tracks = link_tracks(dets, bin_m, lam)
    print(f"[fm  ] {nd} dwells in {time.time()-tw:.0f} s; cancellation median {np.median(canc):.1f} dB; "
          f"{len(dets)} detections -> {len(tracks)} tracks of >= 3 dwells")
    for i, tr in enumerate(tracks[:12]):
        print(f"  track {i}: {tr['n']} dwells, {tr['span_s']:.1f}s, {tr['range0_m']/1e3:.1f}->{tr['range1_m']/1e3:.1f} km, "
              f"fd {tr['doppler_mean_hz']:+.1f} Hz | observed {tr['observed_rate_ms']:+.1f} m/s vs predicted "
              f"{tr['predicted_rate_ms']:+.1f} ({tr['agreement_pct']:.0f}% off), max {max(tr['snr_db']):.1f} dB")
    json.dump(dict(capture=os.path.abspath(a.capture), fs_hz=fs, fc_hz=fc, station=st,
                   dwell_s=a.dwell, hop_s=a.hop, ref=a.ref, taps=a.taps, cancellation_db=[float(c) for c in canc],
                   n_detections=len(dets), tracks=tracks, detections=dets,
                   note="FM illuminator; CANDIDATES until truth_compare says otherwise"),
              open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
