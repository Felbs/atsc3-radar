"""track_hunt.py -- follow EVERY detection across dwells, not just the best.

Why this exists: the 20-frame flattening sweep stored only the strongest
peak per dwell, and two candidate tracks were spotted in that thin slice.
"The best peak moved" is a much weaker statement than "this peak persisted
and moved while others did not" -- and a track that merely DROPS BELOW the
per-frame maximum looks, in a best-only record, exactly like a track that
ended. Every conclusion about persistence needs the full detection list.

It also fixes the span problem the 60-frame attempt hit: find_frames'
bootstrap detector misses ~1/3 of frames over a 14.8 s load (measured: 40
found where 61 were expected), which inflates mean(diff) and trips the
frame-spacing gate. That gate is correct and stays. So this walks the
capture in BATCHES the detector handles, re-anchoring each batch, and
stitches the results on absolute sample index -- never on frame number.

Association is deliberately dumb: nearest-neighbour in (range, Doppler)
with a gate sized from the kinematics, no Kalman, no smoothing. A tracker
that interpolates is a tracker that can invent a track, and the whole
question here is whether the motion is real.

Output: lab/tracks.json -- every detection with its absolute time, plus
the associated tracks and, for each, the kinematic consistency number
that is the only thing making any of this interesting:

    does the observed range walk match the range rate its Doppler implies?

Range comes from correlation delay; Doppler from phase evolution across
batches. Nothing in the processor couples them, so agreement is evidence
and disagreement is disqualifying.

Run:  python tools/track_hunt.py [--batches 3] [--frames-per-batch 18]
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pr_core as PR                                          # noqa: E402
import pr_pilots as PP                                        # noqa: E402
import pilot_cancel_test as PCT                               # noqa: E402
sys.path.insert(0, r"Z:\src\atsc3\lab")
import spec_pilots as P                                       # noqa: E402

FS = PCT.FS
FC = PCT.FC
LAM = 2.998e8 / FC
BIN_M = 2.998e8 / FS
LSH_VOID_DB = 10.0      # pilot LS-fit-one-H below this = illuminator inaudible


def flattened(x, g, amps, signs, d0, guard, taps):
    """The validated flattening path, verbatim from pilot_cancel_test."""
    amps_sp = dict(amps)
    amps_sp["cp_common"] = 0.0
    amps_cp = {k: (amps["cp_common"] if k == "cp_common" else 0.0)
               for k in amps if not k.startswith("_")}
    al = PCT.Aligner(x, g, signs)
    cfo, per, lsh, _ = PCT.refine_alignment(al, x_cfo0, float(g["T"]), 0, amps)
    ahead = d0 + guard
    n = len(x)
    ref_sp = PCT.synth_window(g, amps_sp, signs, d0, ahead, n, cfo, per)
    ref_cp = PCT.synth_window(g, amps_cp, signs, d0, ahead, n, cfo, per)
    w, _ = PCT.fit_groups(x, [ref_sp, ref_cp], taps)
    b_sp = 1.0 - 1.0 / amps["scattered"]
    b_cp = 1.0 - 1.0 / amps["cp_common"]
    return PCT.subtract_groups(x, [ref_sp, ref_cp], taps, w, (b_sp, b_cp)), lsh


def detections(y, min_doppler, save_to=None):
    clean = PR.cancel_direct(y, y, PCT.ECA_TAPS)
    rd, rax, dax = PR.caf(y, clean, FS, PCT.MAX_RANGE_M, PCT.DOPPLER_HZ,
                          PCT.N_DOPPLER)
    det = [d for d in PR.detect(rd, rax, dax, pfa=PCT.PFA)
           if abs(d["doppler_hz"]) >= min_doppler]
    for d in det:
        d.pop("bin_idx", None)
    if save_to:
        # CFAR-normalised dB (map / local noise), which is what the eye
        # should judge - raw power is dominated by the residue pedestal
        noise, _alpha = PR.cfar_threshold_map(rd, pfa=PCT.PFA)
        snr = 10.0 * np.log10(np.maximum(rd, 1e-300)
                              / np.maximum(noise, 1e-300))
        np.save(save_to, snr.astype(np.float16))
        np.save(os.path.join(os.path.dirname(save_to), "axes.npy"),
                np.array([rax.astype(np.float32), None, dax.astype(np.float32)],
                         dtype=object), allow_pickle=True)
    return det


def comb_cells(g, n=8):
    r1 = PP.lattice_range_m(g["nfft"], g["pattern"])
    return [n_ * r1 for n_ in range(1, n + 1)]


def main():
    global x_cfo0
    ap = argparse.ArgumentParser()
    ap.add_argument("capture", nargs="?",
                    default=r"Z:\src\atsc3\data\fox25\rf25_fox_g4_0810_2254.cs16")
    ap.add_argument("--l1", default=os.path.join(os.path.dirname(HERE),
                                                 "lab", "l1_radar_capture.json"))
    ap.add_argument("--t0", type=float, default=30.0)
    ap.add_argument("--batches", type=int, default=3)
    ap.add_argument("--frames-per-batch", type=int, default=18,
                    help="18 keeps find_frames' bootstrap detector honest; "
                         "60 in one load made it miss 1/3 of frames")
    ap.add_argument("--min-doppler", type=float, default=PCT.MIN_DOPPLER_HZ)
    ap.add_argument("--comb-guard-m", type=float, default=130.0,
                    help="ignore detections within this of a comb tooth")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(HERE),
                                                  "lab", "tracks.json"))
    ap.add_argument("--save-maps", action="store_true",
                    help="save each dwell's flattened range-Doppler map "
                         "(dB, float16) under lab/maps/<out-stem>/ for "
                         "pr_panel.py; ~1.4 MB per dwell")
    a = ap.parse_args()
    maps_dir = None
    if a.save_maps:
        stem = os.path.splitext(os.path.basename(a.out))[0]
        maps_dir = os.path.join(os.path.dirname(a.out), "maps", stem)
        os.makedirs(maps_dir, exist_ok=True)

    npass, nfail = PP.selftest(verbose=False)
    print(f"[gate] pr_pilots selftest {npass} pass, {nfail} fail")
    if nfail:
        return 1
    g, _ = PCT.geometry_from_l1(a.l1)
    amps = PP.pilot_amplitudes(g["pattern"], g["boost"])
    noc = P.NOC[(g["nfft"], g["cred"])]
    signs = PP.reference_signs(noc)
    teeth = comb_cells(g)
    print(f"[l1  ] {g['nfft']}/{g['gi']}/{g['pattern']}  "
          f"comb teeth at {', '.join(f'{t/1e3:.2f}' for t in teeth)} km")

    n_dwell = int(0.150 * FS)
    all_det = []
    n_void = 0
    for b in range(a.batches):
        t0 = a.t0 + b * a.frames_per_batch * (g["frame_samples"] / FS)
        try:
            frames, spacing, err = PCT.find_frames(a.capture, t0,
                                                   a.frames_per_batch, g)
        except SystemExit as e:
            print(f"[batch {b}] anchor refused at t0={t0:.2f}s: {e}")
            continue
        print(f"[batch {b}] t0={t0:7.2f}s  {len(frames)} frames, "
              f"spacing err {err/g['T']:+.4f} symbols")
        for fi, fr in enumerate(frames):
            d0 = fr["data0_abs"]
            x = PR.load_cs16(a.capture, n_samples=n_dwell, offset=d0)
            x_cfo0 = fr["cfo_hz"]
            y, lsh = flattened(x, g, amps, signs, d0, PCT.PILOT_GUARD,
                               PCT.PILOT_TAPS)
            if lsh < LSH_VOID_DB:
                # 8/20 discone run: the illuminator was ~20 dB down and
                # the pilot fit never locked (0.1 dB). Flattening was a
                # no-op and the floor got reported as detections. A dwell
                # whose pilots are inaudible is VOID for radar purposes -
                # the illuminator must be audible before the echoes can be.
                n_void += 1
                print(f"   f{fi:2d} t={d0 / FS:8.4f}s lsh {lsh:5.1f} dB  "
                      f"VOID (pilot fit < {LSH_VOID_DB:.0f} dB - illuminator "
                      f"inaudible, not processed)")
                continue
            save_to = (os.path.join(maps_dir, f"rd_{d0 / FS:09.4f}.npy")
                       if maps_dir else None)
            det = detections(y, a.min_doppler, save_to=save_to)
            keep = [d for d in det
                    if min(abs(d["range_m"] - t) for t in teeth) > a.comb_guard_m]
            t_abs = d0 / FS
            for d in keep:
                d["t_s"] = round(t_abs, 4)
                d["batch"] = b
                all_det.append(d)
            print(f"   f{fi:2d} t={t_abs:8.4f}s lsh {lsh:5.1f} dB  "
                  f"det {len(det):2d} -> {len(keep):2d} off-comb"
                  + (f"  best {max(keep, key=lambda d: d['snr_db'])['snr_db']:5.1f} dB "
                     f"@ {max(keep, key=lambda d: d['snr_db'])['range_m']/1e3:7.3f} km "
                     f"{max(keep, key=lambda d: d['snr_db'])['doppler_hz']:+7.1f} Hz"
                     if keep else ""))

    # --- association: nearest neighbour, gate from kinematics, no smoothing
    all_det.sort(key=lambda d: d["t_s"])
    times = sorted({d["t_s"] for d in all_det})
    tracks = []
    for d in all_det:
        d["_used"] = False
    for ti, t in enumerate(times):
        for d in [x for x in all_det if x["t_s"] == t and not x["_used"]]:
            trk = [d]
            d["_used"] = True
            for t2 in times[ti + 1:]:
                last = trk[-1]
                dt = t2 - last["t_s"]
                if dt > 0.75:                     # a gap this long ends it
                    break
                # where a target with THIS Doppler must be by t2
                pred = last["range_m"] - last["doppler_hz"] * LAM * dt
                cands = [x for x in all_det
                         if x["t_s"] == t2 and not x["_used"]
                         and abs(x["range_m"] - pred) < 3 * BIN_M
                         and abs(x["doppler_hz"] - last["doppler_hz"]) < 60.0]
                if not cands:
                    continue
                nxt = min(cands, key=lambda x: abs(x["range_m"] - pred))
                nxt["_used"] = True
                trk.append(nxt)
            if len(trk) >= 3:
                tracks.append(trk)
    for d in all_det:
        d.pop("_used", None)

    print(f"\n[track] {len(all_det)} off-comb detections -> "
          f"{len(tracks)} tracks of >=3 frames"
          + (f"   ({n_void} dwells VOID: pilot fit < {LSH_VOID_DB:.0f} dB)"
             if n_void else ""))
    if n_void and n_void >= 0.5 * a.batches * a.frames_per_batch:
        print("[track] *** MOST DWELLS VOID - this capture cannot test "
              "anything; the illuminator is inaudible on this antenna ***")
    out = []
    for k, trk in enumerate(sorted(tracks, key=len, reverse=True)):
        t = [p["t_s"] for p in trk]
        r = [p["range_m"] for p in trk]
        f = [p["doppler_hz"] for p in trk]
        dt = t[-1] - t[0]
        obs = (r[-1] - r[0]) / dt if dt > 0 else float("nan")
        pred = -float(np.mean(f)) * LAM
        agree = abs(obs - pred) / abs(pred) * 100 if pred else float("nan")
        rec = dict(n=len(trk), t0=t[0], t1=t[-1], span_s=round(dt, 4),
                   range0_m=r[0], range1_m=r[-1],
                   doppler_mean_hz=round(float(np.mean(f)), 1),
                   observed_rate_ms=round(obs, 1),
                   predicted_rate_ms=round(pred, 1),
                   agreement_pct=round(agree, 1),
                   walk_bins=round(abs(r[-1] - r[0]) / BIN_M, 2),
                   snr_db=[round(p["snr_db"], 1) for p in trk],
                   points=[dict(t_s=p["t_s"], range_m=p["range_m"],
                                doppler_hz=p["doppler_hz"],
                                snr_db=round(p["snr_db"], 2)) for p in trk])
        out.append(rec)
        print(f"  track {k}: {len(trk)} frames, {dt:.3f}s, "
              f"{r[0]/1e3:.3f}->{r[-1]/1e3:.3f} km, fd {np.mean(f):+.1f} Hz | "
              f"observed {obs:+.1f} m/s vs predicted {pred:+.1f} "
              f"({agree:.0f}% off), walk {abs(r[-1]-r[0])/BIN_M:.1f} bins")

    meta = dict(capture=os.path.abspath(a.capture), fs_hz=FS, fc_hz=FC,
                batches=a.batches, frames_per_batch=a.frames_per_batch,
                comb_teeth_m=teeth, comb_guard_m=a.comb_guard_m,
                min_doppler_hz=a.min_doppler,
                n_detections=len(all_det), tracks=out,
                detections=all_det,
                note="CANDIDATES. Kinematic consistency is necessary, not "
                     "sufficient. No ADS-B truth exists for this capture.")
    with open(a.out, "w") as fh:
        json.dump(meta, fh, indent=1)
    print(f"\n-> {a.out}")
    print("CANDIDATES ONLY. Kinematic consistency is a necessary test, not a "
          "sufficient one; this capture has no ADS-B truth.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
