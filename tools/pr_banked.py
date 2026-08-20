"""pr_banked.py -- first look: run the proven CAF chain over a banked capture.

v0 of the actual radar. Reference = the received channel itself (the direct
path dominates by 70+ dB, so the self-referenced CAF is the standard
poor-man's PBR arrangement; the re-encoded perfect reference is v2).

Per dwell: ECA-cancel zero-Doppler clutter, CAF, CFAR. Across dwells: keep
every detection with its dwell timestamp so a human (or track_check) can ask
the only question that matters at v0 -- do peaks MOVE COHERENTLY, dwell to
dwell, like something flying, or do they twinkle like residue?

Honesty rails, from FEASIBILITY.md:
  * dwell defaults to 150 ms -- the range-walk limit for a fast target;
    longer dwells silently smear exactly the targets we care about.
  * Doppler span +-1000 Hz -- sized off the 899 Hz radial-airliner number,
    not the 449 Hz sum-path number (the corrected convention).
  * Every output row carries the capture file, gain tag, dwell offset and
    all processing parameters -- a detection that cannot name its
    provenance is not evidence (captures-must-record-their-gain law).
  * This tool DETECTS CANDIDATES. It does not claim aircraft. The claim
    needs the controls: the g2/g4 gain pair, the VOID capture, and
    eventually ADS-B truth. It prints that caveat on every run.

Usage:
  python tools/pr_banked.py <capture.cs16> [--t0 30] [--dwells 20] \
      [--out lab/first_look.json]
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pr_core as PR                                          # noqa: E402

FS = 6.912e6          # Hz -- every atsc3 capture in data/ is this rate
FC = 539.0e6          # Hz -- RF25 centre (WBFF multiplex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--t0", type=float, default=30.0,
                    help="seconds into the capture to start (skip AGC settle)")
    ap.add_argument("--dwells", type=int, default=20)
    ap.add_argument("--dwell-s", type=float, default=0.150,
                    help="range-walk-limited; see FEASIBILITY.md")
    ap.add_argument("--max-range-km", type=float, default=60.0)
    ap.add_argument("--doppler-hz", type=float, default=1000.0)
    ap.add_argument("--n-doppler", type=int, default=512)
    ap.add_argument("--eca-taps", type=int, default=96)
    ap.add_argument("--pfa", type=float, default=1e-6)
    ap.add_argument("--min-doppler-hz", type=float, default=15.0,
                    help="ignore the clutter ridge this close to 0 Hz")
    ap.add_argument("--out", default=str(HERE.parent / "lab" /
                                         "first_look.json"))
    ap.add_argument("--save-maps", action="store_true",
                    help="also dump each rd_map as .npy next to --out")
    a = ap.parse_args()

    n_total = PR.cs16_len(a.capture)
    n_dwell = int(a.dwell_s * FS)
    need = int(a.t0 * FS) + a.dwells * n_dwell
    if need > n_total:
        sys.exit(f"capture holds {n_total / FS:.1f}s; asked for "
                 f"{need / FS:.1f}s")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    res = PR.caf_resolution(FS, a.max_range_km * 1e3, a.doppler_hz,
                            a.n_doppler, n_dwell)
    meta = dict(capture=os.path.abspath(a.capture),
                capture_bytes=os.path.getsize(a.capture),
                fs_hz=FS, fc_hz=FC, t0_s=a.t0, dwell_s=a.dwell_s,
                dwells=a.dwells, max_range_km=a.max_range_km,
                doppler_hz=a.doppler_hz, n_doppler=a.n_doppler,
                eca_taps=a.eca_taps, pfa=a.pfa,
                min_doppler_hz=a.min_doppler_hz, resolution=res,
                note="CANDIDATES ONLY - no aircraft claim without the "
                     "gain pair, the VOID capture, and ADS-B truth")
    print(f"[pr_banked] {os.path.basename(a.capture)}  "
          f"{n_total / FS:.1f}s of RF25 @ {FS / 1e6:.3f} MS/s")
    print(f"[pr_banked] dwell {a.dwell_s * 1e3:.0f} ms x {a.dwells}, "
          f"range to {a.max_range_km:.0f} km, doppler +-{a.doppler_hz:.0f} Hz")
    print(f"[pr_banked] resolution: {res}")

    all_det = []
    t_start = time.time()
    for k in range(a.dwells):
        off = int(a.t0 * FS) + k * n_dwell
        x = PR.load_cs16(a.capture, n_samples=n_dwell, offset=off)
        clean = PR.cancel_direct(x, x, a.eca_taps)
        canc_db = PR.cancellation_db(x, clean)
        rd, rax, dax = PR.caf(x, clean, FS, a.max_range_km * 1e3,
                              a.doppler_hz, a.n_doppler)
        det = PR.detect(rd, rax, dax, pfa=a.pfa)
        det = [d for d in det if abs(d["doppler_hz"]) >= a.min_doppler_hz]
        for d in det:
            d["dwell"] = k
            d["t_s"] = round(a.t0 + k * a.dwell_s, 3)
            d["vel_ms"] = round(PR.bistatic_velocity(d["doppler_hz"], FC), 1)
            d.pop("bin_idx", None)
        all_det += det
        if a.save_maps:
            np.save(os.path.splitext(a.out)[0] + f"_rd{k:03d}.npy", rd)
        print(f"  dwell {k:3d} t={a.t0 + k * a.dwell_s:7.3f}s  "
              f"ECA {canc_db:5.1f} dB  detections {len(det):3d}"
              + (f"  best {max(d['snr_db'] for d in det):5.1f} dB "
                 f"@ {max(det, key=lambda d: d['snr_db'])['range_m'] / 1e3:6.2f} km "
                 f"{max(det, key=lambda d: d['snr_db'])['doppler_hz']:+7.1f} Hz"
                 if det else ""))

    meta["wall_s"] = round(time.time() - t_start, 1)
    meta["detections"] = all_det
    with open(a.out, "w") as f:
        json.dump(meta, f, indent=1)
    print(f"\n[pr_banked] {len(all_det)} candidate detections -> {a.out}")
    print("[pr_banked] REMEMBER: candidates, not aircraft. Next controls: "
          "same dwells on the g2 capture (intermod moves with gain), the "
          "VOID capture (must be empty), coherent dwell-to-dwell motion.")


if __name__ == "__main__":
    sys.exit(main())
