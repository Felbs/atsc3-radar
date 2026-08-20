"""Does pilot-flattening EAT ECHOES? The gate that decides if it is usable.

Flattening subtracts a structured, signal-correlated quantity from the
surveillance channel. That is exactly the shape of an operation that can
remove the thing you are hunting along with the thing you are cancelling
-- and the artifact-killing result (comb tooth 29.1 -> 4.8 dB while the
scrambled control stayed at 28.7) says NOTHING about whether a target
survives it. pilot_cancel_test's author flagged this as unverified. This
is it, and until it passes the comb result must not be used for target
work.

METHOD -- echo injected into REAL air, not into a synthetic signal:
a bistatic echo physically IS a delayed, Doppler-shifted, attenuated copy
of the illuminator, so that is what gets added to a real dwell. Then the
IDENTICAL flattening path runs (same alignment, same fitted weights, same
beta), and we ask whether the injected truth is still at its known cell.

Injection sites probe the specific risk:
  * BETWEEN comb teeth      -- the ordinary case, and the one that must pass
  * ON a comb tooth         -- adversarial: a target hiding exactly where
                               we cancel. Reported, NOT graded: sharing a
                               cell with the cancelled structure is
                               expected to hurt, and pretending otherwise
                               would be the dishonest part.
  * weak, between teeth     -- does flattening cost sensitivity?

PRE-REGISTERED PASS CONDITIONS (written before the first run):
  1. each off-comb echo is recovered within 1 range bin and 1 Doppler bin;
  2. its SNR loss vs the un-flattened map is < 3 dB;
  3. the no-echo control does not INVENT a peak at the probe cell.

Run:  python tools/test_flatten_preserves.py
Exit 0 = flattening preserves targets off the comb.
"""
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

CAP = r"Z:\src\atsc3\data\fox25\rf25_fox_g4_0810_2254.cs16"
L1 = os.path.join(os.path.dirname(HERE), "lab", "l1_radar_capture.json")
FS = PCT.FS
C = PR.C_LIGHT
DWELL_S = 0.150
ECHO_DB = -35.0          # echo level under the direct path, per dwell
WEAK_DB = -45.0

RESULTS = []


def rec(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def inject(x, range_m, doppler_hz, level_db):
    """A bistatic echo IS a delayed + Doppler-shifted + attenuated copy."""
    d = int(round(range_m / C * FS))
    e = np.zeros_like(x)
    e[d:] = x[:len(x) - d]
    t = np.arange(len(x), dtype=np.float64) / FS
    e = e * np.exp(2j * np.pi * doppler_hz * t).astype(x.dtype)
    return (x + (10 ** (level_db / 20.0)) * e).astype(np.complex64)


def cell_snr(y, range_m, doppler_hz, tol=1):
    """CFAR SNR at the truth cell, through pr_banked's exact pipeline."""
    clean = PR.cancel_direct(y, y, PCT.ECA_TAPS)
    rd, rax, dax = PR.caf(y, clean, FS, PCT.MAX_RANGE_M, PCT.DOPPLER_HZ,
                          PCT.N_DOPPLER)
    noise, alpha = PR.cfar_threshold_map(rd, pfa=PCT.PFA)
    snr = 10.0 * np.log10(np.maximum(rd, 1e-300)
                          / np.maximum(noise * alpha, 1e-300))
    i = int(np.argmin(np.abs(rax - range_m)))
    j = int(np.argmin(np.abs(dax - doppler_hz)))
    w = (slice(max(0, i - tol), i + tol + 1),
         slice(max(0, j - tol), j + tol + 1))
    sub = snr[w]
    k = np.unravel_index(int(np.argmax(sub)), sub.shape)
    ii, jj = max(0, i - tol) + k[0], max(0, j - tol) + k[1]
    return dict(snr_db=float(sub.max()), dr=abs(ii - i), dd=abs(jj - j),
                range_m=float(rax[ii]), doppler_hz=float(dax[jj]))


def main():
    print("=" * 70)
    print("flatten-preserves-targets gate -- injected truth in REAL air")
    print("=" * 70)
    npass, nfail = PP.selftest(verbose=False)
    print(f"[gate] pr_pilots selftest {npass} pass, {nfail} fail")
    if nfail:
        return 1

    g, _l1 = PCT.geometry_from_l1(L1)
    amps = PP.pilot_amplitudes(g["pattern"], g["boost"])
    noc = P.NOC[(g["nfft"], g["cred"])]
    signs = PP.reference_signs(noc)
    frames, spacing, err = PCT.find_frames(CAP, 30.0, 1, g)
    fr = frames[0]
    d0 = fr["data0_abs"]
    n = int(DWELL_S * FS)
    x0 = PR.load_cs16(CAP, n_samples=n, offset=d0)

    al = PCT.Aligner(x0, g, signs)
    cfo, per, lsh, _ = PCT.refine_alignment(al, fr["cfo_hz"], float(g["T"]),
                                            0, amps)
    print(f"[align] CFO {cfo:+.3f} Hz, period {per:.4f}, "
          f"LS-fit-one-H {lsh:.2f} dB")

    amps_sp = dict(amps)
    amps_sp["cp_common"] = 0.0
    amps_cp = {k: (amps["cp_common"] if k == "cp_common" else 0.0)
               for k in amps if not k.startswith("_")}
    ahead = d0 + PCT.PILOT_GUARD
    ref_sp = PCT.synth_window(g, amps_sp, signs, d0, ahead, n, cfo, per)
    ref_cp = PCT.synth_window(g, amps_cp, signs, d0, ahead, n, cfo, per)
    b_sp = 1.0 - 1.0 / amps["scattered"]
    b_cp = 1.0 - 1.0 / amps["cp_common"]

    def flatten(y):
        """The identical flattening path: fit on THIS signal, subtract beta."""
        w, _ = PCT.fit_groups(y, [ref_sp, ref_cp], PCT.PILOT_TAPS)
        return PCT.subtract_groups(y, [ref_sp, ref_cp], PCT.PILOT_TAPS, w,
                                   (b_sp, b_cp))

    tooth = PP.lattice_range_m(g["nfft"], g["pattern"])
    print(f"[comb ] tooth spacing {tooth:.1f} m; beta sp {b_sp:.4f} "
          f"cp {b_cp:.4f}\n")

    cases = [
        ("between teeth", 11_100.0, +260.0, ECHO_DB, True),
        ("weak, between teeth", 26_000.0, -430.0, WEAK_DB, True),
        ("ON a comb tooth", tooth, +198.5, ECHO_DB, False),
    ]
    for label, rng, dop, lvl, graded in cases:
        xi = inject(x0, rng, dop, lvl)
        raw = cell_snr(xi, rng, dop)
        fl = cell_snr(flatten(xi), rng, dop)
        loss = raw["snr_db"] - fl["snr_db"]
        detail = (f"{rng / 1e3:.2f} km {dop:+.0f} Hz @ {lvl:.0f} dB | "
                  f"raw {raw['snr_db']:.1f} -> flat {fl['snr_db']:.1f} dB "
                  f"(loss {loss:+.1f}) | offset {fl['dr']}r/{fl['dd']}d")
        if graded:
            rec(f"recovered after flattening: {label}",
                fl["snr_db"] > 0 and fl["dr"] <= 1 and fl["dd"] <= 1
                and loss < 3.0, detail)
        else:
            print(f"  [INFO] {label} (reported, not graded) - {detail}")

    ctl = cell_snr(flatten(x0), 11_100.0, +260.0)
    rec("no-echo control: nothing invented at the probe cell",
        ctl["snr_db"] < 12.0, f"{ctl['snr_db']:.1f} dB at 11.10 km +260 Hz")

    print("=" * 70)
    bad = [n_ for n_, ok in RESULTS if not ok]
    print(f"{len(RESULTS) - len(bad)}/{len(RESULTS)} passed")
    for n_ in bad:
        print(f"  FAILED: {n_}")
    print("\nVERDICT: " + ("flattening preserves injected truth off the comb."
                           if not bad else
                           "flattening is NOT cleared for target work."))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
