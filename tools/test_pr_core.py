#!/usr/bin/env python3
"""test_pr_core.py -- synthetic ground truth for the passive-radar core.

WHY THIS FILE EXISTS BEFORE ANY REAL CAPTURE IS TOUCHED
-------------------------------------------------------
We have 4.7 GB of RF25 already on disk and every aircraft that flew
through those 180 s is in it.  We are not going to look yet.  On real
air there is no truth to check against, so a sign error, an off-by-one
in the delay indexing, or a Doppler axis that runs backwards would not
look like a bug -- it would look like a discovery, and we would write it
up.  docs/FEASIBILITY.md made this control #1 for exactly that reason.

So: inject an echo at a KNOWN delay and a KNOWN Doppler and require the
processor to hand back that answer.  Then the two controls that decide
whether the answer means anything:

  * THE NEGATIVE CONTROL.  Put the echo 60 dB under the direct path and
    show it is MASKED without cancellation and FOUND with it.  Without
    this, a passing test proves only that the test was easy.  The
    masking is not stage-managed: at these numbers the direct path's own
    range-Doppler sidelobe pedestal sits ~3 dB ABOVE the target, which
    is the real physics of the problem and the reason ECA exists.
  * THE FALSE-ALARM CONTROL.  Pure noise in must give zero detections at
    pfa=1e-6.  A detector that finds aircraft in noise is worse than no
    detector, because it produces publishable-looking output forever.

Everything is seeded, so a pass is reproducible and a failure is
debuggable.

STATED IDEALISATIONS (they make the test easier than live air, and the
first two make ECA look better than it will be on the RSPdx):
  * the reference channel is noise-free -- the "re-encode the decoded
    bitstream" case FEASIBILITY calls v2.  A noisy reference limits
    cancellation depth; nothing here measures that limit.
  * the direct path is a static, perfectly linear multipath channel.  No
    front-end compression, no intermodulation -- and E65/E73 already
    caught this rig compressing on this exact signal.
  * the target is a point at constant delay and constant Doppler over
    the dwell.  No range-walk, no acceleration smearing.
"""
from __future__ import annotations

import math
import os
import sys
import tempfile
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import pr_core as PR                                            # noqa: E402

# ---------------------------------------------------------------------------
# the scenario -- one place, so every test quotes the same truth
# ---------------------------------------------------------------------------

FS = 6.912e6                # samples/s, the ATSC 3.0 capture rate we bank at
FC = PR.RF25_HZ             # Hz, RF25 = 539 MHz, our illuminator of record

MAX_RANGE_M = 20_000.0      # m of extra bistatic path to search
DOPPLER_HZ = 400.0          # Hz, half-span (= +-222 m/s of bistatic range-rate)
N_DOPPLER = 64              # Doppler bins

# Dwell chosen so n_batch lands exactly on N_DOPPLER: L = fs/(2*400) = 8640,
# and 64*8640 = 552960 samples = 80.0 ms.  No zero-padding, so the Doppler
# BIN SPACING and the Doppler RESOLUTION are the same 12.5 Hz and "within
# one bin" means what it sounds like.
BATCH_L = int(math.floor(FS / (2 * DOPPLER_HZ)))
N_SAMP = N_DOPPLER * BATCH_L

# The target, in truth units.
TGT_DELAY = 100             # samples of delay  -> 100 * c/fs = 4337.6 m
TGT_DOPPLER = 130.0         # Hz  -> -72.3 m/s of bistatic range-rate
ECHO_DB = -60.0             # dB below the direct path

# Static clutter: the direct ray plus two motionless reflections, all at zero
# Doppler and all inside the ECA tap span.  Delay in samples, complex gain.
CLUTTER = [(0, 1.0 + 0.0j), (5, 0.30 - 0.12j), (17, 0.08 + 0.05j)]
ECA_TAPS = 32               # samples of clutter horizon = 1.39 km of delay
                            # (must stay well below TGT_DELAY or ECA eats the
                            #  target along with the clutter)

NOISE_P = 1e-4              # thermal noise power per sample, relative to a
                            # unit-power direct path -> 40 dB direct-to-noise

SEED = 20260820


def _rng(salt=0):
    return np.random.default_rng(SEED + salt)


def noise_ref(n, rng):
    """An OFDM-like reference: flat spectrum, noise-like, unit power.

    ATSC 3.0 is ~4000 randomly-modulated subcarriers, so over a dwell it is
    statistically a band-limited complex Gaussian -- which is exactly what
    makes its ambiguity function a thumbtack on a 1/(BT) pedestal.  Using
    real ATSC 3.0 here would add pilot structure (periodic, so it puts
    ridges in the CAF) and nothing else this test needs.  That is a real
    difference from live air and it is not modelled.
    """
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2.0)
    return x.astype(np.complex64)


def delayed(x, d):
    """x delayed by d samples, zero-filled.  d >= 0."""
    out = np.zeros_like(x)
    if d == 0:
        return x.copy()
    out[d:] = x[:-d]
    return out


def echo(x, d, fd, amp, fs=FS):
    """A point target: delay d samples, Doppler fd Hz, complex amplitude amp."""
    n = np.arange(x.size, dtype=np.float64)
    ph = np.exp(2j * np.pi * fd * n / fs).astype(np.complex64)
    return (amp * delayed(x, d) * ph).astype(np.complex64)


def build_scene(rng, targets, with_clutter=True, noise_p=NOISE_P, n=N_SAMP):
    """(ref, surv).  targets = [(delay_samples, doppler_hz, amp_dB), ...]."""
    ref = noise_ref(n, rng)
    surv = np.zeros(n, dtype=np.complex64)
    if with_clutter:
        for d, g in CLUTTER:
            surv += (np.complex64(g) * delayed(ref, d)).astype(np.complex64)
    for d, fd, a_db in targets:
        surv += echo(ref, d, fd, 10.0 ** (a_db / 20.0))
    if noise_p > 0:
        w = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        surv += (math.sqrt(noise_p / 2.0) * w).astype(np.complex64)
    return ref, surv


def nearest_det(dets, i_r, i_d, tol_r=3, tol_d=3):
    """The strongest detection inside +-tol of a truth bin, or None."""
    best = None
    for d in dets:
        ir, idd = d["bin_idx"]
        if abs(ir - i_r) <= tol_r and abs(idd - i_d) <= tol_d:
            if best is None or d["snr_db"] > best["snr_db"]:
                best = d
    return best


def truth_bins(range_axis, doppler_axis, delay, fd):
    i_r = int(np.argmin(np.abs(range_axis - delay * PR.C_LIGHT / FS)))
    i_d = int(np.argmin(np.abs(doppler_axis - fd)))
    return i_r, i_d


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def t_resolution_is_what_we_claim():
    """The cell sizes must be the ones the asserts below lean on."""
    res = PR.caf_resolution(FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER, N_SAMP)
    notes = [f"range bin {res['range_bin_m']:.2f} m x {res['n_range']} bins",
             f"doppler bin {res['doppler_bin_hz']:.3f} Hz, "
             f"resolution {res['doppler_res_hz']:.3f} Hz",
             f"dwell {res['dwell_s']*1e3:.2f} ms, {res['n_batch']} batches of "
             f"{res['batch_len']}",
             f"unambiguous doppler +-{res['doppler_span_hz']:.1f} Hz"]
    assert abs(res["range_bin_m"] - PR.C_LIGHT / FS) < 1e-6
    assert res["n_batch"] == N_DOPPLER, "dwell should give exactly N_DOPPLER batches"
    # no zero-padding => spacing == resolution.  If these ever diverge, every
    # "within one bin" claim below silently gets easier.
    assert abs(res["doppler_bin_hz"] - res["doppler_res_hz"]) < 1e-6
    assert res["doppler_span_hz"] >= DOPPLER_HZ - 1e-9
    return "; ".join(notes)


def t_caf_recovers_known_range_and_doppler():
    """The core claim: a KNOWN delay and Doppler come back out, within a bin.

    No direct path here and the echo only 20 dB down -- this isolates the
    CAF's indexing and sign conventions from the cancellation problem.  If
    the delay axis were reversed or the Doppler axis conjugated, this is
    the test that catches it.
    """
    rng = _rng(1)
    ref, surv = build_scene(rng, [(TGT_DELAY, TGT_DOPPLER, -20.0)],
                            with_clutter=False, noise_p=1e-6)
    rd, rax, dax = PR.caf(ref, surv, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)

    assert rd.dtype == np.float32, f"rd_map dtype {rd.dtype}, want float32"
    assert rd.shape == (len(rax), len(dax))

    i_r, i_d = np.unravel_index(int(np.argmax(rd)), rd.shape)
    t_r, t_d = truth_bins(rax, dax, TGT_DELAY, TGT_DOPPLER)
    r_bin = float(rax[1] - rax[0])
    d_bin = float(dax[1] - dax[0])
    err_m = abs(rax[i_r] - TGT_DELAY * PR.C_LIGHT / FS)
    err_hz = abs(dax[i_d] - TGT_DOPPLER)

    assert abs(i_r - t_r) <= 1, f"range bin {i_r}, truth {t_r}"
    assert err_m <= r_bin + 1e-6, f"range error {err_m:.1f} m > one bin {r_bin:.1f} m"
    assert abs(i_d - t_d) <= 1, f"doppler bin {i_d}, truth {t_d}"
    assert err_hz <= d_bin + 1e-9, f"doppler error {err_hz:.2f} Hz > one bin {d_bin:.2f} Hz"

    v = PR.bistatic_velocity(dax[i_d], FC)
    return (f"peak at bin ({i_r},{i_d}) = {rax[i_r]:.1f} m, {dax[i_d]:+.2f} Hz; "
            f"truth {TGT_DELAY * PR.C_LIGHT / FS:.1f} m, {TGT_DOPPLER:+.2f} Hz; "
            f"err {err_m:.1f} m (bin {r_bin:.1f}), {err_hz:.2f} Hz (bin {d_bin:.2f}); "
            f"range-rate {v:+.1f} m/s")


def t_masked_without_cancellation_found_with():
    """THE NEGATIVE CONTROL.  -60 dB echo: invisible raw, detected after ECA.

    This is the test that proves cancellation is doing the work.  If it
    passed with the target visible in the raw map too, every other pass in
    this file would be uninformative.
    """
    rng = _rng(2)
    ref, surv = build_scene(rng, [(TGT_DELAY, TGT_DOPPLER, ECHO_DB)])

    rd0, rax, dax = PR.caf(ref, surv, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    t_r, t_d = truth_bins(rax, dax, TGT_DELAY, TGT_DOPPLER)
    dets0 = PR.detect(rd0, rax, dax, pfa=1e-6)
    hit0 = nearest_det(dets0, t_r, t_d)

    # what the target's cell looks like against its own CFAR threshold
    noise0, alpha0 = PR.cfar_threshold_map(rd0, pfa=1e-6)
    margin0 = 10 * math.log10(float(rd0[t_r, t_d]) / (alpha0 * noise0[t_r, t_d]))

    assert hit0 is None, (f"target detected WITHOUT cancellation at "
                          f"{hit0['range_m']:.0f} m -- the negative control "
                          f"is broken, this test proves nothing")

    surv_c = PR.cancel_direct(surv, ref, n_taps=ECA_TAPS)
    canc = PR.cancellation_db(surv, surv_c)
    rd1, rax1, dax1 = PR.caf(ref, surv_c, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    dets1 = PR.detect(rd1, rax1, dax1, pfa=1e-6)
    hit1 = nearest_det(dets1, t_r, t_d)

    noise1, alpha1 = PR.cfar_threshold_map(rd1, pfa=1e-6)
    margin1 = 10 * math.log10(float(rd1[t_r, t_d]) / (alpha1 * noise1[t_r, t_d]))
    floor_drop = 10 * math.log10(float(np.median(noise0)) / float(np.median(noise1)))

    assert hit1 is not None, (f"target NOT detected after {canc:.1f} dB of "
                              f"cancellation (cell is {margin1:+.1f} dB vs "
                              f"threshold)")
    assert abs(hit1["bin_idx"][0] - t_r) <= 1 and abs(hit1["bin_idx"][1] - t_d) <= 1

    return (f"raw: target cell {margin0:+.1f} dB vs CFAR threshold (masked, "
            f"{len(dets0)} detections, none at the target); "
            f"ECA {canc:.1f} dB on the surveillance channel dropped the map "
            f"floor {floor_drop:.1f} dB; "
            f"then target at {hit1['range_m']:.0f} m {hit1['doppler_hz']:+.1f} Hz, "
            f"SNR {hit1['snr_db']:.1f} dB, {margin1:+.1f} dB over threshold")


def t_pure_noise_gives_zero_detections():
    """THE FALSE-ALARM CONTROL.  No target in, no target out.

    Surveillance is thermal noise only -- not even a direct path -- so
    every cell of the map is an exponential draw and the expected number of
    CFAR crossings is pfa * n_cells, quoted below.  Anything more than that
    means the alpha scaling is wrong, and a detector with the wrong alpha
    reports aircraft forever.
    """
    rng = _rng(3)
    ref = noise_ref(N_SAMP, rng)
    w = (rng.standard_normal(N_SAMP) + 1j * rng.standard_normal(N_SAMP))
    surv = (math.sqrt(NOISE_P / 2.0) * w).astype(np.complex64)

    rd, rax, dax = PR.caf(ref, surv, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    dets = PR.detect(rd, rax, dax, pfa=1e-6)
    expect = 1e-6 * rd.size
    worst = max((d["snr_db"] for d in dets), default=float("-inf"))
    assert len(dets) == 0, (f"{len(dets)} false alarms in pure noise "
                            f"(expected {expect:.3f}); strongest {worst:.1f} dB")

    # ...and the same map at a deliberately loose pfa MUST light up, or the
    # zero above would be proof of a dead detector rather than a clean map.
    loose = PR.detect(rd, rax, dax, pfa=1e-2)
    assert len(loose) > 0, "no detections even at pfa=1e-2: the detector is dead"
    return (f"0 detections at pfa=1e-6 over {rd.size} cells (expected "
            f"{expect:.3f}); {len(loose)} at pfa=1e-2, so the detector is live")


def t_two_targets_both_recovered():
    """Two aircraft, different range AND different Doppler, both found."""
    rng = _rng(4)
    tgts = [(TGT_DELAY, TGT_DOPPLER, ECHO_DB),
            (250, -212.0, ECHO_DB + 3.0)]     # 10842 m, +118 m/s range-rate
    ref, surv = build_scene(rng, tgts)
    surv_c = PR.cancel_direct(surv, ref, n_taps=ECA_TAPS)
    rd, rax, dax = PR.caf(ref, surv_c, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    dets = PR.detect(rd, rax, dax, pfa=1e-6)

    lines = []
    for d, fd, _ in tgts:
        t_r, t_d = truth_bins(rax, dax, d, fd)
        hit = nearest_det(dets, t_r, t_d)
        assert hit is not None, (f"missed target at {d * PR.C_LIGHT / FS:.0f} m "
                                 f"{fd:+.1f} Hz among {len(dets)} detections")
        assert abs(hit["bin_idx"][0] - t_r) <= 1 and abs(hit["bin_idx"][1] - t_d) <= 1
        v = PR.bistatic_velocity(hit["doppler_hz"], FC)
        lines.append(f"{hit['range_m']:.0f} m {hit['doppler_hz']:+.1f} Hz "
                     f"({v:+.1f} m/s) SNR {hit['snr_db']:.1f} dB")
    return f"{len(dets)} detections total; both truths matched: " + " | ".join(lines)


def _cell_loss(targets, delay, fd, salt, taps=ECA_TAPS, noise_p=1e-8):
    """dB the target's own CAF cell loses when ECA is run over the scene."""
    ref, surv = build_scene(_rng(salt), targets, with_clutter=False,
                            noise_p=noise_p)
    rd0, rax, dax = PR.caf(ref, surv, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    rd1, _, _ = PR.caf(ref, PR.cancel_direct(surv, ref, n_taps=taps),
                       FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    i_r, i_d = truth_bins(rax, dax, delay, fd)
    return 10 * math.log10(float(rd1[i_r, i_d]) / float(rd0[i_r, i_d]))


def t_eca_blind_zone_is_doppler_not_delay():
    """Map the ECA blind zone -- and it is NOT the shape the naive story says.

    The naive story is "ECA deletes everything inside n_taps of delay", and
    it is WRONG, which this test exists to keep on the record.  The basis is
    STATIC delayed copies of the reference, so what it spans is delay <
    n_taps AND Doppler ~ 0 (within about 1/T, the dwell's own Doppler
    resolution).  A moving target survives even at short delay, because its
    phase ramp is orthogonal to every column.  Three cases:

      1. moving target OUTSIDE the tap span  -> untouched (this is the case
         that would break if n_taps were ever set to the search range);
      2. STATIC scatterer INSIDE the tap span -> destroyed, which is the
         clutter this whole stage exists to remove;
      3. moving target INSIDE the tap span   -> largely survives.

    Case 3 is the useful, non-obvious one: the short-range blind zone of an
    ECA passive radar is a zero-Doppler NOTCH, not a dead range gate.  A
    crossing aircraft with near-zero bistatic Doppler falls into it at ANY
    range, and that is a real limitation of the instrument (see
    bistatic_velocity's caveat) rather than of the tap count.
    """
    moving_outside = _cell_loss([(TGT_DELAY, TGT_DOPPLER, -20.0)],
                                TGT_DELAY, TGT_DOPPLER, salt=5)
    assert moving_outside > -0.5, (f"ECA cost a target outside its tap span "
                                   f"{moving_outside:.2f} dB -- it is eating targets")

    static_inside = _cell_loss([(8, 0.0, -20.0)], 8, 0.0, salt=6)
    assert static_inside < -30.0, (f"a STATIC scatterer at delay 8 (inside "
                                   f"{ECA_TAPS} taps) only lost "
                                   f"{static_inside:.1f} dB -- ECA is not cancelling")

    moving_inside = _cell_loss([(8, TGT_DOPPLER, -20.0)], 8, TGT_DOPPLER, salt=7)
    assert moving_inside > -3.0, (f"a MOVING target at delay 8 lost "
                                  f"{moving_inside:.1f} dB -- the notch is "
                                  f"wider in Doppler than advertised")

    return (f"moving target outside {ECA_TAPS} taps: {moving_outside:+.3f} dB; "
            f"STATIC scatterer inside the taps: {static_inside:.1f} dB (killed); "
            f"MOVING target inside the taps: {moving_inside:+.2f} dB (survives) "
            f"-- the blind zone is a zero-Doppler notch, not a dead range gate")


def t_doppler_taps_widen_the_notch():
    """n_doppler_taps>0 must cancel clutter that drifts off zero Doppler.

    Wind, a drifting transmitter clock, and our own uncorrected frequency
    offset all smear clutter a bin or two wide.  Plain delay-only ECA models
    a STATIC channel and leaves that residue behind, sitting right where
    slow targets live.

    THE NUMBER THIS TEST PINS DOWN, AND IT IS NOT THE ONE WE ASSUMED
    ---------------------------------------------------------------
    The drift here is 0.6 of a Doppler bin -- deliberately OFF-GRID.  An
    exactly-one-bin drift is exactly one basis column, and ECA then cancels
    it to machine precision; that would be a test of arithmetic, not of
    interpolation.  Off-grid, the basis has to synthesise the drift out of
    integer-bin exponentials, and a rectangular dwell's leakage falls off
    only as 1/(bin offset) in AMPLITUDE.  So the achievable cancellation of
    an off-grid drift converges SLOWLY in the number of taps.  Measured on
    this scene (2026-08-20): delay-only 8.3 dB, +-1 tap 15.6 dB (+7.4),
    +-3 taps 19.7 dB (+11.4).  TRIPLING the tap count bought 4 dB.  That is
    nothing like the "it is in the subspace, so it dies" one might expect.

    That is a budget item for the real campaign, recorded here rather than
    discovered later: Doppler taps are worth a lot but they are not a
    substitute for shortening the dwell until the clutter really is static.
    """
    rng = _rng(11)
    n = N_SAMP
    ref = noise_ref(n, rng)
    # 1/T for this dwell -- one CAF Doppler bin, in Hz
    bin_hz = FS / (BATCH_L * N_DOPPLER)
    fd_drift = 0.6 * bin_hz
    surv = (delayed(ref, 0) + 0.5 * echo(ref, 5, fd_drift, 1.0)).astype(np.complex64)
    w = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    surv = (surv + math.sqrt(1e-8 / 2.0) * w).astype(np.complex64)

    def db(q):
        return PR.cancellation_db(
            surv, PR.cancel_direct(surv, ref, n_taps=ECA_TAPS, n_doppler_taps=q))

    plain, q1, q3 = db(0), db(1), db(3)
    # plain ~= 7 dB is just "the direct died, the drifting half survived"
    assert q1 > plain + 5.0, f"+-1 Doppler tap bought only {q1 - plain:.1f} dB"
    assert q3 > q1 + 2.0, f"+-3 taps did not improve on +-1 ({q1:.1f} -> {q3:.1f})"
    assert q3 > plain + 8.0, f"+-3 Doppler taps bought only {q3 - plain:.1f} dB"
    return (f"clutter at delay 5 drifting {fd_drift:.2f} Hz = 0.6 bin (off-grid): "
            f"delay-only {plain:.1f} dB, +-1 tap {q1:.1f} dB (+{q1-plain:.1f}), "
            f"+-3 taps {q3:.1f} dB (+{q3-plain:.1f}); "
            f"convergence is slow because rectangular-dwell leakage decays as "
            f"1/offset -- budget for it, do not assume it")


def t_eca_solver_paths_agree():
    """The SVD path and the blocked normal-equations path must give one answer.

    cancel_direct() has two solvers because a 3.5e6-sample basis cannot be
    materialised, and the big-data path is the one every real run will take
    while the small-data path is the one every quick experiment takes.  Two
    solvers is two chances to be wrong, and normal equations SQUARE the
    condition number -- so this pins them together on a scene small enough
    for both, rather than trusting that they agree.
    """
    n = 60_000
    rng = _rng(12)
    ref = noise_ref(n, rng)
    surv = np.zeros(n, dtype=np.complex64)
    for d, g in CLUTTER:
        surv += (np.complex64(g) * delayed(ref, d)).astype(np.complex64)
    surv += echo(ref, TGT_DELAY, TGT_DOPPLER, 10.0 ** (ECHO_DB / 20.0))
    w = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    surv = (surv + math.sqrt(NOISE_P / 2.0) * w).astype(np.complex64)

    old = PR.ECA_DIRECT_MAX_ENTRIES
    try:
        PR.ECA_DIRECT_MAX_ENTRIES = n * ECA_TAPS + 1        # force SVD
        a = PR.cancel_direct(surv, ref, n_taps=ECA_TAPS)
        PR.ECA_DIRECT_MAX_ENTRIES = 0                       # force blocked
        b = PR.cancel_direct(surv, ref, n_taps=ECA_TAPS)
    finally:
        PR.ECA_DIRECT_MAX_ENTRIES = old

    da, db_ = PR.cancellation_db(surv, a), PR.cancellation_db(surv, b)
    # relative to the residual, not to the input: agreeing to 1% of a signal
    # that was 40 dB bigger would be no agreement at all.
    rel = float(np.sqrt(np.mean(np.abs(a - b) ** 2) / np.mean(np.abs(a) ** 2)))
    assert da > 30.0, f"SVD path only cancelled {da:.1f} dB"
    assert rel < 1e-4, (f"solver paths differ by {rel:.2e} of the residual "
                        f"({da:.2f} dB vs {db_:.2f} dB)")
    return (f"n={n}, {ECA_TAPS} taps: SVD {da:.2f} dB, blocked normal equations "
            f"{db_:.2f} dB; residuals agree to {rel:.2e} relative")


def t_bistatic_velocity_geometry():
    """Sign and scale of the range-rate conversion, plus its round trip.

    ALSO PINS A FACTOR-OF-TWO THE DOC LEAVES AMBIGUOUS.  FEASIBILITY's
    resolution table reads "2.0 Hz = 1.11 m/s bistatic range-rate" and
    "Airliner (250 m/s radial) 449 Hz".  Both are fd = rate / lambda, i.e.
    the doc's number is the rate of change of the SUM path Rt+Rr, and this
    module matches it.  But the doc calls 250 m/s "radial", and for a
    genuinely 250 m/s-radial target in near-monostatic geometry (tower and
    receiver in the same direction from the target, beta -> 0) the sum path
    changes at 2 x 250 m/s and the Doppler is 899 Hz, not 449.  Both numbers
    are asserted here so the ambiguity is recorded rather than inherited:
    the doc's 449 Hz is CONSERVATIVE by 2x for the fastest geometry, which
    matters when sizing a Doppler search window.
    """
    lam = PR.C_LIGHT / FC                       # m, wavelength at 539 MHz
    v = PR.bistatic_velocity(TGT_DOPPLER, FC)
    assert v < 0, "positive Doppler must mean a CLOSING (shortening) sum path"
    assert abs(abs(float(v)) - TGT_DOPPLER * lam) < 1e-6
    back = PR.range_rate_to_doppler(v, FC)
    assert abs(float(back) - TGT_DOPPLER) < 1e-9, "round trip broken"

    fd_doc = float(PR.range_rate_to_doppler(-250.0, FC))     # sum-path 250 m/s
    assert 445.0 < fd_doc < 453.0, f"sum-path 250 m/s -> {fd_doc:.1f} Hz, doc says 449"
    fd_mono = float(PR.range_rate_to_doppler(-2 * 250.0, FC))  # beta -> 0
    assert abs(fd_mono - 2 * fd_doc) < 1e-6

    # 1.11 m/s per 2.0 Hz, the other doc row
    assert abs(float(PR.bistatic_velocity(-2.0, FC)) - 1.11) < 0.01
    return (f"lambda {lam*100:.2f} cm; {TGT_DOPPLER:+.1f} Hz -> {float(v):+.2f} m/s "
            f"sum-path range-rate; 250 m/s of SUM path -> {fd_doc:.1f} Hz "
            f"(doc says 449); the same target at beta->0 -> {fd_mono:.1f} Hz "
            f"-- size Doppler windows off the larger number")


def t_cs16_io_roundtrip_and_guard():
    """cs16 load: exact values, offset/length honoured, and the size guard bites."""
    n = 4096
    rng = _rng(8)
    i = rng.integers(-32768, 32767, n, dtype=np.int16)
    q = rng.integers(-32768, 32767, n, dtype=np.int16)
    inter = np.empty(2 * n, dtype=np.int16)
    inter[0::2], inter[1::2] = i, q
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "synthetic.cs16")
        inter.tofile(p)
        assert PR.cs16_len(p) == n
        x = PR.load_cs16(p)
        assert x.dtype == np.complex64 and x.size == n
        want = (i.astype(np.float32) + 1j * q.astype(np.float32)) / np.float32(32768.0)
        assert np.max(np.abs(x - want.astype(np.complex64))) == 0.0, "values differ"
        y = PR.load_cs16(p, n_samples=100, offset=1000)
        assert y.size == 100 and np.max(np.abs(y - x[1000:1100])) == 0.0
        blocks = list(PR.iter_cs16(p, block=777))
        assert np.array_equal(np.concatenate(blocks), x), "iter != load"
        # the guard: ask for more than the cap and get a message, not a swap storm
        old = PR.LOAD_LIMIT_SAMPLES
        try:
            PR.LOAD_LIMIT_SAMPLES = 10
            try:
                PR.load_cs16(p)
                raise AssertionError("load_cs16 ignored LOAD_LIMIT_SAMPLES")
            except ValueError as e:
                assert "iter_cs16" in str(e), "guard must name the fix"
        finally:
            PR.LOAD_LIMIT_SAMPLES = old
    return (f"{n} samples exact bit-for-bit; offset/length slice exact; "
            f"iter_cs16 == load_cs16; size guard raises and names iter_cs16")


def t_decimation_agrees_on_the_answer():
    """decim must trade range resolution for speed, not move the target."""
    rng = _rng(9)
    ref, surv = build_scene(rng, [(TGT_DELAY, TGT_DOPPLER, -20.0)],
                            with_clutter=False, noise_p=1e-6)
    rd, rax, dax = PR.caf(ref, surv, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER, decim=4)
    i_r, i_d = np.unravel_index(int(np.argmax(rd)), rd.shape)
    true_m = TGT_DELAY * PR.C_LIGHT / FS
    r_bin = float(rax[1] - rax[0])
    err = abs(float(rax[i_r]) - true_m)
    assert err <= r_bin, f"decim=4 moved the target {err:.1f} m (bin {r_bin:.1f} m)"
    assert abs(float(dax[i_d]) - TGT_DOPPLER) <= abs(dax[1] - dax[0])
    return (f"decim=4: range bin {r_bin:.1f} m (was {PR.C_LIGHT/FS:.1f}), peak at "
            f"{float(rax[i_r]):.1f} m vs truth {true_m:.1f} m, err {err:.1f} m; "
            f"doppler {float(dax[i_d]):+.1f} Hz")


TESTS = [
    ("resolution is what we claim", t_resolution_is_what_we_claim),
    ("CAF recovers a KNOWN range and Doppler", t_caf_recovers_known_range_and_doppler),
    ("NEGATIVE CONTROL: -60 dB echo masked raw, found after ECA",
     t_masked_without_cancellation_found_with),
    ("FALSE-ALARM CONTROL: pure noise -> zero detections",
     t_pure_noise_gives_zero_detections),
    ("two targets both recovered", t_two_targets_both_recovered),
    ("ECA blind zone is a zero-Doppler notch, not a range gate",
     t_eca_blind_zone_is_doppler_not_delay),
    ("Doppler taps widen the clutter notch", t_doppler_taps_widen_the_notch),
    ("ECA solver paths (SVD vs blocked) agree", t_eca_solver_paths_agree),
    ("bistatic range-rate geometry", t_bistatic_velocity_geometry),
    ("cs16 I/O round trip + size guard", t_cs16_io_roundtrip_and_guard),
    ("decimation does not move the target", t_decimation_agrees_on_the_answer),
]


def main():
    print("=" * 78)
    print("passiveTuna pr_core -- synthetic ground truth")
    print(f"  fs {FS/1e6:.3f} MS/s   fc {FC/1e6:.1f} MHz   dwell "
          f"{N_SAMP/FS*1e3:.1f} ms ({N_SAMP} samples)")
    print(f"  truth: delay {TGT_DELAY} samples = {TGT_DELAY*PR.C_LIGHT/FS:.1f} m, "
          f"Doppler {TGT_DOPPLER:+.1f} Hz, echo {ECHO_DB:.0f} dB under direct")
    print(f"  seed {SEED} -- a pass is reproducible, a failure is debuggable")
    print("=" * 78)

    n_pass = n_fail = 0
    for name, fn in TESTS:
        t0 = time.time()
        try:
            note = fn()
            dt = time.time() - t0
            print(f"[PASS] {name}  ({dt:.1f}s)")
            if note:
                for line in str(note).split("; "):
                    print(f"       {line}")
            n_pass += 1
        except Exception:
            dt = time.time() - t0
            print(f"[FAIL] {name}  ({dt:.1f}s)")
            for line in traceback.format_exc().rstrip().split("\n"):
                print(f"       {line}")
            n_fail += 1

    print("=" * 78)
    print(f"{n_pass} passed, {n_fail} failed")
    if n_fail:
        print("RESULT: FAIL -- do not point this at a real capture")
    else:
        print("RESULT: PASS -- processor recovers known truth, is masked without "
              "cancellation,")
        print("        finds nothing in noise. Cleared for banked captures.")
    print("=" * 78)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
