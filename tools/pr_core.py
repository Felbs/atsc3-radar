#!/usr/bin/env python3
"""pr_core.py -- the passive-bistatic-radar DSP core: CAF, ECA, CFAR.

An ATSC 3.0 tower is a half-megawatt floodlight that never blinks.  Our
RSPdx hears two things at once: the DIRECT path from the tower, and --
70 to 90 dB weaker -- echoes off aircraft, delayed by the extra path
length and Doppler-shifted by the target's motion along the ellipse.
This module turns the second thing into a range-Doppler map.

WHAT THIS FILE EXISTS TO PREVENT
--------------------------------
docs/FEASIBILITY.md did the link budget and reached a conclusion that
sets the whole design: thermal noise is NOT the wall.  A 6.912 MHz
signal integrated for 0.5 s is a time-bandwidth product of 3.5 million
(+65 dB), which puts an airliner at 20/20 km some 60 dB above thermal.
The wall is the direct path, which is 73-89 dB stronger than that echo
and arrives at the same antenna.  So:

  * `caf()` is the processing gain.  It must be FFT-based, because the
    naive delay-Doppler double loop over 3.5 million samples x hundreds
    of range bins x hundreds of Doppler bins is not a slow program, it
    is a program that never finishes.
  * `cancel_direct()` is the experiment.  Everything the project can
    ever claim rests on suppressing the direct path and its zero-Doppler
    clutter far enough that an echo is visible underneath it.
  * `detect()` is the honesty gate.  A correlation peak is not an
    aircraft.  CFAR fixes the false-alarm rate FIRST and only then asks
    what crosses it, so "we saw something" has a stated probability of
    being noise instead of an eyeball's opinion.

NOTHING HERE OPENS THE RADIO and nothing here is validated against live
air.  It is validated against synthetic ground truth in
`test_pr_core.py`, which injects an echo at a KNOWN range and a KNOWN
Doppler and requires this code to hand back that answer -- including the
negative controls (masked without cancellation, zero detections in pure
noise).  That ordering is deliberate: on real air we have no truth to
check against yet, so a bug would look exactly like a discovery.

UNITS: every physics constant below carries its units in the comment.
Sample indices are `n`, delays in samples are `tau`, delays in metres of
EXTRA path are `range_m`.  Frequencies in Hz, velocities in m/s.
"""
from __future__ import annotations

import math
import os

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import fft as SFFT
from scipy.ndimage import maximum_filter, uniform_filter
from scipy.signal import resample_poly

# ---------------------------------------------------------------------------
# physics constants -- units in the comment, always
# ---------------------------------------------------------------------------

C_LIGHT = 299_792_458.0        # m/s, speed of light in vacuum (exact, SI)

# Our illuminator of record.  RF25 = 539.0 MHz centre, 6.912 MS/s captures
# already banked (docs/FEASIBILITY.md).  Only used as a default for
# bistatic_velocity(); every function takes its own fs/fc explicitly.
RF25_HZ = 539.0e6              # Hz, ATSC 3.0 RF channel 25 centre frequency
FS_ATSC3 = 6.912e6             # samples/s, the STVT/atsc3 capture rate

# load_cs16() will refuse to materialise more than this many complex samples
# in one call.  4.7 GB of cs16 is 1.175e9 samples; as complex64 that is
# 9.4 GB of RAM, i.e. the machine, not the file, is what breaks.  The cap
# exists so the failure is a message that names the fix instead of a swap
# storm.  512 Mi complex64 = 4 GiB.
LOAD_LIMIT_SAMPLES = 1 << 29

# cancel_direct() solves by SVD while the basis has at most this many entries
# (16 B each in complex128, so 4e6 entries = 64 MiB), and by blocked normal
# equations above it.  Module-level so the test can flip it and prove the two
# paths agree -- an untested second solver is a second chance to be wrong.
ECA_DIRECT_MAX_ENTRIES = 4_000_000


# ---------------------------------------------------------------------------
# 1 -- capture I/O.  The file is 4.7 GB; it is never read whole.
# ---------------------------------------------------------------------------

def cs16_len(path) -> int:
    """Number of complex samples in an interleaved-int16 capture."""
    nbytes = os.path.getsize(path)
    if nbytes % 4:
        # 4 bytes per complex sample (int16 I + int16 Q).  A remainder means
        # a truncated capture -- report it rather than silently dropping a
        # half sample, because capture-integrity is a standing lab law.
        raise ValueError(f"{path}: {nbytes} bytes is not a whole number of "
                         f"cs16 samples (remainder {nbytes % 4})")
    return nbytes // 4


def iter_cs16(path, block=1 << 22, offset=0, n_samples=None):
    """Yield complex64 blocks of a cs16 capture.  Bounded memory, always.

    This is the streaming entry point and the one to reach for on the real
    captures.  `block` is in COMPLEX samples; the default 4 Mi samples is
    32 MiB of complex64 per block.
    """
    total = cs16_len(path)
    if offset < 0 or offset > total:
        raise ValueError(f"offset {offset} outside capture of {total} samples")
    want = total - offset if n_samples is None else min(n_samples, total - offset)
    if want <= 0:
        return
    # memmap, not read(): the OS pages in only what we touch, so a 4.7 GB
    # file costs 4.7 GB of address space and ~32 MiB of resident memory.
    raw = np.memmap(path, dtype=np.int16, mode="r")
    done = 0
    while done < want:
        take = min(block, want - done)
        i0 = 2 * (offset + done)
        chunk = np.asarray(raw[i0:i0 + 2 * take], dtype=np.float32)
        out = chunk[0::2] + 1j * chunk[1::2]
        # /32768 puts full-scale int16 at |x| ~ 1.0.  32768 = 2**15, the
        # magnitude of the most negative int16, so the scaling is symmetric
        # and nothing can land outside [-1, 1).
        yield (out / np.float32(32768.0)).astype(np.complex64)
        done += take
    del raw


def load_cs16(path, n_samples=None, offset=0) -> np.ndarray:
    """Load a slice of an interleaved-int16 capture as complex64, scaled /32768.

    The file is memory-mapped and converted in bounded chunks -- it is never
    read into memory in one piece, and no float64 intermediate of the whole
    capture is ever built.  What IS bounded by RAM is the RETURN value: a
    complex64 array is 8 bytes per sample against the file's 4, so asking
    for a whole 4.7 GB capture asks for 9.4 GB of RAM.  Past
    LOAD_LIMIT_SAMPLES this raises rather than thrashes, and the message
    names `iter_cs16` as the fix.

    n_samples=None means "to end of file", subject to that cap.
    """
    total = cs16_len(path)
    want = total - offset if n_samples is None else min(int(n_samples), total - offset)
    if want > LOAD_LIMIT_SAMPLES:
        raise ValueError(
            f"{path}: refusing to materialise {want} complex samples "
            f"({want * 8 / 2**30:.1f} GiB of complex64); cap is "
            f"{LOAD_LIMIT_SAMPLES}. Pass n_samples=, or stream with "
            f"iter_cs16(path).")
    out = np.empty(max(want, 0), dtype=np.complex64)
    pos = 0
    for blk in iter_cs16(path, offset=offset, n_samples=want):
        out[pos:pos + blk.size] = blk
        pos += blk.size
    return out[:pos]


# ---------------------------------------------------------------------------
# 2 -- the cross-ambiguity function
# ---------------------------------------------------------------------------

def caf(ref, surv, fs, max_range_m, doppler_hz, n_doppler, decim=1):
    """Cross-ambiguity function by the batches (correlate-then-FFT) algorithm.

    Returns `(rd_map, range_axis_m, doppler_axis_hz)` with `rd_map` of shape
    `(n_range, n_doppler)`, dtype float32, magnitude-SQUARED (power).

    WHAT IS BEING COMPUTED
    ----------------------
    The CAF is

        chi(tau, fd) = sum_n  surv[n] * conj(ref[n - tau]) * exp(-j2 pi fd n / fs)

    Done literally that is n_range x n_doppler x N multiply-accumulates.
    For our banked capture (N = 3.5e6 at T = 0.5 s, ~1400 range bins for
    60 km, 512 Doppler bins) that is 2.5e12 operations per dwell.  Not slow:
    infeasible.

    THE ALGORITHM AND WHY THIS ONE
    ------------------------------
    Split the dwell into `n_batch` batches of `L` samples.  Inside one
    batch the target's Doppler phase is essentially constant, so the batch
    only needs a plain cross-correlation over delay -- which is three FFTs.
    The Doppler then lives entirely in how that correlation's PHASE rotates
    from batch to batch, so one FFT ACROSS the batch index resolves it, per
    range bin.  Cost drops to ~n_batch FFTs of length L plus n_range FFTs of
    length n_doppler.  For the numbers above that is ~1e9 operations: four
    orders of magnitude, and it is exact for a constant-Doppler target, not
    an approximation to it.

    The alternative (Doppler-shift the reference, correlate, repeat per
    Doppler bin) costs n_doppler full correlations and buys nothing here.
    The batch method's one real cost is the batch-rate sampling of Doppler:
    the unambiguous Doppler span is +-fs/(2L), which is why L is derived
    FROM `doppler_hz` below rather than chosen for convenience.

    THE CORRELATION IS LINEAR, NOT CIRCULAR
    ---------------------------------------
    Each batch takes its L surveillance samples but L + n_lag - 1 reference
    samples, reaching back into the previous batch, and transforms at
    M >= L + n_lag - 1.  A circular correlation of length L (the cheaper and
    more common shortcut) wraps the first `tau` reference samples around and
    contaminates exactly the long-delay bins -- which is where the distant
    targets are.  Batch 0 alone has no previous batch to reach into and is
    zero-padded there; that costs batch 0 a little correlation energy at
    long lags and nothing else.

    ARGUMENTS
    ---------
    ref, surv    complex arrays, same length, same time origin.  `ref` is
                 the direct-path reference channel; `surv` is the
                 surveillance channel (ideally after cancel_direct()).
    fs           samples/s of BOTH inputs, before `decim`.
    max_range_m  metres of EXTRA (bistatic minus baseline) path to search.
                 Delay bin 0 is the direct path.  Negative delays are not
                 searched: an echo cannot arrive before the direct signal.
    doppler_hz   half-span of the Doppler search, Hz.  Sets L.
    n_doppler    number of Doppler bins returned (the across-batch FFT
                 length).  If fewer batches than bins exist the transform
                 is zero-padded -- that INTERPOLATES the axis, it does not
                 add resolution, which stays 1/T_dwell.
    decim        integer decimation applied to both channels first, with an
                 anti-alias filter.  Trades range resolution (c/fs, 43 m at
                 6.912 MS/s) for speed.  Doppler is untouched by it.

    NORMALISATION: the map is divided by (L * n_batch), so a coherent target
    peaks at roughly its per-sample power and the noise-like pedestal sits
    at power/(L*n_batch).  CFAR is scale-invariant, so this is purely so the
    numbers are readable.
    """
    ref = np.asarray(ref)
    surv = np.asarray(surv)
    if ref.ndim != 1 or surv.ndim != 1:
        raise ValueError("ref and surv must be 1-D")
    n = min(ref.size, surv.size)
    ref, surv = ref[:n], surv[:n]
    decim = int(decim)
    if decim > 1:
        # resample_poly, not slicing: plain slicing aliases every out-of-band
        # interferer straight onto our range bins.
        ref = resample_poly(ref, 1, decim).astype(np.complex64)
        surv = resample_poly(surv, 1, decim).astype(np.complex64)
        fs = fs / decim
        n = min(ref.size, surv.size)
        ref, surv = ref[:n], surv[:n]

    # range bins: tau samples <-> c*tau/fs metres of extra path
    n_lag = int(math.floor(max_range_m * fs / C_LIGHT)) + 1
    if n_lag < 1:
        raise ValueError("max_range_m too small: no range bins")

    # L from the Doppler span.  Batch rate is fs/L; the across-batch FFT is
    # therefore sampled at fs/L and unambiguous over +-fs/(2L).  We need
    # fs/(2L) >= doppler_hz.
    if doppler_hz <= 0:
        raise ValueError("doppler_hz must be > 0")
    L = int(math.floor(fs / (2.0 * doppler_hz)))
    if L < 2:
        raise ValueError(f"doppler_hz={doppler_hz} too wide for fs={fs}")
    if L <= n_lag:
        raise ValueError(
            f"batch length L={L} <= n_lag={n_lag}: the Doppler span "
            f"({doppler_hz} Hz) and the range span ({max_range_m} m) cannot "
            f"both be had at fs={fs}. Widen decim or narrow one of them.")
    n_batch = n // L
    if n_batch < 1:
        raise ValueError(f"need at least {L} samples for one batch, have {n}")
    n_doppler = int(n_doppler)
    if n_batch > n_doppler:
        # More batches than Doppler bins would alias the across-batch FFT.
        # Truncate the dwell instead: honest, and the message says how to
        # use the rest of the data.
        n_batch = n_doppler

    # ref padded with n_lag-1 leading zeros so batch 0 can be indexed the
    # same way as every other batch (see docstring).
    ref_pad = np.concatenate([np.zeros(n_lag - 1, dtype=ref.dtype), ref])
    M = SFFT.next_fast_len(L + n_lag - 1)

    xcorr = np.empty((n_batch, n_lag), dtype=np.complex128)
    # Chunk the batches so the (chunk, M) work array stays bounded even for a
    # multi-second dwell.
    per = max(1, int(64 * (1 << 20) // (16 * M)))
    for b0 in range(0, n_batch, per):
        b1 = min(b0 + per, n_batch)
        nb = b1 - b0
        a = np.zeros((nb, M), dtype=np.complex128)
        r = np.zeros((nb, M), dtype=np.complex128)
        for k in range(nb):
            s = (b0 + k) * L
            a[k, :L] = surv[s:s + L]
            r[k, :L + n_lag - 1] = ref_pad[s:s + L + n_lag - 1]
        A = SFFT.fft(a, axis=-1)
        R = SFFT.fft(r, axis=-1)
        # corr[k] = sum_m a[m] conj(r[m+k])  for k = 0 .. n_lag-1.
        # ifft(conj(A)*R)[k] = sum_m conj(a[m]) r[m+k]; conjugate it.
        corr = np.conj(SFFT.ifft(np.conj(A) * R, axis=-1))[:, :n_lag]
        # r was built starting n_lag-1 samples EARLY, so k and tau run
        # opposite ways: tau = n_lag-1-k.  Reverse.
        xcorr[b0:b1] = corr[:, ::-1]
        del a, r, A, R, corr

    # Doppler: a target at fd advances the batch-to-batch phase by
    # exp(+j2 pi fd L / fs), so a forward FFT over the batch axis puts it at
    # fftfreq(n_doppler, d=L/fs) == fd.  No sign flip anywhere.
    Y = SFFT.fft(xcorr, n=n_doppler, axis=0)
    Y = SFFT.fftshift(Y, axes=0)
    Y /= float(L * n_batch)

    rd_map = np.ascontiguousarray((Y.real**2 + Y.imag**2).T.astype(np.float32))
    range_axis = (np.arange(n_lag, dtype=np.float64) * C_LIGHT / fs)   # metres
    doppler_axis = SFFT.fftshift(SFFT.fftfreq(n_doppler, d=L / fs))    # Hz
    return rd_map, range_axis, doppler_axis


def caf_resolution(fs, max_range_m, doppler_hz, n_doppler, n_samples, decim=1):
    """The cell sizes `caf()` will actually deliver, without running it.

    Reported separately because "within one bin" is a claim that needs the
    bin size stated, and because a zero-padded Doppler axis has a BIN
    SPACING finer than its RESOLUTION -- conflating the two is how an
    interpolated peak gets sold as precision it does not have.
    """
    fs_e = fs / decim
    n = n_samples // decim
    n_lag = int(math.floor(max_range_m * fs_e / C_LIGHT)) + 1
    L = int(math.floor(fs_e / (2.0 * doppler_hz)))
    n_batch = min(n // L, int(n_doppler))
    dwell_s = n_batch * L / fs_e
    return {
        "range_bin_m": C_LIGHT / fs_e,          # metres per delay sample
        "n_range": n_lag,
        "batch_len": L,
        "n_batch": n_batch,
        "dwell_s": dwell_s,
        "doppler_bin_hz": (fs_e / L) / n_doppler,   # axis spacing
        "doppler_res_hz": 1.0 / dwell_s,            # true resolution, 1/T
        "doppler_span_hz": fs_e / (2.0 * L),        # unambiguous half-span
    }


# ---------------------------------------------------------------------------
# 3 -- ECA: kill the direct path and everything that stands still
# ---------------------------------------------------------------------------

def _eca_block(ref_win, dop_phase, r0, r1):
    """One row-block of the ECA basis, shape (r1-r0, n_taps*len(dop_phase))."""
    base = ref_win[r0:r1]                       # (rows, n_taps), delay 0..n_taps-1
    if dop_phase is None:
        return np.ascontiguousarray(base)
    cols = [base * ph[r0:r1, None] for ph in dop_phase]
    return np.concatenate(cols, axis=1)


def cancel_direct(surv, ref, n_taps, n_doppler_taps=0, ridge=1e-9):
    """Extensive Cancellation Algorithm: project out the direct path + clutter.

    Returns a cleaned copy of `surv`, same dtype/length.

    WHAT IS BEING KILLED, AND WHY IT IS THE WHOLE PROJECT
    ----------------------------------------------------
    The tower is 73-89 dB above the echo we want (docs/FEASIBILITY.md).
    But the direct path and every static scatterer -- terrain, buildings,
    the neighbour's roof -- share one property an aircraft does not have:
    they do not move, so they sit at ZERO Doppler.  Zero-Doppler energy
    would be harmless if the CAF were perfect, but it is not: a noise-like
    OFDM reference has a range-Doppler sidelobe pedestal about 1/(BT) below
    its own peak, spread over the entire map.  With the direct 80 dB up and
    only ~65 dB of processing gain, that pedestal BURIES every target
    everywhere, including range-Doppler cells the direct path never
    occupies.  Cancelling zero-Doppler is not about clearing bin (0,0); it
    is about lowering the floor of the whole map back to thermal.

    THE METHOD
    ----------
    Build a basis of delayed copies of the reference, `ref[n-i]` for
    i = 0..n_taps-1, optionally times Doppler shifts exp(j2 pi q n / N) for
    q = -n_doppler_taps..+n_doppler_taps.  Any linear time-invariant
    multipath channel from the tower to us -- direct ray plus every static
    reflection inside n_taps samples of delay -- is a linear combination of
    those columns.  Least-squares project `surv` onto that subspace and
    subtract:  clean = surv - X (X^H X)^-1 X^H surv.

    The Doppler taps exist because clutter is not perfectly static: wind in
    trees, a slowly drifting transmitter clock, and our own uncorrected
    frequency offset all smear clutter a few 1/T-wide bins off zero.  Each
    q costs another n_taps columns, so they are off by default.

    A TARGET IS NOT CANCELLED because it is at a delay beyond n_taps (set
    n_taps to cover only the clutter horizon, NOT the target search range)
    and because its Doppler puts it out of the subspace even if it were not.
    The residual self-cancellation is ~ P/N in power, negligible while the
    number of columns P stays far below the dwell length N.  Choosing
    n_taps larger than the clutter horizon is not free caution -- it starts
    eating targets.

    NUMERICS
    --------
    Two paths.  Small problems go through `lstsq` (SVD), which is the
    stable one and squares nothing.  Large problems accumulate the normal
    equations X^H X and X^H y in complex128 over row blocks, because
    materialising X for a 3.5e6-sample dwell is tens of gigabytes.  Normal
    equations square the condition number, so that path adds a ridge term
    `ridge * trace(G)/P` on the diagonal -- scale-relative, so it means the
    same thing whatever units the capture is in -- and falls back to
    `lstsq` on the (P x P) system if the solve is singular.
    """
    surv = np.asarray(surv)
    ref = np.asarray(ref)
    n = min(surv.size, ref.size)
    surv_c = surv[:n]
    ref = ref[:n]
    n_taps = int(n_taps)
    n_doppler_taps = int(n_doppler_taps)
    if n_taps < 1:
        raise ValueError("n_taps must be >= 1")
    if n_taps >= n:
        raise ValueError(f"n_taps={n_taps} needs a dwell longer than {n}")

    # ref_win[m, k] = ref[m - (n_taps-1-k)]; reversing k gives delay 0..n_taps-1.
    # sliding_window_view is a VIEW, so this costs no memory.
    ref_pad = np.concatenate([np.zeros(n_taps - 1, dtype=ref.dtype), ref])
    ref_win = sliding_window_view(ref_pad, n_taps)[:, ::-1]

    if n_doppler_taps > 0:
        qs = np.arange(-n_doppler_taps, n_doppler_taps + 1)
        idx = np.arange(n)
        # q bins of 1/T each: exp(j2 pi q n / N) is exactly one DFT bin of
        # Doppler per q, matching the CAF's own Doppler resolution.
        dop_phase = [np.exp(2j * np.pi * q * idx / n).astype(np.complex128)
                     for q in qs]
    else:
        dop_phase = None

    n_col = n_taps * (2 * n_doppler_taps + 1)
    if n_col >= n:
        raise ValueError(f"ECA basis has {n_col} columns for {n} samples: "
                         "underdetermined, the 'clean' output would be zero")

    y = surv_c.astype(np.complex128)

    # Direct SVD path while the basis fits comfortably (16 B/entry).
    if n * n_col <= ECA_DIRECT_MAX_ENTRIES:
        X = _eca_block(ref_win, dop_phase, 0, n).astype(np.complex128)
        w, *_ = np.linalg.lstsq(X, y, rcond=None)
        clean = y - X @ w
        return clean.astype(surv.dtype, copy=False)

    rows = max(1 << 12, int(ECA_DIRECT_MAX_ENTRIES // max(n_col, 1)))
    G = np.zeros((n_col, n_col), dtype=np.complex128)
    b = np.zeros(n_col, dtype=np.complex128)
    for r0 in range(0, n, rows):
        r1 = min(r0 + rows, n)
        B = _eca_block(ref_win, dop_phase, r0, r1).astype(np.complex128)
        G += B.conj().T @ B
        b += B.conj().T @ y[r0:r1]
        del B
    lam = ridge * (np.trace(G).real / n_col)
    G[np.diag_indices_from(G)] += lam
    try:
        w = np.linalg.solve(G, b)
    except np.linalg.LinAlgError:
        w, *_ = np.linalg.lstsq(G, b, rcond=None)

    clean = y.copy()
    for r0 in range(0, n, rows):
        r1 = min(r0 + rows, n)
        B = _eca_block(ref_win, dop_phase, r0, r1).astype(np.complex128)
        clean[r0:r1] -= B @ w
        del B
    return clean.astype(surv.dtype, copy=False)


def cancellation_db(before, after) -> float:
    """Power ratio before/after, in dB.  The one number that says if ECA worked."""
    pb = float(np.mean(np.abs(before) ** 2))
    pa = float(np.mean(np.abs(after) ** 2))
    if pa <= 0:
        return float("inf")
    return 10.0 * math.log10(pb / pa)


# ---------------------------------------------------------------------------
# 4 -- CFAR: fix the false-alarm rate before looking
# ---------------------------------------------------------------------------

def cfar_threshold_map(rd_map, pfa=1e-6, guard=(2, 2), train=(8, 8)):
    """Return `(noise_est, alpha)` for 2-D cell-averaging CFAR.

    `noise_est[i,j]` is the mean power of the training ring around cell
    (i,j): the full (2*(g+t)+1) window minus the (2*g+1) guard box.  Guard
    cells exist so a strong target does not raise its own threshold by
    leaking into the estimate -- self-masking, which is how a real aircraft
    gets thresholded away.
    """
    gr, gd = int(guard[0]), int(guard[1])
    tr, td = int(train[0]), int(train[1])
    so_r, so_d = 2 * (gr + tr) + 1, 2 * (gd + td) + 1
    sg_r, sg_d = 2 * gr + 1, 2 * gd + 1
    n_out, n_gd = so_r * so_d, sg_r * sg_d
    n_train = n_out - n_gd
    if n_train < 1:
        raise ValueError("train window must add cells beyond the guard box")

    # float64: the direct-path peak can be 1e9 x the floor, and a float32
    # box sum next to it loses the floor entirely.
    m = np.asarray(rd_map, dtype=np.float64)
    # mode='nearest' clamps at the edges.  Not 'reflect': reflecting folds a
    # target back into its own training ring at the map edge.
    s_out = uniform_filter(m, size=(so_r, so_d), mode="nearest") * n_out
    s_gd = uniform_filter(m, size=(sg_r, sg_d), mode="nearest") * n_gd
    noise = np.maximum((s_out - s_gd) / n_train, 0.0)

    # CA-CFAR scaling for exponentially distributed power (the |.|^2 of a
    # complex-Gaussian cell): P_fa = (1 + alpha/n_train)^-n_train, so
    #   alpha = n_train * (pfa^(-1/n_train) - 1).
    # This is per-cell, not per-map.
    alpha = n_train * (pfa ** (-1.0 / n_train) - 1.0)
    return noise, float(alpha)


def detect(rd_map, range_axis, doppler_axis, pfa=1e-6, guard=(2, 2),
           train=(8, 8), min_separation=None):
    """CA-CFAR detections on a range-Doppler map.

    Returns a list of dicts sorted by descending SNR:
        {"range_m", "doppler_hz", "snr_db", "power", "threshold",
         "bin_idx": (i_range, i_doppler)}

    WHY CFAR AND NOT A FIXED THRESHOLD
    ----------------------------------
    The floor of a passive-radar map is not flat.  It rises near zero
    Doppler (clutter residue), near zero range (direct-path residue), and
    with any interferer that drifts through the band.  A fixed threshold
    tuned on quiet cells fires continuously in the loud ones; tuned on the
    loud ones it goes deaf everywhere else.  CFAR estimates the local floor
    from a ring around each cell and thresholds RELATIVE to it, so the
    false-alarm probability is `pfa` per cell wherever you stand.

    Expected false alarms in a clean map are pfa * rd_map.size -- state
    that number, do not assume it is zero.  At pfa=1e-6 on a 461x64 map it
    is 0.03; on a 1400x512 map it is 0.7.

    Only local maxima are returned: a real target lights several adjacent
    cells, and reporting each as its own "aircraft" would inflate the count
    without adding information.  `min_separation` (bins, per axis) defaults
    to the guard box.
    """
    rd_map = np.asarray(rd_map)
    if rd_map.shape != (len(range_axis), len(doppler_axis)):
        raise ValueError(f"rd_map {rd_map.shape} does not match axes "
                         f"({len(range_axis)}, {len(doppler_axis)})")
    noise, alpha = cfar_threshold_map(rd_map, pfa=pfa, guard=guard, train=train)
    thresh = alpha * noise
    m = np.asarray(rd_map, dtype=np.float64)

    over = m > thresh
    if not over.any():
        return []
    sep = (2 * int(guard[0]) + 1, 2 * int(guard[1]) + 1) if min_separation is None \
        else (2 * int(min_separation[0]) + 1, 2 * int(min_separation[1]) + 1)
    peak = m >= maximum_filter(m, size=sep, mode="nearest")
    hits = np.argwhere(over & peak)

    out = []
    for i, j in hits:
        nz = noise[i, j]
        snr = 10.0 * math.log10(m[i, j] / nz) if nz > 0 else float("inf")
        out.append({
            "range_m": float(range_axis[i]),
            "doppler_hz": float(doppler_axis[j]),
            "snr_db": float(snr),
            "power": float(m[i, j]),
            "threshold": float(thresh[i, j]),
            "bin_idx": (int(i), int(j)),
        })
    out.sort(key=lambda d: -d["snr_db"])
    return out


# ---------------------------------------------------------------------------
# 5 -- geometry
# ---------------------------------------------------------------------------

def bistatic_velocity(doppler_hz, f_carrier_hz=RF25_HZ):
    """Bistatic range-RATE in m/s from a measured Doppler shift in Hz.

        fd = -(1/lambda) * d/dt (Rt + Rr)     [Hz]
        =>  d/dt (Rt + Rr) = -fd * c / fc     [m/s]

    Sign convention: a target CLOSING the total path (tower->target->us
    getting shorter) gives a POSITIVE Doppler and a NEGATIVE range-rate.

    THE CAVEAT, WHICH MATTERS MORE THAN THE FORMULA
    -----------------------------------------------
    This is the rate of change of the SUM path Rt + Rr, not the target's
    ground speed and not a radial velocity toward us.  Converting it to
    aircraft speed needs the bistatic angle beta (the angle at the target
    between the tower and us):

        d/dt(Rt + Rr) = 2 * v * cos(delta) * cos(beta/2)

    where delta is the angle between the velocity vector and the bisector
    of beta.  We measure neither beta nor delta with one omnidirectional
    receiver and one illuminator -- the measurement is an ellipsoid, not a
    position (docs/FEASIBILITY.md, "Honest limits").  So:

      * a target crossing the baseline (beta -> 180 deg) shows near-ZERO
        bistatic Doppler however fast it is flying: cos(beta/2) -> 0.  A
        fast aeroplane can be Doppler-invisible.  That is geometry, not a
        broken detector.
      * quoting this number as "aircraft speed" is wrong by a factor of up
        to 2 in one direction and unboundedly in the other.  It is a
        range-rate.  Say range-rate.

    Fusing with ADS-B truth (the plan in FEASIBILITY step 3) works the
    other way round: predict fd from a known state vector and check our
    peak, never invert one measurement into a position.
    """
    return -np.asarray(doppler_hz, dtype=np.float64) * C_LIGHT / float(f_carrier_hz)


def range_rate_to_doppler(range_rate_ms, f_carrier_hz=RF25_HZ):
    """Inverse of bistatic_velocity: bistatic range-rate (m/s) -> Doppler (Hz)."""
    return -np.asarray(range_rate_ms, dtype=np.float64) * float(f_carrier_hz) / C_LIGHT


def bistatic_range_m(delay_s):
    """Extra path length (m) for an echo delay (s) measured against the direct path."""
    return np.asarray(delay_s, dtype=np.float64) * C_LIGHT


__all__ = [
    "C_LIGHT", "RF25_HZ", "FS_ATSC3", "LOAD_LIMIT_SAMPLES",
    "ECA_DIRECT_MAX_ENTRIES",
    "cs16_len", "iter_cs16", "load_cs16",
    "caf", "caf_resolution",
    "cancel_direct", "cancellation_db",
    "cfar_threshold_map", "detect",
    "bistatic_velocity", "range_rate_to_doppler", "bistatic_range_m",
]
