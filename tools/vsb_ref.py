"""vsb_ref.py -- ATSC 1.0 (8-VSB) as a passive-radar illuminator, with a
reference rebuilt from our own symbol decisions.

8-VSB is one carrier carrying 8-level PAM at Rs = 4.5 MHz x 684/286 =
10.762 Msym/s through a vestigial-sideband filter, plus a pilot tone at
the lower band edge. With one receiver there is no clean copy of the
transmitted signal, so we make one: lock to the pilot, recover the
transmitter's symbol clock from the segment syncs, equalize, slice to the
8 levels, and re-modulate. The decisions are the transmitted symbols (a
symbol error rate of 1e-3 costs the reference about -40 dB), so the
rebuilt signal carries no echoes, no multipath and no noise.

Pipeline, per ~1 s block (blocks are independent):
  1. resample the capture to 2 samples/symbol of the NOMINAL symbol rate
  2. find the pilot (FFT peak near -Rs/4), mix it to DC, track its phase
     with a 1 kHz low-pass and remove it
  3. segment syncs (+5 -5 -5 +5 every 832 symbols): matched-filter, fold,
     track, and fit position vs segment index -> the transmitter's clock
     (offset in ppm) and symbol phase
  4. windowed-sinc resample the block onto the TRANSMITTER's 2-sps grid,
     so the reference and the surveillance share one time base and the
     ECA never sees clock drift
  5. least-squares fractionally-spaced equalizer, trained first on the
     sync symbols, then decision-directed (batch), slice to 8 levels
  6. rebuild: symbols + pilot -> complex VSB pulse -> reference.
     A second, BLANKED reference (segment syncs, field-sync segments and
     the pilot removed) is what the CAF correlates against, because the
     periodic parts of 8-VSB make false targets: segment syncs every
     77.3 us -> range ghosts every 23.2 km; field syncs every 24.2 ms ->
     Doppler ghosts every 41.3 Hz; the pilot -> a zero-Doppler ridge.
     The ECA still cancels with the FULL reference.
  7. shift to channel centre, decimate by 3 -> 7.17 MS/s (41.8 m range
     cells), ECA, CAF, CFAR, link tracks (track_hunt's JSON shape, so
     truth_compare / truth_trajectory / unmatched_tracks run unchanged
     with --station=).

Modes:
  python tools/vsb_ref.py selftest
  python tools/vsb_ref.py survey lab/vsb_survey_rf*.cs16
  python tools/vsb_ref.py check <capture.cs16> [--t0 S]          (one block: lock, SER, cancellation)
  python tools/vsb_ref.py radar <capture.cs16> --station=<json> [--secs S] [--save-maps] [--out lab/tracks_V.json]
"""
import json
import math
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.signal import resample_poly, fftconvolve, oaconvolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pr_core as P                                            # noqa: E402

C = 2.998e8
RS_FR = Fraction(4500000 * 684, 286)          # exact symbol rate, Hz
RS = float(RS_FR)
FS2 = 2 * RS                                  # working rate: 2 samples/symbol
SEG = 832                                     # symbols per segment
FIELD = 313                                   # segments per field
SYNC = np.array([5.0, -5.0, -5.0, 5.0])
LEVELS = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=np.float32)
PILOT = 1.25
ALPHA = 0.1152
DECIM3 = 3                                    # FS2 / 3 = 7.1748 MS/s
FS3 = FS2 / DECIM3


# ------------------------------------------------------------------ pulses
def rrc(sps, alpha, span):
    """Root-raised-cosine, unit energy, `sps` samples per symbol, +-span symbols."""
    t = np.arange(-span * sps, span * sps + 1) / sps
    h = np.zeros_like(t)
    for i, ti in enumerate(t):
        if abs(ti) < 1e-12:
            h[i] = 1 - alpha + 4 * alpha / np.pi
        elif abs(abs(ti) - 1 / (4 * alpha)) < 1e-9:
            h[i] = (alpha / np.sqrt(2)) * ((1 + 2 / np.pi) * np.sin(np.pi / (4 * alpha))
                                           + (1 - 2 / np.pi) * np.cos(np.pi / (4 * alpha)))
        else:
            h[i] = (np.sin(np.pi * ti * (1 - alpha)) + 4 * alpha * ti * np.cos(np.pi * ti * (1 + alpha))) \
                / (np.pi * ti * (1 - (4 * alpha * ti) ** 2))
    return h / np.sqrt((h ** 2).sum())


def vsb_pulse(span=24):
    """Complex VSB pulse at FS2 (2 sps of Rs): an RRC whose Nyquist edge is
    Rs/4, shifted up by Rs/4, so the band runs from the pilot (DC) to Rs/2."""
    h = rrc(4, ALPHA, span)                        # 'symbol rate' Rs/2 sampled at 2Rs -> 4 sps
    n = np.arange(len(h)) - (len(h) - 1) / 2
    return (h * np.exp(2j * np.pi * n / 8)).astype(np.complex64)   # Rs/4 = FS2/8


H_VSB = vsb_pulse()


def modulate(sym2, h=H_VSB):
    """Symbols (one per symbol) -> complex baseband at FS2, pilot-at-DC domain."""
    up = np.zeros(2 * len(sym2), np.complex64)
    up[0::2] = sym2
    return oaconvolve(up, h, mode="same").astype(np.complex64)


# --------------------------------------------------------------- resampling
try:
    import numba

    @numba.njit(parallel=True, fastmath=True, cache=True)
    def _frac_interp_nb(xr, xi, tau, half, beta, i0b):
        n = tau.shape[0]
        outr = np.zeros(n, np.float32)
        outi = np.zeros(n, np.float32)
        nx = xr.shape[0]
        for m in numba.prange(n):
            t = tau[m]
            i0 = int(np.floor(t))
            mu = t - i0
            ar = 0.0
            ai = 0.0
            for k in range(-half + 1, half + 1):
                d = k - mu
                if abs(d) < 1e-12:
                    sc = 1.0
                else:
                    sc = np.sin(np.pi * d) / (np.pi * d)
                u = 1.0 - (d / half) ** 2
                if u <= 0.0:
                    continue
                # modified Bessel I0 by series
                xx = beta * np.sqrt(u) / 2.0
                term = 1.0
                ssum = 1.0
                for j in range(1, 25):
                    term *= (xx / j) ** 2
                    ssum += term
                w = sc * ssum / i0b
                idx = i0 + k
                if idx < 0:
                    idx = 0
                elif idx >= nx:
                    idx = nx - 1
                ar += w * xr[idx]
                ai += w * xi[idx]
            outr[m] = ar
            outi[m] = ai
        return outr, outi

    def frac_interp(x, tau, half=8, beta=8.0, chunk=None):
        """x evaluated at fractional sample positions tau (windowed sinc, 2*half taps)."""
        xr = np.ascontiguousarray(x.real, dtype=np.float32)
        xi = np.ascontiguousarray(x.imag, dtype=np.float32)
        r, i = _frac_interp_nb(xr, xi, np.ascontiguousarray(tau, dtype=np.float64), half, beta, float(np.i0(beta)))
        return (r + 1j * i).astype(np.complex64)

except ImportError:                                            # pragma: no cover
    def frac_interp(x, tau, half=8, beta=8.0, chunk=1 << 21):
        """x evaluated at fractional sample positions tau (windowed sinc, 2*half taps)."""
        out = np.zeros(len(tau), np.complex64)
        ks = np.arange(-half + 1, half + 1)
        i0b = np.i0(beta)
        for a in range(0, len(tau), chunk):
            t = tau[a:a + chunk]
            i0 = np.floor(t).astype(np.int64)
            mu = (t - i0).astype(np.float64)
            acc = np.zeros(len(t), np.complex128)
            for k in ks:
                d = k - mu
                w = np.sinc(d) * np.i0(beta * np.sqrt(np.clip(1 - (d / half) ** 2, 0, 1))) / i0b
                idx = np.clip(i0 + k, 0, len(x) - 1)
                acc += w * x[idx]
            out[a:a + chunk] = acc
        return out


def to_fs2(x, fs):
    fr = (2 * RS_FR / Fraction(fs)).limit_denominator(4000)
    return resample_poly(x, fr.numerator, fr.denominator).astype(np.complex64), float(fr)


# ------------------------------------------------------------------ stages
def pilot_lock(z, fs2=FS2, search_hz=60e3, lp_s=1e-3):
    """Mix the pilot (near -Rs/4) to DC and remove its phase (1 kHz tracking)."""
    n = len(z)
    nfft = 1 << int(math.ceil(math.log2(min(n, 1 << 22))))
    Z = np.fft.fft(z[:nfft])
    f = np.fft.fftfreq(nfft, 1 / fs2)
    sel = np.abs(f + RS / 4) < search_hz
    fp = f[sel][np.argmax(np.abs(Z[sel]))]
    # refine with a parabola on the peak
    k = np.where(sel)[0][np.argmax(np.abs(Z[sel]))]
    a, b, c = np.abs(Z[k - 1]), np.abs(Z[k]), np.abs(Z[k + 1])
    fp += 0.5 * (a - c) / (a - 2 * b + c) * (fs2 / nfft)
    m = np.exp(-2j * np.pi * fp * np.arange(n) / fs2).astype(np.complex64)
    zz = z * m
    L = max(3, int(lp_s * fs2))
    cs = np.concatenate([[0], np.cumsum(zz.astype(np.complex128))])
    h2 = L // 2
    lo = np.clip(np.arange(n) - h2, 0, n)
    hi = np.clip(np.arange(n) + (L - h2), 0, n)
    lp = ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.complex64)
    pil_db = 10 * np.log10(np.mean(np.abs(lp) ** 2) / max(np.mean(np.abs(zz) ** 2), 1e-30))
    zz = zz * np.conj(lp) / np.maximum(np.abs(lp), 1e-12)
    return zz.astype(np.complex64), float(fp), float(pil_db)


def segment_sync_fit(zz, group=16):
    """Fit the transmitter's segment-sync positions in a pilot-locked FS2-nominal
    block. One segment's sync stands only ~2 sigma above random data (sync
    metric 100 vs a data spread of ~46), so positions are measured on GROUPS
    of `group` segments folded along the current clock estimate, then fitted
    by a line, iterated. Returns (p0, period_samples, n_groups,
    rms_residual_samples, metric_db)."""
    y = np.real(oaconvolve(zz, np.conj(H_VSB[::-1]), mode="same"))
    c = (y[:-6] - y[2:-4] - y[4:-2] + y[6:]) * 5
    P0 = 2 * SEG
    nseg = len(c) // P0 - 2
    if nseg < 4 * group:
        raise ValueError("block too short for a sync fit")
    fold = np.zeros(P0)
    for s_ in range(min(40, nseg)):
        fold += c[s_ * P0:(s_ + 1) * P0]
    q, per = float(np.argmax(fold)), float(P0)
    offs = np.arange(-6, 7, dtype=float)

    def group_peak(g0, q, per):
        pos = q + np.arange(g0, g0 + group)[:, None] * per + offs[None, :]
        ok = (pos[:, 0] >= 0) & (pos[:, -1] < len(c) - 1)
        pos = pos[ok]
        b = np.floor(pos).astype(np.int64)
        f = pos - b
        acc = (c[b] * (1 - f) + c[b + 1] * f).sum(axis=0)
        j = int(np.argmax(acc))
        frac = 0.0
        if 0 < j < len(acc) - 1:
            a, bb, d = acc[j - 1], acc[j], acc[j + 1]
            den = a - 2 * bb + d
            frac = 0.5 * (a - d) / den if den != 0 else 0.0
        return offs[j] + frac, acc[j] / group

    # pass 1: track forward, refitting the line every group, so drift never
    # leaves the +-6 sample window (15 ppm = 80 samples over 0.25 s)
    gc, gp, gpk = [], [], []
    for g0 in range(0, nseg - group, group):
        if len(gc) >= 3:
            k1, k0 = np.polyfit(np.array(gc), np.array(gp), 1)
            q_eff, per_eff = k0, k1
        else:
            q_eff, per_eff = q, per
        d, pk = group_peak(g0, q_eff, per_eff)
        center = g0 + (group - 1) / 2
        gc.append(center); gp.append(q_eff + center * per_eff + d); gpk.append(pk)
    gc, gp = np.array(gc), np.array(gp)
    good = np.ones(len(gc), bool)
    for _ in range(3):
        k1, k0 = np.polyfit(gc[good], gp[good], 1)
        r = gp - (k0 + k1 * gc)
        good = np.abs(r) < max(3 * np.std(r[good]), 0.02)
    q, per = float(k0), float(k1)
    # pass 2: re-measure every group on the final line
    gd = np.array([group_peak(g0, q, per)[0] for g0 in range(0, nseg - group, group)])
    good = np.ones(len(gc), bool)
    for _ in range(3):
        k1, k0 = np.polyfit(gc[good], gd[good], 1)
        r = gd - (k0 + k1 * gc)
        good = np.abs(r) < max(3 * np.std(r[good]), 0.02)
    q += k0
    per += k1
    rms = float(np.sqrt(np.mean(r[good] ** 2)))
    metric_db = 10 * np.log10(np.median(gpk) ** 2 / max(np.median(c ** 2), 1e-30))
    return float(q), float(per), int(good.sum()), rms, float(metric_db)


def to_tx_grid(zz, p0, period):
    """Resample onto the transmitter's 2-sps grid; output sample 2k = symbol k,
    with symbol 0 = first sync symbol of some segment. Returns (zt, m0) where
    m0 is the output index of the first fitted sync."""
    r = period / (2 * SEG)                      # input samples per output sample
    m0 = 2 * SEG * int(math.ceil((p0 - 16) / period))   # leave room at the front
    m0 = max(m0, 0)
    start = p0 - m0 * r
    while start < 16:
        m0 -= 2 * SEG
        start = p0 - m0 * r
        if m0 < 0:
            m0 += 2 * SEG
            start = p0 - m0 * r
            break
    n_out = int((len(zz) - 16 - start) / r)
    n_out -= n_out % (2 * SEG)
    tau = start + np.arange(n_out) * r
    return frac_interp(zz, tau), m0, start, r


def equalize_and_slice(zt, m0, L=96, n_train=400_000, iters=2):
    """LS fractionally-spaced equalizer on complex zt (2 sps); real output at the
    symbol instants. Trained on sync symbols, then decision-directed (batch)."""
    n_sym = len(zt) // 2
    k = np.arange(n_sym)
    seg_pos = (k - m0 // 2) % SEG
    is_sync = seg_pos < 4
    sync_val = np.where(is_sync, SYNC[np.minimum(seg_pos, 3)], 0.0)
    c = L // 2

    def design(idx):
        m = 2 * idx
        X = np.stack([zt[np.clip(m - c + i, 0, len(zt) - 1)] for i in range(L)], axis=1)
        return np.hstack([X.real, -X.imag, np.ones((len(idx), 1), np.float32)]).astype(np.float64)

    def apply(w):
        wc = (w[:L] + 1j * w[L:2 * L]).astype(np.complex64)
        full = oaconvolve(zt, wc[::-1], mode="full")
        m = 2 * k
        return np.real(full[m - c + L - 1]) + w[-1]

    def ls(X, y):
        G = X.T @ X
        G[np.diag_indices_from(G)] += 1e-9 * np.trace(G) / len(G)
        return np.linalg.solve(G, X.T @ y)

    margin = L
    valid = (k > margin) & (k < n_sym - margin)
    tr = np.where(valid & is_sync)[0]
    w = ls(design(tr), sync_val[tr])
    u = apply(w)
    for _ in range(iters):
        a = LEVELS[np.clip(np.round((u + 7) / 2).astype(int), 0, 7)]
        a = np.where(is_sync, sync_val, a)
        idx = np.where(valid)[0]
        if len(idx) > n_train:
            idx = idx[np.linspace(0, len(idx) - 1, n_train).astype(int)]
        w = ls(design(idx), a[idx])
        u = apply(w)
    a = LEVELS[np.clip(np.round((u + 7) / 2).astype(int), 0, 7)].astype(np.float32)
    a = np.where(is_sync, sync_val, a).astype(np.float32)
    err = u[valid] - a[valid]
    snr_db = 10 * np.log10(np.mean(a[valid] ** 2) / np.mean(err ** 2))
    sync_err = float(np.mean(np.abs(u[valid & is_sync] - sync_val[valid & is_sync]) > 1.0))
    # SER estimate from a Gaussian error model at this SNR (8 levels)
    sigma = np.sqrt(np.mean(err ** 2))
    from math import erfc
    ser = 2 * (7 / 8) * 0.5 * erfc(1 / (sigma * np.sqrt(2)))
    return a, is_sync, dict(eq_snr_db=float(snr_db), sync_err=sync_err, ser_est=float(ser), sigma=float(sigma))


def field_sync_mask(a, is_sync, m0):
    """Symbols belonging to field-sync segments: the segment whose data repeats
    exactly FIELD segments later (PN511 is identical every field)."""
    n_sym = len(a)
    s0 = (m0 // 2) % SEG
    first = s0 if s0 == 0 else s0 - SEG + SEG
    seg_starts = np.arange(first % SEG, n_sym - SEG, SEG)
    mask = np.zeros(n_sym, bool)
    nseg = len(seg_starts)
    if nseg <= FIELD:
        return mask
    score = np.zeros(nseg)
    for i in range(nseg - FIELD):
        x = a[seg_starts[i] + 4:seg_starts[i] + 515]
        y = a[seg_starts[i + FIELD] + 4:seg_starts[i + FIELD] + 515]
        score[i] = np.mean(x == y)
    hits = np.where(score > 0.9)[0]
    for h in hits:
        for j in (h, h + FIELD):
            mask[seg_starts[j]:seg_starts[j] + SEG] = True
    return mask


def build_block(x, fs):
    """One block of capture -> (surv, ref_full, ref_blank) at FS3, plus diagnostics."""
    z, ratio = to_fs2(x, fs)
    zz, fp, pil_db = pilot_lock(z)
    p0, per, nseg, rms, mdb = segment_sync_fit(zz)
    ppm = (per / (2 * SEG) - 1) * 1e6
    zt, m0, start, r = to_tx_grid(zz, p0, per)
    a, is_sync, eq = equalize_and_slice(zt, m0)
    fmask = field_sync_mask(a, is_sync, m0)
    ref_full = modulate(a + PILOT)
    a_blank = np.where(is_sync | fmask, 0.0, a).astype(np.float32)
    ref_blank = modulate(a_blank)
    # to channel centre, decimate to FS3
    sh = np.exp(-2j * np.pi * np.arange(len(zt)) / 8).astype(np.complex64)
    surv3 = resample_poly(zt * sh, 1, DECIM3).astype(np.complex64)
    rf3 = resample_poly(ref_full * sh, 1, DECIM3).astype(np.complex64)
    rb3 = resample_poly(ref_blank * sh, 1, DECIM3).astype(np.complex64)
    t_first = (start / ratio) / fs           # capture time (s) of output sample 0, block-relative
    diag = dict(pilot_hz=fp, pilot_db=pil_db, clock_ppm=ppm, n_sync=nseg, sync_rms_samples=rms,
                sync_metric_db=mdb, field_sync_segments=int(fmask.sum() // SEG), t_first_s=t_first, **eq)
    return surv3, rf3, rb3, diag


# ----------------------------------------------------------------- selftest
def synth(seconds=1.0, fs=8e6, ppm=15.0, cfo=3100.0, seed=1):
    rng = np.random.default_rng(seed)
    n_sym = int(seconds * RS) + 4000
    a = rng.choice(LEVELS, n_sym).astype(np.float32)
    k = np.arange(n_sym)
    a[(k % SEG) < 4] = SYNC[(k % SEG)[(k % SEG) < 4]]
    pn = rng.choice([-5.0, 5.0], SEG - 4).astype(np.float32)
    for s in range(0, n_sym // SEG, FIELD):
        a[s * SEG + 4:(s + 1) * SEG] = pn
    tx = modulate(a + PILOT)                                   # pilot-DC domain, FS2
    tx = tx * np.exp(-2j * np.pi * np.arange(len(tx)) / 8)     # to channel centre
    # channel: direct + 2 static multipaths + a moving target 60 dB down
    def delayed(x, d_samples, dop=0.0):
        tau = np.arange(len(x)) - d_samples
        y = frac_interp(x, np.clip(tau, 0, len(x) - 1))
        if dop:
            y = y * np.exp(2j * np.pi * dop * np.arange(len(x)) / FS2)
        return y
    tgt_km, tgt_hz = 20.0, 150.0
    d_t = tgt_km * 1e3 / C * FS2
    rx = tx + 0.3 * delayed(tx, 7.3) + 0.1 * delayed(tx, 23.6) \
        + 10 ** (-60 / 20) * delayed(tx, d_t, tgt_hz)
    rx = rx * np.exp(2j * np.pi * cfo * np.arange(len(rx)) / FS2)
    # receiver clock: samples at n / (fs (1+ppm)) true time -> tx index
    n_out = int(seconds * fs)
    tau = np.arange(n_out) * FS2 / (fs * (1 + ppm * 1e-6)) + 50.0
    y = frac_interp(rx.astype(np.complex64), tau)
    pw = np.mean(np.abs(y) ** 2)
    y = y + (rng.standard_normal(n_out) + 1j * rng.standard_normal(n_out)).astype(np.complex64) \
        * np.sqrt(pw * 10 ** (-35 / 10) / 2)
    return y.astype(np.complex64), dict(ppm=ppm, cfo=cfo, tgt_km=tgt_km, tgt_hz=tgt_hz)


def selftest():
    print("[self] synthesizing 1.0 s of 8-VSB: +15 ppm clock, +3.1 kHz carrier offset,")
    print("       multipath -10.5 dB @ 0.34 us and -20 dB @ 1.1 us, noise -35 dB, target -60 dB @ 20 km / +150 Hz")
    x, truth = synth()
    t = time.time()
    surv, rf, rb, d = build_block(x, 8e6)
    print(f"[self] build {time.time() - t:.1f} s: pilot {d['pilot_hz']:.0f} Hz (expect {-RS/4 + truth['cfo']:.0f}), "
          f"clock {d['clock_ppm']:+.2f} ppm (expect {truth['ppm']:+.2f}), sync rms {d['sync_rms_samples']:.3f} samples, "
          f"eq SNR {d['eq_snr_db']:.1f} dB, SER~{d['ser_est']:.1e}, sync errors {d['sync_err']:.4f}, "
          f"field-sync segments {d['field_sync_segments']}")
    gates = []
    gates.append(("clock recovered within 0.5 ppm", abs(d["clock_ppm"] - truth["ppm"]) < 0.5))
    f_exp = (-RS / 4 + truth["cfo"]) / (1 + truth["ppm"] * 1e-6)   # a fast receiver clock reads frequencies low
    gates.append(("pilot within 20 Hz (in receiver-clock units)", abs(d["pilot_hz"] - f_exp) < 20))
    gates.append(("equalized SNR > 25 dB", d["eq_snr_db"] > 25))
    n = min(len(surv), 1_200_000)
    s = surv[200:n]; r_full = rf[200:n]; r_bl = rb[200:n]
    t = time.time()
    clean = P.cancel_direct(s, r_full, 64)
    canc = P.cancellation_db(s, clean)
    print(f"[self] ECA 64 taps: {canc:.1f} dB in {time.time() - t:.1f} s")
    gates.append(("cancellation > 30 dB (noise floor is -35 dB)", canc > 30))
    rd, rax, dax = P.caf(r_bl, clean, FS3, 40e3, 1000.0, 512)
    dets = P.detect(rd, rax, dax, pfa=1e-6)
    dets = [q for q in dets if abs(q["doppler_hz"]) > 15]
    best = dets[0] if dets else None
    if best:
        print(f"[self] strongest off-notch detection: {best['range_m']/1e3:.2f} km {best['doppler_hz']:+.1f} Hz "
              f"{best['snr_db']:.1f} dB; {len(dets)} detections total")
    ok_t = bool(best and abs(best["range_m"] / 1e3 - truth["tgt_km"]) < 0.1 and abs(best["doppler_hz"] - truth["tgt_hz"]) < 8)
    gates.append(("target recovered at 20 km / +150 Hz as the strongest detection", ok_t))
    ghosts = [q for q in dets if abs(abs(q["range_m"] / 1e3 - truth["tgt_km"]) - 23.2) < 0.2]
    gates.append(("no 23.2 km segment-sync ghost of the target", not ghosts))
    for g, ok in gates:
        print(f"  [{'PASS' if ok else 'FAIL'}] {g}")
    return 0 if all(ok for _, ok in gates) else 1


# -------------------------------------------------------------------- survey
def survey(paths):
    rows = []
    for p in paths:
        side_p = os.path.splitext(p)[0] + ".json"
        if not os.path.exists(side_p):
            continue
        side = json.load(open(side_p))
        fs = float(side["rate_hz"])
        n = min(os.path.getsize(p) // 4, int(fs))
        raw = np.fromfile(p, dtype=np.int16, count=2 * n)
        x = (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0
        nfft = 1 << 18
        S = np.zeros(nfft)
        for k in range(len(x) // nfft):
            S += np.abs(np.fft.fft(x[k * nfft:(k + 1) * nfft] * np.hanning(nfft))) ** 2
        f = np.fft.fftfreq(nfft, 1 / fs)
        sel = np.abs(f + RS / 4) < 40e3
        pk = S[sel].max()
        loc = np.median(S[sel])
        band = S[np.abs(f) < 2.5e6].mean()
        edge = S[(np.abs(f) > 3.3e6) & (np.abs(f) < 0.48 * fs)].mean()
        rf = os.path.basename(p).split("rf")[-1].split(".")[0]
        rows.append((rf, 10 * np.log10(pk / loc), 10 * np.log10(band / edge),
                     side["headroom"]["peak_dbfs"], side["verdict"]))
    print("RF   pilot over local floor   band over edges   peak dBFS   verdict")
    for r in sorted(rows, key=lambda r: -r[1]):
        print(f"{r[0]:>3}   {r[1]:6.1f} dB               {r[2]:6.1f} dB        {r[3]:6.1f}     {r[4]}")


# ---------------------------------------------------------------------- check
def load(path, fs, t0, secs):
    n0, n = int(t0 * fs), int(secs * fs)
    raw = np.fromfile(path, dtype=np.int16, count=2 * n, offset=4 * n0)
    if len(raw) < 2 * n:
        return None
    return (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0


def check(path, t0=1.0):
    side = json.load(open(os.path.splitext(path)[0] + ".json"))
    fs = float(side["rate_hz"])
    x = load(path, fs, t0, 1.0)
    t = time.time()
    surv, rf, rb, d = build_block(x, fs)
    print(f"[check] {os.path.basename(path)} t0={t0}s: build {time.time()-t:.1f}s")
    for k, v in d.items():
        print(f"   {k:20s} {v}")
    for taps in (32, 64, 128):
        n = min(len(surv), 1_200_000)
        clean = P.cancel_direct(surv[200:n], rf[200:n], taps)
        print(f"   ECA {taps:3d} taps: {P.cancellation_db(surv[200:n], clean):.1f} dB")


# ---------------------------------------------------------------------- radar
def _radar_block(job):
    """One block: build the reference, then ECA + CAF + CFAR per dwell."""
    (path, fs, t_blk, block, dwell_n, taps, max_range_km, doppler, n_doppler, pfa, notch, maps_dir) = job
    x = load(path, fs, t_blk, block + 0.02)
    if x is None:
        return t_blk, None, [], [], "short read"
    try:
        surv, rf, rb, d = build_block(x, fs)
    except Exception as e:                                     # a bad block is VOID, not fatal
        return t_blk, None, [], [], f"{type(e).__name__}: {e}"
    bin_m = C / FS3
    dets, canc = [], []
    start = 200
    k = 0
    while start + dwell_n <= len(surv) - 200 and (start + dwell_n) / FS3 + d["t_first_s"] <= block + 0.01:
        s_ = surv[start:start + dwell_n]
        clean = P.cancel_direct(s_, rf[start:start + dwell_n], taps)
        canc.append(float(P.cancellation_db(s_, clean)))
        rd, rax, dax = P.caf(rb[start:start + dwell_n], clean, FS3, max_range_km * 1e3, doppler, n_doppler)
        q = P.detect(rd, rax, dax, pfa=pfa)
        t_mid = t_blk + d["t_first_s"] + (start + dwell_n / 2) / FS3
        for r in q:
            if abs(r["doppler_hz"]) <= notch or r["range_m"] < 2 * bin_m:
                continue
            dets.append(dict(t_s=round(t_mid, 4), range_m=float(r["range_m"]), doppler_hz=float(r["doppler_hz"]),
                             snr_db=float(r["snr_db"]), power=float(r["power"]), threshold=float(r["threshold"])))
        if maps_dir:
            noise, _ = P.cfar_threshold_map(rd, pfa=pfa)
            np.save(os.path.join(maps_dir, f"rd_{t_mid:09.4f}.npy"),
                    (10 * np.log10(rd + 1e-30) - 10 * np.log10(noise + 1e-30)).astype(np.float16))
            if t_blk == 0.0 and k == 0 or not os.path.exists(os.path.join(maps_dir, "axes.npy")):
                np.save(os.path.join(maps_dir, "axes.npy"), np.array([rax, None, dax], dtype=object), allow_pickle=True)
        start += dwell_n
        k += 1
    d["t_block"] = t_blk
    return t_blk, d, dets, canc, None


def radar(argv):
    import argparse
    from concurrent.futures import ProcessPoolExecutor
    ap = argparse.ArgumentParser()
    ap.add_argument("capture")
    ap.add_argument("--station", required=True)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--secs", type=float, default=None)
    ap.add_argument("--block", type=float, default=1.0)
    ap.add_argument("--dwell-samples", type=int, default=1_180_000)
    ap.add_argument("--taps", type=int, default=160)   # 6.7 km: a 5.7 km static reflector leaked onto the 41.3 Hz field comb at 64
    ap.add_argument("--max-range-km", type=float, default=60.0)
    ap.add_argument("--doppler", type=float, default=1000.0)
    ap.add_argument("--n-doppler", type=int, default=512)
    ap.add_argument("--pfa", type=float, default=1e-6)
    ap.add_argument("--notch-hz", type=float, default=15.0)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--save-maps", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    side = json.load(open(os.path.splitext(a.capture)[0] + ".json"))
    fs = float(side["rate_hz"])
    secs = a.secs or (os.path.getsize(a.capture) / 4 / fs - a.t0 - 0.05)
    out = a.out or os.path.join(os.path.dirname(a.capture), "tracks_V.json")
    maps_dir = None
    if a.save_maps:
        maps_dir = os.path.join(os.path.dirname(out), "maps", os.path.splitext(os.path.basename(out))[0])
        os.makedirs(maps_dir, exist_ok=True)
    bin_m = C / FS3
    lam = C / float(json.load(open(a.station))["fc_hz"])
    print(f"[vsb ] {os.path.basename(a.capture)}  {secs:.0f} s in {a.block:.1f} s blocks x {a.jobs} jobs; "
          f"FS3 {FS3/1e6:.4f} MS/s (bin {bin_m:.1f} m), dwell {a.dwell_samples/FS3*1e3:.0f} ms, ECA {a.taps} taps", flush=True)
    jobs = []
    t_blk = a.t0
    while t_blk + a.block <= a.t0 + secs + 1e-6:
        jobs.append((a.capture, fs, round(t_blk, 6), a.block, a.dwell_samples, a.taps, a.max_range_km,
                     a.doppler, a.n_doppler, a.pfa, a.notch_hz, maps_dir))
        t_blk += a.block
    dets, diags, canc = [], [], []
    tw = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for t_b, d, dd, cc, err in ex.map(_radar_block, jobs):
            if d is None:
                print(f"   block t={t_b:6.2f}s VOID ({err})", flush=True)
                continue
            diags.append(d); dets += dd; canc += cc
            best = max(dd, key=lambda k: k["snr_db"], default=None)
            print(f"   block t={t_b:6.2f}s clock {d['clock_ppm']:+6.2f} ppm eqSNR {d['eq_snr_db']:5.1f} dB "
                  f"SER~{d['ser_est']:.0e} canc {np.median(cc) if cc else 0:5.1f} dB  dets {len(dd):3d}"
                  + (f"  best {best['snr_db']:4.1f} dB @ {best['range_m']/1e3:6.2f} km {best['doppler_hz']:+6.0f} Hz" if best else ""),
                  flush=True)
    import fm_radar as F
    tracks = F.link_tracks(dets, bin_m, lam, max_dr_bins=4, max_df=20.0, min_n=3)
    print(f"[vsb ] {len(canc)} dwells from {len(diags)} blocks in {time.time()-tw:.0f} s; cancellation median "
          f"{np.median(canc):.1f} dB; {len(dets)} detections -> {len(tracks)} tracks")
    for i, tr in enumerate(tracks[:12]):
        print(f"  track {i}: {tr['n']} dwells, {tr['span_s']:.1f}s, {tr['range0_m']/1e3:.2f}->{tr['range1_m']/1e3:.2f} km, "
              f"fd {tr['doppler_mean_hz']:+.1f} Hz | observed {tr['observed_rate_ms']:+.1f} vs predicted "
              f"{tr['predicted_rate_ms']:+.1f} m/s, max {max(tr['snr_db']):.1f} dB")
    json.dump(dict(capture=os.path.abspath(a.capture), fs_hz=FS3, illuminator="ATSC 1.0 8-VSB, decision-rebuilt reference",
                   taps=a.taps, cancellation_db=canc, blocks=diags,
                   n_detections=len(dets), tracks=tracks, detections=dets), open(out, "w"), indent=1)
    print(f"-> {out}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if mode == "selftest":
        sys.exit(selftest())
    elif mode == "survey":
        survey(sys.argv[2:])
    elif mode == "check":
        t0 = float(sys.argv[sys.argv.index("--t0") + 1]) if "--t0" in sys.argv else 1.0
        check(sys.argv[2], t0)
    elif mode == "radar":
        radar(sys.argv[2:])
    else:
        print(__doc__)
