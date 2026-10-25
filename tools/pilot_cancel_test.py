#!/usr/bin/env python3
"""pilot_cancel_test.py -- does a synthesized pilot lattice kill the comb?

THE EXPERIMENT, IN ONE SENTENCE
-------------------------------
Build the ATSC 3.0 pilot lattice from A/322 tables alone (`pr_pilots.py`),
align it to a banked RF25 capture on the capture's own raw clock, subtract
it, and ask what happened to the 7.40 km / +198.5 Hz comb that dominates
every range-Doppler map broadcast-radar has ever made.

Why this is the entry ticket and not a nice-to-have: FIRST_LOOK.md proved
the comb is signal-borne (absent from the antenna-disconnected VOID
capture) and gain-invariant (identical at rfgain 2 and 4), and identified
it to four figures as the transmitter's scattered-pilot lattice.  A passive
radar cannot see past its own illuminator's ambiguity function.  If pilots
alone -- no FEC, no constellation decisions, no LDM -- knock that comb
down, the full decode-remodulate build in REENCODER.md is worth its three
days.  If they do not, we learned that for 300 lines.

WHAT THIS FILE EXISTS TO PREVENT
--------------------------------
1.  **A number without a control.**  Every cancellation figure here is
    reported beside the SAME pipeline run with SCRAMBLED pilot values --
    identical lattice geometry, identical per-class amplitudes,
    uncorrelated signs.  A wrong pilot field that cancels as well as the
    right one means we are fitting noise and the experiment proves
    nothing.  (This is the e69 corrupt-one-bit control, promoted from the
    cell domain to the waveform.)
2.  **A remembered geometry.**  Every parameter comes from
    lab/l1_radar_capture.json -- this capture's OWN decoded L1 -- and the
    frame anchor comes from its own bootstraps.  REENCODER.md Step 0 is
    explicit: never trust a remembered config over the frame's L1.
3.  **A silent alignment failure.**  Alignment is scored before it is
    used: the l mod Dy hypothesis predicted by the frame anchor has to WIN
    against the three wrong ones, and the LS-fit-one-H residual has to
    clear the bar REENCODER.md set for it (pilots alone >> 20 dB).

WHAT IT FOUND (2026-08-20), STATED UP FRONT BECAUSE IT REVERSES THE PLAN
-----------------------------------------------------------------------
Subtracting the pilots is the WRONG operation.  The comb is a
spectral-FLATNESS artifact -- R(nfft/48) is the DFT of the power spectrum
at 1/48, so what combs is pilot carriers holding a DIFFERENT power from
data carriers, bump or hole alike.  Full removal turns a +4.5 dB bump into
a -15 dB hole and buys only 3.6 dB (28.2 -> 24.6 dB of CFAR SNR on the
first tooth) despite suppressing pilot-carrier power by 19.4 dB.
Subtracting the EXCESS instead (A - 1 per class,
`pr_pilots.pilot_excess_amplitudes`) levels the lattice to the data cells
and takes the tooth to 4.8 dB -- 23.4 dB off -- while the scrambled
control moves it 0.3 dB.  Both are run here, side by side, because the
pair is the evidence and the surprising half is the one that matters.

Usage:
  python tools/pilot_cancel_test.py \
      Z:/src/atsc3/data/fox25/rf25_fox_g4_0810_2254.cs16 \
      [--t0 30] [--frames 3] [--out lab/pilot_cancel.json]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pr_core as PR                                           # noqa: E402
import pr_pilots as PP                                         # noqa: E402

# pr_pilots has already put the atsc3 tree on sys.path (read-only).
import spec_pilots as P                                        # noqa: E402
from atsc3 import bootstrap as bs                              # noqa: E402
from m2_stage_a import find_bootstraps                         # noqa: E402
from m2_stage_b import resample_to                             # noqa: E402

FS = PP.FS_ATSC3          # samples/s -- the capture rate, and the decoder's
                          # no-op resample path (REENCODER.md sample-rate note)
FC = PR.RF25_HZ           # Hz -- RF25 centre, 539.0 MHz

# Processing parameters copied from pr_banked.py so BEFORE reproduces the
# first-look run exactly.  Changing one of these changes what "before"
# means, which is why they are named here and echoed into the JSON.
MAX_RANGE_M = 60e3
DOPPLER_HZ = 1000.0
N_DOPPLER = 512
ECA_TAPS = 96             # the self-referenced ECA of the v0 pipeline
PFA = 1e-6
MIN_DOPPLER_HZ = 15.0

PILOT_TAPS = 24           # taps for the PILOT basis.  8 already saturates
                          # (measured: 8 -> 19.47 dB, 128 -> 19.63 dB of
                          # pilot suppression), so this is the clutter
                          # horizon with margin and not a licence to eat
                          # targets.  See PILOT_GUARD for why 24, not 16.

PILOT_GUARD = 8           # samples of NEGATIVE delay the pilot basis must
                          # cover, and the single most expensive lesson of
                          # this build.  `cancel_direct` spans delays
                          # 0..n_taps-1, so a channel whose first arrival
                          # lands EARLIER than the reference is not in the
                          # subspace at all.  The bootstrap anchor is
                          # quantised at 6.144 MS/s = 1.125 samples of our
                          # clock, so the sign of that fractional error is
                          # luck: measured direct-path delays over three
                          # consecutive frames were +1.64, +0.25 and -0.08
                          # samples, and pilot suppression tracked them at
                          # 19.5, 10.9 and 7.1 dB.  Generating the
                          # reference GUARD samples ahead moves the whole
                          # tap window to -GUARD..n_taps-1-GUARD and all
                          # three frames land on 19.2-19.6 dB.  A result
                          # that varies frame to frame at a constant
                          # alignment score is a bookkeeping bug, not
                          # propagation (REENCODER risk 1).


# ---------------------------------------------------------------------------
# 1 -- geometry, from the capture's own decoded L1
# ---------------------------------------------------------------------------

def geometry_from_l1(path):
    """Every downstream number flows from this file, never from memory."""
    with open(path) as f:
        j = json.load(f)
    b = j["l1_basic"]
    if b["L1B_papr_reduction"] != 0:
        raise SystemExit("L1B_papr_reduction != 0: tone reservation makes "
                         "some carriers unknowable; resynthesis is not clean")
    if b["L1B_first_sub_miso"] != 0:
        raise SystemExit("L1B_first_sub_miso != 0: unmodeled precoder")
    fft_size = {0: 8192, 1: 16384, 2: 32768}[b["L1B_first_sub_fft_size"]]
    import spec_l1syntax as LS
    g = dict(
        nfft=fft_size,
        gi=LS.GUARD_INTERVAL_SAMPLES[b["L1B_first_sub_guard_interval"]],
        cred=b["L1B_first_sub_reduced_carriers"],
        pattern=P.PATTERNS[b["L1B_first_sub_scattered_pilot_pattern"]],
        boost=b["L1B_first_sub_scattered_pilot_boost"],
        # BOTH counts are "+1 encoded" fields (m6_cells.py:369-372).  The
        # closed check is below: bootstrap + (np + nsym) symbols must equal
        # the measured bootstrap-to-bootstrap spacing.
        n_sym=b["L1B_first_sub_num_ofdm_symbols"] + 1,
        np_sym=b["L1B_preamble_num_symbols"] + 1,
        sbs_first=bool(b["L1B_first_sub_sbs_first"]),
        sbs_last=bool(b["L1B_first_sub_sbs_last"]),
    )
    g["T"] = g["nfft"] + g["gi"]
    g["frame_samples"] = (PP.BOOTSTRAP_SAMPLES_POST
                          + (g["np_sym"] + g["n_sym"]) * g["T"])
    sbs = []
    if g["sbs_first"]:
        sbs.append(0)
    if g["sbs_last"]:
        sbs.append(g["n_sym"] - 1)
    g["sbs_l"] = tuple(sbs)
    return g, j


# ---------------------------------------------------------------------------
# 2 -- anchoring: where does a frame start, and how fast is their clock?
# ---------------------------------------------------------------------------

def find_frames(capture, t0_s, n_frames, g, span_pad_s=0.32):
    """Absolute sample index of each frame's FIRST DATA SYMBOL, from air.

    The bootstrap is A/321's fully-known waveform and the one thing in an
    ATSC 3.0 frame we can find without knowing anything else, so it is the
    anchor.  Its spacing is also a free measurement of the transmitter's
    clock against ours -- and a closed check on the frame layout, because
    bootstrap + (Np + Nsym) symbols has to equal it.
    """
    need = int((n_frames + 1) * g["frame_samples"] + span_pad_s * FS)
    off = int(t0_s * FS)
    head = PR.load_cs16(capture, n_samples=need, offset=off)
    hits = find_bootstraps(resample_to(head, FS, bs.FS))
    if len(hits) < 2:
        raise SystemExit(f"found {len(hits)} bootstraps near t={t0_s}s; "
                         "cannot anchor")
    pos = [off + h["position"] * FS / bs.FS for h in hits]
    spacing = float(np.mean(np.diff(pos)))
    err = spacing - g["frame_samples"]
    # A frame layout that is off by even one symbol shows up here as an
    # 8704-sample error, so this is a real gate and not decoration.
    if abs(err) > g["T"] / 2:
        raise SystemExit(
            f"measured frame spacing {spacing:.1f} samples but the L1 says "
            f"{g['frame_samples']} (= bootstrap + {g['np_sym']} preamble + "
            f"{g['n_sym']} data symbols): off by {err / g['T']:.2f} symbols")
    frames = []
    for h, p in zip(hits, pos):
        d0 = int(round(p)) + PP.BOOTSTRAP_SAMPLES_POST + g["np_sym"] * g["T"]
        frames.append(dict(bootstrap_abs=int(round(p)), data0_abs=d0,
                           cfo_hz=float(h["fine_cfo_hz"]),
                           peak_ratio=float(h["mean_peak_ratio"])))
    return frames[:n_frames], spacing, err


# ---------------------------------------------------------------------------
# 3 -- alignment scoring.  Fit-free first, then a fitted objective.
# ---------------------------------------------------------------------------

class Aligner:
    """Holds the per-symbol FFTs of one dwell and scores alignments on them."""

    def __init__(self, x, g, signs):
        self.x = x
        self.g = g
        self.noc = P.NOC[(g["nfft"], g["cred"])]
        self.bins = PP.symbol_bins(g["nfft"], g["cred"])
        # signed frequency index per carrier -- the phase ramp a fractional
        # timing offset puts across the band
        self.freq = (np.fft.fftfreq(g["nfft"]) * g["nfft"])[self.bins]
        self.signs = signs
        self.n_sym = min(int(x.size // g["T"]), g["n_sym"])

    def fields(self, amps, l0):
        """pilot field per emitted symbol index, cached by (l mod Dy, sbs)."""
        g = self.g
        out = []
        for m in range(self.n_sym):
            l = l0 + m
            out.append(PP.pilot_field(g["nfft"], g["cred"], g["pattern"], l,
                                      l in g["sbs_l"], amps, self.signs))
        return out

    def _symbol_ffts(self, cfo, period):
        g = self.g
        for m in range(self.n_sym):
            w = int(round(m * period))
            s = w + g["gi"]
            if s + g["nfft"] > self.x.size:
                return
            seg = self.x[s:s + g["nfft"]]
            n = np.arange(s, s + g["nfft"])
            yield m, w, np.fft.fft(seg * np.exp(-2j * np.pi * cfo * n / FS))

    def coherence(self, cfo, period, l0, amps, n_max=16):
        """Fit-free score of the l mod Dy hypothesis (m3_preamble's statistic).

        At a scattered pilot the transmitted value is a known real +-A, so
        Y/p leaves A*H(k), which is SMOOTH along the regularly-spaced
        scattered lattice.  A wrong lattice phase turns three quarters of
        the "pilots" into data cells and the quotient goes rough.  Nothing
        is fitted, so the statistic cannot be inflated by searching -- and
        the hypothesis the frame anchor PREDICTS is the one that has to
        win.
        """
        g = self.g
        dx, dy = P.dxdy(g["pattern"])
        num = den = 0.0
        for m, _w, Y in self._symbol_ffts(cfo, period):
            if m >= n_max:
                break
            l = l0 + m
            k = np.arange(dx * (l % dy), self.noc, dx * dy)
            q = Y[self.bins[k]] / (amps["scattered"] * self.signs[k])
            num += abs(np.sum(q[:-1] * np.conj(q[1:])))
            den += float(np.sum(np.abs(q) ** 2))
        return num / max(den, 1e-30)

    def ls_one_h(self, cfo, period, l0, amps, start_delay=0.0):
        """dB: pilot-carrier power over the residual after fitting ONE H per
        carrier for the whole dwell.

        This is the gate REENCODER.md wrote for this build ("LS-fit-one-H
        residual on air must drop >= N dB; pilots alone >> 20 dB"), and it
        doubles as the alignment objective because it is smooth in both
        CFO and symbol period -- the fractional-timing phase ramp is
        modelled explicitly, exactly as the synthesizer applies it, so a
        sub-sample clock drift moves the score instead of hiding inside an
        integer window start.
        """
        F = self.fields(amps, l0)
        num = np.zeros(self.noc, dtype=np.complex128)
        den = np.zeros(self.noc)
        keep = []
        for m, w, Y in self._symbol_ffts(cfo, period):
            Yc = Y[self.bins]
            p = F[m] * np.exp(-2j * np.pi * self.freq
                              * (m * period - w + start_delay) / self.g["nfft"])
            num += Yc * np.conj(p)
            den += np.abs(p) ** 2
            keep.append((Yc, p))
        if not keep:
            return -np.inf
        H = np.where(den > 0, num / np.maximum(den, 1e-30), 0.0)
        res = tot = 0.0
        for Yc, p in keep:
            msk = np.abs(p) > 0
            res += float(np.sum(np.abs(Yc[msk] - H[msk] * p[msk]) ** 2))
            tot += float(np.sum(np.abs(Yc[msk]) ** 2))
        return 10.0 * math.log10(tot / max(res, 1e-30))

    def direct_delay(self, cfo, period, l0, amps, n_max=8):
        """Direct-path delay in samples, relative to our assumed symbol start.

        Read off the phase ramp along the regularly-spaced scattered
        lattice: a delay of tau samples tilts adjacent pilots (Dx*Dy = 48
        carriers apart) by -2*pi*48*tau/nfft.  Reported because its SIGN
        decides whether the ECA basis can represent the channel at all --
        see PILOT_GUARD.
        """
        g = self.g
        dx, dy = P.dxdy(g["pattern"])
        acc = 0j
        for m, _w, Y in self._symbol_ffts(cfo, period):
            if m >= n_max:
                break
            l = l0 + m
            k = np.arange(dx * (l % dy), self.noc, dx * dy)
            q = Y[self.bins[k]] / (amps["scattered"] * self.signs[k])
            acc += np.sum(q[:-1] * np.conj(q[1:]))
        return float(np.angle(acc) / (2 * np.pi) * g["nfft"] / (dx * dy))

    def class_power(self, y, cfo, period, l0):
        """mean |Y|^2 per pilot class and for data cells, on THIS timebase.

        Two jobs.  Before cancellation it MEASURES the per-class amplitude
        ratio from air, which REENCODER.md risk 3 demands as an independent
        witness on the boost table.  After cancellation it says which class
        actually went away -- and the `data` row is the built-in control:
        if data-cell power moves, we are not cancelling pilots, we are
        cancelling signal.
        """
        g = self.g
        acc = {k: 0.0 for k in ("scattered", "sbs", "cp_common",
                                "cp_additional", "edge", "data")}
        cnt = dict(acc)
        for m in range(self.n_sym):
            w = int(round(m * period))
            s = w + g["gi"]
            if s + g["nfft"] > y.size:
                break
            n = np.arange(s, s + g["nfft"])
            Y = np.fft.fft(y[s:s + g["nfft"]]
                           * np.exp(-2j * np.pi * cfo * n / FS))[self.bins]
            l = l0 + m
            cls, _ = PP.pilot_classes(g["nfft"], g["cred"], g["pattern"], l,
                                      l in g["sbs_l"])
            used = np.zeros(self.noc, bool)
            for name in acc:
                if name == "data":
                    continue
                k = np.array([kk for kk, v in cls.items() if v == name], int)
                if k.size:
                    acc[name] += float(np.sum(np.abs(Y[k]) ** 2))
                    cnt[name] += k.size
                    used[k] = True
            acc["data"] += float(np.sum(np.abs(Y[~used]) ** 2))
            cnt["data"] += int((~used).sum())
        return {k: (acc[k] / cnt[k] if cnt[k] else 0.0) for k in acc}


def refine_alignment(al, cfo0, period0, l0, amps, rounds=9,
                     step_cfo=4.0, step_per=0.02):
    """Coordinate descent on (CFO, symbol period) maximising ls_one_h.

    Two parameters, not one: the capture is not frequency-corrected AND the
    transmitter's sample clock is not ours.  Leaving the second out costs
    real cancellation -- measured on this capture, refining CFO alone got
    the LS-one-H residual to 9.5 dB while refining both reached 21.3 dB,
    because a 0.8 ppm clock walks the symbol timing ~0.8 samples across a
    150 ms dwell and a time-INVARIANT ECA filter cannot follow a drifting
    delay.
    """
    cfo, per = float(cfo0), float(period0)
    best = al.ls_one_h(cfo, per, l0, amps)
    trace = [(cfo, per, best)]
    sc, sp = step_cfo, step_per
    for _ in range(rounds):
        moved = True
        while moved:
            moved = False
            for dc, dp in ((sc, 0.0), (-sc, 0.0), (0.0, sp), (0.0, -sp)):
                v = al.ls_one_h(cfo + dc, per + dp, l0, amps)
                if v > best + 1e-4:
                    best, cfo, per = v, cfo + dc, per + dp
                    moved = True
        sc /= 2.0
        sp /= 2.0
        trace.append((cfo, per, best))
    return cfo, per, best, trace


# ---------------------------------------------------------------------------
# 4 -- synthesis on the capture's index space, and the projections
# ---------------------------------------------------------------------------

def synth_window(g, amps, signs, data0_abs, abs_start, n, cfo, period):
    """Pilot reference covering capture samples [abs_start, abs_start+n).

    Built by emitting whole symbols from the frame's own data-symbol 0 and
    slicing, never by starting a symbol mid-way: the symbol index l decides
    the lattice phase and which symbols are SBS, so "start here" has to be
    answered in symbols before it can be answered in samples.
    """
    m0 = int(math.floor((abs_start - data0_abs) / period)) - 1
    n_sym = int(math.ceil(n / period)) + 3
    start_abs = data0_abs + int(round(m0 * period))
    ref = PP.synth_pilot_stream(n_sym, m0, g["nfft"], g["gi"], g["cred"],
                                g["pattern"], amps, sbs_l=g["sbs_l"],
                                signs=signs, sym_period=period, cfo_hz=cfo,
                                fs=FS, n0=start_abs)
    o = abs_start - start_abs
    if o < 0 or o + n > ref.size:
        raise AssertionError("synth window slice out of range")
    return ref[o:o + n]


def _basis_block(refs, n_taps, r0, r1):
    cols = [sliding_window_view(
        np.concatenate([np.zeros(n_taps - 1, r.dtype), r]), n_taps)[r0:r1, ::-1]
        for r in refs]
    return np.concatenate(cols, axis=1).astype(np.complex128)


def fit_groups(surv, refs, n_taps, ridge=1e-9, rows=1 << 17):
    """Least-squares filter per reference GROUP -- the two-class projection.

    `cancel_direct` fits one filter to one reference and subtracts the
    whole thing, which is the right tool for full removal.  Flattening
    needs the classes separable, because the amount to subtract differs by
    class (beta = 1 - 1/A, and A_SP != A_CP).  So the basis here is the
    delayed copies of EACH group side by side and the solve returns one
    filter per group, which the caller can then scale independently.

    Blocked normal equations for the same reason pr_core uses them: the
    full basis for a 150 ms dwell at 2 groups x 16 taps is 265 MB in
    complex128, and it does not need to exist all at once.
    """
    n = min(surv.size, min(r.size for r in refs))
    refs = [np.asarray(r)[:n] for r in refs]
    y = np.asarray(surv[:n], dtype=np.complex128)
    n_col = n_taps * len(refs)
    G = np.zeros((n_col, n_col), dtype=np.complex128)
    b = np.zeros(n_col, dtype=np.complex128)
    for r0 in range(0, n, rows):
        r1 = min(r0 + rows, n)
        B = _basis_block(refs, n_taps, r0, r1)
        G += B.conj().T @ B
        b += B.conj().T @ y[r0:r1]
        del B
    G[np.diag_indices_from(G)] += ridge * (np.trace(G).real / n_col)
    try:
        w = np.linalg.solve(G, b)
    except np.linalg.LinAlgError:
        w, *_ = np.linalg.lstsq(G, b, rcond=None)
    return w, n


def subtract_groups(surv, refs, n_taps, w, scale, rows=1 << 17):
    """surv - sum_g scale[g] * (w_g convolved with ref_g)."""
    n = min(surv.size, min(r.size for r in refs))
    ws = w.copy()
    for gi_, s in enumerate(scale):
        ws[gi_ * n_taps:(gi_ + 1) * n_taps] *= s
    out = np.array(surv[:n], dtype=np.complex128)
    for r0 in range(0, n, rows):
        r1 = min(r0 + rows, n)
        B = _basis_block([np.asarray(r)[:n] for r in refs], n_taps, r0, r1)
        out[r0:r1] -= B @ ws
        del B
    return out


# ---------------------------------------------------------------------------
# 5 -- the measurement that matters: what happened to the comb
# ---------------------------------------------------------------------------

def comb_report(y, g, tag, save=None):
    """Run the v0 pipeline on `y` and read the comb teeth off the map.

    The pipeline is pr_banked's, unchanged, so BEFORE is literally the
    first-look run: self-referenced ECA, then a CAF of the raw channel
    against the cancelled one, then CFAR.  Only the input changes.
    """
    clean = PR.cancel_direct(y, y, ECA_TAPS)
    rd, rax, dax = PR.caf(y, clean, FS, MAX_RANGE_M, DOPPLER_HZ, N_DOPPLER)
    det = [d for d in PR.detect(rd, rax, dax, pfa=PFA)
           if abs(d["doppler_hz"]) >= MIN_DOPPLER_HZ]
    noise, _alpha = PR.cfar_threshold_map(rd, pfa=PFA)
    snr = 10.0 * np.log10(np.maximum(rd, 1e-300) / np.maximum(noise, 1e-300))

    r1 = PP.lattice_range_m(g["nfft"], g["pattern"])
    fd = PP.lattice_doppler_hz(g["nfft"], g["gi"], g["pattern"])
    teeth = []
    for n in range(1, 9):
        r = n * r1
        if r > rax[-1]:
            break
        i = int(np.argmin(np.abs(rax - r)))
        j = int(np.argmin(np.abs(dax - fd)))
        win = slice(max(0, i - 2), i + 3), slice(max(0, j - 2), j + 3)
        teeth.append(dict(n=n, range_m=round(r, 1),
                          snr_db=round(float(snr[win].max()), 2),
                          power_db=round(float(10 * math.log10(
                              max(rd[win].max(), 1e-300))), 2)))
    out = dict(tag=tag, n_detections=len(det),
               floor_db=round(float(10 * math.log10(
                   np.median(rd[200:, :]))), 2),
               teeth=teeth,
               best=(dict(range_km=round(det[0]["range_m"] / 1e3, 3),
                          doppler_hz=round(det[0]["doppler_hz"], 1),
                          snr_db=round(det[0]["snr_db"], 2)) if det else None))
    if save:
        np.save(save, rd)
    return out


def _fmt_comb(c):
    t = "  ".join(f"n{d['n']}:{d['snr_db']:5.1f}" for d in c["teeth"])
    b = (f"best {c['best']['snr_db']:5.1f} dB @ "
         f"{c['best']['range_km']:6.2f} km {c['best']['doppler_hz']:+7.1f} Hz"
         if c["best"] else "no detections")
    return f"{c['tag']:26s} det {c['n_detections']:3d}  {b}\n      comb {t}"


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("capture", nargs="?",
                    default=r"Z:\src\atsc3\data\fox25\rf25_fox_g4_0810_2254.cs16")
    ap.add_argument("--l1", default=str(HERE.parent / "lab"
                                        / "l1_radar_capture.json"))
    ap.add_argument("--t0", type=float, default=30.0,
                    help="seconds into the capture; the dwell is SNAPPED "
                         "forward to the next frame's first data symbol")
    ap.add_argument("--dwell-s", type=float, default=0.150)
    ap.add_argument("--frames", type=int, default=3,
                    help="one dwell per frame, each re-anchored (REENCODER "
                         "risk 1: the reference is rebuilt per frame)")
    ap.add_argument("--pilot-taps", type=int, default=PILOT_TAPS)
    ap.add_argument("--guard", type=int, default=PILOT_GUARD,
                    help="samples of NEGATIVE delay the pilot basis covers")
    ap.add_argument("--seed", type=int, default=7,
                    help="RNG seed for the scrambled-pilot negative control")
    ap.add_argument("--out", default=str(HERE.parent / "lab"
                                         / "pilot_cancel.json"))
    a = ap.parse_args()

    t_start = time.time()
    print("=" * 72)
    print("pilot_cancel_test -- can a SYNTHESIZED pilot lattice kill the comb?")
    print("=" * 72)

    # --- gates first.  A cancellation number from an ungated table is a
    # --- guess with a decimal point on it.
    npass, nfail = PP.selftest(verbose=False)
    print(f"[gate] pr_pilots selftest: {npass} pass, {nfail} fail")
    if nfail:
        PP.selftest(verbose=True)
        return 1

    g, l1 = geometry_from_l1(a.l1)
    amps = PP.pilot_amplitudes(g["pattern"], g["boost"])
    excess = PP.pilot_excess_amplitudes(amps)
    noc = P.NOC[(g["nfft"], g["cred"])]
    signs = PP.reference_signs(noc)
    signs_scr = PP.reference_signs(noc, "random", seed=a.seed)
    r1 = PP.lattice_range_m(g["nfft"], g["pattern"])
    fd = PP.lattice_doppler_hz(g["nfft"], g["gi"], g["pattern"])
    print(f"[l1  ] FFT {g['nfft']} / GI {g['gi']} / {g['pattern']} / "
          f"Cred {g['cred']} / boost {g['boost']} / {g['np_sym']} preamble + "
          f"{g['n_sym']} data symbols")
    print(f"[amp ] A_SP {amps['scattered']:.4f}  A_CP {amps['cp_common']:.4f} "
          f" flatten beta: sp {1 - 1 / amps['scattered']:.4f} "
          f"cp {1 - 1 / amps['cp_common']:.4f}")
    print(f"[comb] predicted tooth {r1:.1f} m / {fd:.2f} Hz "
          f"(FIRST_LOOK measured 7370 m / +199.2 Hz, bins 43.4 m / 3.9 Hz)")

    frames, spacing, err = find_frames(a.capture, a.t0, a.frames, g)
    print(f"[anch] {len(frames)} frames; bootstrap spacing {spacing:.1f} "
          f"samples vs L1 layout {g['frame_samples']} "
          f"(err {err:+.1f} = {err / g['T']:+.4f} symbols, "
          f"{err / g['frame_samples'] * 1e6:+.2f} ppm clock)")

    n_dwell = int(a.dwell_s * FS)
    if n_dwell > g["n_sym"] * g["T"]:
        raise SystemExit("dwell longer than the Subframe's data symbols")

    meta = dict(capture=os.path.abspath(a.capture),
                l1=os.path.abspath(a.l1), fs_hz=FS, fc_hz=FC,
                t0_requested_s=a.t0, dwell_s=a.dwell_s, frames=a.frames,
                geometry={k: v for k, v in g.items() if k != "sbs_l"},
                sbs_l=list(g["sbs_l"]),
                amplitudes={k: v for k, v in amps.items()},
                bootstrap_spacing_samples=spacing,
                frame_layout_err_symbols=err / g["T"],
                caf=dict(max_range_m=MAX_RANGE_M, doppler_hz=DOPPLER_HZ,
                         n_doppler=N_DOPPLER, eca_taps=ECA_TAPS, pfa=PFA),
                pilot_taps=a.pilot_taps, seed=a.seed,
                comb_predicted=dict(range_m=r1, doppler_hz=fd),
                dwells=[])

    for fi, fr in enumerate(frames):
        d0 = fr["data0_abs"]
        print("\n" + "-" * 72)
        print(f"frame {fi}: bootstrap @ {fr['bootstrap_abs']} "
              f"(ratio {fr['peak_ratio']:.0f}), data symbol 0 @ {d0} "
              f"= t {d0 / FS:.4f}s, coarse CFO {fr['cfo_hz']:+.2f} Hz")
        x = PR.load_cs16(a.capture, n_samples=n_dwell, offset=d0)
        al = Aligner(x, g, signs)

        # --- l mod Dy: the frame anchor PREDICTS 0 (the dwell starts at
        # --- data symbol 0).  Score all Dy hypotheses; the prediction has
        # --- to win, and by how much is the discrimination.
        _dx, dy = P.dxdy(g["pattern"])
        coh = [al.coherence(fr["cfo_hz"], float(g["T"]), l, amps)
               for l in range(dy)]
        won = int(np.argmax(coh))
        print("[align] pilot coherence by l mod Dy: "
              + "  ".join(f"l={i}:{c:.3f}" for i, c in enumerate(coh))
              + f"   -> predicted 0, won {won}"
              + ("  OK" if won == 0 else "  *** MISMATCH ***"))
        if won != 0:
            raise SystemExit("frame anchor and pilot lattice disagree; "
                             "refusing to report a cancellation number")

        cfo, per, lsh, _tr = refine_alignment(al, fr["cfo_hz"],
                                              float(g["T"]), 0, amps)
        lsh0 = al.ls_one_h(fr["cfo_hz"], float(g["T"]), 0, amps)
        print(f"[align] refined CFO {cfo:+.4f} Hz (from {fr['cfo_hz']:+.2f}), "
              f"symbol period {per:.6f} samples "
              f"({(per - g['T']) / g['T'] * 1e6:+.3f} ppm)")
        print(f"[align] LS-fit-one-H residual: {lsh0:.2f} dB at the coarse "
              f"anchor -> {lsh:.2f} dB refined "
              f"(REENCODER gate: pilots alone >> 20 dB)")

        tau0 = al.direct_delay(cfo, per, 0, amps)
        print(f"[align] direct-path delay {tau0:+.3f} samples "
              f"({tau0 * PR.C_LIGHT / FS:+.0f} m); basis covers "
              f"{-a.guard} .. {a.pilot_taps - 1 - a.guard} samples")
        if not (-a.guard + 1 <= tau0 <= a.pilot_taps - 1 - a.guard):
            print("        *** direct path outside the tap window -- "
                  "raise --guard or --pilot-taps ***")

        # --- per-class amplitude, MEASURED from air.  REENCODER risk 3
        # --- asks for this independently of the table; agreement is the
        # --- point, and disagreement would invalidate everything below.
        cp = al.class_power(x, cfo, per, 0)
        meas_sp = 10 * math.log10(cp["scattered"] / cp["data"])
        meas_cp = 10 * math.log10(cp["cp_common"] / cp["data"])
        print(f"[amp ] measured from air: SP/data {meas_sp:+.2f} dB "
              f"(table {20 * math.log10(amps['scattered']):+.2f}), "
              f"CP/data {meas_cp:+.2f} dB "
              f"(table {20 * math.log10(amps['cp_common']):+.2f})")

        # --- synthesize.  Two groups so the classes can be scaled apart.
        amps_sp = dict(amps)
        amps_sp["cp_common"] = 0.0
        amps_cp = {k: (amps["cp_common"] if k == "cp_common" else 0.0)
                   for k in amps if not k.startswith("_")}
        # The reference is generated GUARD samples AHEAD of the dwell, so
        # the tap window spans -GUARD..n_taps-1-GUARD and a direct path
        # that arrives before our anchor is still inside the subspace.
        ahead = d0 + a.guard
        results = {}
        for mode, sg in (("real", signs), ("scrambled", signs_scr)):
            ref_all = synth_window(g, amps, sg, d0, ahead, n_dwell, cfo, per)
            ref_sp = synth_window(g, amps_sp, sg, d0, ahead, n_dwell, cfo, per)
            ref_cp = synth_window(g, amps_cp, sg, d0, ahead, n_dwell, cfo, per)

            # (a) FULL removal, through pr_core.cancel_direct as-is.
            x_full = PR.cancel_direct(x, ref_all, a.pilot_taps)
            # (b) FLATTEN: same fitted channel, but subtract only the
            #     EXCESS over the data-cell level, per class.
            w, _n = fit_groups(x, [ref_sp, ref_cp], a.pilot_taps)
            b_sp = 1.0 - 1.0 / amps["scattered"]
            b_cp = 1.0 - 1.0 / amps["cp_common"]
            x_flat = subtract_groups(x, [ref_sp, ref_cp], a.pilot_taps, w,
                                     (b_sp, b_cp))
            results[mode] = dict(ref=ref_all, full=x_full, flat=x_flat,
                                 canc_full=PR.cancellation_db(x, x_full),
                                 canc_flat=PR.cancellation_db(x, x_flat))
            cf = al.class_power(x_full, cfo, per, 0)
            results[mode]["suppression_db"] = {
                k: (round(10 * math.log10(cp[k] / cf[k]), 2)
                    if cp[k] > 0 and cf[k] > 0 else None) for k in cp}
            print(f"[canc] {mode:9s}: total {results[mode]['canc_full']:.3f} dB "
                  f"(full removal), {results[mode]['canc_flat']:.3f} dB "
                  f"(flatten); pilot-carrier suppression "
                  + ", ".join(f"{k} {v:+.1f}"
                              for k, v in results[mode]["suppression_db"].items()
                              if v is not None and k in ("scattered",
                                                         "cp_common", "data")))

        # --- the measurement that matters
        maps = {}
        maps["before"] = comb_report(x, g, "BEFORE (raw capture)")
        maps["after_full"] = comb_report(results["real"]["full"], g,
                                         "AFTER full pilot removal")
        maps["after_flat"] = comb_report(results["real"]["flat"], g,
                                         "AFTER pilot flattening")
        maps["ctl_full"] = comb_report(results["scrambled"]["full"], g,
                                       "CONTROL scrambled, full")
        maps["ctl_flat"] = comb_report(results["scrambled"]["flat"], g,
                                       "CONTROL scrambled, flatten")
        if fi == 0:
            # positive control: does the SYNTHETIC lattice, alone, make the
            # artifact?  If it does not, nothing downstream is meaningful.
            maps["synthetic_only"] = comb_report(results["real"]["ref"], g,
                                                 "synthetic pilots ALONE")
        print()
        for k in ("before", "after_full", "after_flat", "ctl_full",
                  "ctl_flat", "synthetic_only"):
            if k in maps:
                print(_fmt_comb(maps[k]))

        b1 = maps["before"]["teeth"][0]["snr_db"]
        print(f"\n[VERDICT frame {fi}] tooth n=1 ({r1 / 1e3:.2f} km, "
              f"{fd:+.1f} Hz): {b1:.1f} dB -> "
              f"{maps['after_full']['teeth'][0]['snr_db']:.1f} dB full "
              f"({maps['after_full']['teeth'][0]['snr_db'] - b1:+.1f}), "
              f"{maps['after_flat']['teeth'][0]['snr_db']:.1f} dB flatten "
              f"({maps['after_flat']['teeth'][0]['snr_db'] - b1:+.1f});  "
              f"CONTROL {maps['ctl_flat']['teeth'][0]['snr_db']:.1f} dB "
              f"({maps['ctl_flat']['teeth'][0]['snr_db'] - b1:+.1f})")

        meta["dwells"].append(dict(
            frame=fi, data0_abs=d0, t_s=round(d0 / FS, 6),
            cfo_coarse_hz=fr["cfo_hz"], cfo_refined_hz=cfo,
            symbol_period_samples=per, direct_delay_samples=round(tau0, 4),
            clock_ppm=(per - g["T"]) / g["T"] * 1e6,
            coherence_by_l=[round(c, 4) for c in coh],
            ls_one_h_db=dict(coarse=round(lsh0, 2), refined=round(lsh, 2)),
            measured_amplitude_db=dict(sp_over_data=round(meas_sp, 2),
                                       cp_over_data=round(meas_cp, 2),
                                       table_sp=round(20 * math.log10(
                                           amps["scattered"]), 2),
                                       table_cp=round(20 * math.log10(
                                           amps["cp_common"]), 2)),
            cancellation_db=dict(
                real_full=round(results["real"]["canc_full"], 3),
                real_flat=round(results["real"]["canc_flat"], 3),
                scrambled_full=round(results["scrambled"]["canc_full"], 3),
                scrambled_flat=round(results["scrambled"]["canc_flat"], 3)),
            pilot_suppression_db=dict(
                real=results["real"]["suppression_db"],
                scrambled=results["scrambled"]["suppression_db"]),
            maps={k: v for k, v in maps.items()}))
        del results, x

    meta["wall_s"] = round(time.time() - t_start, 1)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(meta, f, indent=1)
    print(f"\n[out ] {a.out}   ({meta['wall_s']:.0f}s)")

    # --- summary across frames
    d = meta["dwells"]
    def mean(sel):
        return float(np.mean([sel(z) for z in d]))
    print("\n" + "=" * 72)
    print(f"SUMMARY over {len(d)} frames")
    print(f"  LS-fit-one-H residual (refined) {mean(lambda z: z['ls_one_h_db']['refined']):6.2f} dB")
    print(f"  total power cancelled, real     {mean(lambda z: z['cancellation_db']['real_flat']):6.3f} dB "
          f"(flatten) / {mean(lambda z: z['cancellation_db']['real_full']):.3f} dB (full)")
    print(f"  total power cancelled, CONTROL  {mean(lambda z: z['cancellation_db']['scrambled_flat']):6.3f} dB "
          f"(flatten) / {mean(lambda z: z['cancellation_db']['scrambled_full']):.3f} dB (full)")
    print(f"  scattered-pilot suppression     {mean(lambda z: z['pilot_suppression_db']['real']['scattered']):6.2f} dB "
          f"real / {mean(lambda z: z['pilot_suppression_db']['scrambled']['scattered']):.2f} dB CONTROL")
    for tag, key in (("comb tooth n=1, before ", "before"),
                     ("  after full removal   ", "after_full"),
                     ("  after flattening     ", "after_flat"),
                     ("  CONTROL full         ", "ctl_full"),
                     ("  CONTROL flatten      ", "ctl_flat")):
        print(f"  {tag} {mean(lambda z: z['maps'][key]['teeth'][0]['snr_db']):6.2f} dB SNR   "
              f"(detections {mean(lambda z: z['maps'][key]['n_detections']):.1f})")
    print("=" * 72)
    print("CANDIDATES ONLY -- nothing here claims an aircraft.  What it "
          "claims is that the comb is the transmitter's pilot lattice and "
          "that a lattice synthesized from A/322 alone can level it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
