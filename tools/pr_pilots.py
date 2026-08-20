#!/usr/bin/env python3
"""pr_pilots.py -- synthesize the ATSC 3.0 PILOT LATTICE as a radar reference.

The smallest piece of the decode-remodulate build (docs/REENCODER.md), and
the piece that carries the whole premise.  FIRST_LOOK.md identified the
7.4 km / 198.5 Hz comb that dominates every range-Doppler map we have made
as the transmitter's own scattered-pilot lattice: within-symbol pilot
spacing Dx*Dy = 48 carriers -> 8192/48 = 170.67 samples = 7.402 km of
apparent range, repeating every Dy = 4 symbols -> 794.118/4 = 198.53 Hz of
apparent Doppler.  Nothing in the sky put it there.

WHAT THIS FILE EXISTS TO PREVENT
--------------------------------
A passive radar whose reference IS its surveillance channel cannot tell a
target from the waveform's own ambiguity function -- that is tell #1 of the
first-look run (175 dB of "cancellation" that was annihilation).  The fix is
a reference we build rather than receive.  Pilots are the cheapest honest
slice of that reference:

  * they are SINGLE-LAYER.  The RF25 signal is LDM (5.0 dB injection,
    256QAM-NUC enhanced layer we have never decoded), which caps data-cell
    cancellation at ~5 dB.  Pilots carry no LDM, so they can cancel fully.
  * they are KNOWN WITHOUT DECODING ANYTHING.  Positions come from the
    A/322 tables, values from the published PRBS, amplitudes from Table
    9.17.  No FEC, no constellation decisions, no re-encode.
  * they are EXACTLY the structure the CAF is currently seeing.

So if synthesizing pilots alone kills the comb, the re-encoder premise is
validated for the price of ~300 lines instead of ~1000.  And if it does not,
we learn that before spending the three days.

THE ONE THING THAT CAN QUIETLY RUIN THIS (REENCODER.md risk 3)
--------------------------------------------------------------
A wrong A_SP : A_CP : A_edge ratio injects a periodic error on precisely the
lattice we are trying to cancel -- a residual that impersonates the
artifact and cancels "almost" as well as the truth.  Two defences, both
here: `gate_table_9_17()` checks the computed amplitudes against the
amplitude column A/322 prints for itself, and `pilot_cancel_test.py` runs
the whole experiment again with SCRAMBLED pilot signs.  A wrong pilot field
that cancels as well as the right one means we were fitting noise.

UNITS.  Carrier indices are RELATIVE (0 .. NoC-1) unless named `abs`.
Amplitudes are linear voltage ratios; the spec's boost tables are in dB of
POWER, so A = 10^(dB/20).  Delays are in samples at fs = 6.912 MS/s
(1 sample = 43.4 m of bistatic range).  `l` is the index of a symbol within
its Subframe's data symbols and only matters modulo Dy.

NOTHING HERE OPENS THE RADIO and nothing here writes to Z:/src/atsc3 -- the
A/322 tables are imported read-only from that tree, which is their home and
where they are gated against both printed editions of the standard.
"""
from __future__ import annotations

import math
import os
import re
import sys

import numpy as np

# The A/322 tables live in the atsc3 tree and stay there.  spec_pilots does
# `import m3_spec as S` internally, so the lab directory itself has to be
# importable, not just its parent.
ATSC3_ROOT = os.environ.get("ATSC3_ROOT", r"Z:\src\atsc3")
ATSC3_LAB = os.path.join(ATSC3_ROOT, "lab")
for _p in (ATSC3_ROOT, ATSC3_LAB):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import m3_spec as S                                            # noqa: E402
import spec_pilots as P                                        # noqa: E402
import spec_l1syntax as LS                                     # noqa: E402
from m3_preamble import fft_bins                               # noqa: E402


# ---------------------------------------------------------------------------
# frame geometry constants -- units in the comment, always
# ---------------------------------------------------------------------------

FS_ATSC3 = 6.912e6            # samples/s, post-bootstrap rate for bsr_coeff 2
BOOTSTRAP_SAMPLES_POST = 13824  # samples, 4 x 3072 @ 6.144 MS/s scaled to 6.912

# A/322 8.1.4.2 Table 8.5: common continual pilot boost, all FFT sizes.  The
# dB value is the exact one; 2.67 is the spec's own rounded linear print.
A_CP_DB = S.A_CP_DB           # 8.52 dB of power -> 2.6669 in voltage

# Where the printed Table 9.17 lives, for the amplitude gate.  Absent is a
# reported SKIP, never a silent pass.
A322_TEXT = os.path.join(ATSC3_LAB, "spec", "A322_2026.txt")


def symbol_rate_hz(nfft: int, gi: int, fs: float = FS_ATSC3) -> float:
    """OFDM symbols per second.  8192+512 @ 6.912 MS/s -> 794.118 Hz."""
    return fs / float(nfft + gi)


def lattice_range_m(nfft: int, pattern: str, fs: float = FS_ATSC3) -> float:
    """Delay, in metres of bistatic range, of the pilot comb's first tooth.

    A symbol whose only non-zero carriers sit every Dx*Dy bins is periodic
    in time with period nfft/(Dx*Dy) samples, so its autocorrelation combs
    at multiples of that delay.  For 8K / SP12_4: 8192/48 = 170.67 samples
    = 7402 m.  FIRST_LOOK.md measured 7.37 km in a 43.4 m bin.
    """
    dx, dy = P.dxdy(pattern)
    from pr_core import C_LIGHT
    return (nfft / float(dx * dy)) * C_LIGHT / fs


def lattice_doppler_hz(nfft: int, gi: int, pattern: str,
                       fs: float = FS_ATSC3) -> float:
    """Apparent Doppler of the comb: the scattered lattice repeats every Dy
    symbols, so it beats at symbol_rate/Dy.  794.118/4 = 198.53 Hz."""
    _dx, dy = P.dxdy(pattern)
    return symbol_rate_hz(nfft, gi, fs) / float(dy)


# ---------------------------------------------------------------------------
# 1 -- pilot amplitudes.  The table the whole tree was missing.
# ---------------------------------------------------------------------------

def pilot_amplitudes(pattern: str, boost_idx: int) -> dict:
    """Linear voltage amplitude per pilot class for one Subframe config.

    Returns a dict with keys `scattered`, `cp_common`, `cp_additional`,
    `edge`, `sbs`, plus `_asp_db` for provenance.

    Every demodulator in the atsc3 tree passes amp = 1.0 because it DIVIDES
    the received pilot by the reference to get H -- the boost cancels and
    never has to be known.  A synthesizer has no such luxury: it has to emit
    the pilots at the right level RELATIVE to each other or it paints a
    periodic error on the lattice (REENCODER.md risk 3).  Hence this
    function, and hence `gate_table_9_17()`.

    Sources, one per class:
      scattered      A/322 Table 9.16 (dB) -> SCATTERED_PILOT_BOOST_DB
      cp_common      A/322 Table 8.5, 8.52 dB for every FFT size
      cp_additional  A/322 8.1.4.2 last line: "For additional continual
                     pilots, the modulation amplitude of the scattered
                     pilot pattern (ASP) shall be used."
      edge           A/322 8.1.5: "exactly the same as for the scattered
                     pilots"
      sbs            A/322 8.1.7.2: "boosted by the same factor as the
                     scattered pilots, ASP"

    So the ONLY independent ratio in the whole map is A_CP / A_SP.  For
    SP12_4 boost 1 that is 2.6669 / 1.6788 = 1.5886 (+4.02 dB of power),
    and it is the one number a wrong reading would get wrong.
    """
    if pattern not in LS.SCATTERED_PILOT_BOOST_DB:
        raise KeyError(f"unknown scattered pilot pattern {pattern!r}")
    tbl = LS.SCATTERED_PILOT_BOOST_DB[pattern]
    boost_idx = int(boost_idx)
    if boost_idx not in tbl:
        raise KeyError(f"L1D_scattered_pilot_boost={boost_idx} is RFU for "
                       f"{pattern} (A/322 Table 9.16 defines 0..4)")
    a_sp_db = float(tbl[boost_idx])
    a_sp = 10.0 ** (a_sp_db / 20.0)          # dB of POWER -> voltage ratio
    a_cp = 10.0 ** (A_CP_DB / 20.0)
    return {
        "scattered": a_sp,
        "cp_common": a_cp,
        "cp_additional": a_sp,
        "edge": a_sp,
        "sbs": a_sp,
        "_asp_db": a_sp_db,
        "_acp_db": A_CP_DB,
    }


def pilot_excess_amplitudes(amps: dict) -> dict:
    """The amplitudes to SUBTRACT to flatten the spectrum: A - 1 per class.

    WHY THIS EXISTS -- measured 2026-08-20, and it overturns the obvious
    plan.  The CAF's comb tooth at delay tau = nfft/(Dx*Dy) is

        R(tau) = sum_k |X_k|^2 exp(j2*pi*k*Dx*Dy/nfft)

    -- the DFT of the POWER SPECTRUM at 1/48, not a correlation of pilot
    amplitudes.  What combs is therefore the fact that pilot carriers hold
    DIFFERENT power from data carriers, and the sign of the tooth is the
    sign of that difference.  Cancelling the pilots outright does not
    flatten the spectrum: it replaces a bump of +4.5 dB with a hole of
    -15 dB, and a hole combs about as hard as a bump.  Measured on the
    radar capture: 19.5 dB of pilot-carrier suppression bought only 5.9 dB
    off the comb, because the residual tooth came back in ANTI-phase.

    Subtracting A - 1 instead leaves each pilot carrier at amplitude 1 --
    the data-cell level -- so |X_k|^2 is flat across the lattice and the
    tooth nulls.  The scale factor per class is beta = 1 - 1/A, and it is
    PREDICTED by the amplitude table with nothing fitted: for SP12_4 boost
    1 that is 0.4043 (scattered) and 0.6250 (common CP).  Measured effect
    at those exact values, mean of three consecutive frames: the tooth at
    7.40 km / +198.5 Hz falls from 28.2 dB to 4.8 dB of CFAR SNR, against
    24.6 dB for full removal -- and 0.8x and 1.2x of those betas are both
    visibly worse, so the table's value sits on the null it did not fit.
    That is the strongest evidence the amplitude map is right that this
    project has.

    Data cells are taken to have unit RMS, which is what A/322's
    constellations are normalised to and what the air agrees with here
    (measured SP/data +4.42 dB against A_SP^2 = +4.50 dB).
    """
    out = {k: (v - 1.0 if not k.startswith("_") else v)
           for k, v in amps.items()}
    for k, v in out.items():
        if not k.startswith("_") and v < 0:
            raise ValueError(f"class {k} has amplitude below the data-cell "
                             f"level ({amps[k]:.3f}); flattening it would "
                             f"mean ADDING energy, which this path does not do")
    return out


def _parse_table_9_17(path=A322_TEXT):
    """{pattern: [printed amplitudes...]} from the pdftotext dump of A/322.

    Table 9.17 is the amplitude restatement of Table 9.16's dB values and is
    therefore a FREE independent witness on our arithmetic.  Its layout does
    not survive pdftotext intact: the SPx_2 rows print only two columns and
    the leftover SPx_4 columns spill into unlabelled continuation lines at
    the bottom of the page.  So this parser claims only what survives
    unambiguously -- the label and the first two amplitude columns
    (boost 0 and boost 1), which is the column our capture actually uses.
    Claiming more from a mangled extraction is how a gate turns into a
    rubber stamp.
    """
    if not os.path.exists(path):
        return None
    row = re.compile(r"^(SP\d+_\d)\s*/\s*MP\d+_\d\s+((?:[\d.]+\s+){1,}).*$")
    out = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = row.match(line.strip())
            if not m:
                continue
            vals = [float(v) for v in m.group(2).split()]
            # The dB table (9.16) prints in the same shape; keep only the
            # rows whose first column is the amplitude 1.000, not 0.00 dB.
            if not vals or abs(vals[0] - 1.0) > 1e-9:
                continue
            out.setdefault(m.group(1), vals)
    return out or None


def gate_table_9_17(path=A322_TEXT, tol=6e-4):
    """Check `pilot_amplitudes` against the amplitude column A/322 prints.

    Returns (n_pass, n_fail, n_skip, notes).  A missing spec dump is a SKIP
    with a note, never a pass: "the file was not there" and "the numbers
    agreed" must not produce the same result.

    tol defaults to 6e-4 because the table prints three decimals, so a
    correct value can differ from ours by up to half a unit in the last
    place (5e-4) plus float slop.
    """
    printed = _parse_table_9_17(path)
    if printed is None:
        return 0, 0, 1, [f"SKIP: no printed Table 9.17 at {path}"]
    npass = nfail = 0
    notes = []
    for pat, vals in sorted(printed.items()):
        for b in (0, 1):                     # see _parse_table_9_17
            if b >= len(vals):
                continue
            try:
                got = pilot_amplitudes(pat, b)["scattered"]
            except KeyError:
                continue
            if abs(got - vals[b]) <= tol:
                npass += 1
            else:
                nfail += 1
                notes.append(f"FAIL {pat} boost {b}: computed {got:.4f} "
                             f"vs printed {vals[b]:.3f}")
    return npass, nfail, 0, notes


# ---------------------------------------------------------------------------
# 2 -- the pilot field of one symbol
# ---------------------------------------------------------------------------

def reference_signs(noc: int, mode: str = "spec", shift: int = 0, seed: int = 0):
    """+-1 per relative carrier index k -- the sign half of Re{c} = 2A(1/2-r_k).

    `mode="spec"` is A/322 8.1.2's PRBS (G(x) = 1 + X^9 + X^10 + X^12 + X^13,
    reloaded every symbol), gated in m3_spec against the standard's own
    printed 24-bit test vector.

    The other two modes exist ONLY as negative controls and are the reason
    this function is separate from `pilot_field`:
      "shift"   the right PRBS read from the wrong phase -- the failure a
                one-carrier origin error would actually produce;
      "random"  independent random signs at the right POSITIONS with the
                right AMPLITUDES.

    A control that keeps the lattice geometry and destroys only the values
    is the sharp one: if cancellation survives it, we are fitting the
    lattice's ENVELOPE (or noise), not the transmitted pilots.
    """
    if mode == "spec":
        return np.array(S.pilot_values(noc, 1.0), dtype=np.float64)
    if mode == "shift":
        base = np.array(S.pilot_values(noc + abs(shift) + 1, 1.0),
                        dtype=np.float64)
        return np.roll(base, shift)[:noc]
    if mode == "random":
        rng = np.random.default_rng(seed)
        return rng.choice(np.array([-1.0, 1.0]), size=noc)
    raise ValueError(f"unknown sign mode {mode!r}")


def pilot_classes(nfft: int, cred: int, pattern: str, l: int, sbs: bool):
    """Relative carrier index -> class name, for every pilot-bearing carrier.

    Resolution order matters and follows A/322:
      * common CPs are chosen so they NEVER land on scattered pilot
        locations (8.1.4.1), so a collision there would mean our carrier
        origin is wrong -- it is counted and returned, not swallowed;
      * additional CPs MAY land on scattered locations, and A/322's MIMO
        note settles the reading: such a cell "shall be treated as
        scattered".  Both classes carry A_SP anyway, so the choice is
        invisible in SISO -- which is exactly why it must be written down
        rather than left to whichever dict updated last.
    """
    noc = P.NOC[(nfft, cred)]
    dx, dy = P.dxdy(pattern)
    lo, _hi = S.carrier_abs_range(nfft, cred)

    sp = (set(range(0, noc, dx)) if sbs
          else set(range(dx * (l % dy), noc, dx * dy)))
    cp_c = {c - lo for c in S.common_cps(nfft, cred)}
    cp_a = {k for k in P.additional_cps(nfft, pattern, cred) if 0 <= k < noc}
    edge = {0, noc - 1}

    collisions = len(sp & cp_c)
    cls = {}
    for k in cp_c:
        cls[k] = "cp_common"
    for k in cp_a:
        cls[k] = "cp_additional"
    for k in edge:
        cls[k] = "edge"
    for k in sp:                    # scattered / SBS wins every tie
        cls[k] = "sbs" if sbs else "scattered"
    return cls, collisions


def pilot_field(nfft: int, cred: int, pattern: str, l: int, sbs: bool,
                amps: dict, signs=None) -> np.ndarray:
    """complex128[NoC]: the transmitted value of every relative carrier.

    Pilots carry their class amplitude times their PRBS sign; every data
    carrier is ZERO.  That is the whole point of a pilot-only reference --
    the data cells are the part we cannot describe (LDM enhanced layer), so
    we describe nothing there rather than describe it wrongly.

    Cross-checked against `spec_pilots.pilot_carriers`, which is the gated
    path the demodulator uses: the set of non-zero carriers here must equal
    it exactly, or the two halves of the codebase disagree about what a
    pilot is.
    """
    noc = P.NOC[(nfft, cred)]
    cls, _coll = pilot_classes(nfft, cred, pattern, l, sbs)
    want = P.pilot_carriers(nfft, cred, pattern, l, sbs)
    if set(cls) != set(want):
        raise AssertionError(
            f"pilot set disagrees with spec_pilots.pilot_carriers: "
            f"{len(cls)} vs {len(want)} carriers "
            f"(nfft={nfft} cred={cred} {pattern} l={l} sbs={sbs})")
    if signs is None:
        signs = reference_signs(noc)
    signs = np.asarray(signs, dtype=np.float64)
    if signs.size != noc:
        raise ValueError(f"signs has {signs.size} entries, need NoC={noc}")

    field = np.zeros(noc, dtype=np.complex128)
    for k, name in cls.items():
        field[k] = amps[name] * signs[k]      # Im{c} = 0 for every pilot class
    return field


# ---------------------------------------------------------------------------
# 3 -- one OFDM symbol of pilots, in the time domain
# ---------------------------------------------------------------------------

def symbol_bins(nfft: int, cred: int):
    """Indices into an UNSHIFTED length-nfft spectrum for relative carriers.

    `m3_preamble.fft_bins` returns positions in an fftshift()ed spectrum
    (that is the convention every demod path in the atsc3 tree uses).  A
    synthesizer wants the plain ifft order, so the shift is folded in here
    once, as a reindex, instead of an ifftshift of a mostly-zero array.
    """
    noc = P.NOC[(nfft, cred)]
    lo, _hi = S.carrier_abs_range(nfft, cred)
    b = fft_bins(nfft, lo + np.arange(noc))
    if b.min() < 0 or b.max() >= nfft:
        raise ValueError("carrier range falls outside the FFT -- wrong nfft/cred")
    return (b - nfft // 2) % nfft          # fftshifted -> natural order


def synth_pilot_symbol(nfft: int, gi: int, cred: int, pattern: str, l: int,
                       sbs: bool, amps: dict, signs=None,
                       frac_delay: float = 0.0) -> np.ndarray:
    """complex128[nfft+gi]: cyclic prefix + body of one pilot-only symbol.

    `frac_delay` (samples, may be negative) delays the symbol by a
    FRACTIONAL number of samples.  It is applied as a phase ramp across
    carriers rather than by resampling, which is exact here and not an
    approximation: over the CP+body span the transmitted symbol IS periodic
    with period nfft, so a sub-CP delay is a circular shift.

    It is not a nicety.  The transmitter's clock and ours differ by about
    half a part per million (measured: the frame period lands 1676287.1
    samples apart against a nominal 1676288), which walks the symbol timing
    by ~0.6 samples across a 150 ms dwell.  An ECA basis fits ONE filter for
    the whole dwell and cannot track a drifting delay, so an un-drifted
    reference would put a floor under the cancellation and we would blame
    the pilot model.
    """
    field = pilot_field(nfft, cred, pattern, l, sbs, amps, signs)
    X = np.zeros(nfft, dtype=np.complex128)
    X[symbol_bins(nfft, cred)] = field
    if frac_delay:
        # exp(-j2*pi*f*d/N) with f the SIGNED frequency index: a delay of d
        # samples.  fftfreq(N)*N gives [0,1,..,N/2-1,-N/2,..,-1] exactly.
        f = np.fft.fftfreq(nfft) * nfft
        X *= np.exp(-2j * np.pi * f * float(frac_delay) / nfft)
    body = np.fft.ifft(X) * math.sqrt(nfft)   # scale is arbitrary (LS absorbs it)
    return np.concatenate([body[nfft - gi:], body])


# ---------------------------------------------------------------------------
# 4 -- a stream of them, on the CAPTURE's clock
# ---------------------------------------------------------------------------

def synth_pilot_stream(n_symbols: int, l0: int, nfft: int, gi: int, cred: int,
                       pattern: str, amps: dict, sbs_l=(), signs=None,
                       sym_period: float | None = None, cfo_hz: float = 0.0,
                       fs: float = FS_ATSC3, n0: int = 0,
                       start_delay: float = 0.0) -> np.ndarray:
    """The pilot-only reference, sample-aligned to a raw capture.

    n_symbols   how many consecutive OFDM symbols to emit
    l0          data-symbol index of the FIRST emitted symbol (only l mod Dy
                changes the lattice, but the absolute l is what decides
                which symbols are Subframe Boundary Symbols)
    sbs_l       the set of absolute l values that are SBS (first/last data
                symbol of the Subframe)
    sym_period  samples per symbol on the TRANSMITTER's clock as measured on
                OUR clock -- a float, nominally nfft+gi.  See
                `synth_pilot_symbol` on why this is not an integer.
    cfo_hz      carrier frequency offset to IMPOSE on the reference.  The
                capture is not frequency-corrected and we refuse to correct
                it: the surveillance channel must stay the raw evidence, so
                the model moves to meet it, not the other way round.
    n0          absolute sample index of output[0] within the capture, so
                the CFO phase is continuous with whatever came before.
    start_delay fractional-sample delay of the whole stream.

    REENCODER.md risk 1 is raw-timebase bookkeeping, and this signature is
    the answer to it: everything the decode chain would have done to the
    capture (resample, de-rotate, notch, re-anchor) is instead done to the
    reference, in the opposite direction.
    """
    T = nfft + gi
    period = float(T if sym_period is None else sym_period)
    sbs_set = set(int(v) for v in sbs_l)
    n_out = int(math.ceil((n_symbols - 1) * period)) + T
    out = np.zeros(n_out, dtype=np.complex128)

    starts = [int(round(m * period)) for m in range(n_symbols + 1)]
    for m in range(n_symbols):
        l = l0 + m
        frac = m * period - starts[m] + start_delay
        sym = synth_pilot_symbol(nfft, gi, cred, pattern, l,
                                 l in sbs_set, amps, signs, frac_delay=frac)
        s0 = starts[m]
        # Write only up to the next symbol's start.  With a non-integer
        # period consecutive starts differ by T or T-1, and letting symbols
        # overwrite each other's tails would put a one-sample discontinuity
        # in the reference every couple of symbols -- broadband, and exactly
        # the kind of self-inflicted structure a CAF loves.
        s1 = min(s0 + T, starts[m + 1] if m + 1 < n_symbols else s0 + T, n_out)
        out[s0:s1] = sym[:s1 - s0]

    if cfo_hz:
        n = np.arange(n_out, dtype=np.float64) + float(n0)
        out *= np.exp(2j * np.pi * float(cfo_hz) * n / float(fs))
    return out


# ---------------------------------------------------------------------------
# 5 -- self-check
# ---------------------------------------------------------------------------

def selftest(verbose=True):
    """Gates that must hold before any cancellation number means anything."""
    npass = nfail = 0
    lines = []

    def chk(name, cond, detail=""):
        nonlocal npass, nfail
        if cond:
            npass += 1
        else:
            nfail += 1
            lines.append(f"  FAIL {name}  {detail}")
        if verbose and cond:
            lines.append(f"  ok   {name}  {detail}")

    # G1 -- the amplitude table against A/322's own printed amplitudes.
    p, f, s, notes = gate_table_9_17()
    chk("G1 Table 9.17 amplitudes", f == 0 and p > 0,
        f"{p} checked, {f} failed, {s} skipped" + ("; " + "; ".join(notes)
                                                   if notes else ""))

    # G2 -- the one ratio nothing else pins.
    a = pilot_amplitudes("SP12_4", 1)
    chk("G2 SP12_4 boost1 A_SP", abs(a["scattered"] - 1.6788) < 5e-4,
        f"{a['scattered']:.4f} (Table 9.17 prints 1.679)")
    chk("G2 A_CP", abs(a["cp_common"] - 2.6669) < 5e-4,
        f"{a['cp_common']:.4f} (Table 8.5 prints 2.67)")

    # G3 -- common CPs must never collide with scattered pilots (A/322
    # 8.1.4.1 says the CP32 set was chosen so they do not).  A collision
    # would mean our carrier origin is off, and it would be invisible in
    # any demodulator because they all divide the pilot out.
    coll = max(pilot_classes(8192, 3, "SP12_4", l, False)[1] for l in range(4))
    chk("G3 no common-CP / SP collision", coll == 0, f"{coll} collisions")

    # G4 -- pilot count per symbol must equal NoC - Table 7.3/7.4 data cells,
    # the closed identity spec_pilots.verify() is built on.
    noc = P.NOC[(8192, 3)]
    for l in range(4):
        n_pilot = len(pilot_classes(8192, 3, "SP12_4", l, False)[0])
        chk(f"G4 pilots(l={l}) closes Table 7.4",
            noc - n_pilot == P.avail_data_cells(8192, 3, "SP12_4"),
            f"{n_pilot} pilots, {noc - n_pilot} data vs "
            f"{P.avail_data_cells(8192, 3, 'SP12_4')}")

    # G5 -- the artifact's two observables must fall out of the geometry.
    r = lattice_range_m(8192, "SP12_4")
    d = lattice_doppler_hz(8192, 512, "SP12_4")
    chk("G5 comb range == 7.40 km", abs(r - 7402.0) < 2.0, f"{r:.1f} m")
    chk("G5 comb doppler == 198.53 Hz", abs(d - 198.53) < 0.02, f"{d:.2f} Hz")

    # G6 -- the comb, from first principles.  Carriers spaced Dx*Dy = 48
    # apart make a body that repeats every nfft/48 = 170.667 samples, so the
    # first INTEGER repeat is 512 samples = 3 x 170.667 -- and 512 samples
    # is 22.2 km, the third tooth of the comb FIRST_LOOK.md measured.
    amps = pilot_amplitudes("SP12_4", 1)
    noc = P.NOC[(8192, 3)]

    def selfcorr(v, lag):
        a1, a2 = v[:-lag], v[lag:]
        return abs(np.vdot(a1, a2)) / math.sqrt(
            float(np.vdot(a1, a1).real * np.vdot(a2, a2).real))

    # G6a: the scattered sub-lattice alone, which is the periodic part.
    cls, _ = pilot_classes(8192, 3, "SP12_4", 0, False)
    sp_only = np.zeros(noc, dtype=np.complex128)
    sgn = reference_signs(noc)
    for k, name in cls.items():
        if name == "scattered":
            sp_only[k] = amps["scattered"] * sgn[k]
    Xsp = np.zeros(8192, dtype=np.complex128)
    Xsp[symbol_bins(8192, 3)] = sp_only
    body_sp = np.fft.ifft(Xsp)
    chk("G6a scattered-only body repeats @ 512",
        selfcorr(body_sp, 512) > 0.99, f"{selfcorr(body_sp, 512):.3f}")
    chk("G6a control lag 256 does not",
        selfcorr(body_sp, 256) < 0.2, f"{selfcorr(body_sp, 256):.3f}")

    # G6b: in the FULL pilot symbol the continual pilots -- 46 carriers at
    # A_CP = 2.667, i.e. 327 of the 716 units of pilot power -- sit off the
    # 48-carrier lattice and break strict periodicity.  So the full symbol
    # correlates at ~the scattered power fraction, not at 1.  Stating that
    # number here is the point: it is also the ceiling on how much of the
    # comb a scattered-only reference could ever cancel.
    body = synth_pilot_symbol(8192, 512, 3, "SP12_4", 0, False, amps)[512:]
    c512, c256 = selfcorr(body, 512), selfcorr(body, 256)
    chk("G6b full body: comb lag beats control", c512 > 3.0 * c256,
        f"lag 512 {c512:.3f} vs lag 256 {c256:.3f}")

    # G7 -- fractional delay is a pure delay: delaying by an integer must
    # reproduce a circular shift of the undelayed body, to numerical noise.
    b0 = synth_pilot_symbol(8192, 512, 3, "SP12_4", 1, False, amps)[512:]
    b7 = synth_pilot_symbol(8192, 512, 3, "SP12_4", 1, False, amps,
                            frac_delay=7.0)[512:]
    err = np.max(np.abs(b7 - np.roll(b0, 7))) / np.max(np.abs(b0))
    chk("G7 frac_delay(7) == roll(7)", err < 1e-9, f"max rel err {err:.2e}")

    # G8 -- the negative-control signs must keep the LATTICE and destroy the
    # VALUES.  Same carriers, same amplitudes, uncorrelated field.
    f_true = pilot_field(8192, 3, "SP12_4", 0, False, amps)
    f_scr = pilot_field(8192, 3, "SP12_4", 0, False, amps,
                        signs=reference_signs(P.NOC[(8192, 3)], "random",
                                              seed=1))
    same_support = np.array_equal(f_true != 0, f_scr != 0)
    same_power = abs(float(np.vdot(f_true, f_true).real
                           - np.vdot(f_scr, f_scr).real)) < 1e-6
    rho = abs(np.vdot(f_true, f_scr)) / float(np.vdot(f_true, f_true).real)
    chk("G8 control keeps support+power", same_support and same_power,
        f"support {same_support}, power {same_power}")
    chk("G8 control decorrelates", rho < 0.25, f"|rho| = {rho:.3f}")

    # G9 -- the flattening amplitudes.  beta = 1 - 1/A is the null point of
    # the comb; getting A wrong moves the null and the comb comes back.
    ex = pilot_excess_amplitudes(amps)
    chk("G9 excess A_SP", abs(ex["scattered"] - (amps["scattered"] - 1)) < 1e-12,
        f"beta_sp = {1 - 1 / amps['scattered']:.4f}")
    chk("G9 excess A_CP", abs(ex["cp_common"] - (amps["cp_common"] - 1)) < 1e-12,
        f"beta_cp = {1 - 1 / amps['cp_common']:.4f}")

    if verbose:
        print("\n".join(lines))
        print(f"[pr_pilots] selftest {npass} pass, {nfail} fail")
    return npass, nfail


__all__ = [
    "FS_ATSC3", "BOOTSTRAP_SAMPLES_POST", "A_CP_DB",
    "symbol_rate_hz", "lattice_range_m", "lattice_doppler_hz",
    "pilot_amplitudes", "pilot_excess_amplitudes", "gate_table_9_17",
    "reference_signs", "pilot_classes", "pilot_field",
    "symbol_bins", "synth_pilot_symbol", "synth_pilot_stream",
    "selftest",
]


if __name__ == "__main__":
    sys.exit(0 if selftest()[1] == 0 else 1)
