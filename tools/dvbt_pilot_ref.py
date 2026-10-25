"""dvbt_pilot_ref.py -- a DVB-T radar reference synthesized from the pilots alone, from ONE channel.

The single-receiver trick of this repo (pr_pilots.py does it for ATSC 3.0)
applied to DVB-T: nothing is decoded. The receiver measures the OFDM
symbol timing, sample-clock drift and carrier offset from the signal's own
cyclic prefix, then synthesizes exactly what the transmitter put on the
pilot carriers -- values fixed by the standard -- and nothing else. That
waveform is the reference; the same channel, unmodified, is the
surveillance. Works far below the SNR a decode needs, because the pilots
are known rather than decided.

Pilots (EN 300 744 4.5), 8K mode, carriers k = 0..6816:
  continual  177 fixed carriers (dvbt_decode.CONT_K)
  scattered  k = 3 (l mod 4) + 12 p in symbol l
  value      4/3 x 2 (1/2 - w_k), w_k from the PRBS x^11 + x^2 + 1, all-ones start
TPS carriers are left out (68 of 6817; their bits would need decoding).
The pilots carry ~10 % of the carriers at 16/9 the power of a data cell:
about 9 dB less echo energy than a full reference, and they repeat every
4 symbols, which puts range ghosts at multiples of Tu/12 (22.4 km) and
Doppler ghosts at multiples of the symbol rate / 4 (223 Hz) -- far
outside a drone at under a kilometre.

`check` tests the PRBS and the scattered-pilot phase against the received
symbols before anything is built: adjacent pilots 12 carriers apart see
almost the same channel, so the right pilot values make the ratio
received/expected agree between neighbours; wrong values scramble it.

Usage:
  python tools/dvbt_pilot_ref.py check <file.npy> --fs 8e6 [--channel 1] [--t0 5]
  python tools/dvbt_pilot_ref.py build <file.npy> --fs 8e6 [--channel 1] --t0 0 --secs 20 --out lab/x/pilot_ref.npy
     -> a (2, N) complex64 .npy at 64/7 MS/s: row 0 the synthesized reference, row 1 the same channel
        resampled; then tools/two_channel_radar.py <out> --fs 9142857.142857143
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.signal import resample_poly

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FS = 64e6 / 7
NFFT, NCAR, KMIN = 8192, 6817, 8192 // 2 - 3408       # carrier k sits at FFT bin k + KMIN (DC = k 3408)
GUARD = {"1/4": 2048, "1/8": 1024, "1/16": 512, "1/32": 256}


def cont_k():
    import importlib.util
    spec = importlib.util.spec_from_file_location("_dd", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                      "dvbt_decode.py"))
    # dvbt_decode imports gnuradio; read only its CONT_K literal
    src = open(spec.origin, encoding="utf-8").read()
    a = src.index("CONT_K = np.array(")
    b = src.index("])", a) + 2
    ns = {"np": np}
    exec(src[a:b], ns)
    return ns["CONT_K"]


CONT = cont_k()


def prbs(variant):
    """w_k for k = 0..6816. variant 'a': taps 11 and 9 of the register as drawn in EN 300 744
    fig. 9 (x^11 + x^2 + 1); 'b': the mirror reading. check decides which is right."""
    reg = [1] * 11
    w = np.zeros(NCAR, np.int8)
    for k in range(NCAR):
        w[k] = reg[10]
        nb = reg[10] ^ (reg[8] if variant == "a" else reg[1])
        reg = [nb] + reg[:10]
    return w


def pilot_value(w):
    return (4.0 / 3.0) * 2.0 * (0.5 - w.astype(np.float64))


def scattered(l, phase):
    kp = 3 * ((l + phase) % 4) + 12 * np.arange(569)
    return kp[kp < NCAR]


def load_channel(path, fs, ch, t0, secs):
    X = np.load(path, mmap_mode="r", allow_pickle=False)
    o, n = int(t0 * fs), int(secs * fs)
    x = np.asarray(X[ch, o:o + n], np.complex64)
    if abs(fs - FS) > 1:
        from fractions import Fraction
        r = Fraction(FS / fs).limit_denominator(1000)
        x = resample_poly(x, r.numerator, r.denominator).astype(np.complex64)
    return x


def sync(x, cp, n_sym=None):
    """Symbol starts (fractional samples) and carrier offset from the cyclic prefix.
    Returns (starts, cfo_hz, metric_peak)."""
    L = NFFT + cp
    p = x[:-NFFT] * np.conj(x[NFFT:])
    m = np.convolve(p, np.ones(cp, np.complex64), "valid")
    e = np.convolve(np.abs(x[:-NFFT]) ** 2, np.ones(cp, np.float32), "valid")
    r = np.abs(m) / (e + 1e-12)
    n = n_sym or (len(r) // L - 1)
    fold = r[:n * L].reshape(n, L).mean(0)
    # 1. drift first: the folded peak of each ~1 s block, unwrapped, fitted on a line. A fixed
    #    search window around ONE offset loses the peaks once a 1 ppm sample clock has walked
    #    ~180 samples over 20 s; the fit was then biased and the synthesized reference slid
    #    against the channel at ~9 samples/s (peaks walking +295 m/s at 0 Hz).
    blk = max(1, int(round(FS / L)))
    cb = []
    for b0 in range(0, n - blk + 1, blk):
        f = r[b0 * L:(b0 + blk) * L].reshape(blk, L).mean(0)
        cb.append((b0 + blk / 2, float(np.argmax(f))))
    cb = np.array(cb)
    off = np.unwrap(cb[:, 1] * 2 * np.pi / L) * L / (2 * np.pi)
    if len(cb) >= 2:
        d_slope, d_icpt = np.polyfit(cb[:, 0], off, 1)
    else:
        d_slope, d_icpt = 0.0, float(np.argmax(fold))
    # 2. per-symbol peak searched around the PREDICTED position, then refit
    pk = []
    for i in range(n):
        pred = i * L + d_icpt + d_slope * i
        lo = int(round(pred)) - 48
        if lo < 0 or lo + 96 >= len(r):
            continue
        pk.append((i, lo + int(np.argmax(r[lo:lo + 96]))))
    pk = np.array(pk, float)
    A = np.vstack([pk[:, 0], np.ones(len(pk))]).T
    (slope, icpt), *_ = np.linalg.lstsq(A, pk[:, 1], rcond=None)
    for _ in range(2):                                  # reject outliers, refit
        res = pk[:, 1] - (slope * pk[:, 0] + icpt)
        keep = np.abs(res) < 3 * max(1.0, np.std(res))
        (slope, icpt), *_ = np.linalg.lstsq(A[keep], pk[keep, 1], rcond=None)
    # the metric window [j, j+cp) is the cyclic prefix: the FFT window starts at j + cp
    starts = icpt + slope * np.arange(n) + cp
    # carrier offset: the phase of the CP correlation at the peaks
    ph = np.angle(np.sum([m[int(round(s - cp))] for s in starts if 0 <= int(round(s - cp)) < len(m)]))
    cfo = -ph * FS / (2 * np.pi * NFFT)
    return starts, float(cfo), float(fold.max()), float((slope / L - 1) * 1e6)


def integer_offset(x, starts, cfo, span=60):
    """The CP phase gives the carrier offset modulo one carrier spacing (1116 Hz). The whole
    carriers are found by sliding the continual-pilot set: only the right shift makes those
    carriers repeat from symbol to symbol. Returns (shift in carriers, its coherence)."""
    Yf = symbols(x, starts[:120], cfo)
    best = (-1.0, 0)
    for sh in range(-span, span + 1):
        Y = Yf[:, KMIN + sh:KMIN + sh + NCAR]
        d = Y[1:, CONT] * np.conj(Y[:-1, CONT])
        c = float(np.abs(d.sum(1)).sum() / (np.abs(d).sum() + 1e-30))
        if c > best[0]:
            best = (c, sh)
    return best[1], best[0]


def full_sync(x, cp):
    """sync() plus the whole-carrier offset folded into the CFO."""
    starts, cfo, mpk, ppm = sync(x, cp)
    sh, coh = integer_offset(x, starts, cfo)
    return starts, cfo + sh * FS / NFFT, mpk, ppm, sh, coh


def symbols(x, starts, cfo):
    """FFT of each symbol at the fitted grid (integer start + fractional correction), CFO removed."""
    n = np.arange(len(x))
    y = x * np.exp(-2j * np.pi * cfo * n / FS).astype(np.complex64)
    out = []
    kb = np.fft.fftshift(np.fft.fftfreq(NFFT, 1 / NFFT))      # signed bin index after fftshift
    for s in starts:
        i0 = int(np.floor(s)); fr = s - i0
        if i0 + NFFT > len(y):
            break
        Y = np.fft.fftshift(np.fft.fft(y[i0:i0 + NFFT]))
        out.append(Y * np.exp(2j * np.pi * kb * fr / NFFT))     # undo the fractional start
    return np.array(out)


def check(a):
    x = load_channel(a.file, a.fs, a.channel, a.t0, 0.5)
    cp = GUARD[a.guard]
    starts, cfo, mpk, ppm, sh, coh = full_sync(x, cp)
    Y = symbols(x, starts, cfo)[:, KMIN:KMIN + NCAR]
    print(f"[check] ch{a.channel}: {len(Y)} symbols, CP metric {mpk:.2f}, CFO {cfo:+.1f} Hz "
          f"({sh:+d} whole carriers), sample clock {ppm:+.2f} ppm")
    # continual pilots: their received phase must be the same from symbol to symbol
    d = Y[1:, CONT] * np.conj(Y[:-1, CONT])
    print(f"[check] continual-pilot symbol-to-symbol coherence {np.abs(d.sum(1)).sum() / np.abs(d).sum():.2f} "
          "(timing and CFO right if high)")
    best = []
    for v in ("a", "b"):
        pv = pilot_value(prbs(v))
        for ph in range(4):
            sc = []
            for l in range(len(Y)):
                kp = scattered(l, ph)
                r = Y[l, kp] / pv[kp]
                sc.append(np.abs(np.sum(r[1:] * np.conj(r[:-1]))) / np.sum(np.abs(r[1:] * r[:-1])))
            best.append((float(np.mean(sc)), v, ph))
            print(f"   PRBS {v}, scattered phase {ph}: neighbour agreement {np.mean(sc):.3f}")
    # continual pilots are also PRBS-signed: they test the PRBS independently of the scattered phase
    for v in ("a", "b"):
        pv = pilot_value(prbs(v))
        kc = np.sort(CONT)
        r = Y[:, kc] / pv[kc]
        g = np.abs(np.sum(r[:, 1:] * np.conj(r[:, :-1]))) / np.sum(np.abs(r[:, 1:] * r[:, :-1]))
        print(f"   PRBS {v}, continual pilots: neighbour agreement {g:.3f}")
    best.sort(reverse=True)
    print(f"[check] best: PRBS {best[0][1]}, scattered phase {best[0][2]} ({best[0][0]:.3f}; next {best[1][0]:.3f})")
    return best[0], starts, cfo


def build(a):
    t = time.time()
    x = load_channel(a.file, a.fs, a.channel, a.t0, a.secs)
    cp = GUARD[a.guard]
    L = NFFT + cp
    starts, cfo, mpk, ppm, sh, coh = full_sync(x, cp)
    Y = symbols(x, starts[:200], cfo)[:, KMIN:KMIN + NCAR]
    pv_by = {v: pilot_value(prbs(v)) for v in ("a", "b")}
    sc = {}
    for v in ("a", "b"):
        for ph in range(4):
            vals = []
            for l in range(len(Y)):
                kp = scattered(l, ph)
                r = Y[l, kp] / pv_by[v][kp]
                vals.append(np.abs(np.sum(r[1:] * np.conj(r[:-1]))) / np.sum(np.abs(r[1:] * r[:-1])))
            sc[(v, ph)] = float(np.mean(vals))
    (v, ph), q = max(sc.items(), key=lambda kv: kv[1])
    pv = pv_by[v]
    print(f"[build] ch{a.channel} {a.t0:.1f}-{a.t0 + a.secs:.1f} s: {len(starts)} symbols, CFO {cfo:+.1f} Hz, "
          f"clock {ppm:+.2f} ppm, PRBS {v} phase {ph} (agreement {q:.3f})", flush=True)
    # synthesize: per symbol, pilots only, placed at the fitted fractional start, CFO re-applied
    ref = np.zeros(len(x), np.complex64)
    kb = np.arange(NCAR) + KMIN - NFFT // 2                     # signed FFT bin of carrier k
    base = np.zeros(NCAR, np.complex128)
    base[CONT] = pv[CONT]
    for l, s in enumerate(starts):
        X = base.copy()
        kp = scattered(l, ph)
        X[kp] = pv[kp]
        i0 = int(np.floor(s)); fr = s - i0
        F = np.zeros(NFFT, np.complex128)
        F[kb % NFFT] = X * np.exp(-2j * np.pi * kb * fr / NFFT)   # delay by the fractional part
        sym = np.fft.ifft(F) * np.sqrt(NFFT)
        full = np.concatenate([sym[-cp:], sym])                    # cyclic prefix
        a0 = i0 - cp
        if a0 < 0 or a0 + L > len(ref):
            continue
        ref[a0:a0 + L] = full.astype(np.complex64)
    n = np.arange(len(ref))
    ref *= np.exp(2j * np.pi * cfo * n / FS).astype(np.complex64)
    # Put the DIRECT PATH at zero lag and zero Doppler, where the canceller works. The CP
    # timing locks to the start of the channel (an early echo, or another transmitter of the
    # single-frequency network), and the CP carrier estimate leaves ~1 Hz; both show up as the
    # strongest peak of the reference-vs-channel map (UCY ch1: 1771 m, +1.0 Hz, 29.8 dB).
    # Cross-correlation over BOTH signs of lag: a range-Doppler map only looks at positive
    # range, and on UCY ch1 it latched onto an SFN echo 125 samples (4.1 km) behind a direct
    # path that the CP timing had placed at -24 samples.
    w0 = int(min(2.0, a.secs / 4) * FS)
    m = 1 << 22
    rs, xs = ref[w0:w0 + m], x[w0:w0 + m]
    C = np.fft.ifft(np.fft.fft(xs) * np.conj(np.fft.fft(rs)))          # C[k]: channel = ref delayed by k
    span = 1200
    cc = np.concatenate([C[-span:], C[:span]])
    lags = np.arange(-span, span)
    i = int(np.argmax(np.abs(cc)))
    lag = int(lags[i])
    floor_db = 20 * np.log10(np.abs(cc[i]) / np.median(np.abs(cc)))
    # Doppler from the phase walk of the product at that lag, in 20 ms blocks
    xs_l = np.roll(xs, -lag)
    prod = xs_l * np.conj(rs)
    bl = int(0.02 * FS)
    phw = np.unwrap(np.angle(prod[:len(prod) // bl * bl].reshape(-1, bl).sum(1)))
    dop = float(np.polyfit(np.arange(len(phw)) * bl / FS, phw, 1)[0] / (2 * np.pi))
    ref = np.roll(ref, lag)
    if lag > 0:
        ref[:lag] = 0
    elif lag < 0:
        ref[lag:] = 0
    ref *= np.exp(2j * np.pi * dop * n / FS).astype(np.complex64)
    print(f"[build] direct path found at {lag:+d} samples ({lag * 2.998e8 / FS:+.0f} m), {dop:+.2f} Hz, "
          f"{floor_db:.1f} dB over the correlation floor: moved to zero lag, zero Doppler", flush=True)
    out = np.lib.format.open_memmap(a.out, mode="w+", dtype=np.complex64, shape=(2, len(x)))
    out[0] = ref
    out[1] = x
    out.flush()
    meta = dict(source=os.path.abspath(a.file), channel=a.channel, t0=a.t0, secs=a.secs, fs=FS, guard=a.guard,
                cfo_hz=cfo, clock_ppm=ppm, prbs=v, scattered_phase=ph, pilot_agreement=q,
                direct_path_lag_samples=lag, direct_path_doppler_hz=dop,
                n_symbols=len(starts), note="row 0 = pilot-only synthesized reference, row 1 = the same channel")
    json.dump(meta, open(os.path.splitext(a.out)[0] + ".json", "w", encoding="utf-8"), indent=1, default=float)
    # direct-path check: the synthesized reference against the channel at zero range
    m = min(len(x), int(2 * FS))
    c = np.abs(np.vdot(ref[:m], x[:m])) / np.sqrt(np.vdot(ref[:m], ref[:m]).real * np.vdot(x[:m], x[:m]).real)
    print(f"[build] -> {a.out} ({len(x)/FS:.1f} s at {FS/1e6:.4f} MS/s) in {time.time()-t:.0f} s; "
          f"zero-lag correlation ref/channel {c:.3f} (pilot share of power ~ {np.sqrt(0.12):.2f} if all is right)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("check", "build"))
    ap.add_argument("file")
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--channel", type=int, default=1)
    ap.add_argument("--guard", default="1/4", choices=list(GUARD))
    ap.add_argument("--t0", type=float, default=5.0)
    ap.add_argument("--secs", type=float, default=20.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.cmd == "check":
        check(a)
    else:
        if not a.out:
            raise SystemExit("--out required")
        build(a)


if __name__ == "__main__":
    main()
