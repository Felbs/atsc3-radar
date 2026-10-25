"""dvbt_decode.py -- DVB-T (8K, 8 MHz) from a recording to an MPEG transport stream.

GNU Radio's gr-dtv receiver chain (as in its dvbt_rx_8k example), run
headless on a file: resample to the DVB-T clock (64/7 MS/s), OFDM symbol
acquisition, FFT, pilots and channel estimate, demap, symbol and bit
deinterleave, Viterbi, Forney deinterleave, RS(204,188), energy descramble.
Needs radioconda's python (gnuradio + gr-dtv).

The transmitter's settings are not in the file, but DVB-T broadcasts them:
--tps reads the 68 TPS carriers (DBPSK, one bit per symbol, 68-symbol
frames) and prints mode, guard, constellation, hierarchy and code rates.
gr-dtv's symbol acquisition steps the FFT window by a few samples from
symbol to symbol (+-6 on the 2022 UCY recording), which tilts each
symbol's phase across the band; the continual pilots measure the tilt per
symbol pair and it is removed before the TPS bits are read. --tps also
reports the MER of the equalized data cells: the number that says whether
a picture is possible at all (16-QAM 5/6 needs about 14-15 dB). The guard
interval must be given; measure it from the symbol period (the cyclic
prefix repeats the symbol tail). --scan tries every constellation and
rate and scores the transport stream (PIDs that repeat, PAT, continuity
counters, transport-error flags), for when the TPS will not decode.

Input: a .npy of shape (C, N) complex64 (--channel picks the row; read
through the file offset, no copy) or a raw complex64 file.

Usage:
  radioconda/python tools/dvbt_decode.py <file> --fs 8e6 --guard 1/4 --tps [--secs 0.5]   # settings + MER
  radioconda/python tools/dvbt_decode.py <file> --fs 8e6 --guard 1/4 --scan [--secs 6]
  radioconda/python tools/dvbt_decode.py <file> --fs 8e6 --guard 1/4 --qam 16 --rate 2/3 --out lab/x.ts
"""
import argparse
import collections
import os
import sys
import time
from fractions import Fraction

import numpy as np
from gnuradio import blocks, dtv, fft, filter, gr
from gnuradio.fft import window

FS_DVBT = 64e6 / 7
MODS = {"qpsk": dtv.MOD_QPSK, "4": dtv.MOD_QPSK, "16": dtv.MOD_16QAM, "64": dtv.MOD_64QAM}
RATES = {"1/2": dtv.C1_2, "2/3": dtv.C2_3, "3/4": dtv.C3_4, "5/6": dtv.C5_6, "7/8": dtv.C7_8}
GUARDS = {"1/4": (dtv.GI_1_4, 2048), "1/8": (dtv.GI_1_8, 1024), "1/16": (dtv.GI_1_16, 512), "1/32": (dtv.GI_1_32, 256)}

# 8K carrier indices (EN 300 744 tables 6 and 7)
TPS_K = np.array([34, 50, 209, 346, 413, 569, 595, 688, 790, 901, 1073, 1219, 1262, 1286, 1469, 1594, 1687, 1738, 1754,
                  1913, 2050, 2117, 2273, 2299, 2392, 2494, 2605, 2777, 2923, 2966, 2990, 3173, 3298, 3391, 3442, 3458,
                  3617, 3754, 3821, 3977, 4003, 4096, 4198, 4309, 4481, 4627, 4670, 4694, 4877, 5002, 5095, 5146, 5162,
                  5321, 5458, 5525, 5681, 5707, 5800, 5902, 6013, 6185, 6331, 6374, 6398, 6581, 6706, 6799])
CONT_K = np.array([0, 48, 54, 87, 141, 156, 192, 201, 255, 279, 282, 333, 432, 450, 483, 525, 531, 618, 636, 714, 759,
                   765, 780, 804, 873, 888, 918, 939, 942, 969, 984, 1050, 1101, 1107, 1110, 1137, 1140, 1146, 1206,
                   1269, 1323, 1377, 1491, 1683, 1704, 1752, 1758, 1791, 1845, 1860, 1896, 1905, 1959, 1983, 1986,
                   2037, 2136, 2154, 2187, 2229, 2235, 2322, 2340, 2418, 2463, 2469, 2484, 2508, 2577, 2592, 2622,
                   2643, 2646, 2673, 2688, 2754, 2805, 2811, 2814, 2841, 2844, 2850, 2910, 2973, 3027, 3081, 3195,
                   3387, 3408, 3456, 3462, 3495, 3549, 3564, 3600, 3609, 3663, 3687, 3690, 3741, 3840, 3858, 3891,
                   3933, 3939, 4026, 4044, 4122, 4167, 4173, 4188, 4212, 4281, 4296, 4326, 4347, 4350, 4377, 4392,
                   4458, 4509, 4515, 4518, 4545, 4548, 4554, 4614, 4677, 4731, 4785, 4899, 5091, 5112, 5160, 5166,
                   5199, 5253, 5268, 5304, 5313, 5367, 5391, 5394, 5445, 5544, 5562, 5595, 5637, 5643, 5730, 5748,
                   5826, 5871, 5877, 5892, 5916, 5985, 6000, 6030, 6051, 6054, 6081, 6096, 6162, 6213, 6219, 6222,
                   6249, 6252, 6258, 6318, 6381, 6435, 6489, 6603, 6795, 6816])


def npy_layout(path, channel):
    """(byte offset of the channel's first sample, samples per channel) for a C-order (C, N) complex64 .npy."""
    with open(path, "rb") as f:
        ver = np.lib.format.read_magic(f)
        shape, fortran, dtype = (np.lib.format.read_array_header_1_0(f) if ver == (1, 0)
                                 else np.lib.format.read_array_header_2_0(f))
        hdr = f.tell()
    if dtype != np.complex64 or fortran:
        raise SystemExit(f"need C-order complex64, got {dtype} fortran={fortran}")
    if len(shape) == 1:
        return hdr, shape[0]
    return hdr + channel * shape[1] * 8, shape[1]


class FftTap(gr.top_block):
    """Resample, acquire, FFT, and the gr-dtv equalizer; keep both outputs."""
    def __init__(self, path, fs, guard, channel, t0, secs, snr_db, qam="16", rate="2/3"):
        gr.top_block.__init__(self, "dvbt_fft_tap")
        off_b, n_ch = npy_layout(path, channel) if path.endswith(".npy") else (0, os.path.getsize(path) // 8)
        src = blocks.file_source(gr.sizeof_gr_complex, path, False, off_b // 8 + int(t0 * fs), int(secs * fs))
        ratio = Fraction(FS_DVBT / fs).limit_denominator(1000)
        gi, cp = GUARDS[guard]
        acq = dtv.dvbt_ofdm_sym_acquisition(1, 8192, 6817, cp, snr_db)
        ff = fft.fft_vcc(8192, True, window.rectangular(8192), True, 1)
        ref = dtv.dvbt_demod_reference_signals(gr.sizeof_gr_complex, 8192, 6048, MODS[qam], dtv.NH, RATES[rate],
                                               RATES[rate], gi, dtv.T8k, 1, 0)
        self.sf = blocks.vector_sink_c(8192)
        self.se = blocks.vector_sink_c(6048)
        chain = [src]
        if ratio != 1:
            chain.append(filter.rational_resampler_ccc(
                ratio.numerator, ratio.denominator,
                filter.firdes.low_pass(ratio.numerator, fs * ratio.numerator, 3.9e6, 0.4e6)))
        chain += [acq, ff]
        self.connect(*chain)
        self.connect(ff, self.sf)
        self.connect(ff, ref, self.se)


def read_tps(F):
    """F: (n_sym, 8192) shifted-FFT symbols -> (first-carrier bin, median TPS coherence, decoded TPS frames)."""
    best = None
    for b0 in range(687 - 200, 687 + 201):         # integer carrier offset: continual pilots are constant
        z = F[:, CONT_K + b0]
        d = z[1:] * np.conj(z[:-1])
        coh = np.abs(d.sum(1)).sum() / (np.abs(d).sum() + 1e-30)
        if best is None or coh > best[0]:
            best = (coh, b0)
    b0 = best[1]
    slopes = np.linspace(-2 * np.pi * 64 / 8192, 2 * np.pi * 64 / 8192, 4097)
    E = np.exp(-1j * np.outer(slopes, CONT_K))
    bits, q = [], []
    for i in range(len(F) - 1):
        dc = F[i + 1, CONT_K + b0] * np.conj(F[i, CONT_K + b0])
        j = int(np.argmax(np.abs(E @ dc)))           # this pair's timing step, as a phase tilt
        rot = E[j] @ dc
        dt = F[i + 1, TPS_K + b0] * np.conj(F[i, TPS_K + b0]) * np.exp(-1j * slopes[j] * TPS_K) * np.conj(rot) / abs(rot)
        q.append(abs(np.exp(1j * np.angle(dt)).mean()))
        bits.append(int(np.real(dt.sum()) < 0))      # DBPSK: a phase flip is a 1
    s = "".join(map(str, bits))
    C_ = {"00": "QPSK", "01": "16-QAM", "10": "64-QAM"}
    R_ = {"000": "1/2", "001": "2/3", "010": "3/4", "011": "5/6", "100": "7/8"}
    G_ = {"00": "1/32", "01": "1/16", "10": "1/8", "11": "1/4"}
    T_ = {"00": "2K", "01": "8K", "10": "4K"}
    frames = []
    for i in range(len(s) - 67):
        if s[i:i + 16] in ("0011010111101110", "1100101000010001"):
            b = s[i:i + 67]

            def g(a, n, b=b):                         # TPS bit s_a (s_1 = first sync bit)
                return b[a - 1:a - 1 + n]
            frames.append(dict(symbol=i, frame=int(g(23, 2), 2), constellation=C_.get(g(25, 2), g(25, 2)),
                               hierarchy=g(27, 3), rate_hp=R_.get(g(30, 3), g(30, 3)), rate_lp=R_.get(g(33, 3), g(33, 3)),
                               guard=G_.get(g(36, 2), g(36, 2)), mode=T_.get(g(38, 2), g(38, 2))))
    return b0, float(np.median(q)), frames


def mer_db(z, qam):
    lv = {"qpsk": [-1, 1], "4": [-1, 1], "16": [-3, -1, 1, 3], "64": [-7, -5, -3, -1, 1, 3, 5, 7]}[qam]
    Q = np.array([x + 1j * y for x in lv for y in lv])
    Q = Q / np.sqrt(np.mean(np.abs(Q) ** 2))
    z = z / np.sqrt(np.mean(np.abs(z) ** 2))
    d = Q[np.argmin(np.abs(z[:, None] - Q[None]), 1)]
    return float(10 * np.log10(np.mean(np.abs(d) ** 2) / np.mean(np.abs(z - d) ** 2)))


class Rx(gr.top_block):
    def __init__(self, path, fs, guard, qam, rate, out, t0, secs, channel, snr_db):
        gr.top_block.__init__(self, "dvbt_rx_file")
        if path.endswith(".npy"):
            off_b, n_ch = npy_layout(path, channel)
        else:
            off_b, n_ch = 0, os.path.getsize(path) // 8
        assert off_b % 8 == 0
        start = int(t0 * fs)
        n = n_ch - start if secs is None else min(int(secs * fs), n_ch - start)
        src = blocks.file_source(gr.sizeof_gr_complex, path, False, off_b // 8 + start, n)
        ratio = Fraction(FS_DVBT / fs).limit_denominator(1000)
        if ratio != 1:
            taps = filter.firdes.low_pass(ratio.numerator, fs * ratio.numerator, 3.9e6, 0.4e6)
            rs = filter.rational_resampler_ccc(ratio.numerator, ratio.denominator, taps)
        gi, cp = GUARDS[guard]
        mod, cr = MODS[qam], RATES[rate]
        acq = dtv.dvbt_ofdm_sym_acquisition(1, 8192, 6817, cp, snr_db)
        ff = fft.fft_vcc(8192, True, window.rectangular(8192), True, 1)
        ref = dtv.dvbt_demod_reference_signals(gr.sizeof_gr_complex, 8192, 6048, mod, dtv.NH, cr, cr, gi, dtv.T8k, 1, 0)
        dem = dtv.dvbt_demap(6048, mod, dtv.NH, dtv.T8k, 1)
        sdi = dtv.dvbt_symbol_inner_interleaver(6048, dtv.T8k, 0)
        bdi = dtv.dvbt_bit_inner_deinterleaver(6048, mod, dtv.NH, dtv.T8k)
        v2s = blocks.vector_to_stream(gr.sizeof_char, 6048)
        vit = dtv.dvbt_viterbi_decoder(mod, dtv.NH, cr, 768)
        cdi = dtv.dvbt_convolutional_deinterleaver(136, 12, 17)
        rsd = dtv.dvbt_reed_solomon_dec(2, 8, 0x11d, 255, 239, 8, 51, 8)
        eds = dtv.dvbt_energy_descramble(8)
        snk = blocks.file_sink(gr.sizeof_char, out, False)
        chain = [src] + ([rs] if ratio != 1 else []) + [acq, ff, ref, dem, sdi, bdi, v2s, vit, cdi, rsd, eds, snk]
        self.connect(*chain)


def ts_score(path):
    """Transport-stream health of a decoded file."""
    b = np.fromfile(path, np.uint8)
    if len(b) < 188 * 50:
        return dict(packets=len(b) // 188, score=0.0)
    # align to the sync byte
    best = max(range(188), key=lambda o: np.count_nonzero(b[o::188][:2000] == 0x47))
    p = b[best:len(b) - (len(b) - best) % 188].reshape(-1, 188)
    sync = float(np.mean(p[:, 0] == 0x47))
    tei = float(np.mean(p[:, 1] & 0x80 > 0))
    pid = ((p[:, 1].astype(int) & 0x1F) << 8) | p[:, 2]
    cnt = collections.Counter(pid.tolist())
    top = cnt.most_common(12)
    top_frac = sum(c for _, c in top) / len(pid)
    # continuity counters on the busiest real PID (not null 0x1FFF)
    cc_ok = None
    for q, c in top:
        if q != 0x1FFF and c > 50:
            sel = p[pid == q]
            has_payload = (sel[:, 3] & 0x10) > 0
            cc = (sel[has_payload, 3] & 0x0F).astype(int)
            cc_ok = float(np.mean((np.diff(cc) % 16) == 1)) if len(cc) > 2 else None
            break
    score = sync * (1 - tei) * top_frac * (cc_ok if cc_ok is not None else 0.0)
    return dict(packets=len(p), sync=sync, tei=tei, pat=cnt.get(0, 0), top_pids=[f"0x{q:04x}:{c}" for q, c in top[:6]],
                top_frac=top_frac, cc_ok=cc_ok, score=score)


def run(a, qam, rate, out, secs, idle_s=8.0, verbose=False):
    """Run until the output stops growing. gr-dtv's chain does not always
    drain to completion when a file source ends (the OFDM acquisition block
    holds its last partial symbol), so tb.run() can hang; this stops it."""
    tb = Rx(a.file, a.fs, a.guard, qam, rate, out, a.t0, secs, a.channel, a.snr)
    t = time.time()
    tb.start()
    last, last_t = -1, time.time()
    while True:
        time.sleep(1.0)
        n = os.path.getsize(out) if os.path.exists(out) else 0
        if verbose and int(time.time() - t) % 10 == 0:
            print(f"    {time.time()-t:5.0f} s: {n/1e6:.2f} MB", flush=True)
        if n != last:
            last, last_t = n, time.time()
        elif time.time() - last_t > idle_s:
            break
    tb.stop(); tb.wait()
    return time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--guard", default="1/4", choices=list(GUARDS))
    ap.add_argument("--qam", default="16", choices=list(MODS))
    ap.add_argument("--rate", default="2/3", choices=list(RATES))
    ap.add_argument("--channel", type=int, default=0)
    ap.add_argument("--snr", type=float, default=10.0, help="acquisition SNR hint, dB")
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--secs", type=float, default=None)
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--tps", action="store_true", help="read the transmitter's settings from the TPS carriers, and the MER")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    base = os.path.splitext(a.out or a.file)[0]
    if a.tps:
        secs = a.secs or 0.5
        tb = FftTap(a.file, a.fs, a.guard, a.channel, a.t0, secs, a.snr, a.qam, a.rate)
        tb.start()
        t = time.time()
        while len(tb.sf.data()) // 8192 < int(secs * 880) and time.time() - t < 120:
            time.sleep(1)
        tb.stop(); tb.wait()
        F = np.array(tb.sf.data(), np.complex64).reshape(-1, 8192)
        Z = np.array(tb.se.data(), np.complex64).reshape(-1, 6048)
        b0, q, frames = read_tps(F)
        print(f"{len(F)} OFDM symbols; first carrier at FFT bin {b0} (nominal 687); TPS coherence {q:.2f}")
        for f in frames[:4]:
            print(f"  TPS frame {f['frame']} at symbol {f['symbol']}: {f['mode']} guard {f['guard']} {f['constellation']} "
                  f"hierarchy {f['hierarchy']} code rate HP {f['rate_hp']} LP {f['rate_lp']}")
        if not frames:
            print("  no TPS sync word found (wrong guard, or the signal is too weak)")
        if len(Z) > 20:
            qam = {"QPSK": "qpsk", "16-QAM": "16", "64-QAM": "64"}.get(frames[0]["constellation"], a.qam) if frames else a.qam
            print(f"  MER of the gr-dtv-equalized data cells as {qam}: {mer_db(Z[len(Z) // 4:].ravel()[:400000], qam):.1f} dB "
                  "(needs about: 16-QAM 5/6 14-15 dB, 16-QAM 1/2 10, QPSK 1/2 5)")
        return
    if a.scan:
        secs = a.secs or 6.0
        res = []
        for qam in ("qpsk", "16", "64"):
            for rate in RATES:
                out = f"{base}_scan_{qam}_{rate.replace('/', '')}.ts"
                dt = run(a, qam, rate, out, secs)
                s = ts_score(out)
                res.append((s["score"], qam, rate, s))
                print(f"  {qam:>4s} {rate}: {s['packets']:6d} pkts  sync {s.get('sync', 0):.2f}  TEI {s.get('tei', 0):.2f}  "
                      f"PAT {s.get('pat', 0):4d}  top-12 PIDs {s.get('top_frac', 0)*100:5.1f}%  "
                      f"CC ok {s['cc_ok'] if s.get('cc_ok') is not None else '-'}  score {s['score']:.3f}  ({dt:.0f} s)", flush=True)
                os.remove(out)
        res.sort(key=lambda r: -r[0])
        sc, qam, rate, s = res[0]
        print(f"best: {qam} {rate} score {sc:.3f}; PIDs {s.get('top_pids')}")
        return
    out = a.out or base + ".ts"
    dt = run(a, a.qam, a.rate, out, a.secs, verbose=True)
    s = ts_score(out)
    print(f"{a.qam} {a.rate}: {s['packets']} packets in {dt:.0f} s -> {out}")
    print(f"  sync {s.get('sync', 0):.3f}  TEI {s.get('tei', 0):.3f}  PAT {s.get('pat', 0)}  CC ok {s.get('cc_ok')}  PIDs {s.get('top_pids')}")


if __name__ == "__main__":
    main()
