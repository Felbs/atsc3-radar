"""rxtune_dvbt.py -- gr-rxtune on a DVB-T RECORDING: search the receiver's software knobs against MER.

A recording has no gain left to turn: the RF chain was set when it was
made. What is left is the software receiver, and gr-rxtune's Knob
interface takes any function, so its knobs here are:

  path     ch0 | ch1 | both   which antenna, or both aligned and added
                              coherently (the one knob that can add signal:
                              up to +3 dB with two similar antennas)
  acq_snr  5 | 10 | 20 | 30   gr-dtv symbol-acquisition SNR assumption
  t0       5 points in time   where in the recording to listen (the
                              "hour of day" axis, compressed into 70 s)

Dial: MER of gr-dtv's equalized data cells (tools/dvbt_decode.py), with
the cliff at the code rate's published threshold (16-QAM 5/6: 13.5 dB, Gaussian channel).
Liveness, rxtune's "no recommendation without proof": the full
gr-dtv chain to a transport stream, run only for cells within 2 dB of the
cliff, and alive only if the stream is real (PAT, repeating PIDs,
continuity counters). The verdict is rxtune's own.

Run with radioconda's python (gnuradio + gr-dtv). gr-rxtune is imported
from its source tree (--rxtune, default Z:/src/gr-rxtune/src).

Usage:
  radioconda/python tools/rxtune_dvbt.py <file.npy> --fs 8e6 --guard 1/4 --qam 16 --rate 5/6 [--cliff 13.5]
"""
import argparse
import gc
import json
import shutil
import math
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dvbt_decode as D                                        # noqa: E402

# required C/N for DVB-T, EN 300 744 table A.1, Gaussian channel (dB)
CLIFF = {("qpsk", "1/2"): 3.1, ("qpsk", "2/3"): 4.9, ("qpsk", "3/4"): 5.9, ("qpsk", "5/6"): 6.9, ("qpsk", "7/8"): 7.7,
         ("16", "1/2"): 8.8, ("16", "2/3"): 11.1, ("16", "3/4"): 12.5, ("16", "5/6"): 13.5, ("16", "7/8"): 13.9,
         ("64", "1/2"): 14.4, ("64", "2/3"): 16.5, ("64", "3/4"): 18.0, ("64", "5/6"): 19.3, ("64", "7/8"): 20.1}


def window(path, fs, path_knob, t0, secs, scratch):
    """Write the chosen antenna path for [t0, t0+secs) as raw complex64; return its filename."""
    X = np.load(path, mmap_mode="r", allow_pickle=False)
    o, n = int(t0 * fs), int(secs * fs)
    if path_knob in ("ch0", "ch1"):
        x = np.asarray(X[int(path_knob[-1]), o:o + n], np.complex64)
    else:
        a = np.asarray(X[0, o:o + n], np.complex64)
        b = np.asarray(X[1, o:o + n], np.complex64)
        # align b to a: integer lag from the cross-correlation, then one complex gain
        m = 1 << 20
        C = np.fft.ifft(np.fft.fft(a[:m]) * np.conj(np.fft.fft(b[:m])))
        lag = int(np.argmax(np.abs(C)))
        lag = lag if lag < m // 2 else lag - m
        b = np.roll(b, lag)
        g = np.vdot(b, a) / np.vdot(b, b)
        # maximum-ratio weights for two branches with independent noise: each
        # branch scaled by its own signal-to-noise; equal-SNR approximation:
        x = (a + g * b).astype(np.complex64)
    out = os.path.join(scratch, f"w_{path_knob}_{t0:05.1f}.c64")
    x.tofile(out)
    return out


class Recording:
    """rxtune Frontend over a recording: software knobs only."""

    def __init__(self, a, scratch):
        from rxtune.knob import FunctionKnob, KnobSpec
        self.a, self.scratch = a, scratch
        dur = np.load(a.file, mmap_mode="r", allow_pickle=False).shape[1] / a.fs
        self.state = {"path": "ch0", "acq_snr": 10, "t0": 5.0}
        times = tuple(round(float(v), 1) for v in np.linspace(3.0, dur - a.secs - 1.0, 5))

        def knob(name, choices, role, ordered):
            return FunctionKnob(KnobSpec(name, kind="choice", choices=choices, ordered=ordered, sense="none",
                                         role=role, cost="slow", settle_s=0.0),
                                setter=lambda v, n=name: self.state.__setitem__(n, v),
                                getter=lambda n=name: self.state[n])
        paths = tuple(a.paths.split(","))
        self.state["path"] = paths[0]
        self._knobs = {"path": knob("path", paths, "path", False),
                       "acq_snr": knob("acq_snr", (5, 10, 20, 30), "config", True),
                       "t0": knob("t0", times, "virtual", True)}

    def knobs(self):
        return self._knobs

    def level(self):
        return None


class MerMeasurer:
    """One cell = decode a slice of the recording; the dial is gr-dtv's MER."""

    def __init__(self, fe, spec, a):
        self.fe, self.spec, self.a = fe, spec, a
        self.heartbeat = None
        self.log = []

    def measure(self, extra_settle_s=0.0, window_scale=1.0):
        from rxtune.dial import Reading
        a, st = self.a, dict(self.fe.state)
        t = time.time()
        f = window(a.file, a.fs, st["path"], float(st["t0"]), a.secs, self.fe.scratch)
        tb = D.FftTap(f, a.fs, a.guard, 0, 0.0, a.secs, float(st["acq_snr"]), a.qam, a.rate)
        tb.start()
        # until the equalizer stops producing (the file is short; gr-dtv never
        # signals the end), or 120 s
        last, still = -1, 0
        while time.time() - t < 120 and still < 6:
            time.sleep(0.5)
            n = len(tb.se.data())
            still = still + 1 if n == last and n > 0 else 0
            last = n
        tb.stop(); tb.wait()
        Z = np.array(tb.se.data(), np.complex64).reshape(-1, 6048)
        del tb                                   # GNU Radio's file source holds the file until freed (Windows)
        gc.collect()
        Z = Z[min(40, len(Z) // 4):]
        vals = [D.mer_db(Z[i:i + 20].ravel(), a.qam) for i in range(0, len(Z) - 19, 20)]
        if not vals:
            r = Reading(None, None, 0, note="no equalized symbols")
        else:
            v = float(np.median(vals))
            s = self.spec.score(v)
            alive, note = False, "MER more than 2 dB under the cliff: no decode attempted"
            if v >= self.spec.cliff - 2.0:
                ts = os.path.join(self.fe.scratch, "cell.ts")
                ns = argparse.Namespace(file=f, fs=a.fs, guard=a.guard, channel=0, snr=float(st["acq_snr"]), t0=0.0)
                D.run(ns, a.qam, a.rate, ts, a.secs)
                sc = D.ts_score(ts)
                alive = sc["score"] > 0.5 and sc.get("pat", 0) > 0
                note = f"decode: {sc['packets']} packets, PAT {sc.get('pat', 0)}, CC ok {sc.get('cc_ok')}"
            r = Reading(v, s, len(vals), p10=float(np.percentile(vals, 10)), p90=float(np.percentile(vals, 90)),
                        alive=alive, live_rate=0.0, t=time.time(), note=note)
        gc.collect()
        try:
            os.remove(f)
        except OSError:
            pass                                 # removed with the scratch dir at the end
        self.log.append(dict(setting=st, mer=r.value, p10=r.p10, p90=r.p90, alive=r.alive, note=r.note,
                             secs=round(time.time() - t, 1)))
        print(f"  cell {st}: MER {r.value if r.value is None else round(r.value, 2)} dB "
              f"(p10 {r.p10 if r.p10 is None else round(r.p10, 1)}, p90 {r.p90 if r.p90 is None else round(r.p90, 1)})"
              f"  alive {r.alive}  [{r.note}]  {time.time() - t:.0f} s", flush=True)
        return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--fs", type=float, required=True)
    ap.add_argument("--guard", default="1/4")
    ap.add_argument("--qam", default="16")
    ap.add_argument("--rate", default="5/6")
    ap.add_argument("--cliff", type=float, default=None, help="MER decode threshold (default: EN 300 744 table)")
    ap.add_argument("--secs", type=float, default=1.0, help="seconds of recording per cell")
    ap.add_argument("--paths", default="both,ch0,ch1", help="antenna paths to search, comma list")
    ap.add_argument("--max-cells", type=int, default=70)
    ap.add_argument("--rxtune", default="Z:/src/gr-rxtune/src")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    sys.path.insert(0, a.rxtune)
    from rxtune.dial import DialSpec
    from rxtune.loop import tune

    cliff = a.cliff if a.cliff is not None else CLIFF[(a.qam, a.rate)]
    spec = DialSpec("MER", "dB", cliff=cliff, settle_s=0.0, window_s=a.secs, continuous_below_cliff=True,
                    resolution=0.3, min_samples=3)
    scratch = tempfile.mkdtemp(prefix="rxtune_dvbt_", dir=os.path.join(os.path.dirname(os.path.abspath(a.file))))
    fe = Recording(a, scratch)
    meas = MerMeasurer(fe, spec, a)
    print(f"[rxtune] {os.path.basename(a.file)}: DVB-T {a.qam}-QAM {a.rate}, cliff {cliff} dB MER; "
          f"knobs path {fe.knobs()['path'].spec.choices}, acq_snr {fe.knobs()['acq_snr'].spec.choices}, "
          f"t0 {fe.knobs()['t0'].spec.choices}", flush=True)

    class _Dial:                                   # tune() wants a dial object; the measurer does the reading
        def __init__(self, s):
            self.spec = s

        def read(self):
            return None
    rep = tune(fe, _Dial(spec), None, search=["path", "acq_snr", "t0"], measurer=meas, agc_off=False,
               coarse_points=5, max_cells=a.max_cells)
    v = rep.verdict
    print(f"\nVERDICT  {getattr(v, 'shape', '?')}")
    for line in str(v).splitlines():
        print("  " + line)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.file)), "rxtune_dvbt.json")
    json.dump(dict(file=a.file, qam=a.qam, rate=a.rate, cliff=cliff, cells=meas.log, verdict=str(v),
                   warnings=rep.warnings), open(out, "w", encoding="utf-8"), indent=1, default=str)
    print(f"-> {out}")
    gc.collect()
    shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()
