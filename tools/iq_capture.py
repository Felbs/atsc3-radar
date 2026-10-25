#!/usr/bin/env python3
"""iq_capture.py -- a raw cs16 IQ capture at any centre frequency.

atsc3_capture.py with the RF-channel arithmetic replaced by --fc, so the
same capture laws (gain READ BACK not assumed, headroom measured on the
way past, integrity = samples match the wall clock, clipping = VOID,
sidecar after the bytes are safe) apply to an FM illuminator, or any
other, without touching the atsc3 repo.

Usage: python tools/iq_capture.py --fc 107.5e6 --rate 2e6 --secs 120 --ant "Antenna C" --out lab/x.cs16
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_RX, SOAPY_SDR_CS16

SoapySDR.SoapySDR_setLogLevel(SoapySDR.SOAPY_SDR_FATAL)

ap = argparse.ArgumentParser()
ap.add_argument("--fc", type=float, required=True, help="centre frequency, Hz")
ap.add_argument("--secs", type=float, default=120.0)
ap.add_argument("--rate", type=float, default=2.0e6)
ap.add_argument("--ant", default="Antenna C")
ap.add_argument("--ifgr", type=int, default=40)
ap.add_argument("--rfgain", type=int, default=None,
                help="rfgain_sel (LNA state); higher = MORE attenuation; read back")
ap.add_argument("--biast", default=None, choices=("on", "off"))
ap.add_argument("--allow-clipping", action="store_true")
ap.add_argument("--out", required=True)
ap.add_argument("--no-lock", action="store_true",
                help="caller already holds radio_lock (truth_capture does)")
a = ap.parse_args()

# Single-tenant radio: take the fleet lock unless the caller holds it.
# 10/24: a survey run without this tried to open the RSPdx under a
# priority-100 radiosonde hunt.
_locked = False
if not a.no_lock:
    sys.path.insert(0, r"Z:\src\gr-radiotuna\tools")
    import radio_lock
    if not radio_lock.acquire("iq_capture", f"capture {a.fc/1e6:.3f} MHz", 50, wait_s=5):
        h = radio_lock.status() or {}
        print(f"[lock] radio held by {h.get('owner', '?')} ({h.get('purpose', '')}) - not capturing")
        sys.exit(3)
    _locked = True

os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
sdr = SoapySDR.Device("driver=sdrplay")
sdr.setSampleRate(SOAPY_SDR_RX, 0, a.rate)
time.sleep(0.2)
fs = float(sdr.getSampleRate(SOAPY_SDR_RX, 0))
sdr.setAntenna(SOAPY_SDR_RX, 0, a.ant)
try:
    sdr.setGainMode(SOAPY_SDR_RX, 0, False)
except Exception:
    pass
sdr.setGain(SOAPY_SDR_RX, 0, "IFGR", a.ifgr)
if a.rfgain is not None:
    try:
        sdr.writeSetting("rfgain_sel", str(a.rfgain))
    except Exception as e:
        print(f"WARNING: rfgain_sel write failed: {e}")
if a.biast is not None:
    try:
        sdr.writeSetting("biasT_ctrl", "true" if a.biast == "on" else "false")
    except Exception as e:
        print(f"WARNING: biasT_ctrl write failed: {e}")
time.sleep(0.3)
try:
    rb_bias = str(sdr.readSetting("biasT_ctrl"))
except Exception:
    rb_bias = "<unreadable>"
try:
    rb_rf = str(sdr.readSetting("rfgain_sel"))
except Exception:
    rb_rf = "<unreadable>"
rb_if = float(sdr.getGain(SOAPY_SDR_RX, 0, "IFGR"))
rb_ant = sdr.getAntenna(SOAPY_SDR_RX, 0)
sdr.setFrequency(SOAPY_SDR_RX, 0, a.fc)
rb_fc = float(sdr.getFrequency(SOAPY_SDR_RX, 0))
print(f"gain readback: rfgain_sel={rb_rf} IFGR={rb_if:.0f} ant={rb_ant} "
      f"rate={fs/1e6:.3f} Msps biasT={rb_bias} fc={rb_fc/1e6:.4f} MHz")
st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CS16)
sdr.activateStream(st)
buf = np.empty(2 * 65536, np.int16)
n_want = int(a.secs * fs)
got = peak = n_full = n_over = 0
t0 = time.time()
with open(a.out, "wb") as f:
    while time.time() - t0 < 0.3:
        sdr.readStream(st, [buf], 65536, timeoutUs=300000)
    t_start = time.time()
    while got < n_want:
        r = sdr.readStream(st, [buf], 65536, timeoutUs=500000)
        if r.ret > 0:
            blk = buf[:2 * r.ret]
            f.write(blk.tobytes())
            got += r.ret
            m = int(np.abs(blk).max())
            if m > peak:
                peak = m
            n_full += int(np.count_nonzero(np.abs(blk) >= 32000))
        elif r.ret == -4:
            n_over += 1
        elif r.ret < 0:
            print(f"readStream {r.ret}", flush=True)
wall = time.time() - t_start
sdr.deactivateStream(st); sdr.closeStream(st); del sdr
if _locked:
    radio_lock.release("iq_capture")

frac = 100.0 * n_full / max(1, 2 * got)
dbfs = 20.0 * np.log10(peak / 32768.0) if peak else -999.0
print(f"CAPTURED {got} samples ({got/fs:.1f}s) in {wall:.1f}s wall -> {a.out}")
print(f"  headroom: max|s|={peak} ({dbfs:+.1f} dBFS)  full-scale samples={frac:.4f}%  overflows={n_over}")
ok = True
if abs(got / fs - wall) > 0.05 * wall:
    print(f"  VOID: samples/fs={got/fs:.2f}s != wall={wall:.2f}s (gapped)")
    ok = False
else:
    print(f"  integrity OK: samples/fs={got/fs:.2f}s == wall={wall:.2f}s")
if frac > 0.001 or peak >= 32767:
    print("  CLIPPED -> VOID. Re-run with more attenuation (--rfgain higher and/or --ifgr higher).")
    ok = False
else:
    print("  no clipping")
print(f"  VERDICT: {'USABLE' if ok else 'VOID -- do not bank this'}")

side = os.path.splitext(a.out)[0] + ".json"
quality = {}
try:
    n_q = min(got, int(1.0 * fs))
    raw = np.fromfile(a.out, dtype=np.int16, count=2 * n_q)
    iq = (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0
    quality["in_rms"] = float(np.sqrt((np.abs(iq) ** 2).mean()) * 32768)
    # power within +-100 kHz of centre vs the rest: a cheap "is the
    # illuminator there" proxy, the FM twin of the bootstrap peak ratio
    nfft = 1 << 16
    S = np.abs(np.fft.fftshift(np.fft.fft(iq[:nfft] * np.hanning(nfft)))) ** 2
    fax = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))
    inb = S[np.abs(fax) < 100e3].mean()
    outb = S[np.abs(fax) > 300e3].mean() if (np.abs(fax) > 300e3).any() else inb
    quality["centre_vs_band_db"] = round(float(10 * np.log10(inb / max(outb, 1e-30))), 1)
except Exception as e:
    quality["error"] = f"{type(e).__name__}: {e}"
try:
    json.dump({
        "capture": os.path.basename(a.out),
        "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        "centre_hz": rb_fc, "rate_hz": fs,
        "requested": {"rfgain_sel": a.rfgain, "ifgr": a.ifgr, "biast": a.biast, "ant": a.ant, "fc": a.fc},
        "readback": {"rfgain_sel": rb_rf, "ifgr": rb_if, "biast": rb_bias, "ant": rb_ant, "fc": rb_fc},
        "headroom": {"peak_abs": peak, "peak_dbfs": round(dbfs, 2),
                     "full_scale_pct": round(frac, 6), "overflows": n_over},
        "integrity": {"samples": got, "seconds": round(got / fs, 3), "wall_s": round(wall, 3),
                      "ok": bool(abs(got / fs - wall) <= 0.05 * wall)},
        "quality": quality,
        "verdict": "USABLE" if ok else "VOID",
    }, open(side, "w"), indent=1)
    print(f"  sidecar -> {side}")
    if "centre_vs_band_db" in quality:
        print(f"  quality: centre +-100 kHz is {quality['centre_vs_band_db']:+.1f} dB over the band edges, in_rms {quality['in_rms']:.0f}")
except Exception as _e:
    print(f"  ! sidecar step failed ({type(_e).__name__}: {_e}) -- the capture itself is unaffected")

sys.exit(0 if (ok or a.allow_clipping) else 1)
