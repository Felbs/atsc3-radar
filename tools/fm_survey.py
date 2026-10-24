"""fm_survey.py -- rank the FM band by carrier power from wide captures.

Reads one or more cs16 captures (centre + rate from their sidecars) and
reports the power in every 200 kHz FM channel 88.1-107.9 MHz, so the
illuminator is chosen from what THIS antenna hears tonight, not from a
list. No station names here: identity comes from RDS on the capture
itself (gr-radiotuna/tools/rds.py), and coordinates from the FCC table.

Usage: python tools/fm_survey.py lab/fmsurvey_*.cs16
"""
import json
import os
import sys

import numpy as np


def main():
    chans = {}
    for path in sys.argv[1:]:
        side = json.load(open(os.path.splitext(path)[0] + ".json"))
        fc, fs = float(side["centre_hz"]), float(side["rate_hz"])
        n = min(os.path.getsize(path) // 4, int(2.0 * fs))
        raw = np.fromfile(path, dtype=np.int16, count=2 * n)
        iq = (raw[0::2].astype(np.float32) + 1j * raw[1::2].astype(np.float32)) / 32768.0
        nfft = 1 << 15
        nseg = len(iq) // nfft
        win = np.hanning(nfft).astype(np.float32)
        psd = np.zeros(nfft, np.float64)
        for k in range(nseg):
            psd += np.abs(np.fft.fft(iq[k * nfft:(k + 1) * nfft] * win)) ** 2
        psd = np.fft.fftshift(psd / max(nseg, 1))
        fax = fc + np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))
        usable = np.abs(fax - fc) < 0.4 * fs
        for mhz10 in range(881, 1080, 2):
            f0 = mhz10 / 10 * 1e6
            sel = usable & (np.abs(fax - f0) < 90e3)
            if sel.sum() < 10:
                continue
            p = float(10 * np.log10(psd[sel].mean() + 1e-30))
            chans[mhz10 / 10] = max(chans.get(mhz10 / 10, -999.0), p)
    if not chans:
        print("no channels covered"); return
    floor = np.percentile(list(chans.values()), 20)
    print("FM channels by power (dB over the band's 20th percentile):")
    for f, p in sorted(chans.items(), key=lambda kv: -kv[1])[:20]:
        print(f"  {f:6.1f} MHz  {p - floor:+6.1f} dB")
    json.dump(chans, open(os.path.join(os.path.dirname(sys.argv[1]), "fm_survey.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
