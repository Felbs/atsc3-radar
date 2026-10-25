# Someone else's data: the University of Cyprus DVB-T drone recording

The first test of this code on a recording it was not written for: another
continent, another transmitter standard, another radio, another target.

**The data.** KIOS / University of Cyprus passive radar dataset
(kios.ucy.ac.cy/spa, "Passive Radar - Dataset"): a USRP B210 with two
yagis, DVB-T illuminator, a DJI Matrice M300 drone flown over the campus.
The public folder holds one file, `RX1/12_12_56.npy`, 8.96 GB: shape
(2, 560 000 000) complex64. No sample rate, centre frequency or flight log
came with it. Nothing from it is in this repo.

**What the file says about itself.**
- Sample rate: the reference autocorrelates at a lag of exactly 7168
  samples, the 896 us useful-symbol time of DVB-T 8K mode, so 8.000 MS/s
  (70.0 s). The reference fills 95.2 % of the band, DVB-T's 7.61 of 8 MHz.
- Channel 0 is the reference, channel 1 the surveillance: their
  cross-correlation peaks at 0.51, two samples apart, zero Doppler offset.
- Each channel holds the DVB-T signal only a few dB over its own noise;
  in-band coherence between them is 0.57-0.66. A canceller removes only
  what the channels share, so cancellation is capped at about 4 dB and
  `two_channel_radar.py` gets 3.0-3.2 dB on every dwell. That is this
  recording's ceiling, not a failure, and not comparable with the 25-30 dB
  of a single strong channel here: the residue is noise, which the 0.5 s
  coherent integration handles.

**What the radar found** (`tools/two_channel_radar.py <npy> --fs 8e6
--notch-hz 2.5`, 140 dwells of 0.5 s, 37.5 m cells to 3 km, CFAR 1e-6;
97 detections, expected false alarms about 3):

| track | time | range | Doppler | direction |
|---|---|---|---|---|
| 2 | 0.2-2.2 s | 375 -> 300 m | +58 Hz | inbound |
| 0 | 13.2-19.8 s | 262 -> 487 m | -62 Hz | outbound |
| 1 | 50.8-56.8 s | 675 -> 450 m | +73 Hz | inbound |
| 3 | 67.8-69.2 s | 412 -> 450 m | -63 Hz | outbound |

One target flying out and back across the field, at 260-680 m bistatic.

**The check that needed no truth.** The carrier was never given. A real
echo's range walks at -lambda x Doppler, so every track that walks gives
a wavelength; noise and clutter links give random ones, and the wrong sign
half the time. All four tracks walk the way their Doppler says (none the
wrong way). The two long enough to fit give 0.584 m and 0.563 m, within
4 % of each other: 513 and 533 MHz, a UHF TV channel. A 0.574 m wavelength
and 62-73 Hz put the bistatic range rate at 36-42 m/s; near the receiver
with the transmitter far away that is up to twice the ground speed, so
18-21 m/s, inside the M300's 23 m/s.

**Not established.** The drone's position: no flight log came with the
file. The illuminator's frequency: the UCY pages do not give it, and the
only Cyprus channel list found (Troodos, UHF 33 = 570 MHz) is about 16
years old and predates the 2026 DVB-T2 move. The 523 MHz estimate is the
radar's own and unconfirmed. If the authors publish the frequency or the
flight log, both close.

**What carried over unchanged.** `pr_core`'s canceller, CAF, CFAR, and a
gap-tolerant linker, on a waveform this repo cannot decode (DVB-T COFDM),
at 37.5 m cells instead of 42 m, against a 3.7 kg drone at a few hundred
metres instead of an airliner at 50 km. The reference-rebuilding half of
the repo was not exercised: the recording carries its own.
