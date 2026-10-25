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

**Watching it (no).** `tools/dvbt_decode.py` runs gr-dtv's DVB-T
receiver on the file. The transmitter's settings come off the air, from
the TPS carriers, four frames in a row on both channels: 8K, guard 1/4,
16-QAM, non-hierarchical, code rate 5/6. The guard matches the symbol
period measured independently (8960 samples at 8 MS/s). gr-dtv's symbol
acquisition steps its FFT window by up to +-6 samples between symbols,
which tilts each symbol's phase across the band; the TPS reader removes
that tilt with the continual pilots (TPS coherence 0.83 / 0.91). But the
picture needs 14-15 dB of MER at 16-QAM 5/6, and the recording does not
have it: gr-dtv's equalizer gives 5.6-5.9 dB on the reference channel
(-3 to -1 dB on the surveillance channel, whose channel estimate fails),
a per-symbol pilot equalizer 6.9 dB, and the cyclic-prefix correlation,
which needs no equalizer at all, implies about 10 dB. No equalizer
recovers more than the signal has, so this file is at least 4 dB short of
a transport stream. The antennas were aimed for radar, not reception. The
chain is ready for a stronger recording: `--tps` for settings and MER,
then `--qam 16 --rate 5/6 --out x.ts`, then `ffplay x.ts`.

**gr-rxtune on the recording.** A recording has no gain left to turn, so
`tools/rxtune_dvbt.py` gives gr-rxtune the software receiver's knobs
through its `FunctionKnob` interface: which antenna (ch0, ch1, or both
aligned and summed coherently), gr-dtv's acquisition SNR assumption, and
where in the 70 s to listen (its hour-of-day axis, compressed). The dial
is gr-dtv's MER with the cliff at EN 300 744's 13.5 dB for 16-QAM 5/6;
proof of decode would have been a real transport stream. 62 cells:

| path | MER, every setting and time | |
|---|---|---|
| ch0 (reference antenna) | 6.1-6.3 dB | flat: neither the acquisition setting nor the time moves it |
| ch1 (surveillance antenna) | -2.7 to -1.0 dB | gr-dtv's channel estimate fails on this antenna |
| both, summed | 5.0-5.5 dB | equal weights let ch1's echoes in; per-carrier MRC would at best recover ch0 |

rxtune's verdict, unedited: `BELOW_CLIFF - MER 6.3 dB, 7.2 dB short of the
cliff`, best cell `NOT recommended`, advice "try another antenna port,
channel or time of day". It did what it promises: no knob is recommended
without proof of decode, and the number it gives (7.2 dB short) is the
same gap every other measurement found. The advice is right but cannot be
taken on a file; the antenna was chosen in Cyprus in 2022.

## One antenna, no decode: a DVB-T reference from the pilots

The recording's second antenna made the tracks above easy. This repo's
real claim is that one receiver is enough, so the same file was run again
with **channel 1 alone**: `tools/dvbt_pilot_ref.py` synthesizes what the
transmitter put on its pilot carriers -- values fixed by the standard,
known without decoding anything -- and that waveform is the reference for
the same channel. The TV picture needed 13.5 dB; the pilots need no
decode at all.

Getting the reference right took four corrections, each found by a
check before the radar ran:

1. **Pilot values.** The scattered-pilot pattern and the PRBS sign are
   tested against the received symbols: pilots 12 carriers apart see
   nearly the same channel, so the right values make neighbours agree.
   Right: 0.92 (ch0), 0.74 (ch1); every wrong combination 0.02-0.05.
2. **One whole carrier of frequency offset.** The cyclic prefix gives the
   carrier offset only modulo one carrier spacing (1116 Hz). Sliding the
   continual pilots found the rest: -1 carrier, -926 Hz in all (pilot
   coherence 0.08 -> 0.92). This, not the pilot code, is what broke the
   hand-written equalizer earlier.
3. **Sample-clock drift.** The B210 and the transmitter differ by about
   1.7 ppm, ~180 samples over 20 s. A timing search around one fixed
   offset lost the symbols halfway, and the reference slid against the
   channel at ~9 samples/s: every peak walked +295 m/s at 0 Hz. Timing
   is now fitted block by block, then each symbol is searched around its
   predicted place.
4. **Which peak is the direct path.** A range-Doppler map only sees
   positive range. It latched onto a path 125 samples (4.1 km) behind the
   real direct path, which the cyclic-prefix timing had put at -24
   samples. A two-sided cross-correlation finds the direct path 37 dB
   clear; it is moved to zero lag and zero Doppler. The 4.1 km path,
   nearly as strong, is most likely another transmitter of Cyprus's
   single-frequency network, and the reason gr-dtv's equalizer fails on
   this antenna.

Cancellation is then 0.1-0.4 dB, and that is the expected number: the
canceller can only remove what the reference holds, the pilots, about
15 % of the power. The data carriers stay in the channel.

**What it costs, measured on the drone itself**, its own cell, same
window, same channel, the two references swapped:

| time | recorded reference (ch0) | pilot-only (ch1 alone) |
|---|---|---|
| 14 s | 25.8 dB | 6.3 dB, wrong cell |
| 16 s | 17.0 dB | 11.1 dB, at -64 Hz |
| 18 s | 17.4 dB | 6.0 dB |

6-20 dB weaker, more than the 7.5 dB the pilots' share of power alone
predicts. At 0.5 s dwells and a 1e-6 false-alarm rate: 1 detection in
20 s. At 1 s dwells (+3 dB; the drone crosses about one 33 m cell per
second) and 1e-5:

| | two antennas, recorded reference | one antenna, pilot-only |
|---|---|---|
| outbound pass | 13.2-19.8 s, 262 -> 487 m, -62.0 Hz | 17.0-19.5 s, 295 -> 393 m, -63.8 Hz |
| inbound pass | 50.8-56.8 s, 675 -> 450 m, +73.3 Hz | 50.5-56.5 s, 623 -> 393 m, +73.3 Hz |
| wavelength from walk / Doppler | 0.584, 0.563 m | 0.61, 0.49 m |
| short passes (0.2-2.2 s, 67.8-69.2 s) | found | not found |

Both long passes come back from one antenna, in the right place, moving
the right way, the inbound one at the same Doppler to 0.1 Hz. The two
short ones are lost. Ranges differ by a cell or two because each version
measures from its own direct path.

**Unexplained.** The pilot-only runs, and only they, show a cell near
1574 m at +32-34 Hz in three separate windows (3.5-9.5, 13-14.5 and
66-68 s) that does not move in range: a real target at +33 Hz would walk
four cells in six seconds. It is not one of the reference's own
ambiguities (those sit at 22.4 km and 223 Hz). Not claimed as anything.

So the single-receiver method works on DVB-T too, checked against an
independent reference on the same data -- at a price: on this weak
recording the pilots alone give up the short passes and need twice the
dwell.
