# broadcast-radar

**Receive-only. The broadcasters do the transmitting.**

Passive bistatic radar from one SDR: a broadcast tower floods the sky,
aircraft reflect it, and an RSPdx turns those echoes into range and Doppler
with no second channel, because the reference is rebuilt from the broadcast
itself. Three illuminators so far, same receiver, same code, same ADS-B truth:

| illuminator | reference | range cell | Doppler cell | what it is good at |
|---|---|---|---|---|
| **ATSC 3.0** (OFDM, 6 MHz) | synthesized pilot lattice | 43 m | ~6 Hz | *where*: tens of metres |
| **ATSC 1.0** (8-VSB, 5.4 MHz) | symbols decided, re-modulated | 42 m | ~6 Hz | *where*, more towers, VHF too |
| **FM** (analog, 200 kHz) | discriminated and re-modulated | 1.2 km | 0.6 Hz | *how far and how fast*: 97 km, 1 Hz |

No transmitter, no licence, one receiver. Everything a claim rests on
ships with its control, and the observer position is never written here.

Status (2026-10-24, night): the confirming run, the gain-pair control, a
first pass at non-ADS-B targets, and an FM comparison all happened in one
day. ATSC 3.0 off WNUV: eleven aircraft tracked against their own ADS-B
across windows, 17-85 m range error, no cell survived a gain change. FM
off a 50 kW broadcast station on the same receiver: eleven aircraft, one
tracked for 116 s from 26 to 44 km, matches at 77, 79 and 97 km, 1 Hz
Doppler error, detection-grade under the strictest rule. FM sees farther
and reads speed; ATSC 3.0 sees where. See `docs/TRUTH_RUN_1024.md`,
`docs/GAIN_PAIR_1024.md`, `docs/UNMATCHED_1024.md`, `docs/FM_VS_ATSC3_1024.md`.

ATSC 1.0 joined on 2026-10-25: WETA RF31 locked on every block of its
first capture (equalized SNR 29.7 dB, zero sync errors) and its tracks
walk at the rate their own Doppler predicts to 1-3 %; the ADS-B check
waits on the day's adsb.lol history. A sweep of every broadcast channel
the RSPdx can reach (three FM spans, VHF and UHF TV, WNUV) is in
`tools/sweep_capture.py` / `tools/sweep_run.py`.

The goal that counts: aircraft that are *not* on ADS-B.

## Read in this order
1. `docs/FEASIBILITY.md` — the link budget, and why dynamic range (not noise) is the wall
2. `docs/FIRST_LOOK.md` — the first real air, every control, and the pilot-comb kill
3. `docs/REENCODER.md` — the decode-remodulate reference: scoping + Step 0
4. `docs/CANDIDATE_TRACKS.md` and `docs/TRUTH_RUN_0820.md` — candidates, truth runs, nulls
5. `docs/TRUTH_RUN_1024.md` — the confirming run: four aircraft tracked against their own ADS-B, the rule's flaws
6. `docs/GAIN_PAIR_1024.md` — the gain-pair control: no front-end artifacts, why clock-shift nulls fail in approach traffic, the revised rule
7. `docs/UNMATCHED_1024.md` — first pass at the non-ADS-B goal: 31 of 32 far tracks explained, one candidate held pending the tower-end truth
8. `docs/FM_VS_ATSC3_1024.md` — FM broadcast as the illuminator on the same receiver: farther, 1 Hz Doppler, 1.2 km cells; what HD Radio does to cancellation
9. `docs/FM_GAIN_PAIR_1024.md` — the FM gain pair: no front-end products, one aircraft on all 239 dwells, why "fixed cells" on FM are consecutive arrivals
10. `docs/SWEEP_1025.md` — every broadcast channel in one pass: 6 of 11 8-VSB channels and 32 of 59 FM stations lock, the five bugs the sweep found, and a position from three FM towers at once
11. `docs/EXT_KIOS.md` — the first outside recording (University of Cyprus, DVB-T, a drone): 8 MS/s found from the signal, four tracks of one target, its carrier recovered from the tracks to within 4 %

## Run it where you live

Nothing in the processing is tuned to this site: 8-VSB is the same code in
every North American transmitter, FM is FM, and the canceller, CFAR and
tracker adapt to whatever they are fed. What is site-specific lives in one
gitignored file. The recipe, tested here on an RSPdx; other 12-14-bit
SoapySDR radios should work and are unproven:

1. **Site file.** `cp site.example.json lab/site.json`; put the receiver's
   latitude, longitude and altitude in `observer`, the SoapySDR driver and
   antenna port names in `sdr`, and (if your default python lacks SoapySDR)
   the path of one that has it in `python`. `python tools/radar_site.py`
   shows what it read, position withheld.
2. **Stations.** `python tools/make_stations.py --radius-km 80` writes a
   station file for every licensed FM and full-power TV station within
   80 km, from the FCC's public database (tower position, radiation-centre
   height, ERP). The query centre is rounded to 0.1 degree so the FCC
   learns your city, not your house. Add `--atsc3 N` for channels you know
   run ATSC 3.0 (rabbitears.info); the FCC does not say. US and territories
   only; elsewhere, write the files by hand in the same shape.
3. **Check.** `python tools/site_check.py --capture 99.5e6` (any strong FM
   frequency): packages, radios SoapySDR can see, station counts, and a
   2 s capture with its headroom and integrity verdict.
4. **Capture.** `python tools/sweep_capture.py --only fm` records three
   8 MHz FM spans (120 s each); `--only uhf` / `--only vhf` the TV
   channels; or one channel with `tools/iq_capture.py`. The sweep's TV legs
   are this site's channel list; edit `legs()` for yours (`--rf` picks from
   it). Captures are 3-4 GB each and never leave `lab/`.
5. **Run.** `python tools/sweep_run.py lab/sweep_<stamp>` processes every
   leg for its signal type and writes a table: cancellation, lock, tracks,
   and each longest track's range walk against its own Doppler. A locked
   station cancels 20 dB or more; a track whose walk matches its Doppler
   within a few percent is an aircraft whether or not anything confirms it.
6. **Truth.** The day after: `python tools/adsb_history.py fetch <day>`
   (needs `gh`), then `truth <meta.json ...>` and
   `sweep_run.py <dir> --truth --skip-done`. OpenSky's live API works too
   (`truth_capture.py`) within its anonymous quota.

What to expect, and what will not work:

- **Range to the towers.** 20-100 km worked here. Under about 8 dB sync
  SNR an 8-VSB channel cannot be decoded (the Baltimore channels, 100 km
  through trees); within a few km of a megawatt transmitter the front end
  may run out of dynamic range. Untested.
- **8-bit radios** (RTL-SDR) cannot capture a 6 MHz TV channel at all and
  may cap FM cancellation well under the 25-30 dB seen at 14 bits.
- **Antenna farms.** The three-tower position fix (`fuse_fix.py`) needs
  three stations on different towers in one 8 MHz span. A city whose
  stations share one mast gets ranges, Dopplers and tracks, not positions.
- **Outside North America** there is no 8-VSB or ATSC 3.0; FM should work,
  but 100 kHz European channel spacing and the absence of HD Radio are
  untested against the FM tools.
- The thresholds (lock gates, the HD-sideband test, the trellis switch
  point, the field-comb rule) were set on one site's data. They are
  physically motivated and n = 1. A second site is the test that counts.

## Tools
`pr_core.py` CAF/ECA/CFAR · `pr_pilots.py` pilot-lattice synthesis ·
`pilot_cancel_test.py` flattening experiment · `track_hunt.py` detections → tracks ·
`truth_capture.py` capture + ADS-B truth (`--rfgain 2,3` = a gain pair under one poll) ·
`truth_compare.py` truth vs detections (shift + scramble nulls, far-range verdict) ·
`truth_trajectory.py` matched detections per aircraft across windows · `truth_figure.py` the picture ·
`gain_pair.py` fixed-cell artifact test across two gains ·
`unmatched_tracks.py` the matcher inverted: tracks no ADS-B aircraft explains (ghost + truth-box rules) ·
`fm_radar.py` the FM illuminator (reconstructed reference, 1 s dwells) · `iq_capture.py` capture at any centre ·
`fm_survey.py` rank the band · `stations/*.json` illuminator geometry for `--station=` ·
`pr_panel.py` the range-Doppler scope (localhost :8650) ·
`two_channel_radar.py` passive radar on a recorded reference (any coherent two-channel SDR, any waveform) · `dvbt_decode.py` DVB-T settings (TPS), MER and transport stream via gr-dtv · `dvbt_pilot_ref.py` a DVB-T reference from the pilots alone, one channel, no decode · `radar_site.py` the site file · `make_stations.py` station files from the FCC database · `site_check.py` readiness

Tests: `test_pr_core.py` (synthetic truth, 11 gates) · `test_flatten_preserves.py`
(flattening must not eat echoes).

House rules: every claim ships with its control; the observer position is never
written to this repo; captures and truth files live in gitignored `lab/`.
