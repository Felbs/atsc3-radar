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
10. `docs/SWEEP_1025.md` — every broadcast channel in one pass: 6 of 11 8-VSB channels and 32 of 59 FM stations lock, and the five bugs the sweep found

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
`pr_panel.py` the range-Doppler scope (localhost :8650)

Tests: `test_pr_core.py` (synthetic truth, 11 gates) · `test_flatten_preserves.py`
(flattening must not eat echoes).

House rules: every claim ships with its control; the observer position is never
written to this repo; captures and truth files live in gitignored `lab/`.
