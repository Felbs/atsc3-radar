# atsc3-radar

**Receive-only. The TV towers do the transmitting.**

Passive bistatic radar for a single SDR: a local ATSC 3.0 broadcast tower
floods the sky; aircraft reflect it; an RSPdx and a synthesized copy of the
tower's own waveform turn those echoes into range and Doppler. No transmitter,
no licence, one receiver.

Status (2026-10-24, night): the confirming run, the gain-pair control, a
first pass at non-ADS-B targets, and an FM comparison all happened in one
day. ATSC 3.0 off WNUV: eleven aircraft tracked against their own ADS-B
across windows, 17-85 m range error, no cell survived a gain change. FM
off a 50 kW broadcast station on the same receiver: eleven aircraft, one
tracked for 116 s from 26 to 44 km, matches at 77, 79 and 97 km, 1 Hz
Doppler error, detection-grade under the strictest rule. FM sees farther
and reads speed; ATSC 3.0 sees where. See `docs/TRUTH_RUN_1024.md`,
`docs/GAIN_PAIR_1024.md`, `docs/UNMATCHED_1024.md`, `docs/FM_VS_ATSC3_1024.md`.

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
