# atsc3-radar

**Receive-only. The TV towers do the transmitting.**

Passive bistatic radar for a single SDR: a local ATSC 3.0 broadcast tower
floods the sky; aircraft reflect it; an RSPdx and a synthesized copy of the
tower's own waveform turn those echoes into range and Doppler. No transmitter,
no licence, one receiver.

Status (2026-10-24, evening): the confirming run and the gain-pair control
both happened. Across four clean captures on the yagi with ADS-B truth
logged alongside, eleven aircraft were tracked against their own
telemetry across three windows each (best: 119 dwells over 97 s, walk
rate matched to 1 %), an airliner was held at 29 km with zero mean range
error, and no detection cell survived a change of RF gain, which rules
out front-end artifacts. The pre-registered verdict rule passes in some
windows and is shown to fail in others for a measured reason (consecutive
aircraft on one approach route defeat clock-shift nulls); the labels stand
as written and the rule is revised for the next run only. See
`docs/TRUTH_RUN_1024.md` and `docs/GAIN_PAIR_1024.md`.

The goal that counts: aircraft that are *not* on ADS-B.

## Read in this order
1. `docs/FEASIBILITY.md` — the link budget, and why dynamic range (not noise) is the wall
2. `docs/FIRST_LOOK.md` — the first real air, every control, and the pilot-comb kill
3. `docs/REENCODER.md` — the decode-remodulate reference: scoping + Step 0
4. `docs/CANDIDATE_TRACKS.md` and `docs/TRUTH_RUN_0820.md` — candidates, truth runs, nulls

## Tools
`pr_core.py` CAF/ECA/CFAR · `pr_pilots.py` pilot-lattice synthesis ·
`pilot_cancel_test.py` flattening experiment · `track_hunt.py` detections → tracks ·
`truth_capture.py` capture + ADS-B truth · `truth_compare.py` truth vs detections ·
`pr_panel.py` the range-Doppler scope (localhost :8650)

Tests: `test_pr_core.py` (synthetic truth, 11 gates) · `test_flatten_preserves.py`
(flattening must not eat echoes).

House rules: every claim ships with its control; the observer position is never
written to this repo; captures and truth files live in gitignored `lab/`.
