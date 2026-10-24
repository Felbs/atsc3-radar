# atsc3-radar

*(was radarTuna, and before that passiveTuna; the docs keep the old names where they are quoted.)*

**Receive-only. The TV towers do the transmitting.**

Passive bistatic radar for a single SDR: a local ATSC 3.0 broadcast tower
floods the sky; aircraft reflect it; an RSPdx and a synthesized copy of the
tower's own waveform turn those echoes into range and Doppler. No transmitter,
no licence, one receiver.

Status (2026-10-24): the confirming run happened. On the yagi, with
ADS-B truth logged alongside, four aircraft were tracked against their own
telemetry across windows (one across all three, 6.4 -> 2.8 km over 88 s),
and an airliner at 15 km excess range was held for 55 dwells with 78 m /
1.6 Hz mean error. By the pre-registered far-range criterion (R_b >= 5 km,
nulls on the same cells) that is **DETECTION-GRADE in all three windows**;
on the whole field, crowded close-range cells let the +-60 s nulls veto it.
A gain-pair control is still owed before "confirmed" goes unqualified.
See `docs/TRUTH_RUN_1024.md`.

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
