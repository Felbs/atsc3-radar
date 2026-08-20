# The truth run — 2026-08-20 20:37Z

Felbs: "do the fresh capture with the ADS-B truth." This is the
experiment that can turn CANDIDATE_TRACKS.md's candidates into a
measurement — or into a clean negative. Either is a result.

## What was recorded (both halves, simultaneously)

- **Air:** RF25 / WNUV, 6.912 MS/s, **Antenna A (the Philips)**, rfgain_sel
  4 (readback), IFGR 32, bias-T off (readback). 120.0 s, integrity OK
  (829,483,200 samples in 119.99 s wall), no clipping, −18.2 dBFS headroom,
  0 overflows, bootstrap peak ratio 89.3. `lab/truth_rf25_20260820_203730Z.cs16`
- **Truth:** OpenSky `/states/all` every 10 s over a ±0.6°/±0.75° box
  around the receiver, from 20 s before to 20 s after the capture.
  17 polls, 908 state vectors, 0 errors, **38 distinct aircraft**.
- **Lock:** priority 50, held 20:37:50–20:39:52Z, released clean;
  sid_watch skipped those minutes and resumed.

## Geometry facts established

- Illuminator: **WNUV**, Baltimore's ATSC 3.0 lighthouse, Television Hill,
  **39.336222°N 76.649417°W** (public), ERP **750 kW** (my feasibility
  budget assumed 500 — conservative by 1.8 dB), HAAT 372.8 m. Antenna
  height above sea level is an **estimate** (400 m) — a 50 m error moves a
  prediction by well under one range bin at these ranges.
- **Tower–receiver baseline: 96.96 km.** Farther than any case in the
  feasibility table, which *helps*: the direct path arrives weaker, so a
  target near the receiver has a better echo-to-direct ratio than the
  table's 5–30 km cases.
- Our CAF's delay axis is relative to the direct path, so the observable
  is the **excess** path R_b = |T−P| + |P−R| − |T−R|, not the sum path.

## Pre-registered, BEFORE looking at a single detection

Written 20:52Z with `track_hunt.py` still running on the capture.

**Power check (computed from truth alone):** 19 of the 38 aircraft have
predicted (R_b < 60 km, |f_d| < 1 kHz) inside our map for the whole
track-hunt window (t = 30–43 s). So a null result here would be
*informative*, not a geometry accident. Three of the 19 sit inside the
±15 Hz clutter notch (AAL2839, UAL1678, UAL988 — all nearly tangential)
and are **expected invisible by construction**; the test is on the other
16.

Strongest predicted signatures (excess range, bistatic Doppler):

| aircraft | R_b | f_d | note |
|---|---|---|---|
| AAL2175 | 41.2–44.2 km | +397..+438 Hz | fast, far, closing |
| SWA4047 | 51.2–55.6 km | −614..−606 Hz | fast, far, opening |
| UAL1338 | 47.2–50.4 km | +440..+447 Hz | |
| UAL983 | 22.4–23.5 km | +152..+163 Hz | |
| N517KM | 10.8–12.2 km | +187..+201 Hz | low, close |
| N7CE97 | 9.9–10.8 km | +126..+130 Hz | |
| ASH4072 | 10.6–11.1 km | +75..+77 Hz | |

**Predictions:**
1. If the candidate-track machinery is seeing aircraft, detections will
   match predicted (R_b, f_d) for *some* of these 16 within 3 range bins
   (+ timing slop) and 30 Hz, at the same dwell.
2. The two null controls — truth time-reversed, and truth shifted +60 s —
   keep every aircraft's geometry and speed but break the clock. **They
   must score near zero.** A match count that survives either null is
   coincidence and is reported as such.
3. Verdict rule, fixed now: DETECTION-GRADE only if real matches ≥ 3 AND
   real > 3 × max(null). Anything weaker is "not confirmed."
4. Expect MISSES too — small GA aircraft (σ ≈ 5 m²) at 30+ km may simply
   be below the residue floor, and anything crossing the baseline is in
   the notch. Predicted-but-missed is recorded per aircraft, not hidden.

### Result — Antenna A (Philips): a clean NEGATIVE

```
detections (off-comb, flattened)   49 over 54 dwells = 0.91/dwell
  SNR  min 11.5  median 11.7  max 12.7 dB
  expected CFAR false alarms        38.3
  (VOID disconnected-antenna run    ~1.2/dwell at 11.6-12.4 dB)
tracks >= 3 frames                  0
MATCHES to ADS-B (same dwell)       0 of 49
PREDICTED-BUT-MISSED aircraft      19 of 19
NULL time-reversed                  1
NULL shifted +60 s                  0
VERDICT (pre-registered rule)       NO matches
```

Read it straight: **on the Philips, the radar saw nothing but its own
false-alarm floor while 19 aircraft were demonstrably in the map.** The
detection field is statistically indistinguishable from the disconnected-
antenna control — slightly more detections than the nominal 38 (49), all
at the floor, no persistence, no kinematic structure. The detector is
behaving exactly as a noise floor should; there was simply no echo above
it to find.

**What this does and does not say about the 8/10 candidate tracks:**

- It does **not** confirm them. It could not: the 8/10 capture was taken
  on the roof yagi (port B in that era, with the splitter and LNA); today's
  was on the **indoor Philips on port A**, and the yagi is currently not
  connected to anything (B carries the GPS antenna as of the 8/20 map).
  Different antenna, different sky view, different clutter.
- It does **not** refute them either. The 8/10 tracks reached 22 dB with
  persistent motion and survived a permutation control at p ≤ 0.024;
  today's field never rose above 12.7 dB. Those are different regimes, not
  the same regime with a different answer.
- It **does** establish the instrument's behaviour on a known-busy sky:
  no false positives dressed as aircraft. A radar that reports nothing
  when it cannot see is worth more than one that reports something.

**A hypothesis worth writing down, flagged as speculation:** the two slow,
close 8/10 tracks (2.5–3.1 km, +55–60 m/s of *sum-path* rate) are
consistent with **road traffic** — a car at ~30 m/s on a nearby highway
gives roughly that sum-path rate at that range, and cars are strong
reflectors at 2–3 km. ADS-B truth can never validate a car. The two fast,
far tracks (9.5–11 km, ~180 m/s) remain aircraft-shaped. Untested.

**Next, tonight:** same experiment on **Antenna C (the discone)** —
omnidirectional, real sky view, the closest thing to the yagi's geometry
still plugged in. Launched 21:0xZ.

### Result — Antenna C (discone): INCONCLUSIVE, by instrument failure

Capture USABLE (120.0 s, integrity OK, no clipping) but the illuminator is
~20 dB weaker on the discone at 97 km: bootstrap peak ratio **11.9** vs
89.3 on the Philips, rms 109 vs 1090. And that broke the instrument:

```
LS-fit-one-H residual (pilot alignment)   0.1 dB   <- FAILED
  (Philips run 22 dB; 8/10 yagi capture 21 dB; REENCODER gate says >>20)
detections                                 58 over 54 dwells, 11.5-12.7 dB
tracks                                      0
matches to ADS-B                            2 of 58   (GJS4454, ASH6159, one dwell each)
NULL time-reversed / shifted +60 s          0 / 1
VERDICT (pre-registered rule)              NOT confirmed - real does not beat the nulls
```

**Read this one differently from the Philips run.** The Philips run was a
*valid test* that returned a negative: the pilot fit locked (22 dB),
flattening ran, and the map was clean while aircraft were present. The
discone run is **not a valid test at all**: the pilot fit never locked
(0.1 dB means the synthetic lattice found no coherent pilots to fit), so
flattening was a no-op and the "detections" are the raw self-referenced
residue floor. Its two single-dwell "matches" at the floor, against nulls
of 0 and 1, are exactly the chance rate the pre-registered rule exists to
reject — and it rejected them.

Instrument lesson, now a gate: a dwell whose pilot fit is below ~10 dB is
VOID for radar purposes and must be reported as such, not processed. The
illuminator has to be audible before the echoes can be.

## Where the truth run leaves things, 21:30Z

| antenna | illuminator | pilot lock | test valid? | verdict |
|---|---|---|---|---|
| A — Philips (indoor TV) | strong (89) | 22 dB | **yes** | **clean negative**, 19 aircraft in map |
| C — discone (roof) | weak (12) | 0.1 dB | **no** | inconclusive |
| B — yagi (roof, 8/10's antenna) | — | — | — | **not plugged in** |

The 8/10 candidates got no trial tonight, because the antenna that
produced them is disconnected. Their fair trial is one capture on the
yagi with this exact recipe. Until then: four kinematically consistent
tracks on one antenna, a clean negative on another, and an instrument
that has now been shown to say "nothing" when there is nothing to see.
