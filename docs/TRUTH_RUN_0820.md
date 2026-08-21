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

---

# The yagi run — 2026-08-21 00:21Z

Felbs plugged the yagi back into B. The first launch at 22:51Z stood
down because `sonde_rx` held the radio at priority 100 for the 00Z
Sterling hunt — "let the sonde hunt finish, then run it." A waiter keyed
on both the lock file and the process fired at 00:20:51Z.

**Capture:** Antenna B, bias-T **true** (readback), rfgain 2 (RF25's E76
value on this antenna), 120.0 s, integrity exact, no clipping, −14.5 dBFS,
**bootstrap peak ratio 98.0** — the 8/10 regime. Pilot lock 17–19 dB on
every dwell: a *valid* test. **Truth:** 17 polls, 1037 states,
**61 aircraft**, 31 in-map outside the notch — the late arrival push.

## What the radar saw, window by window

Track hunt run on three windows to test the one hypothesis that could
explain a match without an aircraft: a fixed clutter source at ~3 km in
the yagi's beam (the 8/10 capture had a ~3 km track too).

| window | strongest persistent detection | UAL1213 predicted from ADS-B | fixed clutter would be |
|---|---|---|---|
| t = 5–14 s | **2.73 → 2.86 km, −20 Hz**, 28 frames | 2.73 → 2.85 km, −16..−20 Hz | 3.2 km, −30 Hz |
| t = 30–43 s | **3.08 → 3.30 km, −27..−35 Hz**, 50 dwells, 17–19 dB | 3.03 → 3.24 km, −25..−34 Hz | 3.2 km, −30 Hz |
| t = 90–99 s | **4.81 km, −78 Hz**; *nothing* at 3.2 km | 4.52 → 4.96 km, −75..−83 Hz | 3.2 km, −30 Hz |

**The detection moved with the aircraft**, 2.73 → 3.30 → 4.81 km over
92 s, within ~50–250 m and ~5 Hz of where a 10.9 km-altitude, 217 m/s
overflight's own telemetry said it would be, and the walk rate inside
each window matched its Doppler (+14–16 m/s observed vs +16–18 predicted
in the middle window). The stationary-clutter hypothesis is **falsified**
by the late window: the 3.2 km cell is empty once the aircraft has left
it.

Other matches at the same dwell (|ΔR| < 490 m, |Δf| < 30 Hz):

| aircraft | early window | middle | late | note |
|---|---|---|---|---|
| **UAL1213** | 32 dwells | **50** | 2 | the trajectory above |
| **PAT11** | 22 | — | 20 | 1.26–1.6 km; a kinematic track at 1% agreement in the early window; Doppler crossed −46 → +29 Hz as it flew past |
| FFT2269 | 11 | — | — | 2.0 km, +20 Hz |
| SWA4809 | — | — | 6 | 12.9 km, +16 Hz |
| JBU128 | — | — | 6 | **23.4 km**, +51 Hz vs +62 predicted |
| RPA5689 | — | — | 5 | 4.2 km, +66 Hz |
| SWA830 | 6 | — | — | loose (0.95 vs 0.61 km) |

Totals: real matches 71 / 51 / 40 of 96 / 94 / 64 detections.

## The verdict, by the rule written before the data — and why it says what it says

```
                         early   middle   late
real                       71      51      40
NULL time-reversed         71      51      33
NULL shifted +60 s         17       1       0
VERDICT (pre-registered)  NOT     NOT     NOT confirmed
```

**All three windows read NOT confirmed, and that label stands** — the
rule was fixed before the capture and is not being renegotiated after
it. But the rule is now understood to be broken in one specific way, and
that is a statement about the test, not about the sky:

**The time-reversal null has no power against a quasi-static target.**
Over a 9–13 s window a high-altitude overflight's predicted cell moves
~5 range bins and ~10 Hz — inside the match gate in either time
direction. Reversing the clock therefore maps almost the same
predictions onto almost the same dwells, and the null scores whatever
real scores. It scored 71 against 71, 51 against 51, 33 against 40. A
null that cannot lose cannot discriminate; it was the wrong control for
a short window, and `max(null)` in the verdict rule lets it veto
everything.

**The +60 s null is the informative one.** It keeps every aircraft's
geometry and speed and only breaks the clock by a minute — by which time
UAL1213's predicted cell has moved ~1 km. It scored **17 / 1 / 0**
against real's **71 / 51 / 40**. The early-window 17 is itself
instructive: several slow, low aircraft sat near the baseline during
those seconds (approach traffic), so close-range chance matching is not
negligible and the gate is generous there.

## Addendum, 01:1xZ — the first of the new nulls bit, and it should have

After committing the verdict above, `truth_compare` gained ±60/180/300 s
shift nulls (informational tonight, binding on the next run). Middle
window, same detections:

```
shifted -300 s   0      shifted +60 s    1
shifted -180 s   0      shifted +180 s   0
shifted  -60 s  46      shifted +300 s   0
real            51
```

**−60 s scores 46.** Sixty seconds earlier UAL1213's predicted cell was
~2.6 km, and with a 490 m gate (3 bins + 2 s of launch-clock slop at
180 m/s) that still overlaps the 3.1–3.3 km detections at the edge — and
several slow, low approach-traffic aircraft sat near the baseline, so
close-range signatures are crowded. This is the honest size of the
chance-match problem at ≤3 km and small Doppler: **substantial.** The
"50 of 54 dwells" figure is therefore weaker than it sounds on its own.

What it does *not* touch: the ±180 s and ±300 s nulls score 0, so the
match is clock-locked at the minute scale; and the 92-second trajectory
(2.73 → 3.30 → 4.81 km, 3.2 km cell empty at the end) is not something a
crowded close-range cell can produce. That — not the per-dwell count —
is the evidence. Two consequences for the confirming protocol, added now:
the range gate drops its 2 s launch-clock slop by anchoring capture UTC
to the first bootstrap (the file's own clock), and matches at R_b < 5 km
are reported but **do not count** toward the verdict unless the
aircraft's trajectory is tracked across ≥ 2 windows.

## What can honestly be said

1. The instrument, on the yagi, produced a persistent 17–19 dB detection
   whose position, Doppler, and walk rate followed a specific ADS-B
   aircraft's predicted bistatic trajectory across three windows spanning
   92 s, and the one competing explanation that was pre-registered
   (fixed clutter at 3.2 km) is falsified by direct observation.
2. A second aircraft (PAT11) matched in 42 dwells across two windows
   including a 1 %-agreement kinematic track, and a third (JBU128)
   matched at 23 km.
3. By the pre-registered rule it is **not a confirmed detection**,
   because one of the two nulls is structurally unable to lose. That is
   the correct outcome of a flawed rule and it is recorded as such.
4. It is, as of tonight, **the first truth-consistent aircraft candidate
   this instrument has produced**, and the strongest result of the day.

## Pre-registered protocol for the confirming run

Written now, before any further capture:

- Nulls: clock shifts of **±60, ±180, ±300 s** (all must be beaten), and
  an **aircraft-identity permutation** (assign each aircraft's track to
  another aircraft's timeline) repeated 50×, reporting the fraction of
  permutations reaching the real match count.
- Verdict: DETECTION-GRADE iff real ≥ 5 matches on ≥ 2 distinct
  aircraft, real > 3 × every shift null, and permutation p < 0.01.
- Trajectory test is mandatory, not optional: three windows ≥ 60 s apart;
  the matched aircraft's detection must move with its prediction.
- Same recipe: yagi on B, bias-T on (readback), rfgain 2, 120 s, OpenSky
  polled ±20 s, observer coordinates never written.
