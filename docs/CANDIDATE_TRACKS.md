# Two candidate tracks — 2026-08-20 19:45Z

**Status: CANDIDATES. Not aircraft. Read the whole file before quoting a
number from it.**

After pilot-flattening removed the transmitter's own comb (28.1 → 5.4 dB
over 20 frames), detections per dwell fell from 17.8 to 2.1. Two of the
survivors move in a way an artifact has no reason to.

## What was found

20 consecutive 150 ms dwells, `rf25_fox_g4_0810_2254.cs16`, t = 30.07 →
34.68 s (8/10 evening).

**Track A** — frames 0–8, 1.940 s

```
t=30.0675  r=2429 m  fd=-85.9 Hz     t=31.2801  r=2472 m  fd=-89.8 Hz
t=30.3100  r=2429 m  fd=-85.9 Hz     t=31.5226  r=2472 m  fd=-89.8 Hz
t=30.5525  r=2429 m  fd=-85.9 Hz     t=31.7651  r=2516 m  fd=-89.8 Hz
t=30.7950  r=2472 m  fd=-89.8 Hz     t=32.0076  r=2516 m  fd=-89.8 Hz
t=31.0375  r=2472 m  fd=-89.8 Hz
```
observed range rate **+44.8 m/s**, predicted from Doppler **+49.2 m/s**
— agree to 9 %. Total walk 87 m = 2.0 range bins.

**Track B** — frames 14–17, 0.728 s

```
t=33.4627  r=11147 m  fd=+335.9 Hz
t=33.7052  r=11103 m  fd=+332.0 Hz
t=33.9478  r=11060 m  fd=+332.0 Hz
t=34.1903  r=11017 m  fd=+332.0 Hz
```
observed range rate **−178.7 m/s**, predicted from Doppler **−185.2 m/s**
— agree to 4 %. Total walk 130 m = 3.0 range bins, monotonic.

## Why this is worth a second look

**Range and Doppler are independent measurements.** Range comes from
correlation delay; Doppler comes from phase evolution across batches.
Nothing couples them in the processor. An artifact — a sidelobe, a
residue, a numerical ridge — has no mechanism to make one predict the
other. Kinematic consistency is the standard first test a passive-radar
detection has to pass, and both tracks pass it.

**The obvious null fails.** Uncorrected clock drift (−0.790 ppm over
0.2425 s frames) would move *every* cell **+57.4 m per frame in the same
direction*. Track A moves +10.9 m/frame; Track B moves **−43.3**
m/frame. Opposite signs, both smaller than the drift. A shared anchor
error cannot produce this pair.

**They are not comb teeth.** 2.43–2.52 km and 11.0–11.15 km are not
n × 7.402 km, and ±86–336 Hz are not n × 198.53 Hz.

## Why it is NOT a detection yet

1. **Track lengths are short.** 9 and 4 frames. Track A's "motion" is
   two 43 m bin-steps; with quantization that coarse, a slow drift of
   any origin could imitate it. Track B's 3 monotonic steps are stronger
   but still only 0.73 s.
2. **No ADS-B truth exists for 8/10.** The capture predates any thought
   of radar, so there is no independent record of what was in the sky.
   This is the single biggest gap and it cannot be closed retroactively —
   only by a fresh capture with OpenSky state vectors logged alongside.
3. **The self-referenced CAF is still structurally compromised.**
   `comb_report` runs pr_banked's v0 pipeline: ECA(x, x), which
   annihilates. Flattening removed the dominant comb but the underlying
   reference problem is unfixed until decode-remodulate lands.
4. **Track A sits close to the direct path** (2.4 km), where residue
   from imperfect cancellation is strongest.
5. **The 19.4 dB pilot-suppression floor is unexplained.** Whatever sets
   it also sets the residue field these tracks live in.
6. **n = 2.** Two tracks in 20 dwells, found by looking at what survived.
   That is hypothesis generation, not confirmation.

## Pre-registered test (written before running it)

Extend to **60 frames** (t = 30 → ~44.6 s). If Track B is a real target
at −185 m/s sum-path rate, it must walk **2.7 km ≈ 62 range bins** across
that window and appear as a *continuous* track, its Doppler evolving
smoothly as the bistatic geometry changes. If it is residue, it will not
extend — it will sit, jump, or vanish.

**Prediction:** a genuine target continues monotonically; an artifact
does not. Result appended below when the run lands.

### Result — the 60-frame run REFUSED TO EXECUTE, and it was right to

`find_frames` measured a frame spacing of 2621884 samples against the L1's
1676288 and raised rather than proceed. Diagnosis: that spacing is exactly
what you get from **40 detected bootstraps where 61 were expected** — the
bootstrap correlator misses ~1/3 of frames over a 14.8 s load, and
`mean(diff)` inflates accordingly. Every frame would have been
mis-anchored. The gate caught it.

So the test ran instead as **3 batches of 18 frames**, each independently
anchored (spacing error −0.0002 symbols per batch), stitched on absolute
sample index — `tools/track_hunt.py`. That tool also fixes a hole in the
evidence above: `comb_report` stored only the STRONGEST peak per dwell, so
a track that dipped below the per-frame maximum was indistinguishable from
a track that ended.

## 54 frames, every off-comb detection kept

145 off-comb detections → 11 tracks of ≥3 frames. Four are kinematically
consistent (agreement < 15 % AND ≥ 2 range bins of walk):

| frames | span | range | mean f_d | observed | predicted | off | walk |
|---|---|---|---|---|---|---|---|
| **17** | 4.365 s | 2.819 → 3.079 km | −106.2 Hz | +59.6 m/s | +59.0 | **1.0 %** | 6 bins |
| **11** | 3.153 s | 2.559 → 2.732 km | −97.7 Hz | +55.0 m/s | +54.3 | **1.3 %** | 4 bins |
| **9** | 2.183 s | 10.713 → 10.323 km | +321.2 Hz | −178.8 m/s | −178.6 | **0.1 %** | 9 bins |
| **4** | 0.728 s | 11.147 → 11.017 km | +333.0 Hz | −178.8 m/s | −185.2 | 3.5 % | 3 bins |

And the test **discriminates** — it rejects as well as accepts. Tracks at
2.516 km and 6.853 km show **0.0 bins of walk against a Doppler demanding
motion** (100 % off) and are correctly thrown out. A test that only ever
says yes is not a test.

## Two self-checks, because the first version of this was nearly wrong

**1. Circularity — the association rule used Doppler to predict range.**
That would partly manufacture the very consistency being reported.
Re-ran association **Doppler-blind** (pure range proximity, ±4 bins, no
kinematic hint): **byte-identical output.** Every track, every number.
The tracks form on range proximity alone; the kinematic agreement is an
emergent property, not an imposed one.

**2. Permutation control — could 145 scattered detections do this by
chance?** Same detections, same range and Doppler values, **time labels
shuffled**, 40 trials:

```
tracks >= 3 frames    real  11    shuffled mean  8.8  (max 12)
CONSISTENT tracks     real   4    shuffled mean  0.12 (max 2)
                                  0 of 40 shuffles reached 4
p <= 0.024
```

Shuffling produces about as many *tracks* — association is easy — but it
almost never produces a **consistent** one. The structure lives in the
time ordering, which is exactly where real motion lives and where an
artifact's does not.

## Where this now stands

Four candidate tracks survive: two slow and close (2.5–3.1 km, ~55–60 m/s
opening) and two fast and further (9.5–11.1 km, ~180 m/s closing, and the
9-frame one walks **9 range bins** while agreeing with its Doppler to
0.1 %). Those speeds and ranges are aircraft-plausible.

**It is still not a detection**, and the reasons have not changed:
no ADS-B truth exists for an 8/10 capture; the self-referenced CAF is
still structurally compromised until decode-remodulate lands; the 19.4 dB
pilot-suppression floor is unexplained; and n = 4 found by looking at
what survived is hypothesis generation.

**The one experiment that would settle it** is unchanged and now clearly
worth doing: a fresh capture with OpenSky state vectors logged during the
recording window, then predict each known aircraft's bistatic delay and
Doppler from its own telemetry and check whether our peaks land where
geometry says they must. That needs the radio and roughly ten minutes.

## If it does extend

Still not a claim. The next steps in order:
1. Re-run with the same dwells on **g2** (different gain, 3 min later) —
   a real aircraft will be *somewhere else*, not absent; a gain-dependent
   artifact will change level.
2. **Fresh capture with OpenSky truth** logged during the recording
   window. Predict each known aircraft's bistatic delay and Doppler from
   its state vector and check whether our peaks land where geometry says
   they must. That is the test that would make it a detection.
3. Only then does the word "aircraft" belong in a sentence about this.
