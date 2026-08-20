# First look at real air — 2026-08-20

20 × 150 ms dwells, `rf25_fox_g4_0810_2254.cs16` (8/10 evening, gain 4),
t = 30–33 s. 361 CFAR candidates. **None of them are aircraft, and proving
that is the day's second-best result.**

## The three tells, in the order they fired

1. **ECA reported 175.2 dB of cancellation — identically, every dwell.**
   Physical cancellation is bounded by reference SNR; 175 dB is not a
   measurement, it is annihilation. Cause: the v0 driver passed the SAME
   channel as reference and surveillance, so the lag-0 basis column
   contains the entire signal — echoes included — and least squares nulls
   everything. The residue is numerical, and its structure is the
   signal's own.
2. **The best peak sat at exactly 7.37 km for 3.0 seconds while claiming
   ~±197 Hz of Doppler.** 197 Hz at 55.6 cm is ~110 m/s of sum-path rate
   = ~7.6 range bins of walk over those 3 s. A target that reports motion
   and does not move is not a target.
3. **The numbers match the waveform, not the sky:**
   - Doppler +199.2/+195.3 Hz (adjacent 3.9 Hz bins) brackets 198.53 Hz.
   - Range 7.37 km = 169.9 samples ≈ 170.7.

   The CAF found the OFDM pilot lattice. It is exquisitely stable, it is
   strong, and it is painted on every ATSC 3.0 frame ever transmitted.

   **⚠ Identification corrected, same day.** The first version of this
   note attributed the numbers to "32K FFT + CP 2048" (symbol rate
   198.53 Hz) and "N/Dx = 32768/192". The re-encoder scoping study pulled
   the banked RF25 L1 decode (`lab/m8_l1_rf25.json`) and the actual
   config is **8K FFT, GI 512, SP12_4** — and 192 is not even a legal Dx.
   The correct derivation: SP12_4's within-symbol pilot spacing is
   **Dx·Dy = 48 carriers → 8192/48 = 170.67 samples = 7.40 km**, and the
   scattered lattice repeats every **Dy = 4 symbols → 794.1/4 =
   198.53 Hz**. The arithmetic is identical (4 × 8704 = 34816 samples
   either way), which is exactly how a wrong story survives a numeric
   check — the L1 decode, not the numerology, is the authority. The
   conclusion (pilot-lattice artifact, signal-borne) stands unchanged;
   the mechanism is now the right one.

## Pre-registered prediction (written BEFORE reading the control)

`e65_VOID_antenna_disconnected_rf33.cs16` — same pipeline, a capture with
the antenna physically disconnected — contains no ATSC signal. Therefore:
**the 7.4 km / ±199 Hz structure must be ABSENT from the VOID run.** If it
appears there too, the artifact is processor-borne and pr_core has a bug;
if the VOID is clean of it, the artifact is signal-borne as claimed.

### VOID result — prediction CONFIRMED

Same pipeline, 6 dwells over the disconnected-antenna capture:

- **The 7.37 km / ±199 Hz structure is absent.** Not attenuated — absent.
  The pilot-lattice artifact is signal-borne, exactly as predicted, and
  the processor is exonerated.
- 7 residual detections, all weak (11.6–12.4 dB vs the g4 run's 25–28),
  each in a DIFFERENT range/Doppler cell each dwell — twinkle, not
  motion. Expected CFAR false alarms over these 6 dwells:
  6 × 1384 × 512 cells × pfa 1e-6 ≈ **4.3**; observed 7. The detector is
  running at roughly its nominal false-alarm rate even on the residue
  field (the post-ECA correlation penalty the builder warned about is
  mild here because the "signal" is already near-white noise).

So the scoreboard's two live rows both closed today, in the right
direction: the g4 air shows strong stable structure that is provably the
waveform; the VOID shows nothing but nominal noise. The instrument
behaves. What it lacks is a reference clean enough to see past the
transmitter — the decode-remodulate build.

## What this changes

**The re-encoded reference is not v2 polish. It is the entry ticket.**
Single-channel PBR with the raw signal as its own reference is
structurally void (tell #1). With one RSPdx and no second coherent
channel, the only clean reference is the one we uniquely can build: run
the capture through our own ATSC 3.0 decoder, take the decided
constellation points and known pilots, and re-synthesize the ideal
waveform. Noise-free reference, perfect cancellation basis, and the
pilot lattice becomes *known* structure we subtract rather than detect.

Interim mitigation (cheap, imperfect): mask the known ambiguity ridges —
symbol-rate Doppler harmonics and N/Dx delay multiples — and see what
survives. Anything that walks coherently through range OUTSIDE the mask
is worth a second look. But masking hides, it does not cancel; the
decode-remodulate reference is the real instrument.

## Standing scoreboard

| control | status |
|---|---|
| synthetic truth recovered (range+Doppler) | ✅ 11/11 |
| −60 dB echo masked raw, found after ECA | ✅ |
| pure noise → zero detections | ✅ |
| self-reference annihilation identified | ✅ (this run) |
| pilot-lattice artifacts identified to 4 figures | ✅ (this run) |
| VOID capture control | ✅ pre-registered, then confirmed (below) |
| gain-pair control (g2 vs g4) | ✅ comb gain-invariant to 0.18 dB (below) |
| geometry from the capture's own L1 | ✅ Step 0, see REENCODER.md |
| pilot-only cancellation kills the comb | ✅ 23.4 dB, control 0.3 dB (below) |
| flattening PRESERVES targets off-comb | ✅ 0.0 dB cost, pre-registered (below) |
| ADS-B truth on a fresh capture | pending |
| decode-remodulate reference | scoped; its cheapest part already works |

## Gain-pair control — 8/20 17:50Z

Identical processing on `rf25_fox_g2_0810_2257.cs16`, recorded three
minutes after the g4 capture at a *different front-end gain*. E65/E73
proved this rig compresses on this exact signal, and intermod
impersonates targets — so anything that moves with gain is ours, not the
sky's.

| | g4 (rfgain 4) | g2 (rfgain 2) | Δ |
|---|---|---|---|
| detections / dwell | 18.1 | 16.4 | — |
| peak SNR (mean of 20 dwells) | 27.00 dB | 27.18 dB | **−0.18 dB** |
| peak cell | 7.37 km, +199.2 Hz | 7.37 km, +199.2 Hz | identical |

**The structure is gain-invariant.** A 0.18 dB shift across a gain change
that once moved RF33 from 0.0 % to 100.0 % FEC is no shift at all. The
artifact lives in the transmitted waveform, exactly as the pilot-lattice
identification says.

**And the control did better than confirm — it partitioned the
detections.** The range cells shared by both captures are, to the bin,
the pilot ambiguity comb n × 7.402 km:

```
shared top cells   7.37  14.79  22.16  29.58  37.00  44.37  51.79  59.20
predicted n*N/DxDy 7.40  14.80  22.21  29.61  37.01  44.41  51.81  59.22
```

Eight harmonics, eight matches, all within a 43 m bin. Meanwhile every
cell that appears in only ONE capture is *non*-harmonic (g4-only: 2.43,
2.47, 2.52, 4.47 km; g2-only: 2.65, 4.38, 11.71, 22.21 km). So the
gain-invariant comb is the transmitter and the gain-dependent stragglers
are our own front end — which is precisely the discrimination the control
was built to provide, and it cost nothing because the E76 campaign banked
the pair by luck.

**Scoreboard:** gain-pair control ✅. Remaining before any target claim:
a reference clean enough to cancel the comb (the decode-remodulate
build), then ADS-B truth on a fresh capture.

## Pilot-only cancellation — 8/20 18:00Z. The comb dies; targets don't.

A pilot lattice synthesized from A/322 tables alone — **no FEC, no
constellation decisions, no LDM, nothing fitted from the data** — against
the real capture:

| | real pilots | scrambled control |
|---|---|---|
| comb tooth n=1 (7.40 km / +198.5 Hz) | 29.1 → **4.8 dB** | 27.9 → 28.7 dB |
| detections / dwell | 18.3 → **2.0** | 16.7 |
| scattered-pilot-carrier suppression | **+19.4 dB** | −0.17 dB |
| LS-fit-one-H residual (REENCODER's own gate) | **21.2 dB** | — |

**The plan reverses: flatten, don't remove.** Full pilot *removal* only
took 3.6 dB off the tooth, and that is physics. The comb is the DFT of
the power spectrum at 1/48 — it exists because pilot carriers hold a
*different* power from data carriers. Delete them and a +4.5 dB bump
becomes a −15 dB hole, which combs almost as hard in anti-phase. Subtract
only the **excess**, β = 1 − 1/A per class, and it collapses.

β = **0.4043** (scattered) and **0.6250** (common CP) come straight from
the amplitude table with nothing fitted — and they land on the null;
β×0.8 and β×1.2 both leave the tooth 5–10 dB higher. **A prediction
hitting a null it did not fit** is the strongest confirmation of the
amplitude map this project has, and it agrees with an independent
measurement from air (SP/data +4.42 dB measured vs +4.50 tabled).

Alignment, the part that was expected to be hard: bootstrap anchor →
lattice-phase test (pre-registered `l mod Dy = 0` predicted; scored 0.982
against 0.202–0.239 for the three wrong phases) → joint CFO/symbol-period
refinement. Transmitter clock **−0.799 ppm**, confirmed three independent
ways (bootstrap spacing, period fit, pilot phase-ramp slope). CFO
+431.85 Hz.

### The gate that decides whether any of it is usable

Flattening is a *structured, signal-correlated* subtraction — exactly the
shape of an operation that removes the thing you're hunting along with
the thing you're cancelling. Killing the comb says nothing about that.
So: inject a real bistatic echo (a delayed, Doppler-shifted, attenuated
copy of the illuminator) into real air and run the identical path.
Pass conditions pre-registered before the run.

| injection | raw | flattened | loss |
|---|---|---|---|
| between teeth, 11.10 km +260 Hz | 10.3 dB | **10.3 dB** | **+0.0** |
| weak, 26.00 km −430 Hz | 2.2 dB | **2.2 dB** | **−0.0** |
| ON a comb tooth, 7.40 km +198 Hz | 18.7 dB | 13.6 dB | +5.1 |
| no-echo control at the probe cell | — | −4.8 dB | nothing invented |

**Off the comb, flattening is free — zero measurable cost to a target.**
On a tooth it costs 5.1 dB and the target still survives, which is the
honest expected outcome: a target sharing a cell with the structure being
cancelled must suffer, and that case is reported rather than graded.

**So the cancellation path is cleared for target work**, and the
decode-remodulate build is no longer a leap of faith — its cheapest
component already does the heavy lifting.
