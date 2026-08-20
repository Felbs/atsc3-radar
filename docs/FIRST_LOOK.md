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
   - Doppler +199.2/+195.3 Hz (adjacent 3.9 Hz bins) brackets
     **198.53 Hz = the symbol rate of 32K FFT + 2048 CP** at 6.912 MS/s.
   - Range 7.37 km = **169.9 samples ≈ 32768/192 = N/Dx**, the classic
     scattered-pilot delay ambiguity (pilot spacing Dx = 192).

   The CAF found the OFDM pilot lattice. It is exquisitely stable, it is
   strong, and it is painted on every ATSC 3.0 frame ever transmitted.

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
| VOID capture control | pre-registered, pending |
| gain-pair control (g2 vs g4) | pending |
| ADS-B truth on a fresh capture | pending |
| decode-remodulate reference | **the build that unlocks the rest** |
