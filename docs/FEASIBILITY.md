# passiveTuna — feasibility, from our own numbers

**2026-08-20.** Bistatic passive radar using local ATSC 3.0 transmitters as
illuminators of opportunity, an RSPdx as the receiver, and aeroTuna's ADS-B
as ground truth.

Everything below is computed for **our** carrier (RF25, 539 MHz, 6.912 MS/s
captures already on disk) and **our** receiver (RSPdx + discone on Antenna C,
~0 dBi). Assumptions are stated so they can be attacked.

## The headline: thermal noise is not the wall. Dynamic range is.

Bistatic radar equation, EIRP 500 kW (57 dBW), airliner σ = 100 m² (20 dBsm),
Gr = 0 dBi, NF 5 dB, coherent integration T = 0.5 s (processing gain
B·T = 65.4 dB).

> **⚠ CORRECTED 8/20, same day, by the DSP build.** The T = 0.5 s in this
> table is **not achievable for a fast target** — see *Range walk* below.
> A 250 m/s sum-path target walks 2.9 range bins in 0.5 s and smears
> across them, so the coherent gain never lands. Realistic dwell is
> **~170 ms** (or ~87 ms for the β→0 case), costing **4.6–7.6 dB**.
> Subtract that from every "after CAF" figure below. The margins are so
> large that it changes nothing about feasibility — but the number was
> wrong and is now on the record as wrong.

| Rt (km) | Rr (km) | Echo power | SNR in 6 MHz | **SNR after CAF** |
|--------:|--------:|-----------:|-------------:|------------------:|
| 10 | 10 | −121.1 dBW | +9.5 dB | **+74.9 dB** |
| 20 | 20 | −133.1 dBW | −2.5 dB | **+62.9 dB** |
| 30 | 30 | −140.2 dBW | −9.6 dB | **+55.8 dB** |
| 50 | 50 | −149.0 dBW | −18.4 dB | **+47.0 dB** |
| 80 | 80 | −157.2 dBW | −26.6 dB | **+38.8 dB** |
| 120 | 120 | −164.2 dBW | −33.6 dB | **+31.8 dB** |

Small aircraft (Cessna, 5 m² = 7 dBsm) still gives **+49.9 dB at 20/20 km**
and **+37.8 dB at 40/40 km**.

Those margins are enormous, and they are not a trick: a 6.9 MHz-wide
signal integrated for half a second is a time-bandwidth product of 3.5
million. The megawatt floodlight is already paid for. **On thermal noise
alone, this rig should see airliners well past 100 km.**

## So what actually stops us

The echo is not competing with thermal noise. It is competing with the
**direct path from the tower**, which arrives at the same antenna:

| Tower distance | Direct signal | Echo (20/20 km) | Echo is |
|---:|---:|---:|---:|
| 5 km | −44.1 dBW | −133.1 dBW | **89.1 dB down** |
| 15 km | −53.6 dBW | −133.1 dBW | **79.5 dB down** |
| 30 km | −59.6 dBW | −133.1 dBW | **73.5 dB down** |

Three consequences, in order of how much they should worry us:

**1. Cancellation depth is the whole project.** We must suppress the direct
path and all zero-Doppler clutter by 70–90 dB before a target is visible.
ECA (least-squares projection onto delayed copies of the reference)
typically buys 30–70 dB. That gap is the research question. It is a DSP
problem, which is our home turf — but it is not a solved one, and any
claim of detection has to survive the question "is this a real echo or a
cancellation residue?"

**2. The ADC is close to the wall.** The RSPdx is 14-bit: 6.02·14 + 1.76 =
**86 dB** of theoretical dynamic range, and the echo sits 73–89 dB below
the direct signal. A single sample cannot represent both. This works
anyway — an OFDM direct signal is noise-like and self-dithering, so
quantisation noise spreads white across the band and the CAF's 65 dB of
processing gain digs the echo back out. But it means **gain staging is not
a detail here, it is the experiment**: too much gain and the front end
compresses, too little and the echo falls below the quantisation floor.

**3. We have already measured this rig compressing on this exact signal.**
E65/E73 (8/10) found the front end in overload from strong local 8VSB —
`rfgain_sel=2` gave 0.0 % FEC while `rfgain_sel=4` gave 100 %, and the
splitter's 3.5 dB of attenuation turned out to be *therapeutic*. That
history is now directly relevant: **an overloaded front end generates
intermodulation products that look exactly like targets.** Any detection
campaign must record its gain and include a compression control.

## Resolution we can expect

| quantity | value | note |
|---|---|---|
| Range resolution | **43 m** (c/B) | bistatic; 21.7 m if quoted as c/2B |
| Doppler resolution @ T=0.17 s | **5.9 Hz** = 3.3 m/s | the walk-limited dwell |
| Doppler resolution @ T=0.5 s | 2.0 Hz = 1.11 m/s | **only for a slow/near-tangential target** |
| Airliner, 250 m/s **of sum path** | **449 Hz** | ← what the arithmetic below is |
| Airliner, 250 m/s **radial**, β→0 | **899 Hz** | ← what a real aeroplane does |
| Meteor (40 km/s) | **72 kHz** | needs a deliberately wide Doppler search |

**⚠ Doppler convention, corrected 8/20.** An earlier version of this table
said "250 m/s radial → 449 Hz". That is wrong by a factor of two.
f_d = (1/λ)·d(R_t+R_r)/dt — the rate of the **sum** path. 449 Hz is what
you get if 250 m/s *is* the sum-path rate; a target actually moving 250 m/s
radially in near-monostatic geometry changes the sum path at 500 m/s and
gives **899 Hz**. **Size Doppler search windows off the larger number**, or
fast targets fold outside the unambiguous span and vanish.

## Range walk — the constraint I missed

Coherent integration only pays if the target stays in one range cell for
the whole dwell. Ours are 43.4 m:

| sum-path rate | T = 0.5 s | T = 0.25 s | T = 0.17 s | T = 0.08 s |
|---|---|---|---|---|
| 250 m/s | 2.88 bins | 1.44 | 0.86 | 0.46 |
| 500 m/s (β→0, 250 radial) | 5.76 bins | 2.88 | 1.73 | 0.92 |

Longest dwell before the target smears past one bin: **173 ms** at 250 m/s
sum-path, **87 ms** at 500. So the honest processing gain is **60.8 dB**
or **57.8 dB**, not 65.4 — a **4.6 to 7.6 dB** haircut on the whole table
above.

It does not change the verdict (we are working with 30–75 dB of margin),
but it sets the v0 dwell at ~150 ms and names the upgrade: **keystone /
range-walk migration compensation** buys the longer dwell back. Not built.

## The ECA blind zone is a notch, not a gate

Worth stating because I assumed the opposite and the synthetic test
disproved it: the cancellation basis is *static* delayed copies of the
reference, so what it removes is `delay < n_taps` **AND** `Doppler ≈ 0`.
Measured on synthetic truth: a static scatterer inside the tap span is
killed by 180 dB; a *moving* target at the same delay loses **0.00 dB**.

The real consequence is not a dead range gate — it is that **a target
whose bistatic Doppler is near zero is invisible at any range**. That is
geometry, not a bug: an aircraft crossing the tower–receiver baseline
tangentially has no sum-path rate. Such targets will simply not appear,
and a track that fades as it crosses the baseline is expected physics, not
a failure.

## What this buys us that other people's passive radar does not

We **decode the illuminator**. A full ATSC 3.0 receiver means we can
re-encode the decoded bitstream into a mathematically perfect reference
signal — no receiver noise in the reference at all. Most passive radar
work uses the noisy received direct signal as its reference, which limits
cancellation and raises the noise floor. This is v2, and it is the part
of the project that is genuinely ours.

## Experimental design (the part that makes it science)

Detection claims from correlation peaks are cheap. These are the controls:

1. **Synthetic ground truth first.** Inject an echo at a *known* range and
   Doppler into a synthetic reference and prove the processor recovers it.
   Then prove that with the echo 60 dB down it is *masked* without
   cancellation and *found* with it — so cancellation is shown to be doing
   the work. (Built: `tools/test_pr_core.py`.)
2. **False-alarm control.** Pure noise in must yield zero CFAR detections.
   A detector that finds aircraft in noise is worse than no detector.
3. **ADS-B truth.** Predict each known aircraft's bistatic delay and
   Doppler from its ADS-B state vector and check that our peaks land where
   geometry says they must. Peaks that track real aeroplanes are evidence;
   peaks that do not are residue.
4. **The negative control that matters:** a capture taken with the antenna
   disconnected, or off-carrier, must produce nothing. We already have
   `e65_VOID_antenna_disconnected_rf33.cs16` banked — precisely this kind
   of control, from a previous campaign.
5. **Gain control.** Repeat at two `rfgain_sel` values. A "target" that
   moves or vanishes with front-end gain is intermod, not aluminium.

## Data we already have

180 s of RF25 at 6.912 MS/s (`data/fox25/rf25_fox_g4_0810_2254.cs16`,
4.7 GB) plus a second capture at a different gain 3 minutes later
(`_g2_`, 2.4 GB) — **the gain control is already banked**, by luck, from
the E76 campaign. Every aircraft that flew through those 180 s on the
evening of 8/10 is already in the file. We have never looked.

The one thing those captures lack is contemporaneous ADS-B. So:
**v0 = algorithm + candidate detections on banked data** (self-consistency:
do peaks move coherently, like a target, over successive frames?), then
**v1 = a fresh capture with ADS-B truth pulled live from the OpenSky API**
during the recording window — which needs no second SDR.

## Honest limits

- One receiver, one illuminator gives an **ellipsoid, not a position**: we
  measure bistatic range and range-rate, not azimuth or height. Multiple
  illuminators (RF25, RF33, RF8 are physically different towers) intersect
  to a fix — that is the multistatic v2, and it needs either sequential
  captures of a slow target or a second SDR.
- **Rain is invisible.** UHF scatters almost nothing off precipitation;
  weather radar uses far higher frequencies and big dishes. No weather
  claims from this instrument.
- **The moon is a stretch** — 2.5 s round-trip delay makes it
  unmistakable if seen, but the path loss is another ~40 dB beyond
  anything in the table above. We would predict its delay/Doppler from
  ephemeris first and only then look, never the reverse.
- A discone has no directivity, so we cannot spatially separate the direct
  path from targets. A directional receive antenna pointed *away* from the
  tower is the cheapest possible upgrade and would buy real cancellation
  margin for free.
