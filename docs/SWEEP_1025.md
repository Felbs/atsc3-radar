# Every broadcast channel, one receiver: the 2026-10-25 sweep

One RSPdx, one pass, 02:08-02:28 UTC: three 8 MHz FM spans (88-108 MHz) on
the discone, VHF TV RF 7/9/12 on the discone, UHF 8-VSB RF 15-36 on the yagi,
and WNUV (ATSC 3.0, RF25). 16 legs, 15 usable (RF15 gapped during capture and
was voided by the integrity rule). Every leg went through the radar for its
signal type with no per-channel tuning: `tools/sweep_capture.py` then
`tools/sweep_run.py`.

The point was a stress test, and it was one. The sweep found five bugs that
single-channel runs never hit, and every number below is from after the fixes.

**No ADS-B truth yet.** OpenSky's anonymous quota ran out before the sweep;
truth comes from the adsb.lol daily history for 2026-10-25 (published
~03:30 UTC on the 26th) via `tools/adsb_history.py`. Until then a track's
only check is internal: its range walk against its own Doppler.

## TV

| channel | signal | antenna | blocks locked | sync-symbol SNR | cancellation | tracks | longest track | walk vs Doppler |
|---|---|---|---|---|---|---|---|---|
| RF7 WJLA | 8-VSB | discone | 58/59 | 19.5 dB | 21.3 dB | 85 | 100 dwells | 6.0 % |
| RF9 WUSA | 8-VSB | discone | 0/59 | 7.2 dB | 10.2 dB | 3 | 29 | 9.7 % |
| RF12 WBAL | 8-VSB | discone | 0/59 | 1.7 dB | 5.2 dB | 5 | 7 | 6.1 % |
| RF21 WMPT | 8-VSB | yagi | 59/59 | 28.1 dB | 28.1 dB | 33 | 117 | 5.6 % |
| RF22 WMPB | 8-VSB | yagi | 0/59 | -0.8 dB | 0.7 dB | (noise) | | |
| RF26 WBFF | 8-VSB | yagi | 0/59 | 8.0 dB | 12.2 dB | (noise) | | |
| RF27 WMAR | 8-VSB | yagi | 0/59 | 8.0 dB | 13.5 dB | 1 | 3 | |
| RF31 WETA | 8-VSB | yagi | 59/59 | 29.3 dB | 30.6 dB | 49 | 322 | 1.3 % |
| RF34 WRC | 8-VSB | yagi | 59/59 | 30.5 dB | 32.0 dB | 35 | 213 | 5.6 % |
| RF35 WPXW | 8-VSB | yagi | 59/59 | 26.3 dB | 28.8 dB | 81 | 97 | 3.4 % |
| RF36 WTTG | 8-VSB | yagi | 59/59 | 27.9 dB | 29.6 dB | 39 | 171 | 1.1 % |
| RF25 WNUV | ATSC 3.0 | yagi | | | | 7 | 48 | 2.2 % |

Dwells are 164 ms for 8-VSB, so RF31's longest track is 53 s of a 60 s capture.
RF7 is after field-comb rejection (below): 1029 of its 2613 detections were
the comb, and its track count fell from 129 to 85. Its farthest track, 80
dwells from 57.2 to 53.8 km, walks at -257 m/s against -258 predicted. One
block of the rerun came back VOID while the machine was busy with fusion runs.

Six of the eleven 8-VSB channels work as illuminators from here: the five
DC-area UHF transmitters and RF7. The Baltimore channels (RF12, 26, 27) and
RF22 do not. Their true SNR is about 8 dB, 7 dB under the threshold where
8-VSB can be decoded at all. WNUV works from the same direction because the
ATSC 3.0 reference is synthesized from pilots, which needs far less SNR than
deciding every symbol.

## FM

59 stations across three spans, each mixed to DC out of its 8 MHz capture
and processed with `tools/fm_radar.py` at 250 kS/s.

| | stations |
|---|---|
| locked (cancellation >= 10 dB) | 32 |
| cancellation >= 20 dB | 23 |
| cancellation >= 25 dB | 15 |
| locked, with a >= 50-dwell track walking within 10 % of its Doppler | 10 |

The ten with long, self-consistent tracks:

| station | cancellation | tracks | longest | walk vs Doppler |
|---|---|---|---|---|
| 99.5 WIHT | 24.1 dB | 107 | 206 dwells | 5.7 % |
| 107.7 WWWT-FM | 23.7 dB | 39 | 118 | 4.5 % |
| 99.1 WDCH-FM | 25.4 dB | 93 | 104 | 3.3 % |
| 97.9 WIYY | 10.5 dB | 59 | 104 | 9.6 % |
| 101.1 WWDC | 27.6 dB | 36 | 102 | 1.4 % |
| 90.5 WMUC-FM | 16.4 dB | 38 | 94 | 1.8 % |
| 104.5 WGRX | 25.0 dB | 44 | 93 | 8.2 % |
| 92.5 WAIW | 16.6 dB | 15 | 80 | 6.6 % |
| 100.3 WBIG-FM | 26.2 dB | 116 | 58 | 7.6 % |
| 92.7 WPRS | 24.0 dB | 40 | 52 | 7.0 % |

FM dwells hop 0.5 s, so WIHT's 206 dwells are 103 s. Station identities are
the nearest FCC licensee on the frequency; RDS was not checked.

The 27 that did not lock are not random. Seventeen cancel 0.1-0.2 dB, and
every one sits 200 kHz from a station that cancels 22 dB or more: the strong
neighbour's swing reaches into the weak station's band and captures the
discriminator, so no reference is rebuilt. Those stations still produced
detections -- 24,941 of the sweep's 42,865 FM detections came from stations
that never locked -- which is why the table now has a lock column.

## What the sweep broke

1. **FM decimation.** `fm_radar.py` decimated by a fixed 8, so an 8 MS/s span
   landed at 1 MS/s instead of the station's 250 kS/s. Neighbouring channels
   sat inside the surveillance band and every station cancelled a flat
   7-14 dB. Decimation now follows the station file: WFLS 12.6 -> 28.4 dB.
2. **The upper neighbour's pilot (8-VSB).** The 8-VSB surveillance rate spans
   +-3.59 MHz, the channel only +-2.69 MHz, and the next channel up puts its
   pilot 3.31 MHz above centre. RF35 had RF36's pilot 7.7 dB above its own
   in-band level. A channel filter on surveillance and references:
   RF35 7.0 -> 28.8 dB, RF34 15.2 -> 32.0, RF31 27.8 -> 30.6. RF34's
   shortfall had first been blamed on front-end compression; it was not.
3. **Field syncs on weak channels.** The field-sync finder compared hard
   decisions one field apart and needed 90 % agreement; a weak channel
   agrees with itself only ~43 % there, so its field syncs were never
   blanked from the reference. It now uses the soft correlation (field sync
   0.89, data 0.00 on RF27).
4. **"Locked" meant "returned".** RF22 reported 59 locked blocks at +20.9 ppm
   and 3.6 dB. A block now counts as locked only if the known segment-sync
   symbols come back at >= 10 dB SNR. The decision-directed SNR the radar
   printed before reads 15 dB on channels whose sync symbols say 8.
5. **Silent skips and duplicates.** RF22 had no station file and was skipped
   without a row; five FM stations in the overlap of two spans were written
   to one file twice. Both now report or name per span.

Also built and tested on the way:

- **Trellis decoding** (`tools/vsb_trellis.py`, `vsb_ref.py radar --trellis`).
  The code was identified from the RF31 decisions rather than assumed: one
  relation holds on 100.00 % of 1.8 M symbols, every other candidate on 50 %.
  On synthetic symbols it takes 15 dB from 19 % symbol errors to 0.4 %. The
  Baltimore channels are below its threshold, so it does not rescue them; it
  is for channels between roughly 13 and 20 dB.
- **The 8-VSB field comb** (`vsb_ref.py comb_flags`). RF7 put a ladder of
  detections at 0.96 / 1.92 / 2.80 / 5.68 km for 20 s, at +-1, 3 and 4 x
  41.3 Hz (the field rate) in mirror pairs. It is not static clutter (a
  second cancellation pass against the CAF reference removes 0.01 dB), not
  slow clutter (Doppler-shifted cancellation taps change nothing), and no
  second pilot is visible. Its cause is still open; its signature is not:
  one range cell, one dwell, three or more detections on one f0 + k x 41.3 Hz
  comb. Those, and every comb-line detection in a cell that shows it in five
  or more dwells, are flagged and kept out of tracks. On RF21/31/34/35/36 the
  rule flags 0-3 detections and no point of any track of 50 dwells or more.

## A position from one receiver: three FM towers at once

The 102-108 MHz span holds WBQB 101.5, WLNO 104.1 and WJFK 106.7, three
towers in three directions, recorded at the same instant. Each gives a
bistatic range and a Doppler; together they give a position and a velocity.
This is the multistatic fix that otherwise needs a second receiver.

Matching per-station detections (`tools/fix_position.py`) found nothing a
null did not. So did fusing the maps on range alone (`tools/fuse_fix.py`
v1): that searches about a million range-cell triples, and noise reached an
aircraft's score (165.7 against a null of 162.5).

`tools/fuse_fix.py` v2 puts Doppler inside the score. At each grid point, two
stations' Doppler peaks fix a horizontal velocity, which predicts the third
station's Doppler to within one 0.59 Hz cell. The score is the weakest
station's value there, chosen before the run. The nulls shift one station by
60 s, each in turn. A one-dwell shift is not a null: an aircraft moves 100 m
in 0.5 s, inside one 1.2 km cell.

| | |
|---|---|
| highest null peak (the threshold) | 5.45 |
| real peaks above it | 10.66 at 8.0 s, 5.54 at 69.5 s |
| fix A | 234 m/s, heading 53 deg, 11,084 m |
| fix B | 232 m/s, heading 53 deg, 11,102 m |
| A projected 61.5 s along its own velocity, distance from B | 0.62 km |
| follow support, A / B (dwells in +-30 s with two or more stations' detections on the projected range and Doppler, of 121) | 27 / 42 |
| follow support, the 30 top null peaks | median 0, max 13 |

A and B come from windows with no data in common. They are one aircraft at
cruise, positioned twice from one receiver. WBQB, the weakest of the three at
22 dB of cancellation, sees it only in the fix windows; WLNO and WJFK carry
the follow. The positions are in the lab output for the ADS-B check.

## Follow-ups, 2026-10-25 afternoon

**RF15 again.** Recaptured 15:40Z (-3.4 dBFS, 0 overflows): 59/59 blocks
locked, sync-symbol SNR 27.8 dB, cancellation 28.8 dB, 85 tracks, longest
181 dwells at 1.2 %. Seven of the eleven 8-VSB channels work.

**The RF7 comb was an episode.** RF7 recaptured 13.5 h after the sweep:
3 field-comb detections of 2745 (1029 of 2613 at 02:08Z). Something on
the air that night, not the transmitter's or this receiver's fixed
behaviour. The rule stays.

**The 17 FM stations that never locked are not stations.** Power in the
upper half of the channel (+30..+70 kHz) against the lower half: 7-30 dB
lopsided on all 17, within 4 dB on all 42 others. The loud half faces the
strong neighbour: it is that station's HD Radio sideband (129-198 kHz from
its carrier), with no FM carrier under it. No filter rescues them.
`fm_radar.py` now measures this first and writes `status: no carrier`
instead of running (`--force` overrides).

**Trellis decoding, measured where it matters.** RF31 with noise added, a
-50 dB echo planted at 20 km / +150 Hz, the echo's SNR in the map:

| sync-symbol SNR | slicer | trellis |
|---|---|---|
| 29.2 dB | 39.8 dB | 39.8 dB |
| 20.5 dB | 29.1 dB | 30.8 dB |
| 16.9 dB | 17.0 dB | 27.2 dB |
| 14.4 dB | 9.4 dB | 20.8 dB |

Cancellation reads WORSE with the trellis there (17.4 vs 18.5 dB at
16.9): the slicer's errors follow the noise, so its reference "cancels"
noise too. Cancellation is not the measure of a reference at low SNR.
`vsb_ref.py` now decodes per block with the trellis below 25 dB sync SNR
(`--trellis auto`, the default). RF7 (19.8 dB) on the recapture: 88 ->
123 tracks, 2745 -> 3730 detections; an outbound aircraft 28.7 -> 38.3 km
at 0.5 % now holds as one 211-dwell track (the slicer detected 155 of its
points, median 13.5 dB; trellis 211, 15.4 dB).

**The three-tower fix, checked on towers it did not use.** The fix's
position and velocity predict its range and Doppler on every other
station in the capture. Along its own path (+-15 s), against 100-400
random paths at the same altitude and speed:

| held-out station | t = 8 s | t = 69.5 s | random-path 99th percentile |
|---|---|---|---|
| 105.1 WAVA | 9.3 | 6.1 | 1.6 |
| 107.7 WWWT | 10.2 | 28.3 | 1.4-1.9 |
| 103.5 WTOP | 5.1 | 10.6 | 1.4 |
| 107.3 WLVW | 124.5 | 10.8 | 1.6-1.8 |

`fuse_fix.py --heldout` now does this for every top peak and every null
peak. On the 102-108 span the 30 null peaks are seen by 0 held-out
stations each; 15 of the 20 best real peaks by 2-4 of 4, most of them
BELOW the 60 s-shift threshold: the aircraft above, now from t = 3.5 to
92 s, and a second one at the same altitude, 9-11 km from it, on the
reciprocal heading. The shift null is too strict in a busy sky.

The other spans, with the rule "seen by more held-out stations than any
null peak":

| span | set | null peaks seen by | real peaks |
|---|---|---|---|
| 88-96 | WFLS 93.3, WIAD 94.7, WPRS 92.7; held out WETA, WGTS, WPFW, WKYS | <= 2 | 7 seen by 4 of 4: three aircraft (222 deg / 248 m/s; 163 deg / 295 m/s; 283 deg / 167 m/s low) |
| 96-104 | WBIG 100.3, WBQB 101.5, WDCH 99.1; held out WMZQ, WTVR, WWDC, WIHT | <= 3 | best 3: nothing |

Here null peaks are seen too: a null shifts one station, the other two
still hold the real aircraft, and a null peak can sit on one. The 88-96
held-out stations are two towers, not four (WETA/WGTS and WPFW/WKYS share
sites). All of it waits on ADS-B.

## Open

- ADS-B truth for every leg when the 2026-10-25 history publishes, and the
  fused fix's two positions against it.
- What drives the RF7 field comb.
- Fused fixes on the other FM spans (stations on three different towers).
- ADS-B for the 88-96 MHz fixes and the second 102-108 aircraft.
- A quiet-sky session (scheduled 2026-10-26 08:00Z).
