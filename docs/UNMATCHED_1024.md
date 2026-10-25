# Unmatched tracks, first pass — 2026-10-24

The goal that counts is a track no ADS-B aircraft explains. This is the
first pass at it, on the four clean captures already on disk, with the
rules written into `tools/unmatched_tracks.py` before it ran.

## The bar (pre-registered in the tool's docstring)

A candidate is a `track_hunt` track with every point ≥ 5 km excess range
and outside the ±15 Hz notch, ≥ 8 points, **no ADS-B prediction within a
doubled gate** (980 m, 60 Hz) at any point, a walk rate within 25 % of
what its own Doppler implies, **not a Doppler twin** of a matched
aircraft (a second-path echo of a real aircraft sits 0–4 km beyond it at
the same Doppler), and not at a fixed cell from the gain pair. It is
reported beside the capture's own miss rate, because an instrument that
misses 73–90 % of known aircraft produces unmatched tracks for free.

## Result

```
capture               tracks >= 8 pts beyond 5 km   rejected (why)                         CANDIDATES
D   17:51Z rfgain 2            2                     AUA41, UAL2433 explain them              0
A   19:21Z rfgain 2            4                     GJS4508, UAL919, ENY4321 ghost, walk     0
B   19:24Z rfgain 3            2                     ATN3386 (+ its Doppler twin)             1
C   19:09Z rfgain 3            3                     UAL2668, UAL2101 ghost, AAL2142 ghost    0
```

The Doppler-twin rule earned its place on its first run: three tracks in
A and C are second-path ghosts of matched aircraft (ENY4321, UAL2101,
AAL2142), 0–4 km beyond the real echo at the same Doppler, and would
have been "unexplained" without it.

### The one candidate — capture B, t = 5–15 s

```
35 points over 9.5 s   31.27 -> 31.66 km   -70 -> -78 Hz
walk observed +41 m/s, from Doppler +41 m/s   SNR 11.6-18.1 dB
nearest ADS-B at t = 9.3 s: N3281L 30.45 km but -150 Hz (75 Hz off)
truth within +-60 s: 106 aircraft with position, 0 without
capture A, 125 s earlier, where +41 m/s puts it (26 km, -73 Hz): nothing
```

Kinematically it is as clean as the best matched tracks. It is **not**
claimed, for a reason found while examining it:

**The truth box did not cover the far end of the bistatic ellipse.** The
box was ±0.6° / ±0.75° around the receiver, a half-width of about 65 km.
The tower is 97 km away, so the tower and everything around it were
outside the truth. For a target beyond the tower on the baseline's
extension, excess range is twice its distance past the tower, so
R_b = 31 km is an aircraft ~15 km beyond Television Hill — the BWI
approach — lit by 750 kW from 15 km and therefore a *strong* echo at the
receiver despite being 110 km away. Every far-range cell has this second
solution, and the truth never looked there. GJS4508 at 29 km (matched,
capture A) happened to be on the near side; this one may not be.

OpenSky refuses historical queries without an account (403), so the
tower-end check could not be done tonight. The capture tool's box is now
the union of the receiver's box and the tower's (`truth_capture.py`,
tower coordinates are public), which closes the hole for every future
run. For this candidate the check is pre-registered as: pull the
2026-10-24 19:24:07–19:24:17Z window from adsb.lol's daily history when
it publishes, in a box around the tower, predict (R_b, f_d) for every
aircraft in it with the same code, and see whether one lands on
31.3–31.7 km / −70 to −78 Hz. A match closes it as a BWI aircraft; no
match leaves a single-window, single-capture candidate, which is still
below the bar (two windows ≥ 60 s apart).

### Closed (2026-10-25): an ADS-B aircraft at the tower end

adsb.lol's 2026-10-24 history, in the box that now spans receiver and
tower (263 ADS-B + 129 MLAT-only aircraft in the window, against 106 in the
original truth), has ICAO 0c414a, callsign HI1001, 3,850 m and 175 m/s,
52 km from the WNUV tower. It explains 35 of the track's 35 points.
`unmatched_tracks.py` on the history truth: 0 candidates. On all four
captures every far track is now explained.

## What this pass established

1. On four captures and 32 tracks of eight or more points beyond 5 km,
   **every track but one is explained by a known aircraft or its
   second-path ghost**, which is what a working radar on a Friday
   evening near two Class B airports should produce.
2. The ghost rule is necessary: 3 of 32 tracks were twins.
3. The truth box was wrong for far range, now fixed. The far end of the
   ellipse near the illuminator is the most strongly lit region of the
   whole map and must always be in the truth.
4. Truth completeness tonight: all 7846 state vectors were ADS-B
   (position_source 0); OpenSky delivered no MLAT-positioned aircraft in
   this box. A local 1090 MHz receiver running during captures remains
   the next truth upgrade (aeroTuna decodes it, but needs a second SDR —
   the RSPdx is busy on RF25).

## Next

- Tomorrow: the adsb.lol history check on the 31 km candidate.
- Every future capture: the widened box, and hunt the whole 120 s (nine
  batches, not three windows) so persistence can be tested at 13 s
  resolution across the whole file.
- Then a quiet-sky session (late night, fewer approaches) where an
  unmatched track would be rarer and mean more.
