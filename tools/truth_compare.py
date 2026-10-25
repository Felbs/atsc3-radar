"""truth_compare.py -- do our peaks land where the aeroplanes actually were?

The only test that can turn a kinematically-consistent candidate into a
detection. For every ADS-B aircraft OpenSky saw during the capture, predict
from ITS OWN telemetry where it must appear in our range-Doppler map:

    bistatic excess range  R_b(t) = |T-P(t)| + |P(t)-R| - |T-R|
    bistatic Doppler       f_d(t) = -(1/lambda) * dR_b/dt

T = transmitter (WNUV, Television Hill, public FCC coordinates),
R = receiver (our GPS fix, private, read at runtime and never printed),
P = the aircraft's ADS-B position, interpolated to each dwell time.

Our CAF's delay axis is RELATIVE TO THE DIRECT PATH (lag 0 = the tower's
signal), so the observable is the EXCESS range, not the sum path. The
sign convention matches pr_core.bistatic_velocity: +f_d = closing.

Then compare against track_hunt's detections on the same capture, both
ways -- predicted-but-missed is as informative as detected-and-matched --
and against a NULL: the same truth with its timeline reversed, which keeps
every aircraft's geometry and speed but decouples it from our clock. A
match rate that survives time reversal is coincidence.

Honest limits (printed with every run): ADS-B positions are 10 s apart
(interpolated), altitude is geometric where available else barometric,
capture start is known to ~1 s from the launcher's wall clock (bootstrap
anchoring inside the file is exact; the MAP between file sample 0 and UTC
is not), and the tower antenna height is an estimate.

Usage: python tools/truth_compare.py <meta.json> <tracks.json>
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import radar_site                                              # noqa: E402

# WNUV (RF25, Baltimore ATSC 3.0 lighthouse), Television Hill, Woodberry.
# Public: FCC / Wikipedia infobox, fetched 2026-08-20. ERP 750 kW, HAAT 372.8 m.
TX_LAT, TX_LON = 39.336222, -76.649417
TX_ANT_AMSL_M = 400.0        # ESTIMATE: ~120 m ground + ~280 m up the candelabra
FC = 539.0e6
LAM = 2.998e8 / FC
FS = 6.912e6
BIN_M = 2.998e8 / FS
C = 2.998e8

MAX_RANGE_M = 60e3          # map extent used for "in map" predictions
STATION_NAME = "RF25 / WNUV (97 km)"
NOTCH_HZ = 15.0
R_GATE_OVERRIDE = None       # a station json may set the match gates
F_GATE_OVERRIDE = None


def configure(station_json):
    """Point the geometry at another illuminator (FM etc.). Sets the module
    globals every tool reads: TX_LAT/LON/ANT, FC, LAM, FS, BIN_M, map
    extent and match gates. Call before anything else."""
    global TX_LAT, TX_LON, TX_ANT_AMSL_M, FC, LAM, FS, BIN_M, MAX_RANGE_M
    global R_GATE_OVERRIDE, F_GATE_OVERRIDE, STATION_NAME, NOTCH_HZ
    st = json.load(open(station_json))
    STATION_NAME = st.get("short_name", st.get("name", "station"))
    NOTCH_HZ = float(st.get("notch_hz", 15.0))
    TX_LAT, TX_LON = float(st["tx_lat"]), float(st["tx_lon"])
    TX_ANT_AMSL_M = float(st.get("tx_alt_m", 300.0))
    FC = float(st["fc_hz"]); LAM = 2.998e8 / FC
    FS = float(st["fs_hz"]); BIN_M = 2.998e8 / FS
    MAX_RANGE_M = float(st.get("max_range_m", MAX_RANGE_M))
    R_GATE_OVERRIDE = st.get("r_gate_m"); F_GATE_OVERRIDE = st.get("f_gate_hz")
    return st


def station_arg():
    """--station=<json> anywhere on argv: configure and strip it."""
    keep = []
    for a in sys.argv[1:]:
        if a.startswith("--station="):
            configure(a.split("=", 1)[1])
        else:
            keep.append(a)
    sys.argv[1:] = keep


# WGS84
_A = 6378137.0
_E2 = 6.69437999014e-3


def ecef(lat_deg, lon_deg, h_m):
    la, lo = math.radians(lat_deg), math.radians(lon_deg)
    n = _A / math.sqrt(1 - _E2 * math.sin(la) ** 2)
    return np.array([(n + h_m) * math.cos(la) * math.cos(lo),
                     (n + h_m) * math.cos(la) * math.sin(lo),
                     (n * (1 - _E2) + h_m) * math.sin(la)])


def observer():
    """Receiver position from lab/site.json (radar_site.py); never written anywhere."""
    return radar_site.observer()


def load_truth(path):
    """icao24 -> sorted list of (unix_t, lat, lon, alt_m, vel, track, vrate, callsign)."""
    ac = {}
    for line in open(path):
        row = json.loads(line)
        if "states" not in row:
            continue
        for s in row["states"]:
            icao, cs, _c, tpos, _last, lon, lat, baro, onground, vel, trk, vr, \
                _sens, geo = s[:14]
            if onground or lat is None or lon is None or tpos is None:
                continue
            alt = geo if geo is not None else baro
            if alt is None:
                continue
            ac.setdefault(icao, []).append(
                (float(tpos), float(lat), float(lon), float(alt),
                 vel, trk, vr, (cs or "").strip()))
    for k in ac:
        ac[k] = sorted({r[0]: r for r in ac[k]}.values())
    return ac


def interp(samples, t):
    """Linear ECEF interpolation of an aircraft at unix time t, or None."""
    ts = [s[0] for s in samples]
    if t < ts[0] - 12 or t > ts[-1] + 12:
        return None
    if len(samples) == 1:
        s = samples[0]
        return ecef(s[1], s[2], s[3])
    i = int(np.searchsorted(ts, t))
    i0, i1 = max(0, i - 1), min(len(samples) - 1, i)
    if i0 == i1:
        s = samples[i0]
        return ecef(s[1], s[2], s[3])
    s0, s1 = samples[i0], samples[i1]
    w = 0.0 if s1[0] == s0[0] else (t - s0[0]) / (s1[0] - s0[0])
    w = min(1.0, max(0.0, w))
    return (1 - w) * ecef(s0[1], s0[2], s0[3]) + w * ecef(s1[1], s1[2], s1[3])


def predict(samples, t, T, R, dRT):
    p0 = interp(samples, t - 0.5)
    p1 = interp(samples, t + 0.5)
    if p0 is None or p1 is None:
        return None
    rb0 = np.linalg.norm(T - p0) + np.linalg.norm(p0 - R) - dRT
    rb1 = np.linalg.norm(T - p1) + np.linalg.norm(p1 - R) - dRT
    rb = 0.5 * (rb0 + rb1)
    fd = -(rb1 - rb0) / 1.0 / LAM          # +Hz = closing (sum path shrinking)
    return rb, fd


def main():
    station_arg()
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    meta = json.load(open(sys.argv[1]))
    trk = json.load(open(sys.argv[2]))
    truth_path = meta["truth"]
    if not os.path.exists(truth_path):
        # absolute paths in meta predate the passiveTuna -> atsc3-radar -> broadcast-radar renames
        truth_path = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])),
                                  os.path.basename(truth_path))
    truth = load_truth(truth_path)
    lat, lon, alt = observer()
    R = ecef(lat, lon, alt)
    T = ecef(TX_LAT, TX_LON, TX_ANT_AMSL_M)
    dRT = float(np.linalg.norm(T - R))
    t_launch = float(meta["t_launch_unix"])
    # the capture tool opens the device and starts streaming ~1-2 s after
    # launch; bounded here as a stated uncertainty, not hidden
    T0_UNC_S = 2.0

    det = trk["detections"]
    print("=" * 72)
    print("truth_compare -- ADS-B truth vs our range-Doppler detections")
    print("=" * 72)
    print(f"  aircraft seen in the box during the window: {len(truth)}")
    print(f"  tower-receiver baseline: {dRT / 1e3:.2f} km  (direct path = lag 0)")
    print(f"  detections (off-comb, flattened): {len(det)}  "
          f"tracks: {len(trk['tracks'])}")
    print(f"  capture t0 uncertainty: +-{T0_UNC_S:.0f} s  "
          f"(at 180 m/s that is +-{180 * T0_UNC_S / BIN_M:.0f} range bins)")

    # --- predictions at every detection time, per aircraft
    times = sorted({d["t_s"] for d in det})
    preds = []                      # (icao, callsign, t_s, rb, fd)
    for icao, samples in truth.items():
        for ts in times:
            pr = predict(samples, t_launch + ts, T, R, dRT)
            if pr is None:
                continue
            rb, fd = pr
            if 0 < rb < MAX_RANGE_M and abs(fd) < 1000:
                preds.append((icao, samples[0][7], ts, rb, fd))
    in_window = sorted({p[0] for p in preds})
    print(f"\n  aircraft whose predicted (R_b, f_d) fall INSIDE our map "
          f"(<{MAX_RANGE_M/1e3:.0f} km excess, |f_d|<1 kHz): {len(in_window)}")
    for icao in in_window:
        ps = [p for p in preds if p[0] == icao]
        print(f"    {icao} {ps[0][1]:8s}  R_b {min(p[3] for p in ps) / 1e3:6.2f}"
              f"-{max(p[3] for p in ps) / 1e3:6.2f} km  "
              f"f_d {min(p[4] for p in ps):+7.1f}..{max(p[4] for p in ps):+7.1f} Hz")

    # --- match: detection within gate of a prediction at the same dwell
    R_GATE = R_GATE_OVERRIDE or (3 * BIN_M + 180 * T0_UNC_S)      # bins + timing slop
    F_GATE = F_GATE_OVERRIDE or 30.0

    def match(preds_, dets_):
        hits = []
        for d in dets_:
            for (icao, cs, ts, rb, fd) in preds_:
                if abs(ts - d["t_s"]) < 1e-3 and abs(d["range_m"] - rb) < R_GATE \
                        and abs(d["doppler_hz"] - fd) < F_GATE:
                    hits.append((d, icao, cs, rb, fd))
                    break
        return hits

    hits = match(preds, det)
    print(f"\n  MATCHES (|dR| < {R_GATE:.0f} m, |df| < {F_GATE:.0f} Hz, same dwell): "
          f"{len(hits)} of {len(det)} detections")
    by_ac = {}
    for d, icao, cs, rb, fd in hits:
        by_ac.setdefault(icao, []).append((d, rb, fd))
    for icao, rows in by_ac.items():
        cs = next(p[1] for p in preds if p[0] == icao)
        print(f"    {icao} {cs:8s} matched in {len(rows)} dwells; e.g. "
              f"det {rows[0][0]['range_m'] / 1e3:.2f} km {rows[0][0]['doppler_hz']:+.0f} Hz "
              f"vs pred {rows[0][1] / 1e3:.2f} km {rows[0][2]:+.0f} Hz")

    # --- the other direction: predicted inside the map but never detected
    missed = [a for a in in_window if a not in by_ac]
    print(f"\n  PREDICTED-BUT-MISSED aircraft: {len(missed)} of {len(in_window)}")

    # --- NULL: time-reversed truth keeps geometry, breaks the clock.
    # 8/21 LESSON: over a 10 s window a quasi-static target's prediction
    # moves ~5 bins / ~10 Hz - inside the gate in EITHER direction - so
    # this null scores whatever real scores and cannot lose. It stays in
    # the output (the 8/20-21 verdicts were made under it and are not
    # renegotiated), but the shift nulls below are the ones with power.
    rev = [(icao, cs, times[len(times) - 1 - times.index(ts)], rb, fd)
           for (icao, cs, ts, rb, fd) in preds if ts in times]
    null_hits = match(rev, det)

    def shifted_null(dt_s):
        sh = []
        for icao, samples in truth.items():
            for ts in times:
                pr = predict(samples, t_launch + ts + dt_s, T, R, dRT)
                if pr and 0 < pr[0] < MAX_RANGE_M and abs(pr[1]) < 1000:
                    sh.append((icao, samples[0][7], ts, pr[0], pr[1]))
        return len(match(sh, det))
    null2 = shifted_null(60.0)
    shift_nulls = {dt: shifted_null(dt) for dt in (-300, -180, -60, 60, 180, 300)}
    print(f"\n  NULL CONTROLS (same detections, same aircraft, wrong clock):")
    print(f"    time-reversed truth : {len(null_hits)} matches   "
          f"(POWERLESS for quasi-static targets over short windows - see source)")
    for dt, n in shift_nulls.items():
        print(f"    truth shifted {dt:+4d} s : {n} matches")
    print(f"    real                : {len(hits)} matches")
    worst_shift = max(shift_nulls.values())
    print(f"\n  [next-run protocol, informational only tonight] real beats "
          f"every shift null by >3x: "
          f"{'YES' if len(hits) > 3 * max(worst_shift, 1) else 'NO'} "
          f"(worst shift null {worst_shift}); "
          f"distinct aircraft matched: {len(by_ac)}")

    # --- PERMUTATION NULL (pre-registered 8/21 for the confirming run):
    # each aircraft's geometric track is re-timed onto ANOTHER aircraft's
    # timeline (a random derangement, linear re-timing of its sample times
    # onto the other's span), 50x. Keeps every path and the population of
    # signatures, breaks which path flies when. p = fraction of draws whose
    # match count reaches real. Its power is REPORTED, not assumed: if the
    # null median sits near real, it is as blind as time-reversal was.
    rng = np.random.default_rng(20260821)
    icaos = [k for k in truth if len(truth[k]) >= 2]

    def retimed(samples, onto):
        t0, t1 = samples[0][0], samples[-1][0]
        u0, u1 = onto[0][0], onto[-1][0]
        if t1 <= t0 or u1 <= u0:
            return None
        k = (u1 - u0) / (t1 - t0)
        return [(u0 + (s[0] - t0) * k,) + tuple(s[1:]) for s in samples]

    def permutation_null():
        perm = rng.permutation(len(icaos))
        if any(perm[i] == i for i in range(len(icaos))):
            perm = np.roll(np.arange(len(icaos)), int(rng.integers(1, len(icaos))))
        ph = []
        for i, icao in enumerate(icaos):
            samples = retimed(truth[icao], truth[icaos[perm[i]]])
            if samples is None:
                continue
            for ts in times:
                pr = predict(samples, t_launch + ts, T, R, dRT)
                if pr and 0 < pr[0] < MAX_RANGE_M and abs(pr[1]) < 1000:
                    ph.append((icao, samples[0][7], ts, pr[0], pr[1]))
        return len(match(ph, det))

    # Per-aircraft CLOCK SCRAMBLE: every aircraft keeps its own path and
    # speed but gets its own random clock offset, |dt| uniform in
    # [45, 120] s, random sign. Unlike the global shift nulls, no two
    # aircraft move together. Added 2026-10-24 BEFORE today's detections
    # were examined, because the identity permutation above proved blind
    # on the 8/21 data (median == real). Coverage loss is printed: a
    # shifted aircraft can fall outside the ~160 s of truth polling.
    def scramble_null():
        ph = []
        for icao in icaos:
            dt = float(rng.uniform(45.0, 120.0)) * (1 if rng.random() < 0.5 else -1)
            for ts in times:
                pr = predict(truth[icao], t_launch + ts + dt, T, R, dRT)
                if pr and 0 < pr[0] < MAX_RANGE_M and abs(pr[1]) < 1000:
                    ph.append((icao, truth[icao][0][7], ts, pr[0], pr[1]))
        return len(match(ph, det)), len(ph)

    N_PERM = 50
    perm_counts = ([permutation_null() for _ in range(N_PERM)]
                   if len(icaos) >= 3 else [])
    scr = [scramble_null() for _ in range(N_PERM)] if icaos else []
    scr_counts = [c for c, _ in scr]
    scr_cov = (float(np.mean([n for _, n in scr])) / max(len(preds), 1)) if scr else 0.0
    scr_p = (sum(1 for c in scr_counts if c >= len(hits)) / N_PERM) if scr else float("nan")
    scr_med = float(np.median(scr_counts)) if scr else float("nan")
    perm_p = (sum(1 for c in perm_counts if c >= len(hits)) / N_PERM
              if perm_counts else float("nan"))
    perm_med = float(np.median(perm_counts)) if perm_counts else float("nan")
    print(f"    identity permutation: median {perm_med:.0f}, max "
          f"{max(perm_counts) if perm_counts else 0}, p(>= real) = {perm_p:.2f}"
          f"  ({N_PERM} draws)")
    print(f"    per-aircraft clock scramble: median {scr_med:.0f}, max "
          f"{max(scr_counts) if scr_counts else 0}, p(>= real) = {scr_p:.2f}"
          f"  ({N_PERM} draws, prediction coverage {scr_cov:.0%} of real)")
    if perm_counts and perm_med >= 0.5 * max(len(hits), 1):
        print("    ** permutation null median is >= half of real: this null has "
              "little power here and is reported as such **")

    # --- close-range rule (8/21 addendum): matches at R_b < 5 km are
    # crowded by approach traffic and do NOT count toward the verdict
    # unless the aircraft is tracked across >= 2 windows (judged across
    # the per-window runs, outside this script).
    far_hits = [h for h in hits if h[3] >= 5e3]
    far_ac = sorted({h[1] for h in far_hits})
    # the same rule applied to the nulls: a null scored on the crowded
    # close-range cells would veto a far-range result it never tested.
    det_far = [d for d in det if d["range_m"] >= 5e3]
    preds_far = [p for p in preds if p[3] >= 5e3]

    def shifted_null_far(dt_s):
        sh = []
        for icao, samples in truth.items():
            for ts in times:
                pr = predict(samples, t_launch + ts + dt_s, T, R, dRT)
                if pr and 5e3 <= pr[0] < MAX_RANGE_M and abs(pr[1]) < 1000:
                    sh.append((icao, samples[0][7], ts, pr[0], pr[1]))
        return len(match(sh, det_far))
    shift_far = {dt: shifted_null_far(dt) for dt in (-300, -180, -60, 60, 180, 300)}

    def scramble_far():
        ph = []
        for icao in icaos:
            dt = float(rng.uniform(45.0, 120.0)) * (1 if rng.random() < 0.5 else -1)
            for ts in times:
                pr = predict(truth[icao], t_launch + ts + dt, T, R, dRT)
                if pr and 5e3 <= pr[0] < MAX_RANGE_M and abs(pr[1]) < 1000:
                    ph.append((icao, truth[icao][0][7], ts, pr[0], pr[1]))
        return len(match(ph, det_far))
    scr_far = [scramble_far() for _ in range(N_PERM)] if icaos else []
    scr_far_p = (sum(1 for c in scr_far if c >= len(far_hits)) / N_PERM) if scr_far else float("nan")
    worst_shift_far = max(shift_far.values())
    print(f"    R_b >= 5 km only: detections {len(det_far)}, real matches {len(far_hits)}, "
          f"shift nulls " + " ".join(f"{k:+d}s:{v}" for k, v in shift_far.items())
          + f", scramble median {float(np.median(scr_far)) if scr_far else 0:.0f} "
          f"max {max(scr_far) if scr_far else 0} p(>= real) = {scr_far_p:.2f}")
    print(f"    matches at R_b >= 5 km: {len(far_hits)} on {len(far_ac)} aircraft; "
          f"at R_b < 5 km (count only if tracked across windows): "
          f"{len(hits) - len(far_hits)}")

    # 8/20-21 rule, kept for comparability with those nights
    verdict = ("DETECTION-GRADE agreement" if len(hits) >= 3 and
               len(hits) > 3 * max(len(null_hits), null2, 1)
               else "NOT confirmed - real does not beat the nulls decisively"
               if hits else "NO matches - nothing we see is where the ADS-B "
                            "aircraft are (or nothing squawking was in the beam)")
    # confirming-run rule, pre-registered 8/21 01:1xZ (TRUTH_RUN_0820.md)
    beats_shifts = len(hits) > 3 * max(worst_shift, 1)
    verdict2 = ("DETECTION-GRADE" if len(hits) >= 5 and len(by_ac) >= 2 and
                beats_shifts and scr and scr_p < 0.01
                else "NOT confirmed")
    beats_shifts_far = len(far_hits) > 3 * max(worst_shift_far, 1)
    verdict2_far = ("DETECTION-GRADE" if len(far_hits) >= 5 and len(far_ac) >= 2
                    and beats_shifts_far and scr_far and scr_far_p < 0.01
                    else "NOT confirmed")
    print(f"\n  VERDICT (8/20 rule):        {verdict}")
    print(f"  VERDICT (8/21 confirming):  {verdict2}   "
          f"[real>=5: {len(hits) >= 5}, aircraft>=2: {len(by_ac) >= 2}, "
          f"beats shift nulls 3x: {beats_shifts}, scramble p<0.01: "
          f"{scr_p < 0.01 if scr else 'n/a'} (identity-perm p {perm_p:.2f}, blind)]")
    print(f"  VERDICT (8/21, R_b>=5 km only): {verdict2_far}   "
          f"[real>=5: {len(far_hits) >= 5}, aircraft>=2: {len(far_ac) >= 2}, "
          f"beats far shift nulls 3x: {beats_shifts_far} (worst {worst_shift_far}), "
          f"far scramble p<0.01: {scr_far_p < 0.01 if scr_far else 'n/a'}]")
    out = dict(n_aircraft=len(truth), n_in_map=len(in_window),
               n_detections=len(det), matches=len(hits),
               matched_aircraft={k: len(v) for k, v in by_ac.items()},
               matches_far=len(far_hits), matched_aircraft_far=far_ac,
               missed_aircraft=missed, null_reversed=len(null_hits),
               null_shift60=null2, shift_nulls={str(k): v for k, v in shift_nulls.items()},
               perm_counts=perm_counts, perm_p=perm_p,
               scramble_counts=scr_counts, scramble_p=scr_p, scramble_coverage=scr_cov,
               r_gate_m=R_GATE, f_gate_hz=F_GATE,
               t0_uncertainty_s=T0_UNC_S, tx_antenna_amsl_m_ESTIMATE=TX_ANT_AMSL_M,
               verdict=verdict, verdict_confirming=verdict2,
               verdict_confirming_far=verdict2_far,
               shift_nulls_far={str(k): v for k, v in shift_far.items()},
               scramble_counts_far=scr_far, scramble_p_far=scr_far_p)
    outp = Path(sys.argv[2]).with_name("truth_compare.json")
    json.dump(out, open(outp, "w"), indent=1)
    print(f"  -> {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
