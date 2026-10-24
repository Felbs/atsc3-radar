"""gain_pair.py -- the front-end control: two captures minutes apart at two
RF gains, same antenna, same window.

A front-end product (intermodulation of the direct signal, a pilot comb
leftover, a spur) lives at a FIXED (range, Doppler) cell and does not care
what the sky is doing; an aircraft moves hundreds of metres between two
captures and sits where its ADS-B says. So the test is two-sided:

  1. cells detected in BOTH captures at the same (R_b, f_d) and persisting
     in each are artifact candidates -- listed, and checked against the
     ADS-B-matched detections (an aircraft that happens to sit still in
     bistatic range for minutes is possible; the ADS-B match decides);
  2. the far-range ADS-B verdict must hold in both captures independently
     (truth_compare.py does that; this script just reads its JSON).

Usage: python tools/gain_pair.py <tracks_A.json> <tracks_B.json> [truth_compare_A.json truth_compare_B.json]
"""
import json
import sys
from collections import defaultdict

BIN_M = 2.998e8 / 6.912e6
R_TOL = 1.5 * BIN_M      # one bin either side
F_TOL = 6.0              # ~1.5 Doppler bins at a 170 ms dwell
MIN_PERSIST = 3          # dwells in each capture to count as a persistent cell


def cells(trk):
    """Group a capture's detections into persistent (R, f) cells."""
    groups = []
    for d in sorted(trk["detections"], key=lambda d: -d["snr_db"]):
        for g in groups:
            if abs(g["r"] - d["range_m"]) < R_TOL and abs(g["f"] - d["doppler_hz"]) < F_TOL:
                g["n"] += 1
                g["snr"] = max(g["snr"], d["snr_db"])
                g["ts"].append(d["t_s"])
                break
        else:
            groups.append(dict(r=d["range_m"], f=d["doppler_hz"], n=1,
                               snr=d["snr_db"], ts=[d["t_s"]]))
    return [g for g in groups if g["n"] >= MIN_PERSIST]


def main():
    A = json.load(open(sys.argv[1]))
    B = json.load(open(sys.argv[2]))
    tcA = json.load(open(sys.argv[3])) if len(sys.argv) > 4 else None
    tcB = json.load(open(sys.argv[4])) if len(sys.argv) > 4 else None
    ca, cb = cells(A), cells(B)
    print("=" * 72)
    print(f"gain_pair -- A: {len(A['detections'])} detections, {len(ca)} persistent cells;"
          f"  B: {len(B['detections'])} detections, {len(cb)} persistent cells")
    print(f"  (persistent = >= {MIN_PERSIST} dwells within {R_TOL:.0f} m / {F_TOL:.0f} Hz)")
    common = []
    for a in ca:
        for b in cb:
            if abs(a["r"] - b["r"]) < R_TOL and abs(a["f"] - b["f"]) < F_TOL:
                common.append((a, b))
    print(f"\n  cells present in BOTH captures at the same (R_b, f_d): {len(common)}")
    for a, b in sorted(common, key=lambda ab: ab[0]["r"]):
        print(f"    {a['r'] / 1e3:7.2f} km {a['f']:+7.1f} Hz   A: {a['n']:2d} dwells {a['snr']:4.1f} dB"
              f"   B: {b['n']:2d} dwells {b['snr']:4.1f} dB   <- artifact candidate")
    if not common:
        print("    (none) -- nothing persistent survived a gain change at a fixed cell")
    # anything inside the clutter notch is expected to be fixed; say so
    notch = [ab for ab in common if abs(ab[0]["f"]) < 15]
    if notch:
        print(f"  of which inside the +-15 Hz notch (zero-Doppler clutter, expected): {len(notch)}")
    if tcA and tcB:
        print("\n  far-range ADS-B verdict, each capture on its own:")
        for name, tc in (("A", tcA), ("B", tcB)):
            print(f"    {name}: matches {tc['matches']} ({tc['matches_far']} at >= 5 km on "
                  f"{len(tc['matched_aircraft_far'])} aircraft), worst far shift null "
                  f"{max(tc['shift_nulls_far'].values())}, far scramble p {tc['scramble_p_far']:.2f}"
                  f"  -> {tc['verdict_confirming_far']}")
        both = tcA["verdict_confirming_far"].startswith("DETECTION") and \
            tcB["verdict_confirming_far"].startswith("DETECTION")
        print(f"\n  GAIN-PAIR: far-range verdict holds at both gains: {'YES' if both else 'NO'};"
              f" fixed-cell artifact candidates outside the notch: {len(common) - len(notch)}")
    out = dict(common_cells=[dict(r_m=a["r"], f_hz=a["f"], nA=a["n"], nB=b["n"],
                                  snrA=a["snr"], snrB=b["snr"]) for a, b in common],
               n_cells_A=len(ca), n_cells_B=len(cb))
    json.dump(out, open(sys.argv[1].replace(".json", "_gainpair.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
