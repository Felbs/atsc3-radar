"""sweep_run.py -- run the right radar on every leg of a sweep and tabulate.

For each leg in a sweep manifest (sweep_capture.py):
  fm     every FM station file whose frequency falls inside the capture
         -> fm_radar.py (reconstructed reference, 1 s dwells)
  vsb    -> vsb_ref.py radar (decision-rebuilt 8-VSB reference)
  atsc3  -> track_hunt.py on three windows (pilot-synthesized reference)
and, when a truth file has been attached to the leg's meta (adsb_history.py
truth, or a live truth_capture run), truth_compare / truth_trajectory /
unmatched_tracks with the station's geometry.

The output is one table per sweep: did the reference lock, how deep did
the direct path cancel, how many detections and tracks, how many matched
ADS-B and how far -- the same columns for every illuminator, which is
the point of sweeping.

Usage: python tools/sweep_run.py <sweep_dir> [--jobs 8] [--truth] [--only fm|vsb|atsc3]
"""
import argparse
import glob
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(cmd, log):
    with open(log, "w") as f:
        p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    return p.returncode


def station_files(kind, leg):
    if kind == "fm":
        out = []
        for p in glob.glob(str(HERE / "stations" / "fm" / "*.json")):
            st = json.load(open(p))
            if abs(st["fc_hz"] - leg["fc_hz"]) < 0.45 * leg["rate_hz"]:
                out.append(p)
        return sorted(out, key=lambda p: -json.load(open(p)).get("survey_power_db", 0))
    rf = int(leg["tag"][2:])
    for p in glob.glob(str(HERE / "stations" / "*.json")):
        if p.endswith(f"_{rf}.json"):
            return [p]
    return []


def summarize_tracks(path):
    if not os.path.exists(path):
        return dict(detections=0, tracks=0)
    d = json.load(open(path))
    tr = d.get("tracks", [])
    best = max(tr, key=lambda t: t["n"], default=None)
    canc = d.get("cancellation_db")
    blocks = d.get("blocks")
    out = dict(detections=d.get("n_detections", len(d.get("detections", []))), tracks=len(tr),
               longest_track=best["n"] if best else 0,
               longest_walk_pct=best["agreement_pct"] if best else None,
               cancellation_db=round(float(sorted(canc)[len(canc) // 2]), 1) if canc else None)
    if blocks:
        out["eq_snr_db"] = round(sorted(b["eq_snr_db"] for b in blocks)[len(blocks) // 2], 1)
        out["clock_ppm"] = round(sorted(b["clock_ppm"] for b in blocks)[len(blocks) // 2], 2)
        out["blocks_locked"] = len(blocks)
    return out


def summarize_truth(path):
    if not os.path.exists(path):
        return {}
    t = json.load(open(path))
    return dict(matches=t.get("matches"), matches_far=t.get("matches_far"),
                aircraft=len(t.get("matched_aircraft", {})), in_map=t.get("n_in_map"),
                scramble_p_far=t.get("scramble_p_far"), verdict_far=t.get("verdict_confirming_far"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sweep_dir")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--truth", action="store_true", help="also run the truth chain where a truth file exists")
    ap.add_argument("--only", default=None, choices=("fm", "vsb", "atsc3"))
    ap.add_argument("--skip-done", action="store_true")
    a = ap.parse_args()
    sd = Path(a.sweep_dir)
    man = json.load(open(sd / "manifest.json"))
    rows = []
    for leg in man["legs"]:
        if a.only and leg["kind"] != a.only:
            continue
        if leg["rc"] != 0 or not os.path.exists(leg["capture"]):
            rows.append(dict(tag=leg["tag"], kind=leg["kind"], status="capture failed/refused"))
            continue
        meta = leg["meta"]
        for stp in station_files(leg["kind"], leg):
            st = json.load(open(stp))
            name = Path(stp).stem
            out = sd / f"tracks_{name}.json"
            t0 = time.time()
            if not (a.skip_done and out.exists()):
                if leg["kind"] == "fm":
                    rc = run([PY, str(HERE / "fm_radar.py"), leg["capture"], f"--station={stp}", "--out", str(out)],
                             sd / f"run_{name}.log")
                elif leg["kind"] == "vsb":
                    rc = run([PY, str(HERE / "vsb_ref.py"), "radar", leg["capture"], f"--station={stp}",
                              "--jobs", str(a.jobs), "--out", str(out)], sd / f"run_{name}.log")
                else:
                    rc = 0
                    for w in (5, 30, 50):
                        rc |= run([PY, str(HERE / "track_hunt.py"), leg["capture"], "--t0", str(w),
                                   "--out", str(sd / f"tracks_{name}_t{w}.json")], sd / f"run_{name}_t{w}.log")
                    # merge the windows into one tracks file for the table
                    dets, trs = [], []
                    for w in (5, 30, 50):
                        p = sd / f"tracks_{name}_t{w}.json"
                        if p.exists():
                            d = json.load(open(p)); dets += d["detections"]; trs += d["tracks"]
                    json.dump(dict(n_detections=len(dets), detections=dets, tracks=trs), open(out, "w"))
            else:
                rc = 0
            row = dict(tag=leg["tag"], kind=leg["kind"], station=st.get("short_name", name), fc_mhz=st["fc_hz"] / 1e6,
                       rc=rc, seconds=round(time.time() - t0), **summarize_tracks(str(out)))
            truth_path = json.load(open(meta)).get("truth")
            hist_meta = meta.replace(".meta.json", ".hist.meta.json")
            if a.truth and (truth_path or os.path.exists(hist_meta)):
                m = hist_meta if os.path.exists(hist_meta) else meta
                run([PY, str(HERE / "truth_compare.py"), f"--station={stp}", m, str(out)], sd / f"compare_{name}.txt")
                tc_json = sd / "truth_compare.json"
                if tc_json.exists():
                    os.replace(tc_json, sd / f"truth_compare_{name}.json")
                row.update(summarize_truth(str(sd / f"truth_compare_{name}.json")))
                run([PY, str(HERE / "unmatched_tracks.py"), f"--station={stp}", m, str(out)], sd / f"unmatched_{name}.txt")
            rows.append(row)
            print(json.dumps(row), flush=True)
            json.dump(rows, open(sd / "sweep_summary.json", "w"), indent=1)
    # the table
    cols = ["tag", "station", "fc_mhz", "blocks_locked", "eq_snr_db", "cancellation_db", "detections", "tracks",
            "longest_track", "longest_walk_pct", "matches", "matches_far", "aircraft", "in_map", "verdict_far", "seconds"]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + " |")
    open(sd / "sweep_summary.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
