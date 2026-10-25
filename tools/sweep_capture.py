"""sweep_capture.py -- capture every broadcast band the RSPdx can reach, one
leg after another, so the radar code can be run against all of it.

Legs (2026-10-25 plan):
  FM      3 x 8 MHz on the discone (88-96, 96-104, 102-108 MHz), 120 s each
  VHF TV  RF 7, 9, 12 (8-VSB) on the discone, 60 s each
  UHF TV  RF 15 21 22 26 27 31 34 35 36 (8-VSB) on the yagi, 60 s each
  ATSC 3  RF 25 (WNUV) on the yagi, 60 s

Each leg lands as <tag>.cs16 + <tag>.json (iq_capture's sidecar) +
<tag>.meta.json in truth_capture's shape, with t_launch_unix / t_done_unix
so adsb_history.py can attach ADS-B truth after the day's release. No
OpenSky polling here: the anonymous quota is spent and a 25-minute sweep
would spend it again.

iq_capture takes the fleet radio lock per leg; a leg refused by the lock
is recorded as skipped, not retried, so a radiosonde hunt wins.

Usage: python tools/sweep_capture.py [--only fm|vhf|uhf|atsc3] [--secs-tv 60] [--secs-fm 120]
"""
import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent / "lab"
# the capture needs SoapySDR: radioconda has it, the default python does not
RADIOCONDA = r"C:\Users\emane\radioconda\python.exe"
PY = RADIOCONDA if os.path.exists(RADIOCONDA) else sys.executable
IQ = str(HERE / "iq_capture.py")


def tv_fc(rf):
    lo = (174 + (rf - 7) * 6) if rf < 14 else (470 + (rf - 14) * 6)
    return (lo + 3.0) * 1e6


def legs(a):
    L = []
    if a.only in (None, "fm"):
        for c in (92.0, 100.0, 105.0):
            L.append(dict(kind="fm", tag=f"fm{c:.0f}", fc=c * 1e6, rate=8e6, ant="Antenna C",
                          rfgain=6, ifgr=45, secs=a.secs_fm))
    if a.only in (None, "vhf"):
        for rf in (7, 9, 12):
            L.append(dict(kind="vsb", tag=f"rf{rf}", fc=tv_fc(rf), rate=8e6, ant="Antenna C",
                          rfgain=4, ifgr=40, secs=a.secs_tv, rf=rf))
    if a.only in (None, "uhf"):
        hot = {15, 31, 34, 35, 36}            # at/near full scale at rfgain 4 in the 00:22Z survey
        for rf in (15, 21, 22, 26, 27, 31, 34, 35, 36):
            L.append(dict(kind="vsb", tag=f"rf{rf}", fc=tv_fc(rf), rate=8e6, ant="Antenna B",
                          rfgain=5 if rf in hot else 4, ifgr=32, secs=a.secs_tv, rf=rf))
    if a.only in (None, "atsc3"):
        L.append(dict(kind="atsc3", tag="rf25", fc=tv_fc(25), rate=6.912e6, ant="Antenna B",
                      rfgain=2, ifgr=32, secs=a.secs_tv, rf=25))
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, choices=("fm", "vhf", "uhf", "atsc3"))
    ap.add_argument("--secs-tv", type=float, default=60.0)
    ap.add_argument("--secs-fm", type=float, default=120.0)
    a = ap.parse_args()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    sweep_dir = LAB / f"sweep_{stamp}"
    sweep_dir.mkdir(parents=True, exist_ok=True)
    manifest = dict(stamp=stamp, legs=[])
    man_p = sweep_dir / "manifest.json"
    print(f"[sweep] {len(legs(a))} legs -> {sweep_dir}")
    for leg in legs(a):
        cap = sweep_dir / f"{leg['tag']}.cs16"
        cmd = [PY, IQ, "--fc", str(leg["fc"]), "--rate", str(leg["rate"]), "--secs", str(leg["secs"]),
               "--ant", leg["ant"], "--rfgain", str(leg["rfgain"]), "--ifgr", str(leg["ifgr"]),
               "--biast", "off", "--out", str(cap)]
        t_launch = time.time()
        print(f"[sweep] {datetime.now(timezone.utc):%H:%M:%S}Z {leg['tag']:6s} {leg['fc']/1e6:8.3f} MHz "
              f"{leg['ant']} rfgain {leg['rfgain']} {leg['secs']:.0f} s", flush=True)
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=leg["secs"] + 180)
        t_done = time.time()
        tail = [ln for ln in p.stdout.splitlines() if any(k in ln for k in ("headroom", "VERDICT", "lock", "VOID"))]
        print("        " + " | ".join(ln.strip() for ln in tail)[:200], flush=True)
        meta = dict(leg, capture=str(cap), truth=None, rate_hz=leg["rate"], biast="off",
                    t_launch_unix=t_launch, t_done_unix=t_done,
                    t_launch_utc=datetime.fromtimestamp(t_launch, timezone.utc).isoformat(timespec="milliseconds"),
                    capture_rc=p.returncode, capture_stdout=p.stdout[-3000:], capture_stderr=p.stderr[-1000:],
                    note="observer position deliberately absent; truth attached later by adsb_history.py")
        if "fc" in meta:
            meta["fc_hz"] = meta.pop("fc")
        json.dump(meta, open(sweep_dir / f"{leg['tag']}.meta.json", "w"), indent=1)
        manifest["legs"].append(dict(tag=leg["tag"], kind=leg["kind"], fc_hz=leg["fc"], rate_hz=leg["rate"],
                                     ant=leg["ant"], rfgain=leg["rfgain"], secs=leg["secs"], rc=p.returncode,
                                     capture=str(cap), meta=str(sweep_dir / f"{leg['tag']}.meta.json")))
        json.dump(manifest, open(man_p, "w"), indent=1)
    print(f"[sweep] done -> {man_p}")


if __name__ == "__main__":
    main()
