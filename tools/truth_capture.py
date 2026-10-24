"""truth_capture.py -- one RF25 capture WITH aircraft truth logged alongside.

The candidate tracks in docs/CANDIDATE_TRACKS.md cannot be judged against
an 8/10 capture: nobody recorded the sky that night. This records both at
once: a raw 6.912 MS/s capture of RF25 through the proven capture tool,
and OpenSky ADS-B state vectors for the airspace around the receiver,
polled for the whole window plus margins on each side. Afterwards
truth_compare.py predicts every known aircraft's bistatic delay and
Doppler from its own telemetry and asks whether our peaks land there.

Single-tenant radio: acquired through radio_lock at priority 50 (lab),
heartbeated, released in a finally. sid_watch (priority 20) skips its
minutes while this holds and resumes after; that behaviour is already
proven live (GPS validation, 8/20 15:5x).

Observer position is read from the GPS campaign's own fix at runtime
(lab_local/fix_result_joint.json -> joint5) and is NEVER printed, logged,
or written anywhere except the gitignored lab/ directory's truth file --
and there only as the OpenSky bounding box, which is coarse (+-0.6 deg).

Usage: python tools/truth_capture.py [--secs 120] [--ant "Antenna A"]
       [--rfgain 4] [--pre 20] [--post 20]
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAB = HERE.parent / "lab"
LAB.mkdir(exist_ok=True)
sys.path.insert(0, r"Z:\src\gr-radiotuna\tools")
import radio_lock                                             # noqa: E402

PY = sys.executable
CAPTURE_TOOL = r"Z:\src\atsc3\tools\atsc3_capture.py"
OBS_FIX = r"Z:\src\numpy-gps\lab_local\fix_result_joint.json"
OPENSKY = "https://opensky-network.org/api/states/all"
POLL_S = 10.0          # anonymous OpenSky resolution is 10 s; faster is waste
BOX_LAT = 0.6          # +-deg, ~65 km
BOX_LON = 0.75


def utc():
    return datetime.now(timezone.utc)


def observer():
    """(lat, lon, alt_m) from our own GPS fix. Not printed. Ever."""
    d = json.load(open(OBS_FIX))
    j = d["joint5"]
    return float(j["lat"]), float(j["lon"]), float(j.get("alt_m", 0.0))


class TruthPoller(threading.Thread):
    def __init__(self, box, out_path):
        super().__init__(daemon=True)
        self.box = box
        self.out = open(out_path, "a")
        self.stop = threading.Event()
        self.n_polls = 0
        self.n_states = 0
        self.errors = 0

    def run(self):
        lamin, lamax, lomin, lomax = self.box
        url = (f"{OPENSKY}?lamin={lamin:.3f}&lamax={lamax:.3f}"
               f"&lomin={lomin:.3f}&lomax={lomax:.3f}")
        while not self.stop.is_set():
            t0 = time.time()
            try:
                req = urllib.request.Request(
                    url, headers={"User-Agent": "atsc3-radar-truth/0.1"})
                with urllib.request.urlopen(req, timeout=8) as r:
                    payload = json.load(r)
                states = payload.get("states") or []
                row = dict(poll_utc=utc().isoformat(timespec="milliseconds"),
                           opensky_time=payload.get("time"),
                           n=len(states), states=states)
                self.n_polls += 1
                self.n_states += len(states)
            except Exception as e:
                self.errors += 1
                row = dict(poll_utc=utc().isoformat(timespec="milliseconds"),
                           error=f"{type(e).__name__}: {e}")
            self.out.write(json.dumps(row) + "\n")
            self.out.flush()
            self.stop.wait(max(0.5, POLL_S - (time.time() - t0)))
        self.out.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rf", type=int, default=25)
    ap.add_argument("--secs", type=float, default=120.0)
    ap.add_argument("--ant", default="Antenna A",
                    help="8/20 port map: A = the Philips TV antenna")
    ap.add_argument("--rfgain", type=int, default=4,
                    help="4 = the unmeasured-carrier default (asymmetry of "
                         "harm, E79); RF25 measured best at 2 on the yagi, "
                         "but the Philips is unmeasured at RF25")
    ap.add_argument("--ifgr", type=int, default=32)
    ap.add_argument("--biast", default="off", choices=("on", "off"),
                    help="power the port-B bias-T (the yagi's LNA). E86: "
                         "bootstrap peak ratio 64.6 -> 81.0 with it on. "
                         "Read back by the capture tool either way.")
    ap.add_argument("--pre", type=float, default=20.0,
                    help="seconds of truth polling BEFORE the capture")
    ap.add_argument("--post", type=float, default=20.0)
    a = ap.parse_args()

    stamp = utc().strftime("%Y%m%d_%H%M%SZ")
    cap = LAB / f"truth_rf{a.rf}_{stamp}.cs16"
    truth = LAB / f"truth_rf{a.rf}_{stamp}.adsb.jsonl"
    meta_p = LAB / f"truth_rf{a.rf}_{stamp}.meta.json"

    lat, lon, alt = observer()
    box = (lat - BOX_LAT, lat + BOX_LAT, lon - BOX_LON, lon + BOX_LON)
    print(f"[truth] OpenSky box +-{BOX_LAT}/{BOX_LON} deg around the fix "
          f"(coordinates withheld) -> {truth.name}")

    poller = TruthPoller(box, truth)
    poller.start()
    print(f"[truth] polling {a.pre:.0f} s before the capture ...")
    time.sleep(a.pre)

    meta = dict(rf=a.rf, secs=a.secs, ant=a.ant, rfgain=a.rfgain,
                ifgr=a.ifgr, biast=a.biast, rate_hz=6.912e6, capture=str(cap),
                truth=str(truth), box_deg=dict(lat=BOX_LAT, lon=BOX_LON),
                note="observer position deliberately absent; see "
                     "numpy-gps lab_local (private)")
    if not radio_lock.acquire("radartuna_truth",
                              f"RF{a.rf} capture + ADS-B truth", 50,
                              wait_s=30):
        h = radio_lock.status() or {}
        print(f"[lock] radio held by {h.get('owner', '?')} - aborting")
        poller.stop.set()
        return 2
    hb_stop = threading.Event()

    def hb():
        while not hb_stop.wait(20.0):
            radio_lock.heartbeat()
    threading.Thread(target=hb, daemon=True).start()

    rc = -1
    try:
        cmd = [PY, CAPTURE_TOOL, "--rf", str(a.rf), "--secs", str(a.secs),
               "--rate", "6912000", "--ant", a.ant, "--ifgr", str(a.ifgr),
               "--rfgain", str(a.rfgain), "--biast", a.biast, "--out", str(cap)]
        meta["t_launch_utc"] = utc().isoformat(timespec="milliseconds")
        meta["t_launch_unix"] = time.time()
        print(f"[cap  ] {utc():%H:%M:%S}Z launching {a.secs:.0f} s on {a.ant} "
              f"rfgain {a.rfgain}")
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=a.secs + 120)       # pipe-ok: tiny log
        rc = p.returncode
        meta["t_done_utc"] = utc().isoformat(timespec="milliseconds")
        meta["t_done_unix"] = time.time()
        meta["capture_rc"] = rc
        meta["capture_stdout"] = p.stdout[-4000:]
        meta["capture_stderr"] = p.stderr[-2000:]
        print(p.stdout[-1500:])
        if rc != 0:
            print(f"[cap  ] capture tool exit {rc}\n{p.stderr[-800:]}")
    finally:
        hb_stop.set()
        radio_lock.release("radartuna_truth")
        print(f"[lock ] released")

    print(f"[truth] polling {a.post:.0f} s after the capture ...")
    time.sleep(a.post)
    poller.stop.set()
    poller.join(timeout=15)
    meta["truth_polls"] = poller.n_polls
    meta["truth_states_total"] = poller.n_states
    meta["truth_errors"] = poller.errors
    if cap.exists():
        nbytes = cap.stat().st_size
        meta["capture_bytes"] = nbytes
        meta["capture_seconds_by_size"] = nbytes / 4 / 6.912e6
    json.dump(meta, open(meta_p, "w"), indent=1)
    print(f"[done ] capture {meta.get('capture_seconds_by_size', 0):.1f} s by "
          f"size, truth {poller.n_polls} polls / {poller.n_states} states "
          f"/ {poller.errors} errors -> {meta_p.name}")
    return 0 if rc == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
