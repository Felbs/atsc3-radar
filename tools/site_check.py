"""site_check.py -- is this machine ready to run the radar here?

Checks, in order: the site file and observer position (coordinates are
never printed), python packages, SoapySDR and the radios it can see,
`gh` for adsb.lol history downloads, and the station files. With
--capture it records two seconds at an FM frequency and reports headroom,
which proves the whole chain from antenna to .cs16.

Usage: python tools/site_check.py [--capture 99.5e6]
"""
import argparse
import importlib
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import radar_site                                              # noqa: E402


def check(label, ok, detail=""):
    print(f"  [{'ok' if ok else '--'}] {label}" + (f": {detail}" if detail else ""))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", type=float, default=None, help="centre frequency (Hz) for a 2 s test capture")
    a = ap.parse_args()
    good = True

    print("site")
    radar_site.describe()
    try:
        radar_site.observer()
        good &= check("observer position", True, "set (withheld)")
    except SystemExit as e:
        good &= check("observer position", False, str(e))

    print("python", sys.version.split()[0], sys.executable)
    for mod in ("numpy", "scipy", "numba"):
        try:
            m = importlib.import_module(mod)
            good &= check(mod, True, getattr(m, "__version__", ""))
        except ImportError:
            good &= check(mod, False, f"pip install {mod}")

    print("radio")
    cap_py = radar_site.python()
    probe = ("import SoapySDR, json; print(json.dumps([dict(d) for d in SoapySDR.Device.enumerate()]))")
    try:
        r = subprocess.run([cap_py, "-c", probe], capture_output=True, text=True, timeout=60)
        if r.returncode == 0:
            import json
            devs = json.loads(r.stdout.strip().splitlines()[-1])
            names = [d.get("driver", "?") + (f" {d['label']}" if d.get("label") else "") for d in devs]
            good &= check("SoapySDR", True, f"{len(devs)} device(s): " + ("; ".join(names) if names else "none plugged in"))
            want = radar_site.driver()
            if want:
                drv = want.split("driver=")[-1].split(",")[0]
                check(f"configured driver '{drv}' present", any(d.get("driver") == drv for d in devs))
        else:
            good &= check("SoapySDR", False, f"not importable under {cap_py} -- install SoapySDR python bindings "
                          "(radioconda, or your distro's soapysdr-module-* packages) and point site.json \"python\" at it")
    except Exception as e:
        good &= check("SoapySDR", False, f"{type(e).__name__}: {e}")
    check("gh (GitHub CLI, for adsb.lol history)", shutil.which("gh") is not None,
          "" if shutil.which("gh") else "optional; needed by adsb_history.py fetch")

    print("stations")
    files = radar_site.station_files()
    good &= check("station files", len(files) > 0,
                  f"{len(files)}" if files else "none -- python tools/make_stations.py --radius-km 80")
    fm = [st for _, st in files if st["fs_hz"] < 1e6]
    tv = [st for _, st in files if st["fs_hz"] > 1e6]
    print(f"      FM {len(fm)}  TV {len(tv)}  (strongest FM {max((s['erp_kw'] for s in fm), default=0):g} kW, "
          f"TV {max((s['erp_kw'] for s in tv), default=0):g} kW)")

    if a.capture:
        print("capture")
        out = HERE.parent / "lab" / "site_check.cs16"
        cmd = [cap_py, str(HERE / "iq_capture.py"), "--fc", str(a.capture), "--rate", "2e6", "--secs", "2", "--out", str(out)]
        ant = radar_site.antenna_for_hz(a.capture)
        if ant:
            cmd += ["--ant", ant]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        for ln in r.stdout.splitlines():
            if any(k in ln for k in ("readback", "headroom", "integrity", "VERDICT", "quality", "VOID")):
                print("     ", ln.strip())
        if r.returncode != 0 and r.stderr:
            print("     ", r.stderr.strip().splitlines()[-1])
        good &= check("test capture", r.returncode == 0)

    print("READY" if good else "NOT READY -- fix the [--] lines above")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
