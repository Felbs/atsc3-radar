"""radar_site.py -- where the receiver is and what it is, for every tool.

One file, outside the repo, says what is specific to a site:

    lab/site.json            (lab/ is gitignored; or $BROADCAST_RADAR_SITE)
    {
      "observer": {"lat": 0.0, "lon": 0.0, "alt_m": 0.0},
                   -- or {"from_file": "<json>", "key": "<optional sub-key>"}
                      with lat / lon / alt_m inside, if a GPS fix lives elsewhere
      "sdr": {"driver": "sdrplay",            -- SoapySDR driver, or full "driver=...,serial=..." args
              "antennas": {"fm": "Antenna C", "tv": "Antenna B"}},
      "python": "<python with SoapySDR>",    -- optional; captures run under it
      "radio_lock_dir": "<dir with radio_lock.py>",   -- optional single-tenant radio lock
      "stations_dirs": ["lab/stations"]      -- optional extra station-file dirs
    }

The observer position is read from here at run time and is never written
into the repo, into a capture's sidecar, or into a doc. Tools that print
geometry print "coordinates withheld" in its place.

`python tools/radar_site.py` prints the configuration with the position
withheld. site.example.json at the repo root is a template.
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

_cache = None


def path():
    for p in (os.environ.get("BROADCAST_RADAR_SITE"), ROOT / "lab" / "site.json", ROOT / "site.json"):
        if p and Path(p).exists():
            return Path(p)
    return None


def load():
    global _cache
    if _cache is None:
        p = path()
        _cache = json.load(open(p)) if p else {}
    return _cache


def observer():
    """(lat, lon, alt_m) of the receiver. Exits with instructions if unset."""
    obs = load().get("observer")
    if obs and obs.get("from_file"):
        d = json.load(open(obs["from_file"]))
        if obs.get("key"):
            d = d[obs["key"]]
        return float(d["lat"]), float(d["lon"]), float(d.get("alt_m", 0.0))
    if obs and "lat" in obs and "lon" in obs:
        return float(obs["lat"]), float(obs["lon"]), float(obs.get("alt_m", 0.0))
    raise SystemExit("no observer position: copy site.example.json to lab/site.json and fill in "
                     "\"observer\" (or set BROADCAST_RADAR_SITE to such a file). See tools/radar_site.py.")


def sdr():
    return load().get("sdr", {})


def driver():
    return sdr().get("driver")


def antenna(kind):
    """Antenna port name for 'fm' or 'tv', or None (the device default)."""
    return sdr().get("antennas", {}).get(kind)


def antenna_for_hz(fc):
    return antenna("fm" if fc < 174e6 else "tv")


def python():
    p = load().get("python")
    return p if p and os.path.exists(p) else sys.executable


def radio_lock():
    """The fleet's radio_lock module if configured and present, else a stub that always grants."""
    d = load().get("radio_lock_dir")
    if d and Path(d).exists():
        sys.path.insert(0, str(d))
        try:
            import radio_lock                        # noqa: F401
            return radio_lock
        except ImportError:
            pass

    class _NoLock:
        configured = False

        @staticmethod
        def acquire(*a, **k):
            return True

        @staticmethod
        def release(*a, **k):
            return None

        @staticmethod
        def status():
            return None

        @staticmethod
        def heartbeat():
            return None
    return _NoLock


def stations_dirs():
    out = [HERE / "stations", HERE / "stations" / "fm", ROOT / "lab" / "stations"]
    for d in load().get("stations_dirs", []):
        out.append(Path(d) if os.path.isabs(d) else ROOT / d)
    return [d for d in out if d.exists()]


def station_files():
    seen, out = set(), []
    for d in stations_dirs():
        for p in sorted(d.glob("*.json")):
            if p.name.startswith("_"):
                continue
            try:
                st = json.load(open(p))
            except Exception:
                continue
            if "tx_lat" in st and "fc_hz" in st and p.resolve() not in seen:
                seen.add(p.resolve())
                out.append((p, st))
    return out


def truth_box(lat, lon, dlat=0.6, dlon=0.75):
    """A lat/lon box covering the receiver AND every known tower, each padded.

    The bistatic ellipse for a far echo reaches well past the tower (10/24:
    an aircraft 15 km beyond WNUV, 110 km from here, was a strong 31 km
    echo), so truth must cover the tower end too."""
    lats, lons = [lat], [lon]
    for _, st in station_files():
        lats.append(float(st["tx_lat"])); lons.append(float(st["tx_lon"]))
    return (min(lats) - dlat, max(lats) + dlat, min(lons) - dlon, max(lons) + dlon)


def describe():
    p = path()
    print(f"site file: {p if p else 'NONE (copy site.example.json to lab/site.json)'}")
    try:
        observer()
        src = load().get("observer", {})
        print("observer: " + (f"from {src['from_file']}" if src.get("from_file") else "in site file") + " (coordinates withheld)")
    except SystemExit as e:
        print(f"observer: NOT SET -- {e}")
    print(f"sdr driver: {driver() or '(SoapySDR default device)'}   antennas: {sdr().get('antennas') or '(device default)'}")
    print(f"capture python: {python()}")
    print(f"radio lock: {'configured' if radio_lock().__name__ == 'radio_lock' else 'none (single user)'}")
    files = station_files()
    kinds = {}
    for _, st in files:
        k = "ATSC 3.0" if abs(st["fs_hz"] - 6.912e6) < 1 else ("8-VSB" if st["fs_hz"] > 1e6 else "FM")
        kinds[k] = kinds.get(k, 0) + 1
    print(f"stations: {len(files)} in {', '.join(str(d) for d in stations_dirs())}  " + "  ".join(f"{k} {v}" for k, v in sorted(kinds.items())))


if __name__ == "__main__":
    describe()
