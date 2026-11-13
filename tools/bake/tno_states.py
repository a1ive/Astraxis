"""Initial states for the N-body integration of assets/scenes/kuiper_belt.toml.

Fetches JPL Horizons geometric state vectors (ICRF, TDB, km and km/s) relative to the
Sun's center (500@10) at one epoch and prints them as [[nbody.members]] tables. The
giant planets are system barycenters; the small bodies use Horizons' current orbit
solution, whose number (JPL#n) is printed next to each one.

Run from the repository root and paste the output into the scene:
    python tools/bake/tno_states.py > members.toml
"""

import json
import sys
import urllib.parse
import urllib.request

EPOCH_JD_TDB = "2461041.5"  # 2026-01-01 00:00 TDB

# (scene name, Horizons COMMAND)
TARGETS = [
    ("Jupiter", "5"),
    ("Saturn", "6"),
    ("Uranus", "7"),
    ("Neptune", "8"),
    ("2001 QR322", "612243;"),
    ("2008 LC18", "DES=2008 LC18;"),
    ("1995 DA2", "15836;"),
    ("Pluto system", "9"),
    ("Orcus", "90482;"),
    ("1996 TR66", "20161;"),
    ("2000 JG81", "130391;"),
    ("1999 DE9", "26375;"),
    ("Gonggong", "225088;"),
    ("Arrokoth", "486958;"),
    ("Quaoar", "50000;"),
    ("Makemake", "136472;"),
    ("Haumea", "136108;"),
    ("Eris", "136199;"),
    ("Sedna", "90377;"),
    ("2012 VP113", "DES=2012 VP113;"),
    ("2026 ED138", "DES=2026 ED138;"),
]


def fetch(command):
    query = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "'500@10'",
        "TLIST": EPOCH_JD_TDB,
        "TLIST_TYPE": "JD",
        "VEC_TABLE": "2",
        "REF_PLANE": "FRAME",
        "OUT_UNITS": "KM-S",
        "CSV_FORMAT": "YES",
    }
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(query)
    with urllib.request.urlopen(url, timeout=120) as response:
        result = json.load(response)["result"]
    if "$$SOE" not in result:
        sys.exit(f"Horizons gave no vectors for {command}:\n{result[:2000]}")
    line = result.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()[0]
    values = [float(x) for x in line.split(",")[2:8]]
    target = next(l for l in result.splitlines() if l.startswith("Target body name"))
    source = target.split("{source:")[1].rstrip("}").strip() if "{source:" in target else ""
    name = target.split(":", 1)[1].split("{")[0].strip()
    return values, name, source


def main():
    print(f"# Horizons states at JD {EPOCH_JD_TDB} TDB, relative to the Sun's center (ICRF).")
    for scene_name, command in TARGETS:
        state, name, source = fetch(command)
        print()
        print("[[nbody.members]]")
        print(f'body = "{scene_name}"')
        if scene_name not in ("Jupiter", "Saturn", "Uranus", "Neptune"):
            print("gm_km3_s2 = 0.0                # a test particle")
        print('relative_to = ["Sun"]')
        print(f"[nbody.members.icrf_state]      # [HZN] {name}, {source}")
        print("position_km = [{:.9e}, {:.9e}, {:.9e}]".format(*state[:3]))
        print("velocity_km_s = [{:.12e}, {:.12e}, {:.12e}]".format(*state[3:]))


if __name__ == "__main__":
    main()
