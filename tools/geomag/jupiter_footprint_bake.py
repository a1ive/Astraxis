"""Bake System III satellite footpoints; see assets/geomag/SOURCES.md.

Requires JupiterMag==1.5.1 (offline only). Cartesian entries are unit
planetocentric directions at the 900 km spheroid, before the lead angle.
"""

import argparse
import hashlib
import json
from pathlib import Path

import JupiterMag as jm
import numpy as np

A, B, HEIGHT = 71492.0, 66854.0, 900.0
N = 360
MOONS = ("Io", "Europa", "Ganymede")
RADII = np.array((421800.0, 671100.0, 1070400.0)) / A
CON = dict(mu_i=139.6, i_rho=16.7, r0=7.8, r1=51.4, d=3.6,
           xt=9.3, xp=155.8, Edwards=True, Smooth=False,
           equation_type="analytic", azfunc="connerney")


def trace(r, lat, lon, step=0.02):
    r, lat, lon = np.broadcast_arrays(r, lat, lon)
    lat, lon = np.radians(lat.ravel()), np.radians(lon.ravel())
    r = r.ravel()
    max_len = int(8 * r.max() / step) + 1000
    t = jm.TraceField(r * np.cos(lat) * np.cos(lon),
                      r * np.cos(lat) * np.sin(lon), r * np.sin(lat),
                      IntModel="jrm33", ExtModel="Con2020", MaxLen=max_len,
                      MaxStep=step, InitStep=step, MinStep=1e-6,
                      ErrMax=1e-8, ionosphereShape="spheroid",
                      ai=(A + HEIGHT) / A, bi=(B + HEIGHT) / A)
    # Intersect each trace segment with the requested spheroid ourselves.
    # JupiterMag 1.5.1's footprint interpolation sometimes returns NaNs or
    # non-convergent coordinates even for complete, smooth traces.
    axes = np.array(((A + HEIGHT) / A, (A + HEIGHT) / A, (B + HEIGHT) / A))
    points = np.empty((len(r), 2, 3))
    for i, count in enumerate(t.nstep):
        p = np.column_stack((t.x[i, :count], t.y[i, :count], t.z[i, :count]))
        level = np.sum((p / axes)**2, axis=1) - 1
        indices = np.flatnonzero(level[:-1] * level[1:] <= 0)
        if len(indices) != 2:
            raise RuntimeError(f"expected two shell crossings, got {len(indices)}")
        for k in indices:
            lo, hi = 0.0, 1.0
            for _ in range(40):
                mid = (lo + hi) / 2
                v = p[k] + mid * (p[k+1] - p[k])
                if (np.sum((v / axes)**2)-1) * level[k] > 0:
                    lo = mid
                else:
                    hi = mid
            points[i, 0 if v[2] > 0 else 1] = v
    if not np.all(np.isfinite(points)) or np.any(t.nstep >= max_len):
        raise RuntimeError("incomplete magnetic trace")
    shell = np.sum((points / np.array(((A + HEIGHT) / A,
                                      (A + HEIGHT) / A,
                                      (B + HEIGHT) / A))) ** 2, axis=-1)
    if np.max(np.abs(shell - 1)) > 1e-5:
        raise RuntimeError("footpoint is not on the requested spheroid")
    if np.any(points[:, 0, 2] <= 0) or np.any(points[:, 1, 2] >= 0):
        raise RuntimeError("invalid hemispheres")
    return points / np.linalg.norm(points, axis=-1, keepdims=True)


def interpolate(grid, r, lat, lon, nominal):
    x = np.stack(((r / nominal - 0.97) / 0.03, (lat + 2) / 2,
                  (lon % 360) * N / 360), axis=-1)
    i = np.floor(x).astype(int)
    i[:, :2] = np.minimum(i[:, :2], 1)
    f = x - i
    out = np.zeros((len(r), 2, 3))
    for dr in range(2):
        for dl in range(2):
            for dp in range(2):
                w = ((f[:, 0] if dr else 1-f[:, 0]) *
                     (f[:, 1] if dl else 1-f[:, 1]) *
                     (f[:, 2] if dp else 1-f[:, 2]))
                out += w[:, None, None] * grid[i[:, 0]+dr, i[:, 1]+dl, (i[:, 2]+dp) % N]
    return out / np.linalg.norm(out, axis=-1, keepdims=True)


def main():
    if jm.__version__ != "1.5.1":
        raise RuntimeError("this bake requires JupiterMag 1.5.1")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("assets/geomag/jupiter_footprints.txt"))
    parser.add_argument("--report", type=Path, default=Path("build/jovian_footprints/bake.json"))
    args = parser.parse_args()
    jm.Internal.Config(Model="jrm33", Degree=13)
    jm.Con2020.Config(**CON)
    report = {"JupiterMag": "1.5.1", "degree": 13, "con2020": CON, "moons": {}}
    lines = [f"AXJOVIANFP1 {N} 3 3"]
    rng = np.random.default_rng(20261103)
    for name, nominal in zip(MOONS, RADII):
        grid = np.empty((3, 3, N, 2, 3))
        for ri, r in enumerate(nominal * np.array((0.97, 1, 1.03))):
            for li, lat in enumerate((-2.0, 0.0, 2.0)):
                grid[ri, li] = trace(r, lat, np.arange(N) * 360 / N)
                print(name, ri, li, flush=True)
        # Withheld off-grid traces, four times finer integration, including
        # the longitude seam. This measures interpolation, not model accuracy.
        r = nominal * rng.uniform(.97, 1.03, 120)
        lat = rng.uniform(-2, 2, 120)
        lon = rng.uniform(0, 360, 120)
        lon[:4] = (0.001, 359.999, 0.49, 359.51)
        exact = trace(r, lat, lon, step=.005)
        approx = interpolate(grid, r, lat, lon, nominal)
        error = np.linalg.norm(exact-approx, axis=-1) * A
        fixed = interpolate(grid, np.full_like(r, nominal), np.zeros_like(lat), lon, nominal)
        omitted = np.linalg.norm(exact-fixed, axis=-1) * A
        report["moons"][name] = {"max_lookup_error_km": error.max(axis=0).tolist(),
                                   "max_fixed_orbit_error_km": omitted.max(axis=0).tolist()}
        if np.max(error) > 30:
            raise RuntimeError(f"{name}: interpolation error {np.max(error):.2f} km")
        lines.append(f"{name} {nominal*.97:.12f} {nominal*1.03:.12f} -2 2")
        lines.extend(" ".join(f"{v:.10f}" for v in row) for row in grid.reshape(-1, 6))
        print(report["moons"][name], flush=True)
    payload = ("\n".join(lines) + "\n").encode("ascii")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(payload)
    report["sha256"] = hashlib.sha256(payload).hexdigest()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
