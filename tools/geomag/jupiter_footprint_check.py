"""Independent DOP853 trace checks and CON2020 approximation comparison.

Run after jupiter_footprint_bake.py using the same offline environment.
The field evaluator is shared; the integrator and surface crossing are not.
"""

import json
from pathlib import Path

import JupiterMag as jm
import numpy as np
from scipy.integrate import solve_ivp

from jupiter_footprint_bake import A, B, CON, HEIGHT, MOONS, RADII, trace


def independent(r, lat, lon):
    lat, lon = np.radians((lat, lon))
    p = r * np.array((np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)))
    axes = np.array((A+HEIGHT, A+HEIGHT, B+HEIGHT)) / A

    def surface(s, x):
        return np.sum((x/axes)**2) - 1

    surface.terminal = True
    surface.direction = -1
    out = np.zeros((2, 3))
    for sign in (1, -1):
        def rhs(s, x):
            field = (np.array(jm.Internal.Field(*x, MaxDeg=13)) + np.array(jm.Con2020.Field(*x))).ravel()
            return sign * field / np.linalg.norm(field)
        result = solve_ivp(rhs, (0, 5*r), p, method="DOP853", max_step=.025,
                           rtol=2e-10, atol=1e-11, events=surface)
        if len(result.y_events[0]) != 1:
            raise RuntimeError("independent trace missed the ionosphere")
        q = result.y_events[0][0]
        out[0 if q[2] > 0 else 1] = q / np.linalg.norm(q)
    return out


def main():
    if jm.__version__ != "1.5.1":
        raise RuntimeError("this check requires JupiterMag 1.5.1")
    jm.Internal.Config(Model="jrm33", Degree=13)
    jm.Con2020.Config(**CON)
    # C++ integration-test fixtures: input 0 E, 0.37 N, nominal radius.
    # At the seam the periodic continuation averages the two endpoint fits.
    leads = ((7.1, 2.47),
             ((9.07 + 5.78 + 3.29*np.cos(.99*2*np.pi) + 1.41*np.sin(.99*2*np.pi))/2, 3.69),
             ((18.61 + 12.45 + 6.16*np.cos(.98*2*np.pi) + 3.64*np.sin(.98*2*np.pi))/2,
              (7.58 + 12.97 - 5.39*np.cos(.99*2*np.pi) - 4.50*np.sin(.99*2*np.pi))/2))
    for moon, r, lead in zip(MOONS, RADII, leads):
        print("lead fixture", moon,
              [independent(r, .37, lead[h])[h].tolist() for h in range(2)], flush=True)
    report = []
    for moon, r in zip(MOONS, RADII):
        for lon in (0.0, 93.4, 241.0, 359.9):
            lat = .37
            native = trace(r, lat, lon)[0]
            oracle = independent(r, lat, lon)
            error = np.linalg.norm(native-oracle, axis=1) * A
            if error.max() > 10:
                raise RuntimeError(f"{moon} {lon}: independent trace error {error}")
            # Hybrid is much more expensive; retain one longitude per moon
            # as an explicit estimate of the analytic field approximation.
            hybrid_error = None
            if lon == 0:
                jm.Con2020.Config(equation_type="hybrid")
                hybrid = trace(r, lat, lon, step=.04)[0]
                jm.Con2020.Config(equation_type="analytic")
                hybrid_error = (np.linalg.norm(native-hybrid, axis=1)*A).tolist()
            row = dict(moon=moon, radius_rj=float(r), latitude_deg=lat,
                       east_longitude_deg=lon, directions=oracle.tolist(),
                       independent_error_km=error.tolist(), analytic_hybrid_difference_km=hybrid_error)
            report.append(row)
            print(json.dumps(row), flush=True)
    Path("build/jovian_footprints/check.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
