#!/usr/bin/env python3
"""The opacity and phase function of PDS 70's dust by Mie theory.

Portilla-Revelo et al. 2023 (A&A 677, A76, arXiv:2306.16850, Table 1) give the
grains of their model: 0.001-3000 um, n(a) ~ a^-3.9, by volume 60.8% silicate,
19.2% amorphous carbon, 20% vacuum. Their dust is the DIANA standard of Woitke et
al. 2016 (A&A 586, A103), with Mg0.7Fe0.3SiO3 from Dorschner et al. 1995 and
amorphous carbon BE from Zubko et al. 1996. DIANA uses hollow spheres; this is Mie
theory on compact spheres of the Bruggeman mixture, which keeps the backward glory
hollow spheres smooth away. The wavelength is that of the polarized image
[PR23] fit, 1.25 um.

Prints the numbers of [disk] in assets/scenes/pds70.toml: the mass extinction, the
albedo, the asymmetry parameter and three Henyey-Greenstein functions fitted to the
phase function.

The optical constants are optool's copies (github.com/cdominik/optool, lnk_data):
  curl -O https://raw.githubusercontent.com/cdominik/optool/master/lnk_data/pyr-mg70-Dorschner1995.lnk
  curl -O https://raw.githubusercontent.com/cdominik/optool/master/lnk_data/c-z-Zubko1996.lnk

Needs numpy: python -E -P tools/pds70/mie_opacity.py <directory with the .lnk files>
Takes a few seconds.
"""

import os
import sys

import numpy as np

WAVELENGTH_UM = 1.25
A_MIN_UM, A_MAX_UM, Q = 0.001, 3000.0, 3.9
SHARE_SILICATE, SHARE_CARBON, SHARE_VACUUM = 0.608, 0.192, 0.20
SIZES = 500


def read_lnk(path):
    """lambda, n, k and the bulk density from an optool .lnk file."""
    rows = []
    density = None
    with open(path, encoding="utf8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if density is None:
                density = float(parts[1])
                continue
            rows.append([float(v) for v in parts[:3]])
    table = np.array(rows)
    return table[:, 0], table[:, 1], table[:, 2], density


def refractive_index(path, wavelength_um):
    lam, n, k, density = read_lnk(path)
    order = np.argsort(lam)
    x = np.log(wavelength_um)
    return (
        complex(
            np.interp(x, np.log(lam[order]), n[order]),
            np.exp(np.interp(x, np.log(lam[order]), np.log(np.maximum(k[order], 1e-12)))),
        ),
        density,
    )


def bruggeman(m_values, shares):
    """The Bruggeman mixture's refractive index (sum f (e - E) / (e + 2 E) = 0), by Newton."""
    eps = np.array([m * m for m in m_values])
    f = np.array(shares)
    mean = np.sum(f * eps)
    for _ in range(200):
        g = np.sum(f * (eps - mean) / (eps + 2 * mean))
        dg = np.sum(f * (-(eps + 2 * mean) - 2 * (eps - mean)) / (eps + 2 * mean) ** 2)
        step = g / dg
        mean -= step
        if abs(step) < 1e-13:
            break
    return np.sqrt(mean)


def mie_coefficients(m, x):
    """a_n, b_n of a sphere of index m (positive imaginary part) and size parameter x."""
    nmax = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    mx = m * x
    start = int(max(nmax, abs(mx))) + 16
    d = np.zeros(start + 1, complex)
    for n in range(start, 0, -1):
        d[n - 1] = n / mx - 1.0 / (d[n] + n / mx)
    d = d[: nmax + 1]
    # Riccati-Bessel functions, upward.
    psi = np.zeros(nmax + 1)
    chi = np.zeros(nmax + 1)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    psi[0], chi[0] = psi1, chi1
    for n in range(1, nmax + 1):
        k = (2 * n - 1) / x
        psi2 = k * psi1 - psi0
        chi2 = k * chi1 - chi0
        psi[n], chi[n] = psi2, chi2
        psi0, psi1 = psi1, psi2
        chi0, chi1 = chi1, chi2
    xi = psi - 1j * chi
    n = np.arange(1, nmax + 1)
    a_term = d[1:] / m + n / x
    b_term = d[1:] * m + n / x
    an = (a_term * psi[1:] - psi[:-1]) / (a_term * xi[1:] - xi[:-1])
    bn = (b_term * psi[1:] - psi[:-1]) / (b_term * xi[1:] - xi[:-1])
    return an, bn


def mie(m, x, mu):
    """Qext, Qsca, g and the scattered intensity i(mu) = (|S1|^2 + |S2|^2) / 2."""
    an, bn = mie_coefficients(m, x)
    n = np.arange(1, len(an) + 1)
    weight = 2 * n + 1
    q_ext = 2.0 / x**2 * np.sum(weight * (an + bn).real)
    q_sca = 2.0 / x**2 * np.sum(weight * (np.abs(an) ** 2 + np.abs(bn) ** 2))
    term = n[:-1] * (n[:-1] + 2.0) / (n[:-1] + 1.0) * (np.conj(an[:-1]) * an[1:] + np.conj(bn[:-1]) * bn[1:]).real
    cross = (2 * n + 1) / (n * (n + 1.0)) * (an * np.conj(bn)).real
    g = 4.0 / (x**2 * q_sca) * (np.sum(term) + np.sum(cross))

    a = weight / (n * (n + 1.0)) * an
    b = weight / (n * (n + 1.0)) * bn
    s1 = np.zeros_like(mu, complex)
    s2 = np.zeros_like(mu, complex)
    pi0 = np.zeros_like(mu)
    pi1 = np.ones_like(mu)
    for i in range(len(n)):
        k = n[i]
        tau = k * mu * pi1 - (k + 1) * pi0
        s1 += a[i] * pi1 + b[i] * tau
        s2 += a[i] * tau + b[i] * pi1
        pi0, pi1 = pi1, ((2 * k + 1) * mu * pi1 - (k + 1) * pi0) / k
    return q_ext, q_sca, g, 0.5 * (np.abs(s1) ** 2 + np.abs(s2) ** 2)


def self_check():
    """Published values: Wiscombe 1980 (Appl. Opt. 19, 1505), m = 1.5, nonabsorbing."""
    for x, q in ((10.0, 2.8820), (100.0, 2.0944), (1000.0, 2.0139)):
        got = mie_coefficients(1.5 + 0j, x)
        n = np.arange(1, len(got[0]) + 1)
        value = 2.0 / x**2 * np.sum((2 * n + 1) * (got[0] + got[1]).real)
        if abs(value - q) > 1e-3:
            raise SystemExit(f"Mie check failed at x = {x}: {value} against {q}")
    # Rayleigh's small grains: Qsca = 8/3 x^4 |(m^2-1)/(m^2+2)|^2.
    m = 1.5 + 0.0j
    x = 0.01
    q_ext, q_sca, g, _ = mie(m, x, np.array([1.0]))
    expected = 8.0 / 3.0 * x**4 * abs((m * m - 1) / (m * m + 2)) ** 2
    if abs(q_sca / expected - 1.0) > 1e-3 or abs(g) > 1e-3:
        raise SystemExit(f"Rayleigh check failed: {q_sca} against {expected}, g = {g}")


def henyey_greenstein(g, mu):
    d = 1.0 + g * g - 2.0 * g * mu
    return (1.0 - g * g) / (d * np.sqrt(d))


def fit_phase(mu, phase, solid):
    """Three Henyey-Greenstein functions with positive weights summing to 1, in ln P."""
    best = None
    g1s = np.arange(0.40, 0.996, 0.02)
    g2s = np.arange(-0.2, 0.72, 0.04)
    g3s = np.arange(-0.8, 0.2, 0.05)
    w = np.sqrt(solid)  # relative errors, by solid angle
    for g1 in g1s:
        h1 = henyey_greenstein(g1, mu) / phase
        for g2 in g2s:
            if g2 >= g1:
                continue
            h2 = henyey_greenstein(g2, mu) / phase
            for g3 in g3s:
                if g3 >= g2:
                    continue
                h3 = henyey_greenstein(g3, mu) / phase
                # Weights summing to 1: two unknowns.
                a = np.vstack([(h1 - h3) * w, (h2 - h3) * w]).T
                sol, *_ = np.linalg.lstsq(a, (1.0 - h3) * w, rcond=None)
                w1, w2 = sol
                w3 = 1.0 - w1 - w2
                if w1 < 0 or w2 < 0 or w3 < 0:
                    continue
                model = w1 * h1 + w2 * h2 + w3 * h3
                error = np.sum(solid * (model - 1.0) ** 2)
                if best is None or error < best[0]:
                    best = (error, (g1, g2, g3), (w1, w2, w3))
    return best


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    self_check()
    m_si, rho_si = refractive_index(os.path.join(sys.argv[1], "pyr-mg70-Dorschner1995.lnk"), WAVELENGTH_UM)
    m_c, rho_c = refractive_index(os.path.join(sys.argv[1], "c-z-Zubko1996.lnk"), WAVELENGTH_UM)
    m = bruggeman([m_si, m_c, 1.0 + 0j], [SHARE_SILICATE, SHARE_CARBON, SHARE_VACUUM])
    density = SHARE_SILICATE * rho_si + SHARE_CARBON * rho_c
    print(f"silicate {m_si:.4f}, carbon {m_c:.4f}, mixture {m:.4f}, {density:.3f} g/cm^3")

    # Scattering angles: fine near the forward direction, where large grains diffract.
    degrees = np.concatenate([np.geomspace(0.01, 5.0, 120, endpoint=False), np.arange(5.0, 180.01, 1.0)])
    mu = np.cos(np.radians(degrees))
    edges = np.concatenate([[0.0], 0.5 * (degrees[1:] + degrees[:-1]), [180.0]])
    solid = 0.5 * (np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:])))

    radii = np.geomspace(A_MIN_UM, A_MAX_UM, SIZES)
    d_ln = np.log(radii[1] / radii[0])
    k = 2.0 * np.pi / WAVELENGTH_UM
    ext = sca = g_sum = mass = 0.0
    intensity = np.zeros_like(mu)
    share = np.zeros(SIZES)
    for i, a in enumerate(radii):
        x = k * a
        q_ext, q_sca, g, i_mu = mie(m, x, mu)
        number = a ** (1.0 - Q) * d_ln  # n(a) da, in ln a
        area = np.pi * a * a
        ext += area * q_ext * number
        sca += area * q_sca * number
        g_sum += area * q_sca * g * number
        mass += 4.0 / 3.0 * np.pi * a**3 * number
        share[i] = area * q_ext * number
        intensity += i_mu / (k * k) * number
    kappa_ext = ext / mass / (density * 1e-4)  # um^2 / um^3 / (g/cm^3), to cm^2/g
    albedo = sca / ext
    g = g_sum / sca
    # Mean 1 over the sphere: the program's convention.
    phase = 4.0 * np.pi * intensity / sca
    norm = np.sum(solid * phase)
    print(f"phase function's mean over the sphere before normalizing: {norm:.4f}")
    phase = phase / norm

    total = np.cumsum(share) / np.sum(share)
    print(
        "extinction by size: half from below "
        f"{radii[np.searchsorted(total, 0.5)]:.3f} um, 10-90% from "
        f"{radii[np.searchsorted(total, 0.1)]:.3f} to {radii[np.searchsorted(total, 0.9)]:.2f} um"
    )
    print(f"mass extinction {kappa_ext:.0f} cm^2/g, albedo {albedo:.4f}, g {g:.4f}")

    # The fit leaves out the diffraction spike under 1 degree: no line of sight in
    # the scene's views scatters there.
    use = degrees >= 1.0
    error, gs, ws = fit_phase(mu[use], phase[use], solid[use])
    print(f"phase_g = [{gs[0]:.3f}, {gs[1]:.3f}, {gs[2]:.3f}]")
    print(f"phase_w = [{ws[0]:.4f}, {ws[1]:.4f}, {ws[2]:.4f}]")
    model = sum(w * henyey_greenstein(gv, mu) for gv, w in zip(gs, ws))
    print("angle  Mie phase  fit")
    for angle in (1, 5, 10, 20, 38, 50, 70, 90, 110, 130, 142, 160, 175, 180):
        j = np.argmin(np.abs(degrees - angle))
        print(f"{angle:5.0f} {phase[j]:9.4f} {model[j]:9.4f}")
    print("fit's mean cosine", sum(w * gv for gv, w in zip(gs, ws)))


if __name__ == "__main__":
    main()
