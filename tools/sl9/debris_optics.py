"""What the dark debris of comet Shoemaker-Levy 9 does to Jupiter's light: prints the constants of
shaders/impact_scar.frag.hlsl and src/render/impact_pass.cpp, and checks them against [WEST] Fig. 3.

Usage (needs numpy, a few seconds):
    python -E -P tools/sl9/debris_optics.py

[WEST] R. A. West 1996, in "The Collision of Comet Shoemaker-Levy 9 and Jupiter" (Cambridge), pp. 269-292.
  Fig. 4: the imaginary refractive index k, read off the plot: 0.03 at 265 nm, 0.017 at 337, 0.0125 at 408, 0.009
  at 545, 0.0065 at 943 (these from HST, West et al. 1995). Sec. 3: real part 1.44, mean radius 0.21 to 0.28 um,
  growing over the first month; optical depth near 4 at 889 nm in the core of a fresh site, over 12 at 275 nm.
  Fig. 3: a cut across G's core, about 2 hours old. I/F around it and in it: 955 nm 0.36 and 0.23, 275 nm 0.25 and
  0.14. Around the outlying crescent ("PT"): 955 nm 0.38 and 0.32. (893 nm is in a methane band, not modeled here.)

The layer lies above the clouds: a two-stream (delta-Eddington) layer of the Mie grains over a Lambert surface of
the albedo that gives the I/F around the site. k is scaled, within about West's factor of 2, to give his two I/F in
the core, and interpolated between them. The scene's textures aren't calibrated in I/F, so it uses a factor on the
clouds, for each of red, green, blue, with clouds of I/F 0.5 as Jupiter's at these latitudes in visible light:
    factor = T + rho (1 - T),   T = exp(-a tau m),  m = 1 / mu + 1 / mu0,  tau at 890 nm:
thin debris lets the clouds through, thick debris shows only itself. The constants make that equal the two-stream
layer in G's core and crescent at mu0 = 0.77; away from there it is only roughly right (the layer's own brightness
goes differently with mu0). a is the same for all three, the green's.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "pds70"))
from mie_opacity import mie  # noqa: E402

N_REAL = 1.44
RADIUS_UM = 0.25
K_TABLE = np.array([[265.0, 0.03], [337.0, 0.017], [408.0, 0.0125], [545.0, 0.009], [943.0, 0.0065]])
TAU_REF_NM = 890.0
CHANNELS = {"red": 620.0, "green": 540.0, "blue": 460.0}
CORE_TAU_890 = 3.4   # Sanchez-Lavega et al. 1998 Table IV, as impact.cpp
CRESCENT_TAU_890 = 0.8  # Fits West's crescent at both 955 and 275 nm, printed below
MU = 0.75
MU0 = 0.77


# West gives k to about a factor of 2. These factors bring the layer to his I/F in G's core at 955 and 275 nm
# (found by `calibrate`), and are interpolated in between.
K_SCALE = {}


def k_at(nm):
    k = np.exp(np.interp(np.log(nm), np.log(K_TABLE[:, 0]), np.log(K_TABLE[:, 1])))
    if K_SCALE:
        scale = np.interp(np.log(nm), np.log([275.0, 955.0]), np.log([K_SCALE[275.0], K_SCALE[955.0]]))
        k *= np.exp(scale)
    return k


def grain(nm):
    """Extinction efficiency, single-scattering albedo, asymmetry, over a narrow spread of radii."""
    radii = RADIUS_UM * np.exp(np.linspace(-0.36, 0.36, 9))   # log-normal, sigma 0.18 [SL98]
    weights = np.exp(-0.5 * (np.log(radii / RADIUS_UM) / 0.18) ** 2)
    ext = sca = sg = 0.0
    for r, w in zip(radii, weights):
        x = 2.0 * np.pi * r * 1000.0 / nm
        q_ext, q_sca, g, _ = mie(complex(N_REAL, k_at(nm)), x, np.array([1.0]))
        area = w * r * r
        ext += area * q_ext
        sca += area * q_sca
        sg += area * q_sca * g
    return ext, sca / ext, sg / sca


def two_stream(tau, omega, g, mu0, surface):
    """Reflectance over a Lambert surface: delta-Eddington, Meador & Weaver 1980; the I/F is taken as A mu0."""
    f = g * g
    tau_d = (1.0 - omega * f) * tau
    om = (1.0 - f) * omega / (1.0 - omega * f)
    gd = (g - f) / (1.0 - f)
    g1 = (7.0 - om * (4.0 + 3.0 * gd)) / 4.0
    g2 = -(1.0 - om * (4.0 - 3.0 * gd)) / 4.0
    g3 = (2.0 - 3.0 * gd * mu0) / 4.0
    g4 = 1.0 - g3
    kk = np.sqrt(max(g1 * g1 - g2 * g2, 1e-12))
    alpha1 = g1 * g4 + g2 * g3
    alpha2 = g1 * g3 + g2 * g4
    e_p, e_m = np.exp(kk * tau_d), np.exp(-kk * tau_d)
    e0 = np.exp(-tau_d / mu0)
    den = (1.0 - (kk * mu0) ** 2) * ((kk + g1) * e_p + (kk - g1) * e_m)
    r = om * ((1.0 - kk * mu0) * (alpha2 + kk * g3) * e_p - (1.0 + kk * mu0) * (alpha2 - kk * g3) * e_m
              - 2.0 * kk * (g3 - alpha2 * mu0) * e0) / den
    t_diff = -om * e0 * ((1.0 + kk * mu0) * (alpha1 + kk * g4) * e_p - (1.0 - kk * mu0) * (alpha1 - kk * g4) * e_m
                         - 2.0 * kk * (g4 + alpha1 * mu0) / e0) / den
    t_dir = e0
    # Diffuse transmission and reflection for light coming up from the surface: the conservative form.
    gamma = g2 / (g1 + kk)
    td_up = 2.0 * kk / ((kk + g1) * e_p + (kk - g1) * e_m)
    rd_up = g2 * (e_p - e_m) / ((kk + g1) * e_p + (kk - g1) * e_m)
    _ = gamma
    return r + (t_dir + t_diff) * surface * td_up / (1.0 - surface * rd_up)


def scene_output(nm, ref, tau_890, mu0, cloud_if):
    """The I/F over clouds of this I/F."""
    ext, omega, g = grain(nm)
    tau = tau_890 * ext / ref
    return mu0 * two_stream(tau, omega, g, mu0, cloud_if / mu0)


def calibrate(ref_nm=TAU_REF_NM):
    """The factor on k that gives West's I/F in G's core."""
    for nm, around, measured in ((955.0, 0.36, 0.23), (275.0, 0.25, 0.14)):
        lo, hi = 0.1, 8.0
        for _ in range(30):
            mid = np.sqrt(lo * hi)
            K_SCALE.clear()
            K_SCALE.update({275.0: mid, 955.0: mid})
            ref = grain(ref_nm)[0]
            if scene_output(nm, ref, CORE_TAU_890, MU0, around) > measured:
                lo = mid
            else:
                hi = mid
        yield nm, np.sqrt(lo * hi)


def anchored(nm, ref, a=None):
    """a and rho of the scene's factor that give the layer's core and crescent, over clouds of I/F 0.5, at mu0 = 0.77."""
    core = scene_output(nm, ref, CORE_TAU_890, MU0, 0.5) / 0.5
    thin = scene_output(nm, ref, CRESCENT_TAU_890, MU0, 0.5) / 0.5
    m = 2.0 / MU0
    best = None
    for try_a in ([a] if a is not None else np.arange(0.02, 1.5, 0.005)):
        t_core = np.exp(-try_a * CORE_TAU_890 * m)
        t_thin = np.exp(-try_a * CRESCENT_TAU_890 * m)
        rho = (core - t_core) / (1.0 - t_core)
        err = abs(t_thin + rho * (1.0 - t_thin) - thin)
        if best is None or err < best[0]:
            best = (err, try_a, rho, core, thin)
    return best[1:]


def main():
    scales = dict(calibrate())
    K_SCALE.clear()
    K_SCALE.update(scales)
    print(f"k scaled by {scales[955.0]:.2f} at 955 nm and {scales[275.0]:.2f} at 275 nm (West: uncertain by about 2)")
    ref, _, _ = grain(TAU_REF_NM)
    print("nm     k        tau/tau890  omega   g")
    for nm in (275.0, 410.0, 460.0, 540.0, 620.0, 890.0, 955.0):
        ext, omega, g = grain(nm)
        print(f"{nm:5.0f}  {k_at(nm):.4f}   {ext / ref:6.2f}     {omega:.4f}  {g:.3f}")
    print(f"tau 275 / tau 890 = {grain(275.0)[0] / ref:.2f} (West: 6.4 on 23 July, 3.1 on 30 July, 2.0 on 24 August)")

    print("\nthe core against West Fig. 3 (I/F around it -> in it):")
    for nm, around, measured in ((955.0, 0.36, 0.23), (275.0, 0.25, 0.14)):
        print(f"  {nm:.0f} nm, tau890 {CORE_TAU_890}: {around} -> {scene_output(nm, ref, CORE_TAU_890, MU0, around):.3f}  (measured {measured})")
    print("the crescent, measured at 955 nm 0.38 -> 0.32, at 275 nm about 0.26 -> 0.17:")
    for tau in (0.2, 0.4, 0.6, 1.0, 1.5):
        print(f"  tau890 {tau}: {scene_output(955.0, ref, tau, MU0, 0.38):.3f}, {scene_output(275.0, ref, tau, MU0, 0.26):.3f}")

    print("\nthe ratio to clouds of I/F 0.5 in visible light, at mu0 = 0.77:")
    for name, nm in CHANNELS.items():
        core = scene_output(nm, ref, CORE_TAU_890, MU0, 0.5) / 0.5
        thin = scene_output(nm, ref, 0.8, MU0, 0.5) / 0.5
        print(f"  {name:5s} core {core:.3f}  crescent {thin:.3f}")

    print("\nthe scene's factor, T + rho (1 - T), T = exp(-a tau m), tau at 890 nm:")
    a = anchored(CHANNELS["green"], ref)[0]
    print(f"one a for all three, the green's: a = {a:.3f}")
    for name, nm in CHANNELS.items():
        _, rho, core, thin = anchored(nm, ref, a)
        t_thin = np.exp(-a * CRESCENT_TAU_890 * 2.0 / MU0)
        print(f"  {name:5s} rho {rho:.3f}: core {core:.3f}, crescent {t_thin + rho * (1 - t_thin):.3f} (layer {thin:.3f})")


if __name__ == "__main__":
    main()
