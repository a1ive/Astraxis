# White-light solar corona

The corona in `solar_eclipse.toml` and `parker.toml` is computed by this project's
own code. No solar photograph, magnetic map, paper figure or third-party shader
is included. The streamers and plumes are illustrative geometry.

## Scattering and density

- Electron density: Leblanc, Dulk & Bougeret (1998),
  `n_e(r) = 3.3e5/r^2 + 4.1e6/r^4 + 8.0e7/r^6` cm⁻³, with `r` in solar radii.
  The equation is reproduced in Section 3.1, page 8 of
  [Reginald, Newmark & Rastaetter (2023), Solar Physics 298:73](https://ntrs.nasa.gov/api/citations/20230015141/downloads/paper1.pdf?attachment=true),
  DOI 10.1007/s11207-023-02160-3. The cited comparison covers 1.5–3.5 solar radii;
  extending this profile down to 1 and out to 32 is a modeling choice here.
- Thomson scattering follows the finite solar disk treatment described in
  [Inhester (2016), Thomson Scattering in the Solar Corona](https://arxiv.org/abs/1512.00651).
  `src/scene/corona.cpp` numerically integrates the incident disk into coefficients
  `A(r) + B(r) cos²(chi)`, in mean solar brightness per solar radius of path.
  It uses 256 angular samples, a linear limb-darkening coefficient of 0.6 and
  the [CODATA 2022](https://physics.nist.gov/cuu/pdf/wall_2022.pdf) Thomson cross
  section, 6.6524587051e-29 m². The solar radius is 695,700 km.
- The GPU interpolates 128 radial coefficient pairs and integrates each ray
  with 96 samples in `z = b sinh(t)`. The geometric contrast multiplies the
  radial density; it is not normalized to preserve its spherical average.
  The volume fades from 24 to 32 solar radii. Half-resolution integration is
  recomputed at full resolution near opaque silhouettes.
- Solar axes use the IAU pole (286.13°, 63.87°) and prime meridian
  `84.176° + 14.1844° d`, as in
  [NAIF pck00011.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc).

These are independently implemented scientific formulas and published constants;
the code is GPL-3.0-or-later. The papers retain their authors' and publishers'
rights. Their availability through NASA or arXiv does not make their text or
figures public domain.

## Cycle and appearance

The geometric activity parameter is smoothly interpolated through three dates:
April 2014 maximum, December 2019 minimum, and the October 15, 2024 announcement
that the Sun had reached its maximum phase. Sources:
[NASA, Solar Cycle 25 (2020)](https://www.nasa.gov/missions/sdo/what-will-solar-cycle-25-look-like/)
and [NASA/NOAA maximum-phase announcement (2024)](https://science.nasa.gov/science-research/heliophysics/nasa-noaa-sun-reaches-maximum-phase-in-11-year-solar-cycle/).
The 2024 announcement is not an exact peak date. Values outside the anchor span
are held at their endpoints; this is not an eleven-year cycle forecast.
No SILSO sunspot series or GFZ sunspot column is bundled or read.

Low activity emphasizes a wrinkled equatorial sheet and polar plumes. Higher
activity adds a tilted sheet and stronger warping. The angles, widths, density
contrasts and fine striations are project design choices, not reconstructions
of the 2017 or 2024 corona. No PFSS extrapolation or CME is simulated.

The displayed white-light brightness is enhanced: default `display_gain=8e5`,
`display_gamma=0.65`. Ordinary eclipse views reduce the gain according to the
uncovered, limb-darkened solar disk; the Parker event explicitly enables full
enhancement. This contrast compression is separate from the scattering model
and is not a calibrated exposure or naked-eye brightness prediction. The tint
is the project's 5772 K solar continuum color. F-corona, EUV emission,
polarization and instrument response are omitted. The Parker view is not a
WISPR image reconstruction.

## Eclipse geometry

`tools/corona/eclipse_bake.py` generates the Earth and Moon tables and prints
the two observer offsets used by the scene. See
[ephemeris sources](../ephem/SOURCES.md#solar-eclipses-2017-and-2024) for sampling,
verification errors and data attribution. These viewpoints are 1000 km above
a spherical Earth, fixed in geocentric ICRF at each event, not ground sites or
propagated satellite orbits. The events use UTC; the ephemerides use TDB.
The scene opens the 2024 event; the event menu also offers 2017.
