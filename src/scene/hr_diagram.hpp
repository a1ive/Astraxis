#pragma once

#include <vector>

namespace astraxis {

// Stellar models behind the Hertzsprung-Russell diagram in the info panel, for
// solar metallicity (Z = 0.02, the models' own solar value). Masses, luminosities
// and radii in solar units.
//
// [TO96]  Tout, Pols, Eggleton & Han 1996, MNRAS 281, 257: zero-age main-sequence
//         luminosity and radius, fitted to 0.1-100 M_sun.
// [HPT00] Hurley, Pols & Tout 2000, MNRAS 315, 543 (astro-ph/0001295): the
//         main-sequence lifetime, the luminosity and radius at its end and their
//         evolution in between, with the Appendix's coefficients at Z = 0.02. The
//         authors' SSE code settles two misprints in the paper: alpha_R's
//         denominator is a59 + M^a61, and beta_R' rises from 1.06 over
//         (M - 1) / (a74 - 1).
// [TO97]  Tout, Aarseth, Pols & Eggleton 1997, MNRAS 291, 732:
//         white dwarf radius 0.0115 sqrt((M_Ch/M)^(2/3) - (M/M_Ch)^(2/3)) R_sun,
//         M_Ch = 1.44 M_sun, as [HPT00] give it for their white dwarfs.

struct StellarModelPoint {
    double luminosity = 1.0; // L_sun
    double radius = 1.0;     // R_sun
};

// [TO96] at mass `m` (0.1-100 M_sun).
StellarModelPoint zams_star(double m);

// Time from the zero-age main sequence to its end, in years [HPT00].
double main_sequence_lifetime_yr(double m);

// A star of mass `m` at `age_yr` since the zero-age main sequence, held at the
// end of the main sequence once it gets there [HPT00]. The hook (the brief
// contraction near the end of core hydrogen burning, M > 1.0185 M_sun) is left
// out; at the end itself it has no effect.
StellarModelPoint main_sequence_star(double m, double age_yr);

// Radius of a cold white dwarf of mass `m` [TO97].
double white_dwarf_radius(double m);

// Effective temperature of a star of luminosity `l` and radius `r`, with the
// nominal solar 5,772 K (IAU 2015 B3), as the scenes give theirs.
double effective_temperature_k(double l, double r);

// The age of the Universe, 13.787 Gyr (Planck 2018 VI, A&A 641, A6, Table 2).
inline constexpr double kUniverseAgeYr = 13.787e9;

// A point on the diagram: log10 of T_eff (K) and of L (L_sun).
struct HrPoint {
    double log_t = 0.0;
    double log_l = 0.0;
};

// A region of the diagram as a closed outline: its lower edge, then its upper
// edge back.
struct HrBand {
    std::vector<HrPoint> outline;

    bool contains(HrPoint p) const;
    // Distance from `p` to the region, in dex (0 inside).
    double distance(HrPoint p) const;
};

// Where solar-metallicity stars of 0.1-100 M_sun burn hydrogen in their cores
// before `age_yr`: from the zero-age main sequence up to each mass's end of the
// main sequence, or to where it is at `age_yr` if that comes first (the low-mass
// stars, which have not left the band in the age of the Universe).
HrBand main_sequence_band(double age_yr = kUniverseAgeYr);

// Cold white dwarfs of 0.2-1.2 M_sun between `hot_k` and `cool_k`: the band
// between their radii, along which they cool.
HrBand white_dwarf_band(double hot_k = 150000.0, double cool_k = 3500.0);

} // namespace astraxis
