#pragma once

#include <glm/vec3.hpp>

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace astraxis {

// The planetary Kp index (GFZ; Matzka et al. 2021, Space Weather 19,
// e2020SW002641), one value per 3-hour UT interval, as baked by
// tools/geomag/kp_bake.py into assets/geomag/kp.bin.
class KpRecord {
public:
    bool load(const std::filesystem::path& file, std::string* error);

    // Simulation times of the first interval's start and the last one's end.
    double start_tdb() const;
    double end_tdb() const;
    // Kp of the interval containing t; outside the record, its median.
    double at(double t_tdb) const;
    // Mean Kp over [t0, t1] (each interval weighted by its overlap); the parts
    // outside the record count as its median.
    double mean(double t0_tdb, double t1_tdb) const;
    double median() const { return m_median; }

private:
    // Intervals since the record's start (fractional) at a simulation time.
    double interval(double t_tdb) const;
    // Integral of Kp (in thirds) over intervals [0, x], extended by the median.
    double cumulative(double x) const;

    double m_start_jd_utc = 0.0;
    std::vector<uint8_t> m_thirds; // Kp in thirds of a unit
    std::vector<double> m_prefix;  // m_prefix[k] = sum of m_thirds[0 .. k-1]
    double m_median = 2.0;
};

// Auroral boundaries of Starkov (1994, Geomagn. Aeron. 34, 331), as given by
// Sigernes et al. (2011, Opt. Pura Apl. 44, 599, eq. 1-3 and Table 1): the
// corrected geomagnetic colatitudes (degrees) of the equatorward boundary of
// the diffuse aurora and of the equatorward and poleward boundaries of the
// auroral oval, at magnetic local time mlt_hours, for a planetary Kp. Kp
// enters through |AL| = 18 - 12.3 Kp + 27.2 Kp^2 - 2 Kp^3 nT (their eq. 1,
// after Starkov 1994, Geomagn. Aeron. 34, 101).
struct AuroralBoundaries {
    double diffuse_deg = 0.0;
    double equatorward_deg = 0.0;
    double poleward_deg = 0.0;
};
double starkov_al_nt(double kp);
AuroralBoundaries starkov_boundaries(double kp, double mlt_hours);

// Electron hemispheric power (GW) for a planetary Kp: the fit of Emery et al.
// (2008, JGR 113, A06311) to inter-calibrated NOAA POES and DMSP estimates,
// HPe = 2.5103 + 9.9207 Kp - 2.1825 Kp^2 + 0.3033 Kp^3, as quoted by
// Verkhoglyadova et al. (2016, J. Space Weather Space Clim. 6, A20).
double electron_hemispheric_power_gw(double kp);

// The north geomagnetic pole (where the centered dipole's axis leaves the
// northern hemisphere) as a unit vector in the Earth-fixed frame (z toward
// the north pole, x toward the prime meridian), from the dipole terms g10,
// g11, h11 of IGRF-14 (IAGA Working Group V-MOD), interpolated linearly
// between the 5-year epochs 1900-2025 and extrapolated with the 2025-30
// secular variation (held at the ends).
glm::dvec3 geomagnetic_north_pole(double t_tdb);

// Magnetic local time (hours, 0 .. 24) at unit direction n, for the north
// geomagnetic pole and the sun's direction (all in the same frame): the
// angle around the dipole axis from the anti-sun meridian, increasing
// eastward (noon toward the sun, 18 h at dusk). shaders/aurora.frag.hlsl
// computes it the same way.
double magnetic_local_time(const glm::dvec3& n, const glm::dvec3& pole, const glm::dvec3& sun);

// Column emission (kR) per unit electron energy flux (erg cm^-2 s^-1) for
// monoenergetic electrons of mean_energy_kev: Mathieu, Robert & Lamy (2025,
// J. Space Weather Space Clim. 15, 19), Table 2, interpolated in log energy
// between their 1 and 10 keV columns.
struct AuroraLineYields {
    double violet = 0.0; // N2+ first negative (0-0), 391.4 nm
    double blue = 0.0;   // N2+ first negative (0-1), 427.8 nm
    double green = 0.0;  // O(1S), 557.7 nm
    double red = 0.0;    // O(1D), 630.0-636.4-639.2 nm
};
AuroraLineYields aurora_line_yields(double mean_energy_kev);

// Observed line ratios to the blue 427.8 nm band (Dashkevich, Zverev & Ivanov
// 2006, Geomagn. Aeron. 46, 366, as summarized by Qiu et al. 2017, Chinese J.
// Geophys. 60, 1): I557.7 / I427.8 = 4.2 .. 6.4 for I427.8 of 0.1 .. 1.8 kR
// (their middle, kGreenOverBlue), and I630.0 / I427.8 falling from 2 to 0.4
// as I427.8 rises from 0.1 to 3 kR (a power law between those ends, held
// beyond them): faint aurora is redder.
inline constexpr double kGreenOverBlue = 5.3;
double red_over_blue(double blue_kr);

// The aurora of one hemisphere for a Kp, as the renderer draws it: the
// boundaries around the oval at kOvalSamples magnetic local times (sample k
// at 24 k / kOvalSamples hours), and a vertical column through the oval at
// its mean energy flux: its 427.8 nm brightness and the color of its low
// layer (linear RGB in I/F, before any exposure). The red line, in the high
// layer, depends on the local brightness (red_over_blue), so the renderer
// works it out from the blue at each point, with kRedPerKr.
//
// The hemispheric power falls on the oval (between its boundaries) and, at
// kDiffuseFlux of its flux per area, on the diffuse aurora equatorward of it.
// The precipitating electrons have a single mean energy kMeanEnergyKev, which
// sets the yield of the N2+ bands (aurora_line_yields); the green and red
// lines follow from the blue by the observed ratios. The low layer (green,
// blue and violet) peaks at kLowLayerKm, the mean height of green and blue
// aurora of Whiter et al. (2023, Ann. Geophys. 41, 1); the red one at
// kHighLayerKm, where the 630 nm emission of Mathieu et al.'s 7 keV case
// peaks (their Fig. 7). Each layer is a Chapman profile exp(-z - e^-z) / H,
// z = (h - h_peak) / H, of scale height H (km).
struct AuroraEmission {
    static constexpr int kOvalSamples = 48;
    static constexpr double kDiffuseFlux = 0.3;    // assumed
    static constexpr double kMeanEnergyKev = 3.0;  // assumed
    static constexpr double kLowLayerKm = 115.0;
    static constexpr double kLowScaleKm = 10.0;    // assumed
    static constexpr double kHighLayerKm = 221.0;
    static constexpr double kHighScaleKm = 40.0;   // assumed
    // Where the energy flux is spread (km above the 6371 km mean radius).
    static constexpr double kFluxAltitudeKm = 110.0;

    double kp = 0.0;
    double hemispheric_power_gw = 0.0;
    double oval_flux = 0.0; // erg cm^-2 s^-1 on the oval
    // Colatitudes (radians): x = diffuse, y = equatorward, z = poleward.
    std::array<glm::vec3, kOvalSamples> oval{};
    double blue_kr = 0.0;    // 427.8 nm of a column through the oval
    glm::dvec3 low_rgb{0.0}; // the column's low layer
};
// Linear RGB (I/F) of 1 kR of the red line (630.0 nm).
extern const glm::dvec3 kRedPerKr;
AuroraEmission aurora_emission(double kp);

} // namespace astraxis
