#pragma once

#include "scene/shape_model.hpp"

namespace astraxis {

struct RocheSurface;

// An opaque slab around a gainer, with x toward the donor and z along the orbit normal.
struct AccretionDisk {
    ShapeModel mesh;  // photosphere used for the synthetic light curve
    std::vector<double> vertex_temperature_k;
    double outer_radius_km = 0.0;
    double half_height_km = 0.0;
    double rim_temperature_k = 0.0;
    double inner_temperature_k = 0.0;
    bool illustrative_torus = false; // render the gas volume and the gainer
};

AccretionDisk make_accretion_disk(double radius_km, double half_height_km, double core_radius_km,
                                double rim_temperature_k, double inner_temperature_k, double temperature_power,
                                double rim_heating_k, double rim_heating_radius_km);

// Magnitudes relative to maximum light at equally spaced observer azimuths,
// starting on the donor side (+x) and turning toward +y in the disk frame.
std::vector<float> accretion_binary_curve(const AccretionDisk& disk, const RocheSurface& donor, double separation_km,
                                         double observer_z, double wavelength_nm, double limb_darkening);

} // namespace astraxis
