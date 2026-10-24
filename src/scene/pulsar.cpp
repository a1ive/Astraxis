#include "scene/pulsar.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

glm::dvec3 pulsar_spin_axis(const Scene& scene, const Body& body, double t_tdb)
{
    const auto& p = body.pulsar;
    const auto& prec = p.precession;
    if (prec.enabled) {
        const double phase = std::remainder(prec.rate_rad_s * (t_tdb - prec.epoch_tdb), kTwoPi);
        const double tilt = prec.misalignment_deg * kDegToRad;
        const glm::dvec3 y = glm::cross(prec.axis_icrf, prec.x_icrf);
        return std::cos(tilt) * prec.axis_icrf +
               std::sin(tilt) * (std::cos(phase) * prec.x_icrf + std::sin(phase) * y);
    }
    const State s = scene.bodies[static_cast<size_t>(p.spin_axis_orbit_of)].motion->eval(t_tdb);
    const glm::dvec3 h = glm::cross(s.position, s.velocity);
    return glm::length(h) > 0.0 ? glm::normalize(h) : scene.view.up_icrf;
}

bool pulsar_beams(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dvec3 out[2])
{
    if (!body.pulsar.enabled || body.pulsar.spin_axis_orbit_of < 0) {
        return false;
    }
    const glm::dvec3 spin = pulsar_spin_axis(scene, body, t_tdb);
    // Keep the display spin phase independent of the chosen display frame.
    const glm::dvec3 ref = std::abs(spin.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    const glm::dvec3 e1 = glm::normalize(glm::cross(spin, ref));
    const glm::dvec3 e2 = glm::cross(spin, e1);
    const double phase = kTwoPi * std::remainder(real_time, body.pulsar.display_period_s) / body.pulsar.display_period_s;
    const double tilt = body.pulsar.magnetic_inclination_deg * kDegToRad;
    out[0] = scene.current_transform().direction_to_display(
        std::cos(tilt) * spin + std::sin(tilt) * (std::cos(phase) * e1 + std::sin(phase) * e2));
    out[1] = -out[0];
    return true;
}

double pulsar_beam_flash(const Body::Pulsar& pulsar, const glm::dvec3& axis, const glm::dvec3& to_camera)
{
    const double off = std::acos(std::clamp(glm::dot(axis, to_camera), -1.0, 1.0)) /
                       (pulsar.beam_half_angle_deg * kDegToRad);
    const double edge = pulsar.beam_cutoff ? std::clamp((1.0 - off) / 0.2, 0.0, 1.0) : 1.0;
    return std::exp(-off * off) * edge * edge * (3.0 - 2.0 * edge);
}

} // namespace astraxis
