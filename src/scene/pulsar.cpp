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

bool pulsar_corotating_axes(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dmat3& out)
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
    const glm::dvec3 x = std::cos(phase) * e1 + std::sin(phase) * e2;
    const FrameTransform& to_display = scene.current_transform();
    out = glm::dmat3(to_display.direction_to_display(x), to_display.direction_to_display(glm::cross(spin, x)),
                     to_display.direction_to_display(spin));
    return true;
}

bool pulsar_beams(const Scene& scene, const Body& body, double t_tdb, double real_time, glm::dvec3 out[2])
{
    glm::dmat3 axes;
    if (!pulsar_corotating_axes(scene, body, t_tdb, real_time, axes)) {
        return false;
    }
    const double tilt = body.pulsar.magnetic_inclination_deg * kDegToRad;
    out[0] = std::cos(tilt) * axes[2] + std::sin(tilt) * axes[0];
    out[1] = -out[0];
    return true;
}

double light_cylinder_km(double spin_period_s)
{
    return kSpeedOfLightKmS * spin_period_s / kTwoPi;
}

namespace {

// A point of the dipole line r = L sin^2(theta) at magnetic colatitude theta and azimuth phi,
// in the corotating frame (spin along z, the magnetic axis tilted toward +x).
glm::dvec3 dipole_point(double alpha, double l, double theta, double phi)
{
    const double r = l * std::sin(theta) * std::sin(theta);
    const glm::dvec3 axis(std::sin(alpha), 0.0, std::cos(alpha));
    const glm::dvec3 across(std::cos(alpha), 0.0, -std::sin(alpha)); // in the plane of spin and axis
    const glm::dvec3 side(0.0, 1.0, 0.0);
    return r * (std::sin(theta) * (std::cos(phi) * across + std::sin(phi) * side) + std::cos(theta) * axis);
}

double from_spin_axis(const glm::dvec3& p)
{
    return std::hypot(p.x, p.y);
}

} // namespace

double pulsar_last_closed_l(double inclination_deg, double phi, double light_cylinder)
{
    const double alpha = inclination_deg * kDegToRad;
    // The line's greatest distance from the spin axis, per unit L, by sampling and then a
    // golden-section search around the best sample.
    constexpr int kSamples = 720;
    int best = 0;
    double widest = 0.0;
    for (int i = 0; i <= kSamples; ++i) {
        const double d = from_spin_axis(dipole_point(alpha, 1.0, kPi * i / kSamples, phi));
        if (d > widest) {
            widest = d;
            best = i;
        }
    }
    double a = kPi * std::max(best - 1, 0) / kSamples, b = kPi * std::min(best + 1, kSamples) / kSamples;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    for (int k = 0; k < 60; ++k) {
        const double c = b - g * (b - a), d = a + g * (b - a);
        if (from_spin_axis(dipole_point(alpha, 1.0, c, phi)) > from_spin_axis(dipole_point(alpha, 1.0, d, phi))) {
            b = d;
        } else {
            a = c;
        }
    }
    widest = std::max(widest, from_spin_axis(dipole_point(alpha, 1.0, 0.5 * (a + b), phi)));
    return light_cylinder / widest;
}

std::vector<PulsarFieldLine> pulsar_field_lines(double inclination_deg, double radius_km, double light_cylinder, int azimuths)
{
    const double alpha = inclination_deg * kDegToRad;
    // Closed lines from near the star out to the last; open lines from inside the cap's edge.
    constexpr double kClosedFractions[] = {0.02, 0.08, 0.25, 0.55, 1.0};
    constexpr double kOpenFactors[] = {1.6, 4.0};
    constexpr int kSteps = 160;
    std::vector<PulsarFieldLine> lines;
    const auto fade = [&](const glm::dvec3& p) {
        return std::clamp((from_spin_axis(p) / light_cylinder - 0.5) / 0.5, 0.0, 1.0);
    };
    for (int k = 0; k < azimuths; ++k) {
        const double phi = kTwoPi * k / azimuths;
        const double last = pulsar_last_closed_l(inclination_deg, phi, light_cylinder);
        for (double fraction : kClosedFractions) {
            const double l = fraction * last;
            if (l < 1.5 * radius_km) {
                continue;
            }
            const double foot = std::asin(std::sqrt(radius_km / l));
            PulsarFieldLine line;
            for (int i = 0; i <= kSteps; ++i) {
                const double theta = foot + (kPi - 2.0 * foot) * i / kSteps;
                line.points.emplace_back(dipole_point(alpha, l, theta, phi), 0.0);
            }
            lines.push_back(std::move(line));
        }
        for (double factor : kOpenFactors) {
            const double l = factor * last;
            const double foot = std::asin(std::sqrt(radius_km / l));
            for (int pole = 0; pole < 2; ++pole) {
                // From the cap toward the equator until the line reaches the cylinder.
                PulsarFieldLine line;
                line.open = true;
                const auto at = [&](double s) { return dipole_point(alpha, l, pole == 0 ? s : kPi - s, phi); };
                line.points.emplace_back(at(foot), 0.0);
                double previous = foot;
                for (int i = 1; i <= 4 * kSteps; ++i) {
                    const double s = foot + (0.5 * kPi - foot) * 2.0 * i / (4 * kSteps);
                    const glm::dvec3 p = at(s);
                    if (from_spin_axis(p) >= light_cylinder) {
                        // Where along the step the line meets the cylinder.
                        double a = previous, b = s;
                        for (int halving = 0; halving < 50; ++halving) {
                            const double c = 0.5 * (a + b);
                            (from_spin_axis(at(c)) < light_cylinder ? a : b) = c;
                        }
                        line.points.emplace_back(at(0.5 * (a + b)), 1.0);
                        break;
                    }
                    line.points.emplace_back(p, fade(p));
                    previous = s;
                }
                lines.push_back(std::move(line));
            }
        }
    }
    // The axis itself, out of both poles, to the cylinder (or 2 cylinder radii if aligned).
    const glm::dvec3 axis(std::sin(alpha), 0.0, std::cos(alpha));
    const double reach = std::min(light_cylinder / std::max(std::abs(axis.x), 1e-9), 2.0 * light_cylinder);
    for (int pole = 0; pole < 2; ++pole) {
        PulsarFieldLine line;
        line.open = true;
        const glm::dvec3 direction = pole == 0 ? axis : -axis;
        for (int i = 0; i <= 32; ++i) {
            const glm::dvec3 p = direction * (radius_km + (reach - radius_km) * i / 32.0);
            line.points.emplace_back(p, fade(p));
        }
        lines.push_back(std::move(line));
    }
    return lines;
}

double pulsar_beam_flash(const Body::Pulsar& pulsar, const glm::dvec3& axis, const glm::dvec3& to_camera)
{
    const double off = std::acos(std::clamp(glm::dot(axis, to_camera), -1.0, 1.0)) /
                       (pulsar.beam_half_angle_deg * kDegToRad);
    const double edge = pulsar.beam_cutoff ? std::clamp((1.0 - off) / 0.2, 0.0, 1.0) : 1.0;
    return std::exp(-off * off) * edge * edge * (3.0 - 2.0 * edge);
}

} // namespace astraxis
