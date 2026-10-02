#include "scene/scene.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <algorithm>
#include <cmath>

namespace astraxis {

namespace {

// Obliquity of the ecliptic at J2000: 84381.412 arcsec (JPL SSD astrodynamic
// parameters, https://ssd.jpl.nasa.gov/astro_par.html).
constexpr double kObliquityJ2000 = 84381.412 / 3600.0 * kDegToRad;

// Default history length for orbit-trail bodies shown in a rotating frame.
constexpr double kRotatingTrailDays = 30.0;
// Default history length when a motion has no defined start.
constexpr double kDefaultHistoryDays = 365.0;

} // namespace

double lagrange_gamma(double mass_ratio, int point)
{
    const double mu = mass_ratio;
    double g = std::cbrt(mu / 3.0);
    for (int i = 0; i < 50; ++i) {
        double f = 0.0;
        double df = 0.0;
        if (point == 1) {
            f = ((((g - (3.0 - mu)) * g + (3.0 - 2.0 * mu)) * g - mu) * g + 2.0 * mu) * g - mu;
            df = (((5.0 * g - 4.0 * (3.0 - mu)) * g + 3.0 * (3.0 - 2.0 * mu)) * g - 2.0 * mu) * g + 2.0 * mu;
        } else {
            f = ((((g + (3.0 - mu)) * g + (3.0 - 2.0 * mu)) * g - mu) * g - 2.0 * mu) * g - mu;
            df = (((5.0 * g + 4.0 * (3.0 - mu)) * g + 3.0 * (3.0 - 2.0 * mu)) * g - 2.0 * mu) * g - 2.0 * mu;
        }
        const double step = f / df;
        g -= step;
        if (std::abs(step) < 1e-15) {
            break;
        }
    }
    return g;
}

glm::dvec3 ecliptic_pole_icrf()
{
    return rotation_x(kObliquityJ2000) * glm::dvec3(0.0, 0.0, 1.0);
}

int Scene::find(std::string_view body_name) const
{
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].name == body_name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

int Scene::star_index() const
{
    for (size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].kind == BodyKind::Star) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void Scene::set_active_frame(int index)
{
    if (index >= 0 && index < static_cast<int>(frames.size())) {
        m_active_frame = index;
    }
}

State Scene::icrf_state_at(int body, double t_tdb) const
{
    State s;
    for (int b = body; b >= 0; b = bodies[static_cast<size_t>(b)].parent) {
        if (const MotionSource* m = bodies[static_cast<size_t>(b)].motion.get()) {
            const State rel = m->eval(t_tdb);
            s.position += rel.position;
            s.velocity += rel.velocity;
        }
    }
    return s;
}

State Scene::sun_icrf_state_at(double t_tdb) const
{
    const int star = star_index();
    if (star >= 0) {
        return icrf_state_at(star, t_tdb);
    }
    if (origin_heliocentric) {
        const State helio = origin_heliocentric->eval(t_tdb);
        return {-helio.position, -helio.velocity};
    }
    return {glm::dvec3(kAuKm, 0.0, 0.0), glm::dvec3(0.0)};
}

FrameTransform Scene::frame_transform(double t_tdb) const
{
    const DisplayFrame& f = frame();
    FrameTransform out;
    out.origin = icrf_state_at(f.origin, t_tdb).position;
    if (f.type == DisplayFrame::Type::Rotating) {
        const State primary = f.primary >= 0 ? icrf_state_at(f.primary, t_tdb) : sun_icrf_state_at(t_tdb);
        const State secondary = icrf_state_at(f.secondary, t_tdb);
        const glm::dvec3 r = secondary.position - primary.position;
        const glm::dvec3 v = secondary.velocity - primary.velocity;
        const glm::dvec3 x = glm::normalize(r);
        const glm::dvec3 z = glm::normalize(glm::cross(r, v));
        out.axes = glm::dmat3(x, glm::cross(z, x), z);
    }
    return out;
}

void Scene::update(double t_tdb)
{
    m_time = t_tdb;
    const double days = t_tdb / kSecondsPerDay;
    m_transform = frame_transform(t_tdb);

    double extent = 0.0;
    for (Body& body : bodies) {
        glm::dvec3 pos(0.0);
        bool visible = true;
        if (body.parent >= 0) {
            const Body& parent = bodies[static_cast<size_t>(body.parent)];
            pos = parent.icrf_position;
            visible = parent.visible;
        }
        if (body.motion) {
            pos += body.motion->eval(t_tdb).position;
            visible = visible && body.motion->valid_at(t_tdb);
        }
        body.icrf_position = pos;
        body.visible = visible;
        body.world_position = m_transform.to_display(pos);

        const double w = wrap_two_pi((body.pm_w0_deg + body.pm_rate_deg_per_day * days) * kDegToRad);
        body.orientation = glm::transpose(m_transform.axes) * iau_pole_frame(body.pole_ra, body.pole_dec) * rotation_z(w);

        if (visible && body.kind != BodyKind::Star) {
            extent = std::max(extent, glm::length(pos) + body.equatorial_radius_km);
        }
    }
    m_system_extent_km = std::max(extent, bodies.empty() ? 1.0 : bodies[0].equatorial_radius_km);

    const glm::dvec3 sun_icrf = sun_icrf_state_at(t_tdb).position;
    m_sun_position = m_transform.to_display(sun_icrf);

    for (Marker& marker : markers) {
        const glm::dvec3 primary =
            marker.primary >= 0 ? bodies[static_cast<size_t>(marker.primary)].icrf_position : sun_icrf;
        const glm::dvec3 secondary = bodies[static_cast<size_t>(marker.secondary)].icrf_position;
        const glm::dvec3 d = secondary - primary;
        const double gamma = lagrange_gamma(marker.mass_ratio, marker.point);
        const glm::dvec3 offset = d * gamma;
        marker.world_position = m_transform.to_display(marker.point == 1 ? secondary - offset : secondary + offset);
    }
}

std::vector<double> Scene::periapsis_times(int body, double t0, double t1) const
{
    std::vector<double> times;
    const Body& b = bodies[static_cast<size_t>(body)];
    if (!b.motion || b.parent < 0 || t1 <= t0) {
        return times;
    }
    const MotionSource& motion = *b.motion;
    auto radial = [&](double t) {
        const State s = motion.eval(t);
        return glm::dot(s.position, s.velocity);
    };

    // Sample at a small fraction of the osculating period (bounded so that a
    // long window cannot take unbounded time), then bisect each sign change.
    double step = (t1 - t0) / 64.0;
    const double gm = bodies[static_cast<size_t>(b.parent)].gm_km3_s2;
    if (gm > 0.0 && motion.valid_at(t1)) {
        const State s = motion.eval(t1);
        const double energy = 0.5 * glm::dot(s.velocity, s.velocity) - gm / glm::length(s.position);
        if (energy < 0.0) {
            const double a = -gm / (2.0 * energy);
            step = std::min(step, kTwoPi * std::sqrt(a * a * a / gm) / 64.0);
        }
    }
    step = std::max(step, (t1 - t0) / 20000.0);

    bool have_prev = false;
    double prev_t = t0;
    double prev = 0.0;
    for (double t = t0;; t = std::min(t + step, t1)) {
        if (!motion.valid_at(t)) {
            have_prev = false;
        } else {
            const double cur = radial(t);
            if (have_prev && prev < 0.0 && cur >= 0.0) {
                double lo = prev_t;
                double hi = t;
                for (int k = 0; k < 50 && hi - lo > 1e-3; ++k) {
                    const double mid = 0.5 * (lo + hi);
                    (radial(mid) < 0.0 ? lo : hi) = mid;
                }
                times.push_back(0.5 * (lo + hi));
            }
            have_prev = true;
            prev_t = t;
            prev = cur;
        }
        if (t >= t1) {
            break;
        }
    }
    return times;
}

glm::dvec3 Scene::orbit_normal(int body) const
{
    const Body& b = bodies[static_cast<size_t>(body)];
    if (!b.motion) {
        return up_axis();
    }
    const State s = b.motion->eval(m_time);
    const glm::dvec3 h = glm::cross(s.position, s.velocity);
    const double len = glm::length(h);
    return len > 0.0 ? m_transform.direction_to_display(h / len) : up_axis();
}

bool Scene::osculating_apsides(int body, double* periapsis_km, double* apoapsis_km, double* period_s) const
{
    const Body& b = bodies[static_cast<size_t>(body)];
    if (!b.motion || b.parent < 0) {
        return false;
    }
    const double gm = bodies[static_cast<size_t>(b.parent)].gm_km3_s2;
    const State s = b.motion->eval(m_time);
    const double r = glm::length(s.position);
    const double energy = 0.5 * glm::dot(s.velocity, s.velocity) - gm / r;
    if (gm <= 0.0 || r <= 0.0 || energy >= 0.0) {
        return false;
    }
    const double a = -gm / (2.0 * energy);
    const double h = glm::length(glm::cross(s.position, s.velocity));
    const double e = std::sqrt(std::max(0.0, 1.0 - h * h / (gm * a)));
    *periapsis_km = a * (1.0 - e);
    *apoapsis_km = a * (1.0 + e);
    *period_s = kTwoPi * std::sqrt(a * a * a / gm);
    return true;
}

glm::dvec3 Scene::up_axis() const
{
    if (frame_is_rotating()) {
        return glm::dvec3(0.0, 0.0, 1.0);
    }
    const glm::dvec3 up = glm::length(view.up_icrf) > 0.0 ? view.up_icrf : bodies[0].pole;
    return glm::normalize(m_transform.direction_to_display(up));
}

void Scene::trail(int body_index, double t_tdb, int max_points, std::vector<glm::dvec3>& points,
                  std::vector<float>& fades) const
{
    points.clear();
    fades.clear();
    const Body& body = bodies[static_cast<size_t>(body_index)];
    if (!body.visible || !body.motion || body.parent < 0 || body.trail == TrailMode::None) {
        return;
    }

    TrailMode mode = body.trail;
    double history_days = body.trail_history_days;
    if (frame_is_rotating() && mode == TrailMode::Orbit) {
        mode = TrailMode::History;
        history_days = history_days > 0.0 ? history_days : kRotatingTrailDays;
    }

    if (mode == TrailMode::Orbit) {
        if (body.motion->sample_orbit(t_tdb, max_points, m_scratch_points)) {
            const glm::dvec3& parent = bodies[static_cast<size_t>(body.parent)].icrf_position;
            const size_t n = m_scratch_points.size();
            for (size_t k = 0; k < n; ++k) {
                points.push_back(m_transform.to_display(parent + m_scratch_points[k]));
                fades.push_back(static_cast<float>(k) / static_cast<float>(n - 1));
            }
            return;
        }
        // No closed orbit (e.g. hyperbolic): show where it has been instead.
        mode = TrailMode::History;
    }

    // History trail: earliest time at which the body exists, within the requested length.
    double t0 = t_tdb - (history_days > 0.0 ? history_days : kDefaultHistoryDays) * kSecondsPerDay;
    if (history_days <= 0.0 || !body.motion->valid_at(t0)) {
        // Binary search for the start of validity (motions are valid on one interval).
        double lo = t0;
        double hi = t_tdb;
        if (history_days <= 0.0) {
            lo = t_tdb - 300.0 * kDaysPerJulianYear * kSecondsPerDay;
        }
        if (!body.motion->valid_at(lo)) {
            for (int i = 0; i < 60; ++i) {
                const double mid = 0.5 * (lo + hi);
                (body.motion->valid_at(mid) ? hi : lo) = mid;
            }
            lo = hi;
        }
        t0 = lo;
    }
    if (t0 >= t_tdb) {
        return;
    }

    body.motion->history_times(t0, t_tdb, max_points, m_scratch_times);
    const bool fixed_frame = frame().type == DisplayFrame::Type::Inertial && frame().origin == 0;
    const double span = t_tdb - t0;
    for (auto it = m_scratch_times.rbegin(); it != m_scratch_times.rend(); ++it) {
        const double tk = *it;
        const glm::dvec3 icrf = icrf_state_at(body_index, tk).position;
        points.push_back(fixed_frame ? icrf : frame_transform(tk).to_display(icrf));
        fades.push_back(static_cast<float>((t_tdb - tk) / span));
    }
}

} // namespace astraxis
