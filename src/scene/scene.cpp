#include "scene/scene.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/kepler.hpp"

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
// History trails: a segment is split while its midpoint lies farther than this
// fraction of the chord from it (~4.6 deg of turn), at most kTrailRefineDepth
// times, and the trail grows to at most kTrailRefineBudget times the requested
// point count.
constexpr double kTrailMaxSag = 0.01;
constexpr int kTrailRefineDepth = 6;
constexpr size_t kTrailRefineBudget = 4;

} // namespace

double RingSystem::optical_depth(double r_km) const
{
    if (profile.size() < 2 || r_km < inner_km || r_km > outer_km) {
        return 0.0;
    }
    const double x = (r_km - inner_km) / (outer_km - inner_km) * static_cast<double>(profile.size() - 1);
    const size_t i = std::min(static_cast<size_t>(x), profile.size() - 2);
    const double f = x - static_cast<double>(i);
    return profile[i].x + (profile[i + 1].x - profile[i].x) * f;
}

void rasterize_ring_bands(RingSystem& rings, int samples)
{
    rings.profile.clear();
    if (rings.bands.empty() || samples < 2) {
        return;
    }
    double inner = rings.bands[0].inner_km;
    double outer = rings.bands[0].outer_km;
    for (const RingBand& b : rings.bands) {
        inner = std::min(inner, b.inner_km);
        outer = std::max(outer, b.outer_km);
    }
    // One empty sample beyond each end, so the edges fade within the profile.
    const double margin = (outer - inner) / (samples - 3);
    rings.inner_km = std::max(0.0, inner - margin);
    rings.outer_km = outer + margin;
    const double cell = (rings.outer_km - rings.inner_km) / (samples - 1);
    for (int k = 0; k < samples; ++k) {
        const double r = rings.inner_km + cell * k;
        double tau = 0.0;
        double floor_weighted = 0.0;
        for (const RingBand& b : rings.bands) {
            const double width = b.outer_km - b.inner_km;
            const double overlap = std::max(0.0, std::min(b.outer_km, r + 0.5 * cell) - std::max(b.inner_km, r - 0.5 * cell));
            const double inside = overlap / cell;
            tau += inside * b.optical_depth;
            floor_weighted += inside * b.optical_depth * (b.thickness_km / width);
        }
        // Where bands overlap, their thickness bounds are averaged by optical depth.
        rings.profile.emplace_back(static_cast<float>(tau), tau > 0.0 ? static_cast<float>(floor_weighted / tau) : 0.0f);
    }
}

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

glm::dvec3 belt_position(const SceneBelt& belt, size_t index, double t_tdb)
{
    const float* el = &belt.elements[index * SceneBelt::kStride];
    const double a = el[0] * kAuKm;
    const double e = el[1];
    const double n = std::sqrt(kSunGmKm3S2 / (a * a * a)); // rad/s
    KeplerElements k;
    k.a = a;
    k.e = e;
    k.i = el[2];
    k.node = el[3];
    k.arg_peri = el[4];
    k.mean_anomaly = wrap_two_pi(el[5] + n * (t_tdb - belt.epoch_tdb));
    const glm::dvec3 ecliptic = kepler_state(k, n).position;
    // J2000 ecliptic -> ICRF (SBDB uses the IAU 1976 obliquity, 84381.448").
    const double eps = 84381.448 / 3600.0 * kDegToRad;
    return {ecliptic.x, std::cos(eps) * ecliptic.y - std::sin(eps) * ecliptic.z,
            std::sin(eps) * ecliptic.y + std::cos(eps) * ecliptic.z};
}

int Scene::satellite_host(int body) const
{
    const int parent = bodies[static_cast<size_t>(body)].parent;
    if (parent < 0 || bodies[static_cast<size_t>(parent)].kind != BodyKind::Barycenter) {
        return parent;
    }
    for (int k = parent + 1; k < body; ++k) {
        const Body& sibling = bodies[static_cast<size_t>(k)];
        if (sibling.parent == parent && sibling.kind != BodyKind::Barycenter) {
            return k; // the primary
        }
    }
    return satellite_host(parent); // this is the primary
}

void satellite_fades(const Scene& scene, const glm::dvec3& camera, double px_per_radian, double hide_px,
                     double show_px, std::vector<float>& fades)
{
    fades.assign(scene.bodies.size(), 1.0f);
    // Parents and primaries precede their satellites, so host fades are ready.
    for (size_t i = 0; i < scene.bodies.size(); ++i) {
        const int host = scene.satellite_host(static_cast<int>(i));
        if (host < 0) {
            continue;
        }
        const glm::dvec3& host_pos = scene.bodies[static_cast<size_t>(host)].world_position;
        const double r = glm::length(scene.bodies[i].world_position - host_pos);
        const double d = glm::length(host_pos - camera);
        const double px = d > r ? r / d * px_per_radian : show_px;
        const double own = std::clamp((px - hide_px) / (show_px - hide_px), 0.0, 1.0);
        fades[i] = static_cast<float>(own * own * (3.0 - 2.0 * own)) * fades[static_cast<size_t>(host)];
    }
}

int Scene::lighting_stars(int body, StarLight* out, int max_lights) const
{
    const glm::dvec3& at = bodies[static_cast<size_t>(body)].world_position;
    int count = 0;
    for (size_t k = 0; k < bodies.size(); ++k) {
        const Body& star = bodies[k];
        if (star.kind != BodyKind::Star || !star.visible || static_cast<int>(k) == body) {
            continue;
        }
        const glm::dvec3 d = star.world_position - at;
        const StarLight light{static_cast<int>(k), star.luminosity_solar / std::max(glm::dot(d, d), 1.0)};
        // Insertion into the list kept sorted by flux, brightest first.
        int slot = count < max_lights ? count++ : max_lights;
        while (slot > 0 && out[slot - 1].relative_flux < light.relative_flux) {
            if (slot < max_lights) {
                out[slot] = out[slot - 1];
            }
            --slot;
        }
        if (slot < max_lights) {
            out[slot] = light;
        }
    }
    const double brightest = count > 0 ? out[0].relative_flux : 1.0;
    int kept = 0;
    for (int k = 0; k < count; ++k) {
        const double relative = out[k].relative_flux / brightest;
        if (relative >= kMinRelativeStarFlux) {
            out[kept] = out[k];
            out[kept++].relative_flux = relative;
        }
    }
    return kept;
}

int Scene::lighting_star(int body) const
{
    StarLight light;
    return lighting_stars(body, &light, 1) > 0 ? light.star : -1;
}

glm::dvec3 Scene::light_position(int body) const
{
    const int star = lighting_star(body);
    return star >= 0 ? bodies[static_cast<size_t>(star)].world_position : m_sun_position;
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

        double ra = body.pole_ra;
        double dec = body.pole_dec;
        double w_deg = body.pm_w0_deg + body.pm_rate_deg_per_day * days;
        if (body.nut_prec_source >= 0 || body.pole_rate_deg_per_century != glm::dvec2(0.0)) {
            const double centuries = days / kDaysPerJulianCentury;
            ra += body.pole_rate_deg_per_century.x * centuries * kDegToRad;
            dec += body.pole_rate_deg_per_century.y * centuries * kDegToRad;
            if (body.nut_prec_source >= 0) {
                const auto& angles = bodies[static_cast<size_t>(body.nut_prec_source)].nut_prec_angles;
                for (size_t k = 0; k < angles.size(); ++k) {
                    const double theta = (angles[k].x + angles[k].y * centuries) * kDegToRad;
                    if (k < body.nut_prec_ra.size()) {
                        ra += body.nut_prec_ra[k] * std::sin(theta) * kDegToRad;
                    }
                    if (k < body.nut_prec_dec.size()) {
                        dec += body.nut_prec_dec[k] * std::cos(theta) * kDegToRad;
                    }
                    if (k < body.nut_prec_pm.size()) {
                        w_deg += body.nut_prec_pm[k] * std::sin(theta);
                    }
                }
            }
        }
        const double w = wrap_two_pi(w_deg * kDegToRad);
        body.orientation = glm::transpose(m_transform.axes) * iau_pole_frame(ra, dec) * rotation_z(w);

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

bool osculating_elements(const State& relative, double gm, OrbitElements* out)
{
    const double r = glm::length(relative.position);
    if (gm <= 0.0 || r <= 0.0) {
        return false;
    }
    const glm::dvec3 h = glm::cross(relative.position, relative.velocity);
    const double h2 = glm::dot(h, h);
    const double energy = 0.5 * glm::dot(relative.velocity, relative.velocity) - gm / r;
    const double e = std::sqrt(std::max(0.0, 1.0 + 2.0 * energy * h2 / (gm * gm)));
    *out = OrbitElements{};
    out->eccentricity = e;
    out->periapsis_km = h2 / (gm * (1.0 + e)); // = a (1 - e), also for unbound orbits
    if (energy < 0.0) {
        const double a = -gm / (2.0 * energy);
        out->apoapsis_km = a * (1.0 + e);
        out->period_s = kTwoPi * std::sqrt(a * a * a / gm);
    }
    return true;
}

bool Scene::osculating_apsides(int body, double* periapsis_km, double* apoapsis_km, double* period_s) const
{
    const Body& b = bodies[static_cast<size_t>(body)];
    if (!b.motion || b.parent < 0) {
        return false;
    }
    OrbitElements el;
    if (!osculating_elements(b.motion->eval(m_time), bodies[static_cast<size_t>(b.parent)].gm_km3_s2, &el) ||
        el.period_s <= 0.0) {
        return false;
    }
    *periapsis_km = el.periapsis_km;
    *apoapsis_km = el.apoapsis_km;
    *period_s = el.period_s;
    return true;
}

Scene::OrbitCenter Scene::orbit_center(int body) const
{
    const int host = satellite_host(body);
    if (host < 0) {
        return {};
    }
    const Body& b = bodies[static_cast<size_t>(body)];
    const double own = b.body_gm_km3_s2;
    const double host_gm = bodies[static_cast<size_t>(host)].body_gm_km3_s2;
    const int bary = b.parent;
    if (bodies[static_cast<size_t>(bary)].kind != BodyKind::Barycenter) {
        return {host, host_gm > 0.0 ? host_gm + own : 0.0};
    }
    const double system_gm = bodies[static_cast<size_t>(bary)].body_gm_km3_s2;
    if (bodies[static_cast<size_t>(host)].parent != bary) {
        // The primary: the whole system orbits the barycenter's host.
        return {host, host_gm > 0.0 ? host_gm + (system_gm > 0.0 ? system_gm : own) : 0.0};
    }

    // Beside the primary: a circumbinary orbit if another massive member is
    // closer to the barycenter; its GM is then that of everything inside.
    const glm::dvec3& center = bodies[static_cast<size_t>(bary)].icrf_position;
    const double r = glm::length(b.icrf_position - center);
    double inner_gm = own;
    bool circumbinary = false;
    for (size_t k = static_cast<size_t>(bary) + 1; k < bodies.size(); ++k) {
        const Body& s = bodies[k];
        if (s.parent != bary || static_cast<int>(k) == body || !s.visible ||
            glm::length(s.icrf_position - center) >= r) {
            continue;
        }
        inner_gm += s.body_gm_km3_s2;
        if (static_cast<int>(k) != host && s.body_gm_km3_s2 > 0.0 && s.body_gm_km3_s2 >= 0.01 * host_gm) {
            circumbinary = true;
        }
    }
    if (circumbinary) {
        return {bary, inner_gm};
    }
    if (host_gm > 0.0 && own > 0.0) {
        return {host, host_gm + own};
    }
    return {host, system_gm > 0.0 ? system_gm : host_gm};
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
    if (!body.motion || body.parent < 0 || body.trail == TrailMode::None) {
        return;
    }

    // A finished trajectory (e.g. a spacecraft after splashdown) can linger:
    // its history trail ends where the motion ended and fades out over
    // trail_linger_days.
    double t_head = t_tdb; // newest point of the trail
    double linger = 0.0;   // 0 .. 1: how far a finished trail has faded
    if (!body.visible) {
        const double window = body.trail_linger_days * kSecondsPerDay;
        if (body.trail != TrailMode::History || window <= 0.0 ||
            !bodies[static_cast<size_t>(body.parent)].visible) {
            return;
        }
        // Motions are valid on one interval: find a valid time in the window,
        // then bisect for the end of validity.
        double valid = 0.0;
        bool found = false;
        for (int i = 1; i <= 64 && !found; ++i) {
            valid = t_tdb - window * i / 64.0;
            found = body.motion->valid_at(valid);
        }
        if (!found) {
            return;
        }
        double lo = valid;
        double hi = t_tdb;
        for (int i = 0; i < 60; ++i) {
            const double mid = 0.5 * (lo + hi);
            (body.motion->valid_at(mid) ? lo : hi) = mid;
        }
        t_head = lo;
        linger = (t_tdb - t_head) / window;
    }

    TrailMode mode = body.trail;
    double history_days = body.trail_history_days;
    if (frame_is_rotating() && mode == TrailMode::Orbit) {
        mode = TrailMode::History;
        history_days = history_days > 0.0 ? history_days : kRotatingTrailDays;
    }
    if (const double cap = frame().trail_history_days; cap > 0.0 && mode == TrailMode::History) {
        history_days = history_days > 0.0 ? std::min(history_days, cap) : cap;
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
    double t0 = t_head - (history_days > 0.0 ? history_days : kDefaultHistoryDays) * kSecondsPerDay;
    if (history_days <= 0.0 || !body.motion->valid_at(t0)) {
        // Binary search for the start of validity (motions are valid on one interval).
        double lo = t0;
        double hi = t_head;
        if (history_days <= 0.0) {
            lo = t_head - 300.0 * kDaysPerJulianYear * kSecondsPerDay;
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
    if (t0 >= t_head) {
        return;
    }

    body.motion->history_times(t0, t_head, max_points, m_scratch_times);
    const bool fixed_frame = frame().type == DisplayFrame::Type::Inertial && frame().origin == 0;
    const double span = t_head - t0;
    auto display_at = [&](double tk) {
        const glm::dvec3 icrf = icrf_state_at(body_index, tk).position;
        return fixed_frame ? icrf : frame_transform(tk).to_display(icrf);
    };
    auto emit = [&](double tk, const glm::dvec3& p) {
        points.push_back(p);
        const double fade = (t_head - tk) / span;
        fades.push_back(static_cast<float>(fade + (1.0 - fade) * linger)); // lingering: all toward the tail
    };

    // The samples are ephemeris knots. Where they are sparse (Kepler-relative
    // tables put only a few on each orbit) or a moving frame bends the path
    // between them (a rotating frame turns a gentle aphelion arc into a loop),
    // a segment is split while its midpoint sags off the chord.
    const size_t max_refined = static_cast<size_t>(std::max(max_points, 0)) * kTrailRefineBudget;
    auto refine = [&](auto&& self, double ta, const glm::dvec3& pa, double tb, const glm::dvec3& pb,
                      int depth) -> void {
        if (depth == 0 || points.size() >= max_refined) {
            return;
        }
        const double tm = 0.5 * (ta + tb);
        const glm::dvec3 pm = display_at(tm);
        if (glm::length(pm - 0.5 * (pa + pb)) <= kTrailMaxSag * glm::length(pb - pa)) {
            return;
        }
        self(self, ta, pa, tm, pm, depth - 1);
        emit(tm, pm);
        self(self, tm, pm, tb, pb, depth - 1);
    };

    double prev_t = 0.0;
    glm::dvec3 prev_p(0.0);
    for (auto it = m_scratch_times.rbegin(); it != m_scratch_times.rend(); ++it) {
        const double tk = *it;
        const glm::dvec3 p = display_at(tk);
        if (it != m_scratch_times.rbegin()) {
            refine(refine, prev_t, prev_p, tk, p, kTrailRefineDepth);
        }
        emit(tk, p);
        prev_t = tk;
        prev_p = p;
    }
}

} // namespace astraxis
