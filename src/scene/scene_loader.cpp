#include "scene/scene_loader.hpp"

#include "core/math.hpp"
#include "core/time.hpp"
#include "ephem/ephemeris.hpp"
#include "ephem/jpl_approx_planet.hpp"
#include "ephem/kerr_orbit.hpp"
#include "ephem/mean_element_orbit.hpp"
#include "ephem/nbody.hpp"
#include "ephem/visual_orbit.hpp"

#include <glm/vector_relational.hpp>
#include <toml++/toml.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <unordered_map>

namespace astraxis {

namespace {

// Samples of a ring profile rasterized from bands (Jupiter's: ~23 km each, Uranus': ~1.2 km).
constexpr int kBandProfileSamples = 8192;

// Thrown on schema errors; caught at the top level and turned into a message.
struct SchemaError {
    std::string message;
};

[[noreturn]] void fail(const std::string& context, const std::string& what)
{
    throw SchemaError{context + ": " + what};
}

double get_double(const toml::table& t, std::string_view key, const std::string& ctx)
{
    const auto v = t[key].value<double>();
    if (!v) {
        fail(ctx, "missing or non-numeric '" + std::string(key) + "'");
    }
    return *v;
}

double get_double_or(const toml::table& t, std::string_view key, double fallback)
{
    return t[key].value<double>().value_or(fallback);
}

std::string get_string(const toml::table& t, std::string_view key, const std::string& ctx)
{
    const auto v = t[key].value<std::string>();
    if (!v) {
        fail(ctx, "missing or non-string '" + std::string(key) + "'");
    }
    return *v;
}

std::string get_string_or(const toml::table& t, std::string_view key, const std::string& fallback)
{
    return t[key].value<std::string>().value_or(fallback);
}

// Numeric array with between min_n and n elements.
void get_array(const toml::table& t, std::string_view key, const std::string& ctx, double* out, size_t min_n,
               size_t n)
{
    const toml::array* arr = t[key].as_array();
    if (!arr || arr->size() < min_n || arr->size() > n) {
        fail(ctx, "'" + std::string(key) + "' must be an array of " + std::to_string(min_n) +
                      (min_n == n ? "" : "-" + std::to_string(n)) + " numbers");
    }
    for (size_t i = 0; i < arr->size(); ++i) {
        const auto v = (*arr)[i].value<double>();
        if (!v) {
            fail(ctx, "'" + std::string(key) + "' contains a non-number");
        }
        out[i] = *v;
    }
}

// A number (the same for all channels) or an [r, g, b] array.
glm::dvec3 get_rgb_or(const toml::table& t, std::string_view key, const std::string& ctx, glm::dvec3 fallback)
{
    if (const auto v = t[key].value<double>()) {
        return glm::dvec3(*v);
    }
    if (!t.contains(key)) {
        return fallback;
    }
    double rgb[3];
    get_array(t, key, ctx, rgb, 3, 3);
    return glm::dvec3(rgb[0], rgb[1], rgb[2]);
}

// [bodies.atmosphere]: optical depths at 550 nm, scaled to the display
// channels as lambda^-4 (gas) and lambda^-haze_angstrom (haze).
Atmosphere parse_atmosphere(const toml::table& t, const std::string& ctx)
{
    Atmosphere a;
    a.enabled = true;
    a.height_km = get_double(t, "height_km", ctx);
    a.deck_altitude_km = get_double_or(t, "deck_altitude_km", 0.0);
    const double rayleigh = get_double_or(t, "rayleigh_depth", 0.0);
    const double haze = get_double_or(t, "haze_depth", 0.0);
    const double angstrom = get_double_or(t, "haze_angstrom", 0.0);
    for (int c = 0; c < 3; ++c) {
        const double x = kChannelWavelengthsNm[c] / 550.0;
        a.rayleigh_depth[c] = rayleigh * std::pow(x, -4.0);
        a.haze_depth[c] = haze * std::pow(x, -angstrom);
    }
    a.rayleigh_scale_height_km = rayleigh > 0.0 ? get_double(t, "rayleigh_scale_height_km", ctx) : 1.0;
    a.haze_scale_height_km = haze > 0.0 ? get_double(t, "haze_scale_height_km", ctx) : 1.0;
    a.haze_albedo = get_rgb_or(t, "haze_albedo", ctx, glm::dvec3(1.0));
    a.haze_g = get_double_or(t, "haze_g", 0.0);
    a.gain = get_double_or(t, "gain", 1.0);
    if (!(a.height_km > 0.0) || a.deck_altitude_km < 0.0 || rayleigh < 0.0 || haze < 0.0 ||
        !(a.rayleigh_scale_height_km > 0.0) || !(a.haze_scale_height_km > 0.0) || a.gain < 0.0) {
        fail(ctx, "height_km and scale heights must be positive; depths, deck_altitude_km and gain non-negative");
    }
    if (std::abs(a.haze_g) >= 1.0 || glm::any(glm::lessThan(a.haze_albedo, glm::dvec3(0.0))) ||
        glm::any(glm::greaterThan(a.haze_albedo, glm::dvec3(1.0)))) {
        fail(ctx, "|haze_g| must be < 1 and haze_albedo within 0..1");
    }
    return a;
}

// [[bodies.plumes]]
Plume parse_plume(const toml::table& t, const std::string& ctx)
{
    Plume p;
    p.name = get_string(t, "name", ctx);
    const std::string c = ctx + " plume '" + p.name + "'";
    const std::string type = get_string(t, "type", c);
    if (type == "umbrella") {
        p.type = Plume::Type::Umbrella;
    } else if (type == "jets") {
        p.type = Plume::Type::Jets;
    } else if (type == "geyser") {
        p.type = Plume::Type::Geyser;
    } else {
        fail(c, "type must be \"umbrella\", \"jets\" or \"geyser\"");
    }
    double ll[2];
    get_array(t, "lat_lon_deg", c, ll, 2, 2);
    p.lat_lon_deg = glm::dvec2(ll[0], ll[1]);
    p.height_km = get_double(t, "height_km", c);
    p.optical_depth = get_double(t, "optical_depth", c);
    p.albedo = get_rgb_or(t, "albedo", c, glm::dvec3(1.0));
    p.g = get_double_or(t, "g", 0.0);
    switch (p.type) {
    case Plume::Type::Umbrella:
        p.radius_km = get_double_or(t, "radius_km", 2.0 * p.height_km); // ballistic reach
        break;
    case Plume::Type::Jets:
        get_array(t, "end_lat_lon_deg", c, ll, 2, 2);
        p.end_lat_lon_deg = glm::dvec2(ll[0], ll[1]);
        p.count = static_cast<int>(get_double(t, "count", c));
        p.spread_deg = get_double(t, "spread_deg", c);
        if (p.count < 1 || p.count > kMaxPlumeJets || !(p.spread_deg > 0.0 && p.spread_deg < 60.0)) {
            fail(c, "count must be 1.." + std::to_string(kMaxPlumeJets) + " and spread_deg within (0, 60)");
        }
        break;
    case Plume::Type::Geyser:
        p.radius_km = get_double(t, "radius_km", c);
        p.tail_km = get_double(t, "tail_km", c);
        p.tail_azimuth_deg = get_double(t, "tail_azimuth_deg", c);
        break;
    }
    if (!(p.height_km > 0.0) || !(p.radius_km > 0.0) || p.optical_depth < 0.0 || p.tail_km < 0.0 ||
        std::abs(p.g) >= 1.0 || glm::any(glm::lessThan(p.albedo, glm::dvec3(0.0))) ||
        glm::any(glm::greaterThan(p.albedo, glm::dvec3(1.0)))) {
        fail(c, "height and radius must be positive, optical_depth and tail_km non-negative, |g| < 1 and "
                "albedo within 0..1");
    }
    return p;
}

glm::vec3 parse_color(const toml::table& t, std::string_view key, const std::string& ctx, glm::vec3 fallback)
{
    const auto v = t[key].value<std::string>();
    if (!v) {
        return fallback;
    }
    const std::string& s = *v;
    if (s.size() != 7 || s[0] != '#') {
        fail(ctx, "'" + std::string(key) + "' must be \"#rrggbb\"");
    }
    char* end = nullptr;
    const unsigned long rgb = std::strtoul(s.c_str() + 1, &end, 16);
    if (*end != '\0') {
        fail(ctx, "'" + std::string(key) + "' must be \"#rrggbb\"");
    }
    return glm::vec3(static_cast<float>((rgb >> 16) & 0xff), static_cast<float>((rgb >> 8) & 0xff),
                     static_cast<float>(rgb & 0xff)) /
           255.0f;
}

double parse_time(const toml::table& t, std::string_view key, const std::string& ctx)
{
    const std::string text = get_string(t, key, ctx);
    double t_tdb = 0.0;
    if (!parse_utc(text.c_str(), &t_tdb)) {
        fail(ctx, "'" + std::string(key) + "' must be \"YYYY-MM-DD[ HH:MM[:SS]]\" (UTC)");
    }
    return t_tdb;
}

// Campbell elements on the sky. Size: a_au, or a_arcsec with parallax_mas, or
// a_mas with distance_pc, or (when the total `gm` is known) period_days through
// Kepler's third law. Epoch of periastron as a Julian year (t_peri_year) or a
// Julian date (t_peri_jd_tdb).
VisualOrbit parse_visual_orbit(const toml::table& t, const std::string& ctx, double gm = 0.0)
{
    VisualOrbit o;
    if (gm > 0.0 && t.contains("period_days")) {
        const double n = kTwoPi / (get_double(t, "period_days", ctx) * kSecondsPerDay);
        o.a_km = std::cbrt(gm / (n * n));
    } else if (t.contains("a_au")) {
        o.a_km = get_double(t, "a_au", ctx) * kAuKm;
    } else if (t.contains("a_arcsec")) {
        o.a_km = get_double(t, "a_arcsec", ctx) / (get_double(t, "parallax_mas", ctx) / 1000.0) * kAuKm;
    } else if (t.contains("a_mas")) {
        o.a_km = get_double(t, "a_mas", ctx) / 1000.0 * get_double(t, "distance_pc", ctx) * kAuKm;
    } else {
        fail(ctx, "needs a_au, a_arcsec + parallax_mas, or a_mas + distance_pc");
    }
    o.e = get_double(t, "e", ctx);
    o.i_deg = get_double(t, "i_deg", ctx);
    o.node_deg = get_double(t, "node_deg", ctx);
    o.arg_peri_deg = get_double(t, "arg_peri_deg", ctx);
    o.t_peri_tdb = t.contains("t_peri_jd_tdb") ? (get_double(t, "t_peri_jd_tdb", ctx) - kJ2000Jd) * kSecondsPerDay
                                               : tdb_from_julian_year(get_double(t, "t_peri_year", ctx));
    double sky[2];
    get_array(t, "sky_ra_dec_deg", ctx, sky, 2, 2);
    o.ra_deg = sky[0];
    o.dec_deg = sky[1];
    if (o.e < 0.0 || o.e >= 1.0 || o.a_km <= 0.0) {
        fail(ctx, "invalid elements (need 0 <= e < 1, a > 0)");
    }
    return o;
}

// Transit-timing elements (see TransitOrbit). The eccentricity is given either
// as e and arg_peri_deg or as the fitted e_cos_w and e_sin_w.
TransitOrbit parse_transit_orbit(const toml::table& t, const std::string& ctx)
{
    TransitOrbit o;
    o.period_days = get_double(t, "period_days", ctx);
    o.t_transit_tdb = (get_double(t, "t_transit_bjd_tdb", ctx) - kJ2000Jd) * kSecondsPerDay;
    if (t.contains("e_cos_w")) {
        const double ec = get_double(t, "e_cos_w", ctx);
        const double es = get_double(t, "e_sin_w", ctx);
        o.e = std::hypot(ec, es);
        o.arg_peri_deg = std::atan2(es, ec) * kRadToDeg;
    } else {
        o.e = get_double(t, "e", ctx);
        o.arg_peri_deg = get_double(t, "arg_peri_deg", ctx);
    }
    o.i_deg = get_double(t, "i_deg", ctx);
    o.node_deg = get_double_or(t, "node_deg", 0.0);
    o.transit_u_deg = get_double(t, "transit_u_deg", ctx); // required: conventions differ between codes
    double sky[2];
    get_array(t, "sky_ra_dec_deg", ctx, sky, 2, 2);
    o.ra_deg = sky[0];
    o.dec_deg = sky[1];
    if (o.e < 0.0 || o.e >= 1.0 || o.period_days <= 0.0) {
        fail(ctx, "invalid elements (need 0 <= e < 1, period > 0)");
    }
    return o;
}

// Event times: "time" (UTC) or "time_tdb" (TDB, e.g. straight from the bake tool).
double parse_event_time(const toml::table& t, const std::string& ctx)
{
    if (t.contains("time_tdb")) {
        const std::string text = get_string(t, "time_tdb", ctx);
        CalendarDateTime c;
        const int fields = std::sscanf(text.c_str(), "%d-%d-%d %d:%d:%d", &c.year, &c.month, &c.day, &c.hour,
                                       &c.minute, &c.second);
        if (fields != 3 && fields != 5 && fields != 6) {
            fail(ctx, "'time_tdb' must be \"YYYY-MM-DD[ HH:MM[:SS]]\"");
        }
        return (jd_from_calendar(c) - kJ2000Jd) * kSecondsPerDay;
    }
    return parse_time(t, "time", ctx);
}

class Loader {
public:
    explicit Loader(std::filesystem::path asset_root)
        : m_asset_root(std::move(asset_root))
    {
    }

    void parse(const toml::table& root, Scene& out);

private:
    std::unique_ptr<MotionSource> parse_motion(const toml::table& t, const std::string& ctx, const Body* parent);
    void parse_nbody(const toml::table& t, Scene& out);
    std::shared_ptr<const EphemerisTable> load_table(const std::string& file, const std::string& ctx);
    void load_ring_profile(const std::string& file, const std::string& ctx, RingSystem& rings) const;
    void load_belt(const std::string& file, const std::string& ctx, SceneBelt& belt) const;
    int body_ref(const Scene& scene, const toml::table& t, std::string_view key, const std::string& ctx,
                 bool allow_sun) const;

    std::filesystem::path m_asset_root;
    Scene* m_scene = nullptr; // the scene being parsed
    std::map<std::string, std::shared_ptr<const EphemerisTable>> m_tables;
    std::vector<int> m_nbody_bodies; // bodies whose motion comes from [nbody]
    std::vector<std::pair<int, std::string>> m_pulsar_axes; // pulsar, body whose orbit normal is its spin axis
    int m_current_body = -1;
};

std::shared_ptr<const EphemerisTable> Loader::load_table(const std::string& file, const std::string& ctx)
{
    if (const auto it = m_tables.find(file); it != m_tables.end()) {
        return it->second;
    }
    auto table = std::make_shared<EphemerisTable>();
    std::string error;
    if (!table->load(m_asset_root / file, &error)) {
        fail(ctx, error);
    }
    m_tables[file] = table;
    return table;
}

// Belt file (tools/belts/make_belts.py), little endian:
//   char[8] "AXBELT1\0"; uint32 count; uint32 reserved; float64 epoch (JD TDB);
//   float32 a (au), e, i, node, arg_peri, mean anomaly (rad), H per object.
void Loader::load_belt(const std::string& file, const std::string& ctx, SceneBelt& belt) const
{
    std::ifstream in(m_asset_root / file, std::ios::binary);
    char magic[8];
    uint32_t counts[2];
    double epoch_jd = 0.0;
    if (!in || !in.read(magic, 8) || std::memcmp(magic, "AXBELT1", 8) != 0 ||
        !in.read(reinterpret_cast<char*>(counts), sizeof(counts)) ||
        !in.read(reinterpret_cast<char*>(&epoch_jd), sizeof(epoch_jd)) || counts[0] == 0) {
        fail(ctx, "cannot read belt '" + file + "'");
    }
    belt.epoch_tdb = (epoch_jd - kJ2000Jd) * kSecondsPerDay;
    belt.elements.resize(static_cast<size_t>(counts[0]) * SceneBelt::kStride);
    if (!in.read(reinterpret_cast<char*>(belt.elements.data()),
                 static_cast<std::streamsize>(belt.elements.size() * sizeof(float)))) {
        fail(ctx, "truncated belt '" + file + "'");
    }
}

// Ring profile file (tools/rings/make_saturn_rings.py), little endian:
//   char[8] "AXRING1\0"; float64 inner, outer radius (km); uint32 count; uint32 reserved;
//   float32 normal optical depth[count], evenly from inner to outer.
void Loader::load_ring_profile(const std::string& file, const std::string& ctx, RingSystem& rings) const
{
    std::ifstream in(m_asset_root / file, std::ios::binary);
    char magic[8];
    double radii[2];
    uint32_t counts[2];
    if (!in || !in.read(magic, 8) || std::memcmp(magic, "AXRING1", 8) != 0 ||
        !in.read(reinterpret_cast<char*>(radii), sizeof(radii)) ||
        !in.read(reinterpret_cast<char*>(counts), sizeof(counts)) || counts[0] < 2 || !(radii[1] > radii[0])) {
        fail(ctx, "cannot read ring profile '" + file + "'");
    }
    std::vector<float> tau(counts[0]);
    if (!in.read(reinterpret_cast<char*>(tau.data()), static_cast<std::streamsize>(tau.size() * sizeof(float)))) {
        fail(ctx, "truncated ring profile '" + file + "'");
    }
    rings.inner_km = radii[0];
    rings.outer_km = radii[1];
    rings.profile.clear();
    for (float t : tau) {
        rings.profile.emplace_back(std::max(t, 0.0f), 0.0f); // a thin sheet
    }
}

std::unique_ptr<MotionSource> Loader::parse_motion(const toml::table& t, const std::string& ctx, const Body* parent)
{
    const std::string type = get_string(t, "type", ctx);
    const double parent_gm = parent ? parent->gm_km3_s2 : 0.0;

    if (type == "mean_elements") {
        MeanElements el;
        el.epoch_tdb = (get_double_or(t, "epoch_jd_tdb", kJ2000Jd) - kJ2000Jd) * kSecondsPerDay;
        el.a_km = get_double(t, "a_km", ctx);
        el.e = get_double(t, "e", ctx);
        el.arg_peri_deg = get_double(t, "arg_peri_deg", ctx);
        el.mean_anomaly_deg = get_double(t, "mean_anomaly_deg", ctx);
        el.i_deg = get_double(t, "i_deg", ctx);
        el.node_deg = get_double(t, "node_deg", ctx);
        el.period_days = get_double(t, "period_days", ctx);
        el.apsis_period_years = get_double_or(t, "apsis_period_years", 0.0);
        el.node_period_years = get_double_or(t, "node_period_years", 0.0);
        el.mean_anomaly_accel_deg_per_day2 = get_double_or(t, "mean_anomaly_accel_deg_per_day2", 0.0);
        double pole[2] = {0.0, 90.0};
        if (t.contains("laplace_pole_ra_dec_deg")) {
            get_array(t, "laplace_pole_ra_dec_deg", ctx, pole, 2, 2);
        }
        el.pole_ra_deg = pole[0];
        el.pole_dec_deg = pole[1];
        if (el.e < 0.0 || el.e >= 1.0 || el.a_km <= 0.0 || el.period_days <= 0.0) {
            fail(ctx, "invalid orbital elements (need 0 <= e < 1, a > 0, period > 0)");
        }
        return std::make_unique<MeanElementOrbit>(el);
    }

    if (type == "jpl_approx_planet") {
        ApproxPlanetElements el;
        double pair[2];
        get_array(t, "a_au", ctx, pair, 2, 2);
        el.a_au = pair[0];
        el.a_rate = pair[1];
        get_array(t, "e", ctx, pair, 2, 2);
        el.e = pair[0];
        el.e_rate = pair[1];
        get_array(t, "i_deg", ctx, pair, 2, 2);
        el.i_deg = pair[0];
        el.i_rate = pair[1];
        get_array(t, "mean_longitude_deg", ctx, pair, 2, 2);
        el.mean_longitude_deg = pair[0];
        el.mean_longitude_rate = pair[1];
        get_array(t, "lon_peri_deg", ctx, pair, 2, 2);
        el.lon_peri_deg = pair[0];
        el.lon_peri_rate = pair[1];
        get_array(t, "node_deg", ctx, pair, 2, 2);
        el.node_deg = pair[0];
        el.node_rate = pair[1];
        el.b = get_double_or(t, "b", 0.0);
        el.c = get_double_or(t, "c", 0.0);
        el.s = get_double_or(t, "s", 0.0);
        el.f = get_double_or(t, "f", 0.0);
        return std::make_unique<JplApproxPlanetOrbit>(el);
    }

    if (type == "ephemeris") {
        auto table = load_table(get_string(t, "file", ctx), ctx);
        std::unique_ptr<MotionSource> fallback;
        if (const toml::table* fb = t["fallback"].as_table()) {
            fallback = parse_motion(*fb, ctx + " fallback", parent);
        }
        const std::string extrapolate = get_string_or(t, "extrapolate", "none");
        Extrapolation mode = Extrapolation::None;
        if (extrapolate == "linear") {
            mode = Extrapolation::Linear;
        } else if (extrapolate != "none") {
            fail(ctx, "extrapolate must be \"none\" or \"linear\"");
        }
        const double blend_days = get_double_or(t, "fallback_blend_days", 0.0);
        if (blend_days < 0.0 || (blend_days > 0.0 && !fallback)) {
            fail(ctx, "fallback_blend_days needs a fallback and must not be negative");
        }
        std::unique_ptr<MotionSource> motion = std::make_unique<EphemerisMotion>(
            std::move(table), std::move(fallback), mode, parent_gm, blend_days * kSecondsPerDay);
        if (t.contains("relative_to")) {
            // Data relative to a body defined earlier under the same parent, possibly further
            // down (e.g. Mercury for an orbiter, or the Earth under the Earth-Moon barycenter).
            const std::string anchor_name = get_string(t, "relative_to", ctx);
            std::vector<const MotionSource*> chain;
            int b = m_scene->find(anchor_name);
            while (b >= 0 && &m_scene->bodies[static_cast<size_t>(b)] != parent) {
                const Body& body = m_scene->bodies[static_cast<size_t>(b)];
                if (!body.motion) {
                    b = -1;
                    break;
                }
                chain.push_back(body.motion.get());
                b = body.parent;
            }
            if (b < 0 || chain.empty()) {
                fail(ctx, "relative_to needs a body defined earlier below the same parent, with orbits up to it");
            }
            motion = std::make_unique<OffsetMotion>(std::move(motion), std::move(chain));
        }
        return motion;
    }

    if (type == "fixed") {
        double v[3];
        get_array(t, "ra_dec_distance_pc", ctx, v, 3, 3);
        return std::make_unique<FixedMotion>(unit_from_ra_dec(v[0] * kDegToRad, v[1] * kDegToRad) * (v[2] * kParsecKm));
    }

    if (type == "kerr_geodesic") {
        if (!parent || parent->kind != BodyKind::BlackHole || parent->gm_km3_s2 <= 0.0) {
            fail(ctx, "kerr_geodesic needs a black-hole parent with gm_km3_s2");
        }
        const toml::table* el = t["elements"].as_table();
        if (!el) {
            fail(ctx, "kerr_geodesic needs [elements]");
        }
        KerrOrbitSetup setup;
        setup.gm_km3_s2 = parent->gm_km3_s2;
        setup.spin = get_double_or(t, "spin_override", parent->spin);
        setup.elements = parse_visual_orbit(*el, ctx + " elements");
        setup.step_eta = get_double_or(t, "step_eta", setup.step_eta);
        const double a = setup.elements.a_km;
        const double period = kTwoPi * std::sqrt(a * a * a / setup.gm_km3_s2);
        const std::string epoch = get_string_or(t, "elements_epoch", "pericentre");
        if (epoch == "pericentre") {
            setup.elements_epoch_tdb = setup.elements.t_peri_tdb;
        } else if (epoch == "previous_apocentre") {
            setup.elements_epoch_tdb = setup.elements.t_peri_tdb - 0.5 * period;
        } else {
            fail(ctx, "elements_epoch must be \"pericentre\" or \"previous_apocentre\"");
        }
        if (t.contains("span_years")) {
            double years[2];
            get_array(t, "span_years", ctx, years, 2, 2);
            return std::make_unique<KerrOrbitMotion>(setup, tdb_from_julian_year(years[0]),
                                                     tdb_from_julian_year(years[1]));
        }
        return KerrOrbitMotion::rolling(setup, get_double(t, "rolling_window_days", ctx));
    }

    // Fixed two-body orbits on the sky (the companion relative to the primary,
    // for the total gm_km3_s2); `scale` splits them around a barycenter.
    if (type == "visual_orbit" || type == "transit_orbit") {
        const double gm = get_double(t, "gm_km3_s2", ctx);
        if (gm <= 0.0) {
            fail(ctx, "gm_km3_s2 must be positive");
        }
        const VisualOrbit orbit = type == "visual_orbit"
                                      ? parse_visual_orbit(t, ctx, gm)
                                      : visual_orbit_from_transit(parse_transit_orbit(t, ctx), gm);
        return std::make_unique<VisualOrbitMotion>(orbit, gm);
    }

    if (type == "nbody") {
        m_nbody_bodies.push_back(m_current_body);
        return nullptr; // filled in by parse_nbody
    }

    fail(ctx, "unknown motion type '" + type + "'");
}

int Loader::body_ref(const Scene& scene, const toml::table& t, std::string_view key, const std::string& ctx,
                     bool allow_sun) const
{
    const std::string name = get_string(t, key, ctx);
    if (allow_sun && name == "sun") {
        const int star = scene.star_index();
        return star; // -1 when the sun is not a body: handled as "the sun"
    }
    const int index = scene.find(name);
    if (index < 0) {
        fail(ctx, "unknown body '" + name + "' in '" + std::string(key) + "'");
    }
    return index;
}

void Loader::parse(const toml::table& root, Scene& out)
{
    out = Scene{};
    m_scene = &out;
    out.name = get_string(root, "name", "scene");

    if (const toml::table* helio = root["origin_heliocentric"].as_table()) {
        out.origin_heliocentric = parse_motion(*helio, "origin_heliocentric", nullptr);
    }

    const toml::array* bodies = root["bodies"].as_array();
    if (!bodies || bodies->empty()) {
        fail("scene", "needs at least one [[bodies]] entry");
    }

    std::unordered_map<std::string, int> index_of;
    for (size_t i = 0; i < bodies->size(); ++i) {
        const toml::table* t = (*bodies)[i].as_table();
        if (!t) {
            fail("bodies", "entries must be tables");
        }
        Body body;
        body.name = get_string(*t, "name", "bodies[" + std::to_string(i) + "]");
        const std::string ctx = "body '" + body.name + "'";
        if (index_of.count(body.name)) {
            fail(ctx, "duplicate name");
        }

        if (const auto parent = (*t)["parent"].value<std::string>()) {
            const auto it = index_of.find(*parent);
            if (it == index_of.end()) {
                fail(ctx, "parent '" + *parent + "' must be defined earlier");
            }
            body.parent = it->second;
        } else if (i != 0) {
            fail(ctx, "only the first body may have no parent");
        }

        const std::string kind = get_string_or(*t, "kind", "planet");
        if (kind == "planet") {
            body.kind = BodyKind::Planet;
        } else if (kind == "star") {
            body.kind = BodyKind::Star;
        } else if (kind == "spacecraft") {
            body.kind = BodyKind::Spacecraft;
        } else if (kind == "black_hole") {
            body.kind = BodyKind::BlackHole;
        } else if (kind == "barycenter") {
            body.kind = BodyKind::Barycenter;
        } else {
            fail(ctx, "kind must be planet, star, spacecraft, black_hole or barycenter");
        }

        body.gm_km3_s2 = get_double_or(*t, "gm_km3_s2", 0.0);
        body.body_gm_km3_s2 = get_double_or(*t, "body_gm_km3_s2", body.gm_km3_s2);
        body.spin = get_double_or(*t, "spin", 0.0);
        body.temperature_k = get_double_or(*t, "temperature_k", body.temperature_k);
        body.disk_outer_m = get_double_or(*t, "disk_outer_m", 0.0);
        body.disk_temperature_k = get_double_or(*t, "disk_temperature_k", body.disk_temperature_k);
        body.disk_brightness = get_double_or(*t, "disk_brightness", body.disk_brightness);
        body.luminosity_solar = get_double_or(*t, "luminosity_solar", body.luminosity_solar);
        if (body.spin < 0.0 || body.spin >= 1.0) {
            fail(ctx, "spin must be in [0, 1)");
        }

        if (body.kind == BodyKind::BlackHole) {
            // Event horizon r+ = M + sqrt(M^2 - a^2) in Boyer-Lindquist coordinates.
            if (body.gm_km3_s2 <= 0.0) {
                fail(ctx, "black holes need gm_km3_s2");
            }
            const double m = body.gm_km3_s2 / (kSpeedOfLightKmS * kSpeedOfLightKmS);
            body.equatorial_radius_km = body.polar_radius_km = m + m * std::sqrt(1.0 - body.spin * body.spin);
            body.equatorial_radius_b_km = body.equatorial_radius_km;
        } else if (body.kind == BodyKind::Barycenter) {
            body.equatorial_radius_km = body.equatorial_radius_b_km = body.polar_radius_km = 1.0;
        } else {
            // [mean], [equatorial, polar] or [a, b, c] (a along the prime meridian).
            double radii[3] = {0.0, 0.0, 0.0};
            get_array(*t, "radii_km", ctx, radii, 1, 3);
            const bool triaxial = radii[2] > 0.0;
            body.equatorial_radius_km = radii[0];
            body.equatorial_radius_b_km = triaxial ? radii[1] : radii[0];
            body.polar_radius_km = triaxial ? radii[2] : (radii[1] > 0.0 ? radii[1] : radii[0]);
            if (body.equatorial_radius_b_km <= 0.0 || body.polar_radius_km <= 0.0 ||
                body.equatorial_radius_km < body.equatorial_radius_b_km) {
                fail(ctx, "radii must be positive (and a >= b for [a, b, c])");
            }
        }

        double pole[2] = {0.0, 90.0};
        if (t->contains("pole_ra_dec_deg")) {
            get_array(*t, "pole_ra_dec_deg", ctx, pole, 2, 2);
        }
        body.pole_ra = pole[0] * kDegToRad;
        body.pole_dec = pole[1] * kDegToRad;
        body.pole = unit_from_ra_dec(body.pole_ra, body.pole_dec);

        double pm[2] = {0.0, 0.0};
        if (t->contains("prime_meridian_deg")) {
            get_array(*t, "prime_meridian_deg", ctx, pm, 2, 2);
        }
        body.pm_w0_deg = pm[0];
        body.pm_rate_deg_per_day = pm[1];
        if (t->contains("pole_rate_deg_per_century")) {
            double rate[2];
            get_array(*t, "pole_rate_deg_per_century", ctx, rate, 2, 2);
            body.pole_rate_deg_per_century = glm::dvec2(rate[0], rate[1]);
        }
        if (const toml::array* angles = (*t)["nut_prec_angles"].as_array()) {
            for (const auto& node : *angles) {
                const toml::array* pair = node.as_array();
                if (!pair || pair->size() != 2 || !(*pair)[0].is_number() || !(*pair)[1].is_number()) {
                    fail(ctx, "nut_prec_angles must be [deg, deg per century] pairs");
                }
                body.nut_prec_angles.emplace_back((*pair)[0].value<double>().value_or(0.0),
                                                  (*pair)[1].value<double>().value_or(0.0));
            }
        }
        for (const auto& [key, terms] : {std::pair<const char*, std::vector<double>*>{"nut_prec_ra", &body.nut_prec_ra},
                                         {"nut_prec_dec", &body.nut_prec_dec},
                                         {"nut_prec_pm", &body.nut_prec_pm}}) {
            if (const toml::array* values = (*t)[key].as_array()) {
                for (const auto& node : *values) {
                    if (!node.is_number()) {
                        fail(ctx, std::string(key) + " must be numbers (degrees)");
                    }
                    terms->push_back(node.value<double>().value_or(0.0));
                }
            }
        }

        body.color = parse_color(*t, "color", ctx, glm::vec3(0.8f));
        body.orbit_color = parse_color(*t, "orbit_color", ctx, glm::vec3(0.6f));

        if (const toml::table* p = (*t)["pulsar"].as_table()) {
            if (body.kind != BodyKind::Star) {
                fail(ctx, "only stars can be pulsars");
            }
            body.pulsar.enabled = true;
            body.pulsar.display_period_s = get_double_or(*p, "display_period_s", body.pulsar.display_period_s);
            body.pulsar.magnetic_inclination_deg =
                get_double_or(*p, "magnetic_inclination_deg", body.pulsar.magnetic_inclination_deg);
            body.pulsar.beam_half_angle_deg = get_double_or(*p, "beam_half_angle_deg", body.pulsar.beam_half_angle_deg);
            body.pulsar.beam_color = parse_color(*p, "beam_color", ctx + " pulsar", body.pulsar.beam_color);
            if (body.pulsar.display_period_s <= 0.0 || body.pulsar.beam_half_angle_deg <= 0.0) {
                fail(ctx, "pulsar display_period_s and beam_half_angle_deg must be positive");
            }
            m_pulsar_axes.emplace_back(static_cast<int>(out.bodies.size()),
                                       get_string(*p, "spin_axis_orbit_normal", ctx + " pulsar"));
        }

        const std::string style = get_string_or(*t, "style", "solid");
        if (style == "solid") {
            body.style = SurfaceStyle::Solid;
        } else if (style == "gas_giant") {
            body.style = SurfaceStyle::GasGiantBands;
        } else if (style == "saturn") {
            body.style = SurfaceStyle::Saturn;
        } else if (style == "death_star") {
            body.style = SurfaceStyle::DeathStar;
        } else {
            fail(ctx, "unknown style '" + style + "'");
        }

        const std::string trail =
            get_string_or(*t, "trail", body.kind == BodyKind::Spacecraft ? "history" : "orbit");
        if (trail == "orbit") {
            body.trail = TrailMode::Orbit;
        } else if (trail == "history") {
            body.trail = TrailMode::History;
        } else if (trail == "none") {
            body.trail = TrailMode::None;
        } else {
            fail(ctx, "trail must be \"orbit\", \"history\" or \"none\"");
        }
        body.trail_history_days = get_double_or(*t, "trail_history_days", 0.0);
        body.trail_linger_days = get_double_or(*t, "trail_linger_days", 0.0);
        body.mark_periapsides = (*t)["mark_periapsides"].value<bool>().value_or(false);
        if (body.mark_periapsides && body.trail_history_days <= 0.0) {
            fail(ctx, "mark_periapsides needs trail_history_days (the span searched for periapsides)");
        }

        body.texture = get_string_or(*t, "texture", "");
        body.texture_left_lon_deg = get_double_or(*t, "texture_left_lon_deg", -180.0);
        const std::string direction = get_string_or(*t, "texture_lon_direction", "east");
        if (direction != "east" && direction != "west") {
            fail(ctx, "texture_lon_direction must be \"east\" or \"west\"");
        }
        body.texture_west_positive = direction == "west";

        if (const std::string shape = get_string_or(*t, "shape", ""); !shape.empty()) {
            auto model = std::make_shared<ShapeModel>();
            std::string error;
            if (!load_shape_model(m_asset_root / shape, *model, &error)) {
                fail(ctx, error);
            }
            if (!body.texture.empty() && model->map_u.empty()) {
                fail(ctx, "this shape model has no map coordinates for a texture");
            }
            body.shape = std::move(model);
        }

        if (const toml::table* rings = (*t)["rings"].as_table()) {
            const std::string rctx = ctx + " rings";
            body.rings.color = parse_color(*rings, "color", rctx, glm::vec3(1.0f));
            body.rings.gain = get_double_or(*rings, "gain", 1.0);
            body.rings.phase_g = get_double_or(*rings, "phase_g", 0.0);
            if (body.rings.gain < 0.0 || std::abs(body.rings.phase_g) >= 1.0) {
                fail(rctx, "gain must not be negative and |phase_g| must be < 1");
            }
            const toml::array* bands = (*rings)["bands"].as_array();
            const bool has_profile = rings->contains("profile");
            if ((!bands || bands->empty()) == !has_profile) {
                fail(rctx, "needs either a profile file or [[bodies.rings.bands]]");
            }
            if (has_profile) {
                load_ring_profile(get_string(*rings, "profile", rctx), rctx, body.rings);
                if (body.rings.inner_km < body.equatorial_radius_km) {
                    fail(rctx, "the profile starts inside the planet");
                }
            }
            for (size_t k = 0; bands && k < bands->size(); ++k) {
                const toml::node& node = (*bands)[k];
                const toml::table* b = node.as_table();
                if (!b) {
                    fail(rctx, "bands must be tables");
                }
                RingBand band;
                band.name = get_string_or(*b, "name", "");
                band.inner_km = get_double(*b, "inner_km", rctx);
                band.outer_km = get_double(*b, "outer_km", rctx);
                band.optical_depth = get_double(*b, "optical_depth", rctx);
                band.thickness_km = get_double_or(*b, "thickness_km", 0.0);
                if (!(band.inner_km >= body.equatorial_radius_km && band.outer_km > band.inner_km) ||
                    band.optical_depth < 0.0 || band.thickness_km < 0.0) {
                    fail(rctx, "band '" + band.name +
                                   "' needs radius <= inner_km < outer_km and non-negative depth and thickness");
                }
                body.rings.bands.push_back(band);
            }
            if (!body.rings.bands.empty()) {
                rasterize_ring_bands(body.rings, kBandProfileSamples);
            }
        }

        if (const toml::table* atmosphere = (*t)["atmosphere"].as_table()) {
            if (body.kind != BodyKind::Planet || body.shape) {
                fail(ctx, "only ellipsoidal planets and moons can have an atmosphere");
            }
            body.atmosphere = parse_atmosphere(*atmosphere, ctx + " atmosphere");
        }
        if (const toml::array* plumes = (*t)["plumes"].as_array()) {
            if (body.kind != BodyKind::Planet || body.shape) {
                fail(ctx, "only ellipsoidal planets and moons can have plumes");
            }
            for (const toml::node& node : *plumes) {
                const toml::table* p = node.as_table();
                if (!p) {
                    fail(ctx, "plumes must be tables");
                }
                body.plumes.push_back(parse_plume(*p, ctx));
            }
        }

        if (const toml::table* orbit = (*t)["orbit"].as_table()) {
            if (body.parent < 0) {
                fail(ctx, "the root body cannot have an orbit");
            }
            m_current_body = static_cast<int>(out.bodies.size());
            body.motion = parse_motion(*orbit, ctx + " orbit", &out.bodies[static_cast<size_t>(body.parent)]);
            // `scale` applies to the whole motion (fallback included).
            if (orbit->contains("scale")) {
                const double scale = get_double(*orbit, "scale", ctx);
                body.motion = std::make_unique<ScaledMotion>(std::move(body.motion), scale);
                // A pair split around its barycenter: the scale is -m_other / m_total
                // (or m_other / m_total), so the orbit's total GM gives the body's own.
                const double pair_gm = get_double_or(*orbit, "gm_km3_s2", 0.0);
                if (body.body_gm_km3_s2 <= 0.0 && pair_gm > 0.0 && std::abs(scale) < 1.0) {
                    body.body_gm_km3_s2 = pair_gm * (1.0 - std::abs(scale));
                }
            }
        } else if (body.parent >= 0) {
            fail(ctx, "bodies with a parent need an [orbit]");
        }

        if (!body.nut_prec_ra.empty() || !body.nut_prec_dec.empty() || !body.nut_prec_pm.empty()) {
            // The angles come from the nearest ancestor that defines them (or the body itself).
            const size_t count = std::max({body.nut_prec_ra.size(), body.nut_prec_dec.size(), body.nut_prec_pm.size()});
            int source = body.nut_prec_angles.empty() ? body.parent : static_cast<int>(out.bodies.size());
            while (source >= 0 && source != static_cast<int>(out.bodies.size()) &&
                   out.bodies[static_cast<size_t>(source)].nut_prec_angles.empty()) {
                source = out.bodies[static_cast<size_t>(source)].parent;
            }
            const size_t available = source < 0 ? 0
                                     : source == static_cast<int>(out.bodies.size())
                                         ? body.nut_prec_angles.size()
                                         : out.bodies[static_cast<size_t>(source)].nut_prec_angles.size();
            if (count > available) {
                fail(ctx, "nut_prec terms need as many nut_prec_angles (on the body or an ancestor)");
            }
            body.nut_prec_source = source;
        }

        index_of[body.name] = static_cast<int>(out.bodies.size());
        out.bodies.push_back(std::move(body));
    }

    for (const auto& [pulsar, axis_body] : m_pulsar_axes) {
        const int k = out.find(axis_body);
        if (k < 0 || !out.bodies[static_cast<size_t>(k)].motion) {
            fail("body '" + out.bodies[static_cast<size_t>(pulsar)].name + "' pulsar",
                 "spin_axis_orbit_normal needs a body with an orbit");
        }
        out.bodies[static_cast<size_t>(pulsar)].pulsar.spin_axis_orbit_of = k;
    }

    if (const toml::table* nbody = root["nbody"].as_table()) {
        parse_nbody(*nbody, out);
    } else if (!m_nbody_bodies.empty()) {
        fail("scene", "bodies with orbit type \"nbody\" need an [nbody] section");
    }

    // Frames: frames[0] is always inertial around the root.
    out.frames.push_back({out.bodies[0].name + "-centered (ICRF)", DisplayFrame::Type::Inertial, 0, -1, -1});
    if (const toml::array* frames = root["frames"].as_array()) {
        for (const auto& node : *frames) {
            const toml::table* t = node.as_table();
            if (!t) {
                fail("frames", "entries must be tables");
            }
            DisplayFrame f;
            f.name = get_string(*t, "name", "frames");
            const std::string ctx = "frame '" + f.name + "'";
            const std::string type = get_string(*t, "type", ctx);
            if (type == "inertial") {
                f.type = DisplayFrame::Type::Inertial;
                f.origin = body_ref(out, *t, "origin", ctx, false);
            } else if (type == "rotating") {
                f.type = DisplayFrame::Type::Rotating;
                f.primary = body_ref(out, *t, "primary", ctx, true);
                f.secondary = body_ref(out, *t, "secondary", ctx, false);
                f.origin = t->contains("origin") ? body_ref(out, *t, "origin", ctx, false) : f.secondary;
            } else {
                fail(ctx, "type must be \"inertial\" or \"rotating\"");
            }
            f.trail_history_days = get_double_or(*t, "trail_history_days", 0.0);
            if (f.trail_history_days < 0.0) {
                fail(ctx, "trail_history_days must not be negative");
            }
            out.frames.push_back(f);
        }
    }

    if (const toml::array* markers = root["markers"].as_array()) {
        for (const auto& node : *markers) {
            const toml::table* t = node.as_table();
            if (!t) {
                fail("markers", "entries must be tables");
            }
            Marker m;
            m.name = get_string(*t, "name", "markers");
            const std::string ctx = "marker '" + m.name + "'";
            if (get_string(*t, "type", ctx) != "lagrange") {
                fail(ctx, "type must be \"lagrange\"");
            }
            m.point = static_cast<int>(get_double(*t, "point", ctx));
            if (m.point != 1 && m.point != 2) {
                fail(ctx, "point must be 1 or 2");
            }
            m.primary = body_ref(out, *t, "primary", ctx, true);
            m.secondary = body_ref(out, *t, "secondary", ctx, false);
            m.mass_ratio = get_double(*t, "mass_ratio", ctx);
            out.markers.push_back(m);
        }
    }

    auto frame_ref = [&](const toml::table& t, const std::string& ctx) {
        const std::string name = get_string(t, "frame", ctx);
        for (size_t i = 0; i < out.frames.size(); ++i) {
            if (out.frames[i].name == name) {
                return static_cast<int>(i);
            }
        }
        fail(ctx, "unknown frame '" + name + "'");
    };

    if (const toml::array* events = root["events"].as_array()) {
        for (const auto& node : *events) {
            const toml::table* t = node.as_table();
            if (!t) {
                fail("events", "entries must be tables");
            }
            SceneEvent e;
            e.name = get_string(*t, "name", "events");
            const std::string ctx = "event '" + e.name + "'";
            e.t_tdb = parse_event_time(*t, ctx);
            e.focus = t->contains("focus") ? body_ref(out, *t, "focus", ctx, false) : -1;
            e.warp = get_double_or(*t, "warp", 0.0);
            e.distance_km = get_double_or(*t, "distance_km", 0.0);
            e.frame = t->contains("frame") ? frame_ref(*t, ctx) : -1;
            if (t->contains("from_orbit_normal")) {
                e.from_orbit_normal = body_ref(out, *t, "from_orbit_normal", ctx, false);
                if (!out.bodies[static_cast<size_t>(e.from_orbit_normal)].motion) {
                    fail(ctx, "from_orbit_normal needs a body with an orbit");
                }
            }
            if (t->contains("from_body")) {
                e.from_body = body_ref(out, *t, "from_body", ctx, false);
                if (e.focus < 0 || e.from_body == e.focus || e.from_orbit_normal >= 0) {
                    fail(ctx, "from_body needs a focus (other than itself) and no from_orbit_normal");
                }
            }
            e.phase_deg = get_double_or(*t, "phase_deg", 0.0);
            if (e.phase_deg != 0.0 && e.from_body < 0) {
                fail(ctx, "phase_deg needs from_body");
            }
            out.events.push_back(e);
        }
    }

    if (const toml::array* belts = root["belts"].as_array()) {
        for (size_t k = 0; k < belts->size(); ++k) {
            const toml::table* t = (*belts)[k].as_table();
            const std::string ctx = "belts[" + std::to_string(k) + "]";
            if (!t) {
                fail(ctx, "must be a table");
            }
            SceneBelt belt;
            belt.name = get_string(*t, "name", ctx);
            belt.color = parse_color(*t, "color", ctx, glm::vec3(1.0f));
            belt.brightness = get_double_or(*t, "brightness", 1.0);
            belt.point_size_px = get_double_or(*t, "point_size_px", 1.5);
            if (belt.brightness < 0.0 || belt.point_size_px <= 0.0) {
                fail(ctx, "brightness must not be negative, point_size_px must be positive");
            }
            load_belt(get_string(*t, "file", ctx), ctx, belt);
            out.belts.push_back(std::move(belt));
        }
    }

    if (const toml::table* sky = root["sky"].as_table()) {
        out.sky.milky_way = get_string_or(*sky, "milky_way", "");
        out.sky.milky_way_brightness = get_double_or(*sky, "milky_way_brightness", 1.0);
        if (out.sky.milky_way_brightness < 0.0) {
            fail("sky", "milky_way_brightness must not be negative");
        }
        if (const toml::table* c = (*sky)["cluster"].as_table()) {
            const std::string ctx = "sky.cluster";
            ClusterView view;
            double pair[2];
            get_array(*c, "center_ra_dec_deg", ctx, pair, 2, 2);
            view.center_ra_deg = pair[0];
            view.center_dec_deg = pair[1];
            get_array(*c, "viewer_ra_dec_deg", ctx, pair, 2, 2);
            view.viewer_ra_deg = pair[0];
            view.viewer_dec_deg = pair[1];
            view.distance_pc = get_double(*c, "distance_pc", ctx);
            view.viewer_depth_pc = get_double_or(*c, "viewer_depth_pc", 0.0);
            std::string error;
            if (!load_cluster_stars(m_asset_root / get_string(*c, "file", ctx), view, out.sky.stars, &error)) {
                fail(ctx, error);
            }
        }
    }

    if (const toml::table* view = root["view"].as_table()) {
        const std::string ctx = "view";
        if (view->contains("focus")) {
            out.view.focus = body_ref(out, *view, "focus", ctx, false);
        }
        out.view.distance_km = get_double_or(*view, "distance_km", 0.0);
        out.view.pitch_rad = get_double_or(*view, "pitch_deg", out.view.pitch_rad * kRadToDeg) * kDegToRad;
        out.view.yaw_rad = get_double_or(*view, "yaw_deg", out.view.yaw_rad * kRadToDeg) * kDegToRad;
        if (view->contains("frame")) {
            out.view.frame = frame_ref(*view, ctx);
        }
        const std::string start = get_string_or(*view, "start", "now");
        out.view.start_now = start == "now";
        if (!out.view.start_now) {
            out.view.start_tdb = parse_time(*view, "start", ctx);
        }
        out.view.warp = get_double_or(*view, "warp", out.view.warp);
        const std::string up = get_string_or(*view, "up", "");
        const std::string kOrbitNormal = "orbit_normal:";
        if (up == "ecliptic") {
            out.view.up_icrf = ecliptic_pole_icrf();
        } else if (up.rfind(kOrbitNormal, 0) == 0) {
            const int b = out.find(up.substr(kOrbitNormal.size()));
            if (b < 0 || !out.bodies[static_cast<size_t>(b)].motion) {
                fail(ctx, "orbit_normal: needs a body with an orbit");
            }
            const double t = out.view.start_now ? 0.0 : out.view.start_tdb;
            const State s = out.bodies[static_cast<size_t>(b)].motion->eval(t);
            out.view.up_icrf = glm::normalize(glm::cross(s.position, s.velocity));
        } else if (!up.empty()) {
            const int b = out.find(up);
            if (b < 0) {
                fail(ctx, "up must be \"ecliptic\" or a body name");
            }
            out.view.up_icrf = out.bodies[static_cast<size_t>(b)].pole;
        }
    }
    out.set_active_frame(out.view.frame);
}

void Loader::parse_nbody(const toml::table& t, Scene& out)
{
    const std::string ctx = "nbody";
    const double epoch = t.contains("epoch_jd_tdb")
                             ? (get_double(t, "epoch_jd_tdb", ctx) - kJ2000Jd) * kSecondsPerDay
                             : tdb_from_julian_year(get_double(t, "epoch_year", ctx));
    double span[2];
    get_array(t, "span_years", ctx, span, 2, 2);
    const double step = get_double(t, "step_days", ctx) * kSecondsPerDay;
    const double knots_per_orbit = get_double_or(t, "knots_per_orbit", 0.0);

    const toml::array* members = t["members"].as_array();
    if (!members || members->empty()) {
        fail(ctx, "needs [[nbody.members]]");
    }

    std::vector<int> body_of;
    std::vector<NBodyParticle> particles;
    std::vector<double> period_of; // initial osculating period around the relative_to centre (0: none)
    auto member_index = [&](const std::string& name, const std::string& mctx) {
        const int b = out.find(name);
        for (size_t k = 0; k < body_of.size(); ++k) {
            if (body_of[k] == b) {
                return k;
            }
        }
        fail(mctx, "'" + name + "' must be an earlier nbody member");
    };
    auto barycenter = [&](const std::vector<size_t>& idx) {
        State s;
        double gm = 0.0;
        for (size_t k : idx) {
            s.position += particles[k].position * particles[k].gm;
            s.velocity += particles[k].velocity * particles[k].gm;
            gm += particles[k].gm;
        }
        return std::make_pair(State{s.position / gm, s.velocity / gm}, gm);
    };

    for (const auto& node : *members) {
        const toml::table* m = node.as_table();
        if (!m) {
            fail(ctx, "members must be tables");
        }
        const std::string name = get_string(*m, "body", ctx + " member");
        const std::string mctx = "nbody member '" + name + "'";
        const int b = out.find(name);
        // The root may take part as the fixed origin (e.g. a planet host star);
        // every other member is a child of the root with orbit type "nbody".
        const bool is_root = b == 0;
        if (b < 0 ||
            (!is_root && std::find(m_nbody_bodies.begin(), m_nbody_bodies.end(), b) == m_nbody_bodies.end())) {
            fail(mctx, "must be the root or a body with orbit type \"nbody\"");
        }
        if ((!is_root && out.bodies[static_cast<size_t>(b)].parent != 0) || out.bodies[0].motion) {
            fail(mctx, "nbody bodies must be children of a fixed root");
        }

        NBodyParticle p;
        const double body_gm = out.bodies[static_cast<size_t>(b)].gm_km3_s2;
        p.gm = m->contains("gm_km3_s2") || body_gm <= 0.0 ? get_double(*m, "gm_km3_s2", mctx) : body_gm;
        if (out.bodies[static_cast<size_t>(b)].body_gm_km3_s2 <= 0.0) {
            out.bodies[static_cast<size_t>(b)].body_gm_km3_s2 = p.gm;
        }
        double period = 0.0;
        auto osculating_period = [](const State& relative, double gm_pair) {
            const double energy =
                0.5 * glm::dot(relative.velocity, relative.velocity) - gm_pair / glm::length(relative.position);
            if (energy >= 0.0) {
                return 0.0;
            }
            const double a = -gm_pair / (2.0 * energy);
            return kTwoPi * std::sqrt(a * a * a / gm_pair);
        };
        std::vector<size_t> rel;
        if (const toml::array* r = (*m)["relative_to"].as_array()) {
            for (const auto& n : *r) {
                rel.push_back(member_index(n.value_or(std::string()), mctx));
            }
        }
        if (const toml::table* ss = (*m)["sky_state"].as_table()) {
            // Barycentric Cartesian state as photodynamical codes (e.g. ELC) publish
            // it: x and y on the sky, z toward the observer; x points to position
            // angle node_deg (north through east), y 90 deg further. au, au/day.
            if (!rel.empty() || is_root) {
                fail(mctx, "sky_state is barycentric: no relative_to, and not for the root");
            }
            double pos[3];
            double vel[3];
            double sky[2];
            get_array(*ss, "position_au", mctx, pos, 3, 3);
            get_array(*ss, "velocity_au_per_day", mctx, vel, 3, 3);
            get_array(*ss, "sky_ra_dec_deg", mctx, sky, 2, 2);
            const double ra = sky[0] * kDegToRad;
            const double dec = sky[1] * kDegToRad;
            const double pa = get_double_or(*ss, "node_deg", 0.0) * kDegToRad;
            const glm::dvec3 north(-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec));
            const glm::dvec3 east(-std::sin(ra), std::cos(ra), 0.0);
            const glm::dvec3 x = north * std::cos(pa) + east * std::sin(pa);
            const glm::dvec3 y = east * std::cos(pa) - north * std::sin(pa);
            const glm::dmat3 to_icrf(x, y, glm::cross(x, y)); // x cross y points at the observer
            p.position = to_icrf * glm::dvec3(pos[0], pos[1], pos[2]) * kAuKm;
            p.velocity = to_icrf * glm::dvec3(vel[0], vel[1], vel[2]) * (kAuKm / kSecondsPerDay);
            // For knot thinning: the orbit around the earlier members (Jacobi-like).
            if (!particles.empty()) {
                std::vector<size_t> all(particles.size());
                for (size_t k = 0; k < all.size(); ++k) {
                    all[k] = k;
                }
                const auto [center, gm_center] = barycenter(all);
                period = osculating_period({p.position - center.position, p.velocity - center.velocity},
                                           gm_center + p.gm);
            }
        } else if (!rel.empty()) {
            const auto [center, gm_center] = barycenter(rel);
            const double gm_pair = gm_center + p.gm;
            State relative;
            if (const toml::table* vo = (*m)["visual_orbit"].as_table()) {
                relative = visual_orbit_state(parse_visual_orbit(*vo, mctx), gm_pair, epoch);
            } else if (const toml::table* to = (*m)["transit_orbit"].as_table()) {
                const VisualOrbit orbit = visual_orbit_from_transit(parse_transit_orbit(*to, mctx), gm_pair);
                relative = visual_orbit_state(orbit, gm_pair, epoch);
            } else if (const toml::table* gs = (*m)["galactic_state"].as_table()) {
                double pos[3];
                double vel[3];
                get_array(*gs, "position_pc", mctx, pos, 3, 3);
                get_array(*gs, "velocity_km_s", mctx, vel, 3, 3);
                const glm::dmat3 g = galactic_to_icrf();
                relative.position = g * glm::dvec3(pos[0], pos[1], pos[2]) * kParsecKm;
                relative.velocity = g * glm::dvec3(vel[0], vel[1], vel[2]);
            } else {
                fail(mctx, "relative members need [visual_orbit], [transit_orbit] or [galactic_state]");
            }
            p.position = center.position + relative.position;
            p.velocity = center.velocity + relative.velocity;
            period = osculating_period(relative, gm_pair);
        } else if (is_root && !body_of.empty()) {
            fail(mctx, "the root must be the first member");
        }
        body_of.push_back(b);
        particles.push_back(p);
        period_of.push_back(period);
    }
    const bool root_member = body_of.front() == 0;
    if (body_of.size() != m_nbody_bodies.size() + (root_member ? 1 : 0)) {
        fail(ctx, "every body with orbit type \"nbody\" must be listed as a member");
    }

    // Positions are reported relative to the barycentre of the `origin` members.
    NBodyOptions options;
    options.order = static_cast<int>(get_double_or(t, "order", 4.0));
    if (options.order != 4 && options.order != 6) {
        fail(ctx, "order must be 4 or 6");
    }
    if (const toml::array* o = t["origin"].as_array()) {
        for (const auto& n : *o) {
            options.origin.push_back(member_index(n.value_or(std::string()), ctx + " origin"));
        }
    }
    if (options.origin.empty()) {
        for (size_t k = 0; k < particles.size(); ++k) {
            options.origin.push_back(k);
        }
    }
    // The root does not move, so it must be the origin.
    const bool root_is_origin = root_member && options.origin.size() == 1 && options.origin[0] == 0;
    if (root_member && !root_is_origin) {
        fail(ctx, "when the root is a member, origin must be [root]");
    }

    // Thin the output to about knots_per_orbit knots per revolution of each
    // member; members without an orbit of their own get the densest stride.
    if (knots_per_orbit > 0.0) {
        auto stride_for = [&](double period) {
            return std::max(1, static_cast<int>(period / (knots_per_orbit * step)));
        };
        int densest = 0;
        for (double period : period_of) {
            if (period > 0.0) {
                densest = densest == 0 ? stride_for(period) : std::min(densest, stride_for(period));
            }
        }
        for (double period : period_of) {
            options.stride.push_back(period > 0.0 ? stride_for(period) : std::max(1, densest));
        }
    }

    auto knots = integrate_nbody(particles, epoch, tdb_from_julian_year(span[0]), tdb_from_julian_year(span[1]),
                                 step, options);
    for (size_t i = 0; i < particles.size(); ++i) {
        if (body_of[i] == 0) {
            continue; // the root stays at the origin
        }
        auto table = std::make_shared<EphemerisTable>();
        table->set_knots(std::move(knots[i]));
        // Around a single central star, the osculating orbit makes a smooth trail.
        const double parent_gm = root_is_origin ? particles[0].gm + particles[i].gm : 0.0;
        Body& body = out.bodies[static_cast<size_t>(body_of[i])];
        body.motion = std::make_unique<EphemerisMotion>(std::move(table), nullptr, Extrapolation::None, parent_gm);
    }
}

} // namespace

bool load_scene_string(std::string_view text, std::string_view source_name, Scene& out, std::string* error,
                       const std::filesystem::path& asset_root)
{
    try {
        const toml::table root = toml::parse(text, source_name);
        Loader(asset_root).parse(root, out);
        return true;
    } catch (const toml::parse_error& e) {
        if (error) {
            std::ostringstream msg;
            msg << source_name << ":" << e.source().begin.line << ": " << e.description();
            *error = msg.str();
        }
    } catch (const SchemaError& e) {
        if (error) {
            *error = std::string(source_name) + ": " + e.message;
        }
    }
    return false;
}

bool load_scene_file(const std::filesystem::path& path, Scene& out, std::string* error)
{
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        if (error) {
            *error = "cannot open " + path.string();
        }
        return false;
    }
    std::ostringstream ss;
    ss << file.rdbuf();
    // Scenes live in <assets>/scenes/; other asset paths are relative to <assets>.
    return load_scene_string(ss.str(), path.string(), out, error, path.parent_path().parent_path());
}

} // namespace astraxis
