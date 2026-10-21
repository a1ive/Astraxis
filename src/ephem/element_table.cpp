#include "ephem/element_table.hpp"

#include "core/math.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>

namespace astraxis {

namespace {

constexpr char kMagic[8] = {'A', 'X', 'E', 'L', 'M', '1', '\0', '\0'};

} // namespace

bool ElementTable::load(const std::filesystem::path& path, std::string* error)
{
    auto fail = [&](const std::string& what) {
        if (error) {
            *error = path.string() + ": " + what;
        }
        return false;
    };

    std::ifstream file(path, std::ios::binary);
    if (!file) {
        return fail("cannot open");
    }
    char magic[8];
    int32_t ids[2];
    uint32_t counts[2];
    double header[4];
    if (!file.read(magic, 8) || std::memcmp(magic, kMagic, 8) != 0) {
        return fail("not an AXELM1 file");
    }
    if (!file.read(reinterpret_cast<char*>(ids), sizeof(ids)) ||
        !file.read(reinterpret_cast<char*>(counts), sizeof(counts)) || counts[0] < 2 ||
        !file.read(reinterpret_cast<char*>(header), sizeof(header))) {
        return fail("truncated header");
    }
    std::vector<double> raw(static_cast<size_t>(counts[0]) * 8);
    if (!file.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size() * sizeof(double)))) {
        return fail("truncated knot data");
    }

    m_target_id = ids[0];
    m_center_id = ids[1];
    m_plane_to_icrf = frame_from_pole(unit_from_ra_dec(header[0] * kDegToRad, header[1] * kDegToRad));
    m_apsis_rate = header[2];
    m_node_rate = header[3];
    m_knots.resize(counts[0]);
    for (size_t i = 0; i < m_knots.size(); ++i) {
        const double* k = &raw[i * 8];
        m_knots[i] = {k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]};
        if (i > 0 && !(m_knots[i].t > m_knots[i - 1].t)) {
            return fail("knot times are not increasing");
        }
        if (!(m_knots[i].a > 0.0) || !(m_knots[i].e >= 0.0 && m_knots[i].e < 1.0)) {
            return fail("bad elements");
        }
    }
    return true;
}

ElementTable::Knot ElementTable::at(double t, double* lambda_rate) const
{
    if (t <= start() || t >= end()) {
        Knot k = t <= start() ? m_knots.front() : m_knots.back();
        const double dt = t - k.t;
        k.lambda += k.n * dt;
        k.varpi += m_apsis_rate * dt;
        k.node += m_node_rate * dt;
        k.t = t;
        if (lambda_rate) {
            *lambda_rate = k.n;
        }
        return k;
    }
    const auto it =
        std::upper_bound(m_knots.begin(), m_knots.end(), t, [](double value, const Knot& k) { return value < k.t; });
    const Knot& k0 = *(it - 1);
    const Knot& k1 = *it;
    const double h = k1.t - k0.t;
    const double u = (t - k0.t) / h;
    auto lerp = [u](double a, double b) { return a + u * (b - a); };

    Knot out;
    out.t = t;
    out.a = lerp(k0.a, k1.a);
    out.e = lerp(k0.e, k1.e);
    out.varpi = lerp(k0.varpi, k1.varpi);
    out.i = lerp(k0.i, k1.i);
    out.node = lerp(k0.node, k1.node);
    out.n = lerp(k0.n, k1.n);
    // Cubic Hermite through (lambda, n) at both knots.
    const double u2 = u * u;
    const double u3 = u2 * u;
    out.lambda = (2.0 * u3 - 3.0 * u2 + 1.0) * k0.lambda + (u3 - 2.0 * u2 + u) * h * k0.n +
                 (-2.0 * u3 + 3.0 * u2) * k1.lambda + (u3 - u2) * h * k1.n;
    if (lambda_rate) {
        *lambda_rate = (6.0 * u2 - 6.0 * u) / h * k0.lambda + (3.0 * u2 - 4.0 * u + 1.0) * k0.n +
                       (-6.0 * u2 + 6.0 * u) / h * k1.lambda + (3.0 * u2 - 2.0 * u) * k1.n;
    }
    return out;
}

KeplerElements ElementTableOrbit::elements_at(double t_tdb, double* mean_anomaly_rate) const
{
    double lambda_rate = 0.0;
    const ElementTable::Knot k = m_table->at(t_tdb, &lambda_rate);
    KeplerElements el;
    el.a = k.a;
    el.e = k.e;
    el.i = k.i;
    el.node = wrap_two_pi(k.node);
    el.arg_peri = wrap_two_pi(k.varpi - k.node);
    el.mean_anomaly = wrap_two_pi(k.lambda - k.varpi);
    if (mean_anomaly_rate) {
        *mean_anomaly_rate = lambda_rate - m_table->apsis_rate();
    }
    return el;
}

State ElementTableOrbit::eval(double t_tdb) const
{
    double mean_anomaly_rate = 0.0;
    const KeplerElements el = elements_at(t_tdb, &mean_anomaly_rate);
    State s = kepler_state(el, mean_anomaly_rate);
    // The ellipse turns about the plane's pole as the periapsis advances (the nodal
    // regression of these nearly equatorial orbits barely moves them).
    s.velocity += m_table->apsis_rate() * glm::dvec3(-s.position.y, s.position.x, 0.0);
    const glm::dmat3& r = m_table->plane_to_icrf();
    return {r * s.position, r * s.velocity};
}

bool ElementTableOrbit::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    const KeplerElements el = elements_at(t_tdb);
    sample_ellipse(el, solve_kepler(el.mean_anomaly, el.e), count, out);
    for (glm::dvec3& p : out) {
        p = m_table->plane_to_icrf() * p;
    }
    return true;
}

} // namespace astraxis
