#include "ephem/ephemeris.hpp"

#include "ephem/kepler.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>

namespace astraxis {

namespace {

constexpr char kMagic[8] = {'A', 'X', 'E', 'P', 'H', '1', '\0', '\0'};

State hermite(const EphemerisTable::Knot& k0, const EphemerisTable::Knot& k1, double t)
{
    const double h = k1.t - k0.t;
    const double s = (t - k0.t) / h;
    const double s2 = s * s;
    const double s3 = s2 * s;
    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + s;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;
    // Derivatives with respect to t.
    const double d00 = (6.0 * s2 - 6.0 * s) / h;
    const double d10 = 3.0 * s2 - 4.0 * s + 1.0;
    const double d01 = (-6.0 * s2 + 6.0 * s) / h;
    const double d11 = 3.0 * s2 - 2.0 * s;

    State out;
    out.position = h00 * k0.position + h10 * h * k0.velocity + h01 * k1.position + h11 * h * k1.velocity;
    out.velocity = d00 * k0.position + d10 * k0.velocity + d01 * k1.position + d11 * k1.velocity;
    return out;
}

} // namespace

bool EphemerisTable::load(const std::filesystem::path& path, std::string* error)
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
    if (!file.read(magic, 8) || std::memcmp(magic, kMagic, 8) != 0) {
        return fail("not an AXEPH1 file");
    }
    if (!file.read(reinterpret_cast<char*>(ids), sizeof(ids)) ||
        !file.read(reinterpret_cast<char*>(counts), sizeof(counts)) || counts[0] < 2) {
        return fail("truncated header");
    }

    std::vector<double> raw(static_cast<size_t>(counts[0]) * 7);
    if (!file.read(reinterpret_cast<char*>(raw.data()), static_cast<std::streamsize>(raw.size() * sizeof(double)))) {
        return fail("truncated knot data");
    }

    m_target_id = ids[0];
    m_center_id = ids[1];
    m_knots.resize(counts[0]);
    for (size_t i = 0; i < m_knots.size(); ++i) {
        const double* k = &raw[i * 7];
        m_knots[i] = {k[0], {k[1], k[2], k[3]}, {k[4], k[5], k[6]}};
        if (i > 0 && !(m_knots[i].t > m_knots[i - 1].t)) {
            return fail("knot times are not increasing");
        }
    }
    return true;
}

State EphemerisTable::eval(double t) const
{
    if (t <= m_knots.front().t) {
        return {m_knots.front().position, m_knots.front().velocity};
    }
    if (t >= m_knots.back().t) {
        return {m_knots.back().position, m_knots.back().velocity};
    }
    const auto it = std::upper_bound(m_knots.begin(), m_knots.end(), t,
                                     [](double value, const Knot& k) { return value < k.t; });
    const size_t i = static_cast<size_t>(it - m_knots.begin());
    return hermite(m_knots[i - 1], m_knots[i], t);
}

EphemerisMotion::EphemerisMotion(std::shared_ptr<const EphemerisTable> table, std::unique_ptr<MotionSource> fallback,
                                 Extrapolation extrapolation, double parent_gm)
    : m_table(std::move(table))
    , m_fallback(std::move(fallback))
    , m_extrapolation(extrapolation)
    , m_parent_gm(parent_gm)
{
}

State EphemerisMotion::eval(double t_tdb) const
{
    if (m_table->covers(t_tdb)) {
        return m_table->eval(t_tdb);
    }
    if (m_fallback) {
        return m_fallback->eval(t_tdb);
    }
    if (m_extrapolation == Extrapolation::Linear && t_tdb > m_table->end()) {
        const EphemerisTable::Knot& k = m_table->knots().back();
        return {k.position + k.velocity * (t_tdb - k.t), k.velocity};
    }
    return m_table->eval(t_tdb); // clamped
}

bool EphemerisMotion::valid_at(double t_tdb) const
{
    if (m_table->covers(t_tdb) || m_fallback) {
        return true;
    }
    return m_extrapolation == Extrapolation::Linear && t_tdb > m_table->end();
}

bool EphemerisMotion::sample_orbit(double t_tdb, int count, std::vector<glm::dvec3>& out) const
{
    if (m_parent_gm <= 0.0) {
        return false;
    }
    return sample_osculating_ellipse(eval(t_tdb), m_parent_gm, count, out);
}

void EphemerisMotion::history_times(double t0, double t1, int max_points, std::vector<double>& out) const
{
    m_table->history_times(t0, t1, max_points, out);
}

void EphemerisTable::history_times(double t0, double t1, int max_points, std::vector<double>& out) const
{
    out.clear();
    if (t1 <= t0) {
        return;
    }
    const auto& knots = m_knots;
    const auto first = std::upper_bound(knots.begin(), knots.end(), t0,
                                        [](double value, const EphemerisTable::Knot& k) { return value < k.t; });
    const auto last = std::lower_bound(knots.begin(), knots.end(), t1,
                                       [](const EphemerisTable::Knot& k, double value) { return k.t < value; });
    const size_t inner = first < last ? static_cast<size_t>(last - first) : 0;
    const size_t budget = max_points > 2 ? static_cast<size_t>(max_points - 2) : 0;

    out.push_back(t0);
    if (inner > 0 && budget > 0) {
        // Keep every knot if possible; otherwise an even stride (always keeping the
        // knot closest to t1, so the head of the trail stays accurate).
        const size_t stride = (inner + budget - 1) / budget;
        const size_t offset = (inner - 1) % stride;
        for (size_t i = offset; i < inner; i += stride) {
            out.push_back((first + static_cast<std::ptrdiff_t>(i))->t);
        }
    }
    // Between knots beyond the table, fill uniformly (fallback / extrapolation).
    if (inner == 0 && max_points > 2) {
        for (int i = 1; i < max_points - 1; ++i) {
            out.push_back(t0 + (t1 - t0) * static_cast<double>(i) / static_cast<double>(max_points - 1));
        }
    }
    out.push_back(t1);
}

} // namespace astraxis
