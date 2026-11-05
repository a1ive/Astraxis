#include "scene/magnetic_field.hpp"

#include "core/math.hpp"
#include "core/time.hpp"

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <locale>
#include <sstream>

namespace astraxis {

namespace {

constexpr double kStepOfRadius = 0.04; // RK4 step along the line, as a fraction of r
constexpr int kMaxSteps = 4000;

using Point = glm::dvec3;

// Unit field direction, or false near a null.
bool direction_at(const MagneticField& field, const Point& x, double sign, glm::dvec3& out)
{
    const glm::dvec3 b = field.field_nt(x);
    const double length = glm::length(b);
    if (!(length > 1e-6)) {
        return false;
    }
    out = sign * b / length;
    return true;
}

// How far outside the surface a point is: negative inside.
double surface_level(const Point& x, const FieldLineSettings& s)
{
    const double e = s.surface_equatorial;
    const double p = s.surface_polar;
    return (x.x * x.x + x.y * x.y) / (e * e) + x.z * x.z / (p * p) - 1.0;
}

// From the seed along the field until the surface or max_radius; ends on them.
void trace_one_way(const MagneticField& field, const Point& seed, double sign, const FieldLineSettings& s,
                   std::vector<Point>& points)
{
    Point x = seed;
    for (int i = 0; i < kMaxSteps; ++i) {
        const double h = kStepOfRadius * glm::length(x);
        glm::dvec3 k1, k2, k3, k4;
        if (!direction_at(field, x, sign, k1) || !direction_at(field, x + 0.5 * h * k1, sign, k2) ||
            !direction_at(field, x + 0.5 * h * k2, sign, k3) || !direction_at(field, x + h * k3, sign, k4)) {
            return;
        }
        const Point next = x + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        const double level_now = surface_level(x, s);
        const double level_next = surface_level(next, s);
        if (level_next <= 0.0) {
            // Linear in the level is enough over a step this short.
            const double t = level_now / std::max(level_now - level_next, 1e-12);
            points.push_back(x + (next - x) * t);
            return;
        }
        const double r_now = glm::length(x);
        const double r_next = glm::length(next);
        if (r_next >= s.max_radius) {
            const double t = (s.max_radius - r_now) / std::max(r_next - r_now, 1e-12);
            points.push_back(x + (next - x) * t);
            return;
        }
        points.push_back(next);
        x = next;
    }
}

// Ramer-Douglas-Peucker, iterative: the lines have thousands of points at most.
std::vector<Point> simplify(const std::vector<Point>& points, double tolerance)
{
    if (points.size() < 3) {
        return points;
    }
    std::vector<char> keep(points.size(), 0);
    keep.front() = keep.back() = 1;
    std::vector<std::pair<size_t, size_t>> stack{{0, points.size() - 1}};
    while (!stack.empty()) {
        const auto [first, last] = stack.back();
        stack.pop_back();
        const glm::dvec3 chord = points[last] - points[first];
        const double length = glm::length(chord);
        double worst = 0.0;
        size_t at = first;
        for (size_t i = first + 1; i < last; ++i) {
            // A closed line's chord can vanish: then the distance from its end.
            const glm::dvec3 d = points[i] - points[first];
            const double e = length > 1e-12 ? glm::length(glm::cross(d, chord)) / length : glm::length(d);
            if (e > worst) {
                worst = e;
                at = i;
            }
        }
        if (worst > tolerance) {
            keep[at] = 1;
            stack.emplace_back(first, at);
            stack.emplace_back(at, last);
        }
    }
    std::vector<Point> out;
    for (size_t i = 0; i < points.size(); ++i) {
        if (keep[i]) {
            out.push_back(points[i]);
        }
    }
    return out;
}

} // namespace

bool MagneticField::load(const std::filesystem::path& file, std::string* error)
{
    const auto failed = [&](const std::string& what) {
        if (error) {
            *error = file.string() + ": " + what;
        }
        return false;
    };
    std::ifstream in(file);
    if (!in) {
        return failed("cannot open");
    }
    m_degree = 0;
    m_g.clear();
    m_h.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream row(line);
        std::string kind;
        int degree = 0, order = 0;
        double value = 0.0;
        if (!(row >> kind >> degree >> order >> value) || (kind != "g" && kind != "h") || degree < 1 ||
            degree > kMaxDegree || order < 0 || order > degree) {
            return failed("bad line: " + line);
        }
        set_coefficient(kind == "h", degree, order, value);
    }
    if (m_degree < 1 || dipole_nt() <= 0.0) {
        return failed("no dipole");
    }
    return true;
}

void MagneticField::set_coefficient(bool is_h, int degree, int order, double value)
{
    m_degree = std::max(m_degree, degree);
    const size_t size = static_cast<size_t>(index(m_degree, m_degree)) + 1;
    if (m_g.size() < size) {
        m_g.resize(size, 0.0);
        m_h.resize(size, 0.0);
    }
    (is_h ? m_h : m_g)[static_cast<size_t>(index(degree, order))] = value;
}

glm::dvec3 MagneticField::field_nt(const glm::dvec3& x) const
{
    const double r = glm::length(x);
    const double cos_theta = x.z / r;
    const double sin_theta = std::max(std::sqrt(std::max(1.0 - cos_theta * cos_theta, 0.0)), 1e-9);
    const double phi = std::atan2(x.y, x.x);
    const double cos_phi = std::cos(phi);
    const double sin_phi = std::sin(phi);

    // Schmidt functions P[n][m], by rows; cos(m phi), sin(m phi) by recursion.
    std::array<std::array<double, kMaxDegree + 1>, kMaxDegree + 1> p{};
    p[0][0] = 1.0;
    double radial = 0.0, polar = 0.0, azimuthal = 0.0;
    double factor = 1.0 / (r * r); // (a/r)^(n+2) with a = 1
    for (int n = 1; n <= m_degree; ++n) {
        const double dn = n;
        for (int m = 0; m < n; ++m) {
            const double dm = m;
            const double above = n - 1 > m ? std::sqrt((dn - 1) * (dn - 1) - dm * dm) * p[n - 2][m] : 0.0;
            p[n][m] = ((2.0 * dn - 1.0) * cos_theta * p[n - 1][m] - above) / std::sqrt(dn * dn - dm * dm);
        }
        p[n][n] = sin_theta * p[n - 1][n - 1] * (n == 1 ? 1.0 : std::sqrt((2.0 * dn - 1.0) / (2.0 * dn)));

        factor /= r;
        double cos_m = 1.0, sin_m = 0.0;
        for (int m = 0; m <= n; ++m) {
            const double dm = m;
            const double gnm = m_g[static_cast<size_t>(index(n, m))];
            const double hnm = m_h[static_cast<size_t>(index(n, m))];
            const double along = gnm * cos_m + hnm * sin_m;
            const double previous = m < n ? p[n - 1][m] : 0.0;
            const double derivative = (dn * cos_theta * p[n][m] - std::sqrt(dn * dn - dm * dm) * previous) / sin_theta;
            radial += (dn + 1.0) * factor * along * p[n][m];
            polar -= factor * along * derivative;
            azimuthal += dm * factor * (gnm * sin_m - hnm * cos_m) * p[n][m] / sin_theta;
            const double next_cos = cos_m * cos_phi - sin_m * sin_phi;
            sin_m = sin_m * cos_phi + cos_m * sin_phi;
            cos_m = next_cos;
        }
    }
    return {radial * sin_theta * cos_phi + polar * cos_theta * cos_phi - azimuthal * sin_phi,
            radial * sin_theta * sin_phi + polar * cos_theta * sin_phi + azimuthal * cos_phi,
            radial * cos_theta - polar * sin_theta};
}

double MagneticField::dipole_nt() const
{
    return m_g.size() > 2 ? std::sqrt(g(1, 0) * g(1, 0) + g(1, 1) * g(1, 1) + h(1, 1) * h(1, 1)) : 0.0;
}

glm::dvec3 MagneticField::dipole_axis() const
{
    return glm::dvec3(g(1, 1), h(1, 1), g(1, 0)) / dipole_nt();
}

double MagneticField::dipole_tilt_deg() const
{
    return std::acos(std::min(std::abs(g(1, 0)) / dipole_nt(), 1.0)) / kDegToRad;
}

glm::dvec3 MagneticField::dipole_offset() const
{
    if (m_degree < 2) {
        return glm::dvec3(0.0);
    }
    const double g10 = g(1, 0), g11 = g(1, 1), h11 = h(1, 1);
    const double g20 = g(2, 0), g21 = g(2, 1), g22 = g(2, 2), h21 = h(2, 1), h22 = h(2, 2);
    const double m2 = g10 * g10 + g11 * g11 + h11 * h11;
    const double s3 = std::sqrt(3.0);
    const double l0 = 2.0 * g10 * g20 + s3 * (g11 * g21 + h11 * h21);
    const double l1 = -g11 * g20 + s3 * (g10 * g21 + g11 * g22 + h11 * h22);
    const double l2 = -h11 * g20 + s3 * (g10 * h21 - h11 * g22 + g11 * h22);
    const double e = (l0 * g10 + l1 * g11 + l2 * h11) / (4.0 * m2);
    return glm::dvec3(l1 - g11 * e, l2 - h11 * e, l0 - g10 * e) / (3.0 * m2);
}

std::vector<glm::dvec4> trace_field_line(const MagneticField& field, const glm::dvec3& seed,
                                         const FieldLineSettings& s)
{
    std::vector<Point> backward{seed};
    trace_one_way(field, seed, -1.0, s, backward);
    std::vector<Point> points(backward.rbegin(), backward.rend());
    trace_one_way(field, seed, 1.0, s, points);

    const std::vector<Point> kept = simplify(points, s.tolerance);
    std::vector<glm::dvec4> out;
    out.reserve(kept.size());
    // Lines fade over the outer half of the way from the surface to max_radius, also where
    // that is short (Mercury, Ganymede).
    const double fade_start = 0.5 * (s.surface_equatorial + s.max_radius);
    for (const Point& p : kept) {
        const double fade = std::clamp((glm::length(p) - fade_start) / (s.max_radius - fade_start), 0.0, 1.0);
        out.emplace_back(p, fade);
    }
    return out;
}

std::vector<std::vector<glm::dvec4>> trace_field_lines(const MagneticField& field, const FieldLineSettings& s)
{
    const glm::dvec3 axis = field.dipole_axis();
    const glm::dvec3 offset = field.dipole_offset();
    const glm::dvec3 reference = std::abs(axis.y) < 0.9 ? glm::dvec3(0.0, 1.0, 0.0) : glm::dvec3(1.0, 0.0, 0.0);
    const glm::dvec3 e1 = glm::normalize(glm::cross(axis, reference));
    const glm::dvec3 e2 = glm::cross(axis, e1);
    std::vector<std::vector<glm::dvec4>> lines;
    for (size_t shell = 0; shell < s.l_shells.size(); ++shell) {
        for (int k = 0; k < s.lines_per_shell; ++k) {
            // Alternate shells are staggered, so that the lines do not stack up.
            const double angle = kTwoPi * (k + (shell % 2 == 1 ? 0.5 : 0.0)) / s.lines_per_shell;
            const glm::dvec3 seed = offset + s.l_shells[shell] * (std::cos(angle) * e1 + std::sin(angle) * e2);
            if (surface_level(seed, s) <= 0.0 || glm::length(seed) >= s.max_radius) {
                continue;
            }
            lines.push_back(trace_field_line(field, seed, s));
        }
    }
    return lines;
}

bool FieldMap::load(const std::filesystem::path& file, std::string* error)
{
    const auto failed = [&](const std::string& what) {
        if (error) {
            *error = file.string() + ": " + what;
        }
        return false;
    };
    std::ifstream in(file, std::ios::binary);
    if (!in) {
        return failed("cannot open");
    }
    std::string header;
    std::getline(in, header);
    std::istringstream row(header);
    row.imbue(std::locale::classic());
    std::string magic;
    int w = 0, h = 0;
    double low = 0.0, scale = 0.0;
    row >> magic >> w >> h;
    const bool strength = magic == "AXFIELDMAG1";
    if (strength) {
        row >> low;
    }
    if (!(row >> scale) || (magic != "AXFIELDMAP1" && !strength) || w <= 0 || h <= 0 || w > 8192 || h > 4096 ||
        !(scale > low) || low < 0.0) {
        return failed("bad header");
    }
    std::vector<signed char> values(static_cast<size_t>(w) * static_cast<size_t>(h));
    in.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(values.size()));
    if (in.gcount() != static_cast<std::streamsize>(values.size()) || in.peek() != std::char_traits<char>::eof()) {
        return failed("wrong size");
    }
    width = w;
    height = h;
    magnitude = strength;
    scale_nt = scale;
    min_nt = low;
    cells = std::move(values);
    return true;
}

double FieldMap::value_nt(double east_longitude, double latitude) const
{
    if (cells.empty()) {
        return 0.0;
    }
    const double u = east_longitude / kTwoPi + 0.5;
    const int column = static_cast<int>(std::floor((u - std::floor(u)) * width)) % width;
    const int line = std::clamp(static_cast<int>((0.5 - latitude / kPi) * height), 0, height - 1);
    const double value = cells[static_cast<size_t>(line) * static_cast<size_t>(width) + static_cast<size_t>(column)] / 127.0;
    if (magnitude) {
        return min_nt + 0.5 * (value + 1.0) * (scale_nt - min_nt);
    }
    return (value < 0.0 ? -1.0 : 1.0) * value * value * scale_nt;
}

std::vector<unsigned char> FieldMap::texture_rgba() const
{
    std::vector<unsigned char> pixels(cells.size() * 4);
    for (size_t i = 0; i < cells.size(); ++i) {
        const int value = cells[i];
        pixels[4 * i] = static_cast<unsigned char>(std::lround((value + 127) * 255.0 / 254.0));
        pixels[4 * i + 1] = pixels[4 * i + 2] = 0;
        pixels[4 * i + 3] = 255;
    }
    return pixels;
}

double field_line_visibility(double rotation_deg_per_day, double warp)
{
    const double turns_per_second = std::abs(rotation_deg_per_day * warp) / (360.0 * kSecondsPerDay);
    return 1.0 - glm::smoothstep(0.25, 1.0, turns_per_second);
}

double field_line_size_fade(double pixels_per_radius)
{
    return glm::smoothstep(1.5, 6.0, pixels_per_radius);
}

} // namespace astraxis
