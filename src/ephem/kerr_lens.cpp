#include "ephem/kerr_lens.hpp"

#include "ephem/kerr_null.hpp"

#include <glm/geometric.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/matrix.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace astraxis {

namespace {

constexpr int kMaxSteps = 6000;
constexpr double kNearStep = 1.0; // step at most this fraction of the distance to the target
constexpr double kDifferenceStep = 1e-7; // radians, for the finite-difference Jacobian
constexpr double kMaxAngleStep = 0.2;   // radians per Gauss-Newton step
constexpr int kColdIterations = 30;
constexpr int kWarmIterations = 8;
constexpr int kScanSamples = 48;
constexpr double kSameImage = 1e-5; // radians: two solves that found the same image

// Orthonormal tangents of a unit vector d, with e1 x e2 = d.
void tangents(const glm::dvec3& d, glm::dvec3& e1, glm::dvec3& e2)
{
    const glm::dvec3 helper = std::abs(d.z) < 0.9 ? glm::dvec3(0.0, 0.0, 1.0) : glm::dvec3(1.0, 0.0, 0.0);
    e1 = glm::normalize(glm::cross(helper, d));
    e2 = glm::cross(d, e1);
}

struct Probe {
    glm::dvec3 miss{std::numeric_limits<double>::max()};
    glm::dvec3 emission{0.0};
    bool disk_first = false;
    double size() const { return glm::length(miss); }
};

class Solver {
public:
    Solver(const glm::dvec3& camera, const glm::dvec3& target, const LensSettings& settings)
        : m_cam(camera), m_target(target), m_s(settings)
    {
        // Any static frame at the camera: d -> n does not depend on its axes.
        const glm::dvec3 forward = -glm::normalize(camera);
        glm::dvec3 up;
        glm::dvec3 right;
        tangents(forward, right, up);
        m_frame = static_observer_frame(camera, settings.spin, forward, up, right);
        m_to_local = glm::inverse(glm::dmat3(m_frame.right_dir, m_frame.up_dir, m_frame.forward_dir));
        m_r_plus = horizon_radius(settings.spin);
        m_r_stop = 1.2 * std::max(glm::length(camera), glm::length(target)) + 10.0;
    }

    // The backward ray along coordinate direction d: its closest approach to the target.
    Probe probe(const glm::dvec3& d) const
    {
        const double a = m_s.spin;
        Photon ph = frame_photon(m_cam, m_frame, m_to_local * d);
        Probe out;
        double best = std::numeric_limits<double>::max();
        int crossings = 0;
        for (int i = 0; i < kMaxSteps; ++i) {
            const glm::dvec3 prev = ph.x;
            const double near = glm::length(prev - m_target);
            rk4_step(ph, a, std::min(ray_step(prev, a), std::max(kNearStep * near, 1e-9)));
            const glm::dvec3 seg = ph.x - prev;
            const double seg2 = glm::dot(seg, seg);

            double crossing = 2.0; // where along the segment it crosses the disk (> 1: it does not)
            if (m_s.disk_r_out > 0.0 && prev.z * ph.x.z < 0.0) {
                const double t = prev.z / (prev.z - ph.x.z);
                const glm::dvec3 c = prev + seg * t;
                const double rc = kerr_schild_r(glm::dvec3(c.x, c.y, 0.0), a);
                if (rc >= m_s.disk_r_in && rc <= m_s.disk_r_out) {
                    crossing = t;
                }
            }
            const double t = seg2 > 0.0 ? std::clamp(glm::dot(m_target - prev, seg) / seg2, 0.0, 1.0) : 0.0;
            const glm::dvec3 miss = prev + seg * t - m_target;
            const double miss2 = glm::dot(miss, miss);
            if (miss2 < best && seg2 > 0.0) {
                best = miss2;
                out.miss = miss;
                out.emission = seg / std::sqrt(seg2);
                out.disk_first = crossings > 0 || crossing < t;
            }
            if (crossing <= 1.0) {
                ++crossings;
            }
            if (kerr_schild_r(ph.x, a) < m_r_plus + 0.02) {
                break;
            }
            if (glm::dot(ph.x, ph.p) > 0.0 &&
                (glm::dot(ph.x, ph.x) > m_r_stop * m_r_stop || glm::dot(ph.x - m_target, seg) > 0.0)) {
                break; // past the hole (it cannot turn back) and past the target, or gone
            }
        }
        return out;
    }

    // Finite-difference Jacobian at d (miss p there), as a map of tangent vectors.
    glm::dmat3 jacobian(const glm::dvec3& d, const Probe& p) const
    {
        glm::dvec3 e1;
        glm::dvec3 e2;
        tangents(d, e1, e2);
        const glm::dvec3 c1 = (probe(glm::normalize(d + kDifferenceStep * e1)).miss - p.miss) / kDifferenceStep;
        const glm::dvec3 c2 = (probe(glm::normalize(d + kDifferenceStep * e2)).miss - p.miss) / kDifferenceStep;
        return glm::outerProduct(c1, e1) + glm::outerProduct(c2, e2);
    }

    // Gauss-Newton from d. On success fills `out`.
    bool refine(glm::dvec3 d, glm::dmat3 jac, bool have_jacobian, int age, int iterations, LensImage& out) const
    {
        Probe p = probe(d);
        bool fresh = false;
        int updates = 0; // Broyden updates since the Jacobian was measured
        if (!have_jacobian || age >= m_s.jacobian_max_age) {
            jac = jacobian(d, p);
            fresh = true;
            age = 0;
        }
        for (int iter = 0; iter < iterations; ++iter) {
            glm::dvec3 e1;
            glm::dvec3 e2;
            tangents(d, e1, e2);
            const glm::dvec3 c1 = jac * e1;
            const glm::dvec3 c2 = jac * e2;
            const double g11 = glm::dot(c1, c1);
            const double g12 = glm::dot(c1, c2);
            const double g22 = glm::dot(c2, c2);
            const double det = g11 * g22 - g12 * g12;
            if (!(det > 1e-30 * (g11 * g22))) {
                if (fresh) {
                    return false;
                }
                jac = jacobian(d, p);
                fresh = true;
                age = 0;
                updates = 0;
                continue;
            }
            const double b1 = glm::dot(c1, p.miss);
            const double b2 = glm::dot(c2, p.miss);
            double s1 = -(g22 * b1 - g12 * b2) / det;
            double s2 = -(g11 * b2 - g12 * b1) / det;
            const double len = std::hypot(s1, s2);
            if (len < m_s.tolerance_rad) {
                // Broyden's estimate is good enough to steer by, but the image's shape needs
                // the Jacobian itself: measured again after several updates (a cold start)
                // or every few frames.
                if (!fresh && (updates > 1 || age >= m_s.jacobian_max_age)) {
                    jac = jacobian(d, p);
                    age = 0;
                }
                fill(d, jac, age, p, out);
                return true;
            }
            if (len > kMaxAngleStep) {
                s1 *= kMaxAngleStep / len;
                s2 *= kMaxAngleStep / len;
            }
            // Backtracking line search on |miss|.
            bool improved = false;
            glm::dvec3 d1(0.0);
            Probe p1;
            double scale = 1.0;
            for (int k = 0; k < 10 && !improved; ++k, scale *= 0.5) {
                d1 = glm::normalize(d + scale * (s1 * e1 + s2 * e2));
                p1 = probe(d1);
                improved = p1.size() < p.size();
            }
            if (!improved) {
                if (fresh) {
                    return false;
                }
                jac = jacobian(d, p);
                fresh = true;
                age = 0;
                updates = 0;
                continue;
            }
            // Broyden's update along the step taken.
            const glm::dvec3 step = d1 - d;
            jac += glm::outerProduct(p1.miss - p.miss - jac * step, step) / glm::dot(step, step);
            d = d1;
            p = p1;
            fresh = false;
            ++updates;
        }
        return false;
    }

    // Searches one side of the hole (+1: the target's, -1: the other) along the
    // plane through the camera, the hole and the target, then refines the best.
    bool search(double side, LensImage& out) const
    {
        const double dist = glm::length(m_cam);
        const glm::dvec3 u = -m_cam / dist;
        const glm::dvec3 q = m_target - m_cam;
        glm::dvec3 w = q - glm::dot(q, u) * u;
        if (glm::length(w) < 1e-12 * glm::length(q)) {
            glm::dvec3 other;
            tangents(u, w, other);
        }
        w = glm::normalize(w);
        const double target_angle = std::atan2(glm::length(glm::cross(u, q)), glm::dot(u, q));
        const double lo = dist > 2.5 ? std::asin(2.0 / dist) : 0.05;
        const double hi = side > 0.0 ? std::min(std::max(2.0 * target_angle, 4.0 * lo) + 0.05, 3.0)
                                     : std::min(std::max(dist > 40.0 ? std::asin(40.0 / dist) : 1.6, 4.0 * lo), 3.0);
        glm::dvec3 best_d(0.0);
        double best = std::numeric_limits<double>::max();
        for (int i = 0; i < kScanSamples; ++i) {
            const double theta = lo * std::pow(hi / lo, (i + 0.5) / kScanSamples);
            const glm::dvec3 d = std::cos(theta) * u + side * std::sin(theta) * w;
            const double miss = probe(d).size();
            if (miss < best) {
                best = miss;
                best_d = d;
            }
        }
        return best < std::numeric_limits<double>::max() &&
               refine(best_d, glm::dmat3(0.0), false, 0, kColdIterations, out);
    }

    // Which side of the hole a direction is on (> 0: the target's).
    double side_of(const glm::dvec3& d) const
    {
        const glm::dvec3 u = -glm::normalize(m_cam);
        const glm::dvec3 q = m_target - m_cam;
        return glm::dot(d - glm::dot(d, u) * u, q - glm::dot(q, u) * u);
    }

    glm::dvec3 straight() const { return glm::normalize(m_target - m_cam); }

private:
    void fill(const glm::dvec3& d, const glm::dmat3& jac, int age, const Probe& p, LensImage& out) const
    {
        glm::dvec3 e1;
        glm::dvec3 e2;
        tangents(d, e1, e2);
        const glm::dvec3 c1 = jac * e1;
        const glm::dvec3 c2 = jac * e2;
        // Singular values of [c1 c2]: eigenvalues of its 2 x 2 Gram matrix.
        const double g11 = glm::dot(c1, c1);
        const double g12 = glm::dot(c1, c2);
        const double g22 = glm::dot(c2, c2);
        const double mean = 0.5 * (g11 + g22);
        const double spread = std::sqrt(std::max(0.25 * (g11 - g22) * (g11 - g22) + g12 * g12, 0.0));
        const double l1 = mean + spread;
        const double l2 = std::max(mean - spread, 1e-300);
        // Eigenvector of l1 (the image's short axis: the miss grows fastest along it).
        glm::dvec2 v1 = std::abs(g12) > 1e-300 * mean ? glm::dvec2(g12, l1 - g11)
                        : g11 >= g22                    ? glm::dvec2(1.0, 0.0)
                                                        : glm::dvec2(0.0, 1.0);
        if (glm::length(v1) == 0.0) {
            v1 = glm::dvec2(1.0, 0.0);
        }
        v1 = glm::normalize(v1);
        const glm::dvec2 v2(-v1.y, v1.x);

        out.direction = d;
        out.jacobian = jac;
        out.jacobian_age = age + 1;
        out.axis[0] = v1.x * e1 + v1.y * e2;
        out.axis[1] = v2.x * e1 + v2.y * e2;
        out.angle_per_m[0] = 1.0 / std::sqrt(l1);
        out.angle_per_m[1] = 1.0 / std::sqrt(l2);
        const double dist = glm::length(m_target - m_cam);
        const double parity = glm::dot(glm::cross(c1, c2), p.emission) >= 0.0 ? 1.0 : -1.0;
        out.magnification = parity * dist * dist / std::sqrt(l1 * l2);
        out.emission = p.emission;
        out.disk_first = p.disk_first;
        out.miss = p.size();
    }

    glm::dvec3 m_cam;
    glm::dvec3 m_target;
    LensSettings m_s;
    ObserverFrame m_frame;
    glm::dmat3 m_to_local{1.0};
    double m_r_plus = 2.0;
    double m_r_stop = 0.0;
};

} // namespace

int find_lens_images(const glm::dvec3& camera, const glm::dvec3& target, const LensSettings& settings,
                     const LensImage* guesses, int guesses_found, LensImage out[2])
{
    const Solver solver(camera, target, settings);
    int found = 0;
    auto warm = [&](int k) {
        return guesses && (guesses_found & (1 << k)) != 0 &&
               solver.refine(guesses[k].direction, guesses[k].jacobian, true, guesses[k].jacobian_age,
                             kWarmIterations, out[k]);
    };
    if (warm(0) || solver.refine(solver.straight(), glm::dmat3(0.0), false, 0, kColdIterations, out[0]) ||
        solver.search(1.0, out[0])) {
        found |= 1;
    }
    if ((found & 1) && out[0].magnification > 0.0 &&
        out[0].magnification - 1.0 < settings.min_secondary_magnification) {
        return found; // the secondary would be too faint to matter
    }
    if (warm(1)) {
        found |= 2;
    }
    auto same = [&]() {
        return (found & 3) == 3 &&
               glm::length(glm::cross(out[0].direction, out[1].direction)) < kSameImage &&
               glm::dot(out[0].direction, out[1].direction) > 0.0;
    };
    if (same()) {
        found &= 1; // both warm starts slid onto one image (they swap as the target passes behind the hole)
    }
    if (!(found & 2) && settings.search_secondary) {
        if (solver.search(-1.0, out[1])) {
            found |= 2;
        }
        if (same()) {
            found &= 1;
        }
    }
    // The primary is the one on the target's side.
    if ((found & 3) == 3 && solver.side_of(out[1].direction) > solver.side_of(out[0].direction)) {
        std::swap(out[0], out[1]);
    } else if (found == 2 && solver.side_of(out[1].direction) > 0.0) {
        out[0] = out[1];
        found = 1;
    }
    return found;
}

} // namespace astraxis
