#pragma once

#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace astraxis {

// A planet's internal magnetic field: a spherical-harmonic potential with
// Schmidt quasi-normalized Gauss coefficients in nT. Positions are in units of
// the model's reference radius, in the model's own right-handed frame (z along
// its pole, x through its prime meridian). Nothing outside the planet is
// included: no magnetopause, magnetodisc or tail, so it holds near the planet.
class MagneticField {
public:
    static constexpr int kMaxDegree = 30;

    // Lines "g n m value" or "h n m value"; # starts a comment.
    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    void set_coefficient(bool is_h, int degree, int order, double value);

    int degree() const { return m_degree; }
    double g(int degree, int order) const { return m_g[index(degree, order)]; }
    double h(int degree, int order) const { return m_h[index(degree, order)]; }

    glm::dvec3 field_nt(const glm::dvec3& position) const;

    // A uniform field the body stands in (nT, in the model's frame), added to its own:
    // Jupiter's around Ganymede. Lines near the poles then leave for good, and the closed
    // ones stop where the two fields cancel.
    void set_ambient_nt(const glm::dvec3& field) { m_ambient = field; }
    const glm::dvec3& ambient_nt() const { return m_ambient; }

    // The dipole's strength at the equator of the reference sphere (nT), the
    // direction where its field leaves the planet, and its tilt from the model's
    // z axis (the acute angle). The offset is the center of the eccentric dipole
    // that best includes the quadrupole (Schmidt 1934), in radii.
    double dipole_nt() const;
    glm::dvec3 dipole_axis() const;
    double dipole_tilt_deg() const;
    glm::dvec3 dipole_offset() const;

private:
    static int index(int degree, int order) { return degree * (degree + 1) / 2 + order; }

    int m_degree = 0;
    std::vector<double> m_g, m_h;
    glm::dvec3 m_ambient{0.0};
};

// A field model given at epochs (the Earth's IGRF: every 5 years): linear between them,
// then at its rates until `last_year`, held beyond either end. The file has a line
// "epochs y1 y2 ..", a line "rate_until year", and lines "g|h n m v1 v2 .. rate" (nT, and
// nT per year after the last epoch); # starts a comment.
class MagneticFieldSeries {
public:
    bool load(const std::filesystem::path& file, std::string* error = nullptr);

    double first_year() const { return m_epochs.front(); }
    double last_epoch() const { return m_epochs.back(); }
    double last_year() const { return m_rate_until; }
    MagneticField at(double year) const;

private:
    struct Row {
        bool is_h = false;
        int degree = 0;
        int order = 0;
        std::vector<double> values; // one per epoch
        double rate = 0.0;
    };
    std::vector<double> m_epochs;
    double m_rate_until = 0.0;
    std::vector<Row> m_rows;
};

// Where to start lines and where they end, in radii of the field's model.
struct FieldLineSettings {
    std::vector<double> l_shells{2.0, 3.5, 6.0}; // rings around the eccentric dipole, in its equator
    int lines_per_shell = 8;
    double max_radius = 12.0; // lines stop here (the model does not hold beyond)
    double surface_equatorial = 1.0; // the spheroid they end on, z along the model's pole
    double surface_polar = 1.0;
    double tolerance = 0.003; // of the polyline
    // Lines that leave: started on the surface at these magnetic latitudes (degrees from the
    // dipole's equator, in both hemispheres), lines_per_shell around each. For a field with an
    // ambient part; of the rings above only the lines that close on the body are then kept.
    std::vector<double> open_latitudes_deg;
};

// One line through `seed`, from one end to the other: xyz in radii, w = fade
// (0 near the planet .. 1 at max_radius, for the renderer).
std::vector<glm::dvec4> trace_field_line(const MagneticField& field, const glm::dvec3& seed,
                                         const FieldLineSettings& settings);
std::vector<std::vector<glm::dvec4>> trace_field_lines(const MagneticField& field, const FieldLineSettings& settings);

// Geometric approximation to Ganymede's interaction with the corotating plasma. It maps
// dipole-plus-background traces to a tilted open-closed boundary and swept Alfvén wings.
class GanymedeInteraction {
public:
    // Nonzero ambient field (nT) and unit plasma flow direction, in the magnetic model's axes.
    GanymedeInteraction(const glm::dvec3& ambient_nt, const glm::dvec3& downstream);
    glm::dvec3 position(const glm::dvec3& point) const;

private:
    glm::dvec3 m_downstream;
    glm::dvec3 m_field_axis;
    glm::dvec3 m_perpendicular_mach;
    double m_parallel_mach;
};

// A body's field and its lines ([bodies.magnetic_field] in a scene file).
struct BodyMagneticField {
    std::shared_ptr<const MagneticField> model; // null: none
    std::vector<std::vector<glm::dvec4>> lines; // in the model's frame and radii
    // A model that changes over the years: `model` and `lines` are then those of
    // `model_year`, which Scene::update keeps at the whole year nearest the clock (within
    // the series; its last epoch after loading).
    std::shared_ptr<const MagneticFieldSeries> series;
    double model_year = 0.0;
    FieldLineSettings line_settings; // what `lines` were traced with
    // The body stands in its parent's field (Ganymede in Jupiter's, from the parent's table of
    // the field at its moons): Scene::update keeps the model's ambient field and the lines to
    // it, tracing again when it has changed by kAmbientStep of its strength. ambient_nt is
    // what the lines were traced with (model frame); zero until the first update.
    static constexpr double kAmbientStep = 0.03;
    bool ambient_from_parent = false;
    glm::dvec3 ambient_nt{0.0};
    // Set by whoever draws the lines; Scene::update follows the parent's field only for a
    // body whose lines were drawn since the last update (a frame late when they first show),
    // so that a scene with such a moon costs nothing while its lines are off.
    mutable bool lines_in_view = false;
    double reference_radius_km = 1.0;
    glm::vec3 color{1.0f}; // sRGB

    // The model's frame in the body's: first its pole turned over if it is the
    // body's south (Uranus' ULS: north is the IAU's south), then about the body's
    // z by longitude_offset_deg, or else from the body's prime meridian to the
    // model's own, W = pm_w0_deg + pm_rate_deg_per_day * d (Neptune's NLS).
    bool flip_pole = false;
    double longitude_offset_deg = 0.0;
    bool own_prime_meridian = false;
    double pm_w0_deg = 0.0;
    double pm_rate_deg_per_day = 0.0;
};

// A map of the field fixed to a body, in its texture's projection: rows from north
// to south by planetocentric latitude, columns from 180 W eastward. Baked by
// tools/geomag/field_map_bake.py. A radial map holds the square root of the strength,
// value = sign(Br) sqrt(|Br| / scale), saturating at scale; a strength map (`magnitude`)
// is linear, value = 2 (|B| - min) / (max - min) - 1. Either way in -1 .. 1.
struct FieldMap {
    int width = 0;
    int height = 0;
    bool magnitude = false;
    double scale_nt = 0.0; // radial: where it saturates, usually the largest |Br|; strength: the largest |B|
    double min_nt = 0.0;   // strength: the smallest |B|
    std::vector<signed char> cells; // 127 * value

    bool load(const std::filesystem::path& file, std::string* error = nullptr);
    // Br or |B| at a cell (nearest), by east longitude and planetocentric latitude in radians.
    double value_nt(double east_longitude, double latitude) const;
    // Red channel: 0 .. 255 for value -1 .. 1, so that it interpolates through zero.
    std::vector<unsigned char> texture_rgba() const;
};

// A strength map of a model's field on a spheroid (radii of the model, z along its pole), in
// the same cells and range as tools/geomag/field_map_bake.py's: |B| where the planetocentric
// direction of each cell's middle meets the surface.
FieldMap field_strength_map(const MagneticField& field, int width, int height, double surface_equatorial,
                            double surface_polar, double min_nt, double max_nt);

// Which map a view shows.
enum class FieldMapView { None, Surface, Deep, Crust };

// Maps of a body's field ([bodies.magnetic_map] in a scene file): at its surface,
// optionally deeper down, where the dynamo is, and of its crust alone, far weaker than
// the rest. Null: none.
struct BodyMagneticMap {
    std::shared_ptr<const FieldMap> surface;
    std::shared_ptr<const FieldMap> deep;
    std::shared_ptr<const FieldMap> crust;
};

// Lines on screen turn with their planet; too fast to follow, they fade: all
// shown up to 0.25 turns per real second, gone at 1.
double field_line_visibility(double rotation_deg_per_day, double warp);

// Fade with the apparent size of the planet: a tangle of lines around a dot is
// no use. Gone below 1.5 pixels of radius, all shown from 6.
double field_line_size_fade(double pixels_per_radius);

} // namespace astraxis
