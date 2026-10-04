// Plumes by single scattering along the view ray through their bounding box.
// Everything is in the plume's local frame (km): origin at the vent on the
// surface, z up (the local vertical), the planet's center below at u_planet.
// Altitudes above the curved surface use the paraboloid approximation
// alt = z + (x^2 + y^2) / (2 R). Each shape has an analytic density (the
// extinction coefficient, per km):
//   umbrella:      particles on ballistic paths, crowding under their envelope
//                  (the canopy) alt = H (1 - (d / R_c)^2), in filaments;
//   jets:          cones from vents along a row, density ~ exp(-z / H) and
//                  diluted as their width grows, w = z0 + z tan(spread);
//   geyser column: a cylinder up to height H;
//   geyser tail:   along +x at altitude H, widening and fading downwind.
// The plumes are optically thin: sunlight reaching them is cut only by the
// planet's shadow and eclipses by other bodies. Radiance is in I/F units
// (Henyey-Greenstein phase function with a mean of 1, hence the 1/4); alpha is
// the view transmission, so dark plumes (Triton's) dim what lies behind. Slow
// waves run outward along the flow (real time, artistic).

#include "atmosphere.hlsli"

#define MAX_JETS 12
#define MAX_OCCLUDERS 4
#define TYPE_UMBRELLA 0
#define TYPE_JETS 1
#define TYPE_COLUMN 2
#define TYPE_TAIL 3

static const int kSteps = 96;
static const float kTwoPi = 6.2831853;

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_local;  // camera-relative world vector -> local frame (3x3 rotation)
    float4 u_camera;      // xyz = camera (local), w = 1 to start the rays there
    float4 u_sun;         // xyz = unit direction to the sun (local), w = its angular radius
    float4 u_box_min;     // xyz = box corner (local), w = type
    float4 u_box_max;     // xyz = box corner, w = jet count
    float4 u_planet;      // xyz = planet center (local), w = radius (km)
    float4 u_shape;       // umbrella: H, R_c, shell width (in q); jets: H, tan(spread), z0; column: H, radius;
                          // tail: H, radius, length, widening per km
    float4 u_albedo;      // rgb = single-scattering albedo, w = g
    float4 u_params;      // x = time (s), y = occluder count, z = density scale (per km)
    float4 u_occluders[MAX_OCCLUDERS]; // xyz = center (local), w = radius
    float4 u_jet_pos[MAX_JETS];        // vents (local)
    float4 u_jet_axis[MAX_JETS];       // their axes (unit)
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 local    : TEXCOORD1;
};

float gradient_noise(float2 p)
{
    return frac(52.9829189 * frac(dot(p, float2(0.06711056, 0.00583715))));
}

// A wave of period `length` (km) moving outward along the path coordinate s.
float puff(float s, float length)
{
    return 1.0 + 0.3 * sin(kTwoPi * (s / length - u_params.x / 8.0));
}

float altitude(float3 p)
{
    return p.z + dot(p.xy, p.xy) / (2.0 * u_planet.w);
}

// Particles launched from the vent at one speed (reaching H straight up) on
// ballistic paths, horizontal distances scaled so the reach is R_c. Through a
// point at (d, alt) pass two paths, the steeper with launch angle
// tan(theta) = (2 H / d) (1 + q), q = sqrt(1 - alt / H - (d / 2H)^2): q = 0 on
// the canopy (the paths' envelope), where the particles crowd. Filaments are
// ranges of launch angle and azimuth.
float umbrella_density(float3 p)
{
    float h = u_shape.x;
    float alt = altitude(p);
    float d = length(p.xy) * (2.0 * h / u_shape.y);
    float x = d / (2.0 * h);
    float disc = 1.0 - alt / h - x * x;
    if (alt < 0.0 || disc < -0.1) {
        return 0.0;
    }
    float q = sqrt(max(disc, 0.0));
    float shell = disc >= 0.0 ? exp(-q / u_shape.z) : exp(-(disc * disc) / 0.0004);
    float theta = atan2(2.0 * h * (1.0 + q), d);
    float az = atan2(p.y, p.x);
    float filaments = sin(13.0 * theta + 2.0 * sin(3.0 * az) + 1.5 * sin(7.0 * az + 1.3)) * sin(9.0 * az + 0.7);
    filaments = 1.0 + 0.5 * filaments * smoothstep(0.0, 0.2, x); // azimuth is undefined on the axis
    return u_params.z * shell * filaments * puff(d + h * q, 0.3 * h);
}

float jets_density(float3 p)
{
    float sum = 0.0;
    int count = (int)u_box_max.w;
    for (int k = 0; k < count; ++k) {
        float3 q = p - u_jet_pos[k].xyz;
        float z = dot(q, u_jet_axis[k].xyz);
        if (z <= 0.0) {
            continue;
        }
        float rho2 = max(dot(q, q) - z * z, 0.0);
        float w = u_shape.z + z * u_shape.y;
        float dilution = u_shape.z / w;
        sum += exp(-z / u_shape.x - rho2 / (w * w)) * dilution * dilution * puff(z, 0.5 * u_shape.x);
    }
    return u_params.z * sum;
}

float column_density(float3 p)
{
    float top = 1.0 - smoothstep(u_shape.x - 0.5 * u_shape.y, u_shape.x + 0.5 * u_shape.y, p.z);
    return p.z < 0.0 ? 0.0 : u_params.z * exp(-dot(p.xy, p.xy) / (u_shape.y * u_shape.y)) * top;
}

float tail_density(float3 p)
{
    float s = p.x;
    if (s < 0.0) {
        return 0.0;
    }
    float w = u_shape.y + s * u_shape.w;
    float dz = altitude(p) - u_shape.x;
    float r2 = p.y * p.y + dz * dz;
    return u_params.z * (u_shape.y / w) * exp(-r2 / (w * w) - s / u_shape.z) * puff(s, 0.25 * u_shape.z);
}

float density(float3 p)
{
    int type = (int)u_box_min.w;
    if (type == TYPE_UMBRELLA) {
        return umbrella_density(p);
    }
    if (type == TYPE_JETS) {
        return jets_density(p);
    }
    if (type == TYPE_COLUMN) {
        return column_density(p);
    }
    return tail_density(p);
}

// Fraction of the sun's disk seen past a sphere, from p.
float sphere_lit(float3 p, float3 center, float radius)
{
    float3 oc = center - p;
    float along = dot(oc, u_sun.xyz);
    if (along <= 0.0) {
        return 1.0;
    }
    float perp = length(oc - along * u_sun.xyz);
    float penumbra = max(along * u_sun.w, 1e-3 * radius);
    return smoothstep(radius - penumbra, radius + penumbra, perp);
}

float sunlight(float3 p)
{
    float lit = sphere_lit(p, u_planet.xyz, u_planet.w);
    int count = (int)u_params.y;
    for (int k = 0; k < count; ++k) {
        lit *= sphere_lit(p, u_occluders[k].xyz, u_occluders[k].w);
    }
    return lit;
}

// Clips [t0, t1] to the umbrella's canopy, a paraboloid in the local frame:
// z <= top - k (x^2 + y^2), with k = H / R_c^2 + 1 / (2 R) (curvature
// included), a little above the canopy for its soft edge. False if missed.
bool clip_to_canopy(float3 o, float3 dir, inout float t0, inout float t1)
{
    float k = u_shape.x / (u_shape.y * u_shape.y) + 0.5 / u_planet.w;
    float top = 1.05 * u_shape.x;
    float a = k * dot(dir.xy, dir.xy);
    float b = dir.z + 2.0 * k * dot(o.xy, dir.xy);
    float c = o.z - top + k * dot(o.xy, o.xy);
    if (a < 1e-12) {
        if (abs(b) < 1e-12) {
            return c <= 0.0;
        }
        float t = -c / b;
        if (b > 0.0) {
            t1 = min(t1, t);
        } else {
            t0 = max(t0, t);
        }
        return t1 > t0;
    }
    float disc = b * b - 4.0 * a * c;
    if (disc <= 0.0) {
        return false;
    }
    float s = sqrt(disc);
    t0 = max(t0, (-b - s) / (2.0 * a));
    t1 = min(t1, (-b + s) / (2.0 * a));
    return t1 > t0;
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_local, input.rel_pos));
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.local;

    // Slab intersection with the box, then the ground (the planet's sphere).
    float3 inv = 1.0 / select(abs(dir) > 1e-8, dir, 1e-8);
    float3 ta = (u_box_min.xyz - o) * inv;
    float3 tb = (u_box_max.xyz - o) * inv;
    float3 tn = min(ta, tb);
    float3 tf = max(ta, tb);
    float t0 = max(max(max(tn.x, tn.y), tn.z), 0.0);
    float t1 = min(min(tf.x, tf.y), tf.z);
    float3 oc = o - u_planet.xyz;
    float b = dot(oc, dir);
    float3 perp = oc - b * dir;
    float p2 = dot(perp, perp);
    float r2 = u_planet.w * u_planet.w;
    if (p2 < r2) {
        float ground = -b - sqrt(r2 - p2);
        if (ground > 0.0) {
            t1 = min(t1, ground);
        }
    }
    if (t1 <= t0 || ((int)u_box_min.w == TYPE_UMBRELLA && !clip_to_canopy(o, dir, t0, t1))) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }

    float phase = henyey_greenstein_phase(dot(dir, u_sun.xyz), u_albedo.w);
    float dt = (t1 - t0) / kSteps;
    float jitter = gradient_noise(input.position.xy);
    float tau = 0.0;
    float sum = 0.0;
    for (int i = 0; i < kSteps; ++i) {
        float3 p = o + (t0 + (i + jitter) * dt) * dir;
        float sigma = density(p);
        if (sigma > 0.0) {
            sum += sigma * exp(-tau - 0.5 * sigma * dt) * sunlight(p) * dt;
            tau += sigma * dt;
        }
    }
    return float4(0.25 * phase * u_albedo.rgb * sum, exp(-tau));
}
