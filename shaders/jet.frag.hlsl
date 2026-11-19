// A black hole's jet pair in visible light: optically thin synchrotron emission
// integrated along straight view rays (no lensing), added to the scene. The
// formulas are those of src/ephem/jet.cpp (jet_emissivity()), in the jet frame:
// z along the mean axis (+z the approaching jet), x toward the nozzle's tilt, units of M.
//
// Each ray is clipped to quadric envelopes s^2 < A z + B z^2 + C that bound the jet
// (one piece through the tilted nozzle, one beyond it, for each side), and sampled
// uniformly in u = asinh((t - t_c) / b), t_c and b the ray's closest approach to the
// hole: dense near the hole, where the jet is narrow and bright, sparse far away.
// Each span gets the same number of samples.

#include "common.hlsli"

#define MAX_KNOTS 12

cbuffer Uniforms : register(b0, space3)
{
    float4 u_camera;    // xyz = camera, w = shadow radius
    float4 u_right;     // right axis * tan(fov/2) * aspect
    float4 u_up;        // up axis * tan(fov/2)
    float4 u_forward;   // xyz = forward, w = smallest drawn radius per unit distance
    float4 u_shape;     // parabola k, cone slope^2, base, length
    float4 u_speed;     // height of Gamma = 1, Lorentz index, smallest Gamma beta, end of acceleration
    float4 u_speed2;    // deceleration index, Gamma beta at the end of acceleration, spectral index, profile area
    float4 u_tilt;      // x = precession half-angle (rad), y = height where the tilted envelope piece ends
    float4 u_tilt_log;  // ln of the tilt heights
    float4 u_envelope;  // A, B through the tilt; A, B beyond it
    float4 u_envelope2; // x = C
    float4 u_color;     // rgb = color x display gain, w = exposure (inside the compression)
    float4 u_inner;     // x = first knot start, y = jet radius there, z = knot count
    float4 u_knots[MAX_KNOTS]; // start, end, line brightness
};

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

static const int kSamplesPerSpan = 32;
static const float kRidgePsi = 0.82;
static const float kRidgeSigma = 0.15;
static const float kSpineShare = 0.3;
static const float kSpineSigma = 0.2;
static const float kProfileEnd = 1.4;

float jet_radius(float zeta)
{
    return sqrt(u_shape.x * zeta + u_shape.y * zeta * zeta);
}

float gamma_beta(float zeta)
{
    if (zeta <= u_speed.w) {
        float gamma = pow(max(zeta, 1e-6) / u_speed.x, u_speed.y);
        return max(sqrt(max(gamma * gamma - 1.0, 0.0)), u_speed.z);
    }
    return u_speed2.y * pow(zeta / u_speed.w, u_speed2.x);
}

float tilt_fraction(float zeta)
{
    float l = log(max(zeta, 1e-6));
    return smoothstep(u_tilt_log.x, u_tilt_log.y, l) * (1.0 - smoothstep(u_tilt_log.z, u_tilt_log.w, l));
}

float soft_box(float a, float b, float x)
{
    float h = 0.3 * (b - a);
    return smoothstep(a - h, a + h, x) * (1.0 - smoothstep(b - h, b + h, x));
}

float line_brightness(float zeta)
{
    float alpha = u_speed2.z;
    int count = (int)u_inner.z;
    if (count == 0) {
        return pow(jet_radius(zeta), -(1.0 + alpha));
    }
    float z0 = u_inner.x;
    float h = 0.3 * (u_knots[0].y - z0);
    float sum = pow(jet_radius(zeta) / u_inner.y, -(1.0 + alpha)) * (1.0 - smoothstep(z0 - h, z0 + h, zeta)) * u_knots[0].z;
    if (zeta > z0 - h) {
        for (int i = 0; i < count; ++i) {
            sum += u_knots[i].z * soft_box(u_knots[i].x, u_knots[i].y, zeta);
        }
    }
    return sum;
}

float profile(float psi)
{
    float r = (psi - kRidgePsi) / kRidgeSigma;
    float c = psi / kSpineSigma;
    return psi > kProfileEnd ? 0.0 : exp(-0.5 * r * r) + kSpineShare * exp(-0.5 * c * c);
}

// jet_emissivity() of src/ephem/jet.cpp.
float emissivity(float3 x, float3 to_camera, float min_width)
{
    float zeta = abs(x.z);
    if (zeta < u_shape.z || zeta > u_shape.w) {
        return 0.0;
    }
    float side = x.z >= 0.0 ? 1.0 : -1.0;
    float lean = tan(u_tilt.x * tilt_fraction(zeta));
    float2 lateral = float2(x.x - x.z * lean, x.y);
    float s = length(lateral);
    float w = jet_radius(zeta);
    float width2 = w * w + min_width * min_width;
    float psi = s * rsqrt(width2);
    if (psi > kProfileEnd) {
        return 0.0;
    }
    float fade = smoothstep(u_shape.z, 4.0 * u_shape.z, zeta) * (1.0 - smoothstep(0.9 * u_shape.w, u_shape.w, zeta));
    float cross_section = profile(psi) / (u_speed2.w * width2);

    float slope = (u_shape.x + 2.0 * u_shape.y * zeta) / (2.0 * w);
    float3 dir = side * float3(lean, 0.0, 1.0);
    if (s > 0.0) {
        dir += float3(lateral / s, 0.0) * (s / w) * slope;
    }
    dir = normalize(dir);
    float gb = gamma_beta(zeta);
    float doppler = 1.0 / (sqrt(1.0 + gb * gb) - gb * dot(dir, to_camera));
    return fade * line_brightness(zeta) * cross_section * pow(doppler, 2.0 + u_speed2.z);
}

// Where Q(t) = a t^2 + b t + c < 0, within [lo, hi]: up to two intervals (x..y, z..w; empty: x >= y).
float4 quadric_inside(float a, float b, float c, float lo, float hi)
{
    float4 r = float4(1.0, 0.0, 1.0, 0.0);
    // Linear where the quadratic term stays negligible over the range.
    float reach = max(abs(lo), abs(min(hi, 1e30)));
    if (abs(a) * reach * reach < 1e-6 * (abs(b) * reach + abs(c)) + 1e-30) {
        // Linear.
        if (abs(b) < 1e-30) {
            return c < 0.0 ? float4(lo, hi, 1.0, 0.0) : r;
        }
        float t0 = -c / b;
        return b > 0.0 ? float4(lo, min(hi, t0), 1.0, 0.0) : float4(max(lo, t0), hi, 1.0, 0.0);
    }
    float disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        return a < 0.0 ? float4(lo, hi, 1.0, 0.0) : r;
    }
    float q = -0.5 * (b + (b >= 0.0 ? 1.0 : -1.0) * sqrt(disc));
    float t1 = q / a;
    float t2 = abs(q) > 0.0 ? c / q : t1;
    float tmin = min(t1, t2);
    float tmax = max(t1, t2);
    if (a > 0.0) {
        return float4(max(lo, tmin), min(hi, tmax), 1.0, 0.0);
    }
    return float4(lo, min(hi, tmin), max(lo, tmax), hi);
}

// The part of the ray inside one envelope piece (side: +1 approaching jet, -1 the other),
// heights z_lo..z_hi along that side.
float4 piece(float3 o, float3 d, float side, float A, float B, float C, float z_lo, float z_hi, float t_max)
{
    float oz = side * o.z;
    float dz = side * d.z;
    // Heights z_lo..z_hi along the ray.
    float lo = 0.0;
    float hi = t_max;
    if (abs(dz) < 1e-12) {
        if (oz < z_lo || oz > z_hi) {
            return float4(1.0, 0.0, 1.0, 0.0);
        }
    } else {
        float ta = (z_lo - oz) / dz;
        float tb = (z_hi - oz) / dz;
        lo = max(lo, min(ta, tb));
        hi = min(hi, max(ta, tb));
    }
    if (lo >= hi) {
        return float4(1.0, 0.0, 1.0, 0.0);
    }
    float a = d.x * d.x + d.y * d.y - B * dz * dz;
    float b = 2.0 * (o.x * d.x + o.y * d.y) - A * dz - 2.0 * B * oz * dz;
    float c = o.x * o.x + o.y * o.y - A * oz - B * oz * oz - C;
    return quadric_inside(a, b, c, lo, hi);
}

float asinh_(float x)
{
    // Odd, without cancellation for negative x.
    float a = abs(x);
    return sign(x) * log(a + sqrt(a * a + 1.0));
}

float4 main(PSInput input) : SV_Target0
{
    float2 ndc = float2(input.uv.x * 2.0 - 1.0, 1.0 - input.uv.y * 2.0);
    float3 o = u_camera.xyz;
    float3 d = normalize(u_forward.xyz + ndc.x * u_right.xyz + ndc.y * u_up.xyz);

    // Closest approach to the hole; rays passing inside the shadow end there.
    float t_c = -dot(o, d);
    float3 p_c = o + t_c * d;
    float b_c = length(p_c);
    float t_max = 3.0e38;
    if (t_c > 0.0 && b_c < u_camera.w) {
        t_max = t_c;
    }

    // Up to two pieces per side, each up to two intervals.
    float2 spans[8];
    int count = 0;
    float z_base = 0.9 * u_shape.z;
    float z_split = u_tilt.y;
    [unroll]
    for (int k = 0; k < 2; ++k) {
        float side = k == 0 ? 1.0 : -1.0;
        float4 near_part = piece(o, d, side, u_envelope.x, u_envelope.y, u_envelope2.x, z_base, max(z_split, z_base), t_max);
        float4 far_part = piece(o, d, side, u_envelope.z, u_envelope.w, u_envelope2.x, max(z_split, z_base), u_shape.w, t_max);
        if (near_part.x < near_part.y) { spans[count++] = near_part.xy; }
        if (near_part.z < near_part.w) { spans[count++] = near_part.zw; }
        if (far_part.x < far_part.y) { spans[count++] = far_part.xy; }
        if (far_part.z < far_part.w) { spans[count++] = far_part.zw; }
    }
    if (count == 0) {
        discard;
    }

    float bb = max(b_c, 0.5);

    // Interleaved gradient noise: a fixed offset per pixel turns banding into fine grain.
    float jitter = frac(52.9829189 * frac(dot(input.position.xy, float2(0.06711056, 0.00583715))));
    float min_width_per_distance = u_forward.w;
    float3 to_camera = -d;
    float sum = 0.0;
    for (int i = 0; i < count; ++i) {
        float u0 = asinh_((spans[i].x - t_c) / bb);
        float du = (asinh_((spans[i].y - t_c) / bb) - u0) / kSamplesPerSpan;
        for (int j = 0; j < kSamplesPerSpan; ++j) {
            float u = u0 + (j + jitter) * du;
            float dt = bb * sinh(u);
            float3 x = p_c + dt * d;
            float distance = t_c + dt;
            float dl = bb * cosh(u) * du;
            sum += emissivity(x, to_camera, min_width_per_distance * distance) * dl;
        }
    }
    // Compressed like the hot flow's false color: the jet's brightness spans decades within
    // one view (the beamed core seen down the jet, the counter-jet).
    return float4(u_color.rgb * sqrt(u_color.w * sum), 0.0);
}
