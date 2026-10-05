// Comet comae and ion tails: optically thin emission integrated along the view
// ray through a bounding box (vertex shader: plume.vert.hlsl). Everything is in
// the comet's local frame (km): the nucleus at the origin, +x along the ion
// tail. Each emission density is normalized to a volume integral of 1, so the
// radiance is the color times the brightness times the column of the
// normalized density (per km^2), and the total flux is brightness / distance^2.
//   coma:     C2 in the Haser model, parent and daughter scale lengths Lp, Ld:
//             n(r) = (exp(-r / Ld) - exp(-r / Lp)) / (4 pi r^2 (Ld - Lp));
//             along the ray, r = rho / cos(theta) with theta the angle from the
//             point of closest approach, dl / r^2 = d(theta) / rho, so even
//             steps in theta sample the 1 / r^2 exactly;
//   ion tail: a Gaussian cross-section of 1/e radius w = w0 + x tan(opening),
//             fading as exp(-x / L) downstream and within w0 / 2 sunward, with
//             (w0 / w)^2 keeping the flow's flux; rays drift outward. Steps
//             crowd around the ray's closest approach to the axis,
//             t = tc + w(tc) sinh(u) with even steps in u: a few per tail width
//             there, growing geometrically away from it (also when the ray runs
//             along the tail from a camera inside it).
// Rays stop at the occluders (the nucleus, nearby bodies): the box is drawn
// without depth testing when the camera is inside it.

#include "common.hlsli"

#define MAX_OCCLUDERS 4
#define TYPE_COMA 0
#define TYPE_ION_TAIL 1

static const int kComaSteps = 32;
static const int kTailSteps = 64;
static const float kPi = 3.14159265;

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_local;  // camera-relative world vector -> local frame (3x3 rotation)
    float4 u_camera;      // xyz = camera (local), w = 1 to start the rays there
    float4 u_box_min;     // xyz = box corner (local), w = type
    float4 u_box_max;     // xyz = box corner, w = occluder count
    float4 u_shape;       // coma: Lp, Ld, nucleus radius; tail: w0, tan(opening), L; w = brightness
    float4 u_color;       // rgb = color (sRGB), w = time (s)
    float4 u_occluders[MAX_OCCLUDERS]; // xyz = center (local), w = radius
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

float hash3(float3 p)
{
    p = frac(p * 0.1031);
    p += dot(p, p.zyx + 31.32);
    return frac((p.x + p.y) * p.z);
}

// Value noise in [0, 1].
float value_noise(float3 p)
{
    float3 i = floor(p);
    float3 f = p - i;
    f = f * f * (3.0 - 2.0 * f);
    float a = lerp(lerp(hash3(i), hash3(i + float3(1, 0, 0)), f.x),
                   lerp(hash3(i + float3(0, 1, 0)), hash3(i + float3(1, 1, 0)), f.x), f.y);
    float b = lerp(lerp(hash3(i + float3(0, 0, 1)), hash3(i + float3(1, 0, 1)), f.x),
                   lerp(hash3(i + float3(0, 1, 1)), hash3(i + float3(1, 1, 1)), f.x), f.y);
    return lerp(a, b, f.z);
}

// HLSL has no asinh; odd-symmetric form, accurate for large |x|.
float arsinh(float x)
{
    float a = abs(x);
    return sign(x) * log(a + sqrt(a * a + 1.0));
}

// Column of the normalized Haser density along o + t dir, t in [t0, t1].
float coma_column(float3 o, float3 dir, float t0, float t1, float jitter)
{
    float lp = u_shape.x;
    float ld = u_shape.y;
    float tc = -dot(o, dir);
    float rho = max(length(o + tc * dir), u_shape.z);
    float a0 = atan((t0 - tc) / rho);
    float a1 = atan((t1 - tc) / rho);
    float da = (a1 - a0) / kComaSteps;
    float sum = 0.0;
    for (int i = 0; i < kComaSteps; ++i) {
        float r = rho / cos(a0 + (i + jitter) * da);
        sum += exp(-r / ld) - exp(-r / lp);
    }
    return sum * da / (4.0 * kPi * rho * (ld - lp));
}

float tail_density(float3 p)
{
    float w0 = u_shape.x;
    float x = p.x;
    float w = w0 + max(x, 0.0) * u_shape.y;
    float rho2 = dot(p.yz, p.yz);
    float along = x >= 0.0 ? exp(-x / u_shape.z) : exp(2.0 * x / w0);
    // Rays: azimuthal structure, stretched along the tail (log distance), drifting outward.
    float s = log(1.0 + max(x, 0.0) / w0);
    float2 around = p.yz / max(sqrt(rho2), 1e-3 * w0);
    float rays = value_noise(float3(3.0 * around, 2.0 * s - 0.05 * u_color.w));
    float density = (w0 / w) * (w0 / w) * exp(-rho2 / (w * w)) * along * (0.55 + 0.9 * rays * rays);
    return density / (kPi * w0 * w0 * u_shape.z);
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_local, input.rel_pos));
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.local;

    // Slab intersection with the box, then the occluding spheres.
    float3 inv = 1.0 / select(abs(dir) > 1e-8, dir, 1e-8);
    float3 ta = (u_box_min.xyz - o) * inv;
    float3 tb = (u_box_max.xyz - o) * inv;
    float3 tn = min(ta, tb);
    float3 tf = max(ta, tb);
    float t0 = max(max(max(tn.x, tn.y), tn.z), 0.0);
    float t1 = min(min(tf.x, tf.y), tf.z);
    int count = (int)u_box_max.w;
    for (int k = 0; k < count; ++k) {
        float3 oc = o - u_occluders[k].xyz;
        float b = dot(oc, dir);
        float c = dot(oc, oc) - u_occluders[k].w * u_occluders[k].w;
        float disc = b * b - c;
        if (disc > 0.0) {
            float hit = -b - sqrt(disc);
            if (hit > 0.0) {
                t1 = min(t1, hit);
            }
        }
    }
    if (t1 <= t0) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }

    float jitter = gradient_noise(input.position.xy);
    float column = 0.0;
    if ((int)u_box_min.w == TYPE_COMA) {
        column = coma_column(o, dir, t0, t1, jitter);
    } else {
        float across = dot(dir.yz, dir.yz);
        float tc = across > 1e-12 ? clamp(-dot(o.yz, dir.yz) / across, t0, t1) : t0;
        float scale = u_shape.x + max(o.x + tc * dir.x, 0.0) * u_shape.y;
        float u0 = arsinh((t0 - tc) / scale);
        float du = (arsinh((t1 - tc) / scale) - u0) / kTailSteps;
        for (int i = 0; i < kTailSteps; ++i) {
            float u = u0 + (i + jitter) * du;
            column += tail_density(o + (tc + scale * sinh(u)) * dir) * cosh(u);
        }
        column *= scale * du;
    }
    return float4(srgb_to_linear(u_color.rgb) * u_shape.w * column, 1.0);
}
