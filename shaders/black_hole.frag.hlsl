// Black hole ray tracer: light rays traced backward from the camera through the
// Kerr metric (Cartesian Kerr-Schild coordinates, G = c = M = 1, spin along +z),
// with a thin accretion disk in the z = 0 plane and the lensed starfield.
// Same formulas as the CPU reference src/ephem/kerr_null.cpp (tested there).
//
// Far from the hole (impact parameter b > R_s) rays are not integrated: the
// weak-field deflection 4/b + 15 pi / (4 b^2) is applied instead, and the
// integrated rays get the weak-field deflection accumulated outside R_s on
// their way in and out, so the two regions join without a seam. The region
// ends where the deflection drops below about half a pixel (b_max), with a
// soft edge blending into the normally drawn starfield.
//
// Rendered off-screen (possibly at reduced resolution) into a color target and
// a depth target; black_hole_composite.frag.hlsl then writes both into the scene.

#include "common.hlsli"

TextureCube u_sky : register(t0, space2);    // starfield in ICRF directions
Texture2D u_blackbody : register(t1, space2); // log-temperature -> linear color
SamplerState u_sky_sampler : register(s0, space2);
SamplerState u_lut_sampler : register(s1, space2);

cbuffer Uniforms : register(b0, space3)
{
    float4 u_cam;      // xyz = camera position (hole frame, units of M), w = spin a
    float4 u_right;    // xyz = camera right * tan(fov/2) * aspect (hole frame)
    float4 u_up;       // xyz = camera up * tan(fov/2)
    float4 u_forward;  // xyz = camera forward
    float4x4 u_to_sky; // hole frame -> ICRF (rotation)
    float4 u_radii;    // x = horizon r+, y = disk inner radius, z = disk outer radius, w = integration radius R_s
    float4 u_params;   // x = b_max, y = depth of the hole, z = disk brightness, w = disk time (M)
    float4 u_disk;     // x = peak temperature (K), y = 1 / max of the temperature profile, z = max steps, w = edge fade start (fraction of b_max)
};

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

struct PSOutput
{
    float4 color : SV_Target0; // premultiplied rgb, a = coverage (soft edge)
    float4 depth : SV_Target1; // r = scene depth to write (0 = far: lensed sky only)
};

static const float kPi = 3.14159265;

float ks_r(float3 x, float a)
{
    float a2 = a * a;
    float w = dot(x, x) - a2;
    float r2 = 0.5 * (w + sqrt(w * w + 4.0 * a2 * x.z * x.z));
    return sqrt(max(r2, 0.0));
}

// Metric function f, spatial l_i and their derivatives (see kerr_null.cpp).
void kerr_schild(float3 p, float a, out float f, out float3 l, out float3 df, out float3x3 dl)
{
    float x = p.x;
    float y = p.y;
    float z = p.z;
    float a2 = a * a;
    float w = x * x + y * y + z * z - a2;
    float root = sqrt(w * w + 4.0 * a2 * z * z);
    float r2 = 0.5 * (w + root);
    float r = sqrt(r2);
    float3 dr = float3(x * r / root, y * r / root, z * (r2 + a2) / (r * root));

    float num = 2.0 * r2 * r;
    float den = r2 * r2 + a2 * z * z;
    f = num / den;
    float3 dden = 4.0 * r2 * r * dr + float3(0.0, 0.0, 2.0 * a2 * z);
    df = (6.0 * r2 * dr * den - num * dden) / (den * den);

    float s = r2 + a2;
    float nx = r * x + a * y;
    float ny = r * y - a * x;
    l = float3(nx / s, ny / s, z / r);
    float3 dnx = dr * x + float3(r, a, 0.0);
    float3 dny = dr * y + float3(-a, r, 0.0);
    float3 ds = 2.0 * r * dr;
    dl[0] = (dnx * s - nx * ds) / (s * s);
    dl[1] = (dny * s - ny * ds) / (s * s);
    dl[2] = (float3(0.0, 0.0, r) - z * dr) / r2;
}

void derivatives(float3 x, float3 p, float pt, float a, out float3 dx, out float3 dp)
{
    float f;
    float3 l;
    float3 df;
    float3x3 dl;
    kerr_schild(x, a, f, l, df, dl);
    float L = -pt + dot(l, p);
    dx = p - f * L * l;
    // dl[i][j] = d l_i / d x_j: dlp_j = sum_i p_i dl[i][j]
    float3 dlp = p.x * dl[0] + p.y * dl[1] + p.z * dl[2];
    dp = 0.5 * (df * L * L + 2.0 * f * L * dlp);
}

// p_t making the (time-reversed) photon null.
float null_pt(float3 x, float3 p, float a)
{
    float f;
    float3 l;
    float3 df;
    float3x3 dl;
    kerr_schild(x, a, f, l, df, dl);
    float L = dot(l, p);
    float disc = f * f * L * L + (1.0 + f) * (dot(p, p) - f * L * L);
    return (f * L + sqrt(max(disc, 0.0))) / (1.0 + f);
}

// Rotates d toward the hole (the closest-approach point c) by angle alpha.
float3 bend(float3 d, float3 c, float alpha)
{
    float3 n = -normalize(c);
    return normalize(d * cos(alpha) + n * sin(alpha));
}

// Weak-field deflection accumulated beyond radius R on one leg of a ray with
// impact parameter b: (2 / b) (1 - sqrt(R^2 - b^2) / R).
float outer_leg_deflection(float b, float R)
{
    return 2.0 / b * (1.0 - sqrt(max(R * R - b * b, 0.0)) / R);
}

float3 blackbody(float t)
{
    float u = saturate(log(t / 1000.0) / log(40.0)); // 1000 K .. 40000 K
    return u_blackbody.SampleLevel(u_lut_sampler, float2(u, 0.5), 0).rgb;
}

float hash12(float2 q)
{
    float3 p3 = frac(float3(q.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.x + p3.y) * p3.z);
}

// Value noise, periodic in x with integer period `period`.
float periodic_noise(float2 q, float period)
{
    float2 i = floor(q);
    float2 u = q - i;
    u = u * u * (3.0 - 2.0 * u);
    float x0 = i.x - period * floor(i.x / period);
    float x1 = x0 + 1.0 - period * floor((x0 + 1.0) / period);
    float a = hash12(float2(x0, i.y));
    float b = hash12(float2(x1, i.y));
    float c = hash12(float2(x0, i.y + 1.0));
    float d = hash12(float2(x1, i.y + 1.0));
    return lerp(lerp(a, b, u.x), lerp(c, d, u.x), u.y);
}

// Disk texture in [0, 1]: q.x = azimuth in turns (periodic), q.y = log r.
// Cells are ~7x longer in azimuth than in radius at every octave.
float disk_noise(float2 q)
{
    float sum = 0.0;
    float amp = 0.5;
    float period = 8.0;
    for (int k = 0; k < 4; ++k) {
        sum += amp * periodic_noise(float2(q.x * period, q.y * period * 1.2 + 17.0 * k), period);
        period *= 2.0;
        amp *= 0.5;
    }
    return sum / 0.9375;
}

// Emission of the disk at crossing point c by a photon with momentum (p, pt).
float3 disk_emission(float3 c, float3 p, float pt, float a, out float alpha)
{
    float rc = ks_r(float3(c.xy, 0.0), a);
    float r_in = u_radii.y;
    float r_out = u_radii.z;

    // Redshift g = nu_obs / nu_emit for a prograde Keplerian emitter (see disk_redshift()).
    float f;
    float3 l;
    float3 df;
    float3x3 dl;
    kerr_schild(float3(c.xy, 0.0), a, f, l, df, dl);
    float omega = 1.0 / (rc * sqrt(rc) + a);
    float lv = 1.0 + omega * (-c.y * l.x + c.x * l.y);
    float gvv = -1.0 + omega * omega * dot(c.xy, c.xy) + f * lv * lv;
    float ut = 1.0 / sqrt(max(-gvv, 1e-6));
    float g = pt / (ut * (pt - omega * c.y * p.x + omega * c.x * p.y));
    g = clamp(g, 0.0, 5.0);

    // Thin-disk temperature profile T ~ r^-3/4 (1 - sqrt(r_in / r))^1/4.
    float profile = pow(rc, -0.75) * pow(max(1.0 - sqrt(r_in / rc), 0.0), 0.25) * u_disk.y;
    float t_emit = u_disk.x * profile;
    float t_obs = g * t_emit;
    // Bolometric radiance: emitted ~ T^4, observed ~ g^4 (I / nu^3 invariant).
    float radiance = u_params.z * pow(profile, 4.0) * pow(g, 4.0);

    // Turbulent streaks: noise in (log r, phi), stretched along phi and advected
    // at the local Kepler rate, so differential rotation shears it into trailing arcs.
    float phi = atan2(c.y, c.x) - omega * u_params.w;
    float lr = log(rc);
    float n = disk_noise(float2(frac(phi * (1.0 / 6.2831853)), lr)); // frac: keep float precision
    float streaks = 0.35 + 1.3 * n * n;
    alpha = 0.95 * smoothstep(r_in, r_in * 1.15, rc) * (1.0 - smoothstep(r_out * 0.6, r_out, rc));
    return blackbody(max(t_obs, 1000.0)) * radiance * streaks;
}

PSOutput main(PSInput input)
{
    PSOutput output;
    float a = u_cam.w;
    float2 ndc = float2(input.uv.x * 2.0 - 1.0, 1.0 - input.uv.y * 2.0);
    float3 d = normalize(u_forward.xyz + ndc.x * u_right.xyz + ndc.y * u_up.xyz);
    float3 o = u_cam.xyz;

    // Closest approach of the straight ray, and its impact parameter.
    float tc = -dot(o, d);
    float3 c0 = o + max(tc, 0.0) * d;
    float b = length(c0);
    if (b > u_params.x) {
        discard;
    }
    float fade = 1.0 - smoothstep(u_disk.w * u_params.x, u_params.x, b);

    float R = u_radii.w;
    float3 color = 0.0;
    float transmit = 1.0;
    bool opaque = false;
    float3 sky_dir = d;

    if (b > R || (tc < 0.0 && length(o) > R)) {
        // Weak field: rays passing the hole ahead are bent toward it.
        float alpha = tc > 0.0 ? 4.0 / b + 15.0 * kPi / (4.0 * b * b) : 0.0;
        sky_dir = bend(d, c0, alpha);
    } else {
        float3 x = o;
        if (length(o) > R) {
            x = o + (tc - sqrt(R * R - b * b)) * d;
            d = bend(d, c0, outer_leg_deflection(b, R));
        }
        float3 p = d;
        float pt = null_pt(x, p, a);
        float r_plus = u_radii.x;
        bool escaped = false;
        int max_steps = (int)u_disk.z;

        [loop]
        for (int i = 0; i < max_steps; ++i) {
            float r = ks_r(x, a);
            float h = 0.07 * max(r - r_plus, 0.0) + 0.015; // as ray_step() in kerr_null.cpp
            float3 prev = x;

            float3 k1x, k1p, k2x, k2p, k3x, k3p, k4x, k4p;
            derivatives(x, p, pt, a, k1x, k1p);
            derivatives(x + 0.5 * h * k1x, p + 0.5 * h * k1p, pt, a, k2x, k2p);
            derivatives(x + 0.5 * h * k2x, p + 0.5 * h * k2p, pt, a, k3x, k3p);
            derivatives(x + h * k3x, p + h * k3p, pt, a, k4x, k4p);
            x += h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
            p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);

            if (prev.z * x.z < 0.0 && u_params.z > 0.0) {
                float3 c = lerp(prev, x, prev.z / (prev.z - x.z));
                float rc = ks_r(float3(c.xy, 0.0), a);
                if (rc >= u_radii.y && rc <= u_radii.z) {
                    float disk_alpha;
                    float3 e = disk_emission(c, p, pt, a, disk_alpha);
                    color += transmit * disk_alpha * e;
                    transmit *= 1.0 - disk_alpha;
                    opaque = true;
                }
            }
            if (ks_r(x, a) < r_plus + 0.02) {
                opaque = true;
                transmit = 0.0;
                break;
            }
            if (dot(x, x) > R * R * 1.0001 && dot(x, p) > 0.0) {
                float3 dx;
                float3 dp;
                derivatives(x, p, pt, a, dx, dp);
                sky_dir = normalize(dx);
                // Deflection still to come beyond R on the way out.
                float3 cx = x - dot(x, sky_dir) * sky_dir;
                float b_out = length(cx);
                sky_dir = bend(sky_dir, cx, outer_leg_deflection(b_out, length(x)));
                escaped = true;
                break;
            }
            if (transmit < 0.02) {
                break;
            }
        }
        if (!escaped && transmit > 0.0) {
            opaque = true; // ran out of steps: treat as captured
            transmit = 0.0;
        }
    }

    if (transmit > 0.0) {
        float3 sky = mul((float3x3)u_to_sky, sky_dir);
        color += transmit * u_sky.SampleLevel(u_sky_sampler, sky, 0).rgb;
    }

    output.color = float4(color * fade, fade); // premultiplied: bilinear upsampling stays exact
    // Hole and disk sit at the hole's depth; lensed sky only replaces empty background.
    output.depth = float4(opaque ? u_params.y : 0.0, 0.0, 0.0, 1.0);
    return output;
}
