// Black hole ray tracer, shared by the reduced-resolution prepass
// (black_hole.frag.hlsl) and the full-resolution composite, which traces again
// where the prepass cannot be interpolated (black_hole_composite.frag.hlsl).
// Light rays are traced backward from the camera through the Kerr metric
// (Cartesian Kerr-Schild coordinates, G = c = M = 1, spin along +z), with a
// thin accretion disk in the z = 0 plane and the lensed starfield.
// Same formulas as the CPU reference src/ephem/kerr_null.cpp (tested there).
//
// Pixel directions are angles in the rest frame of a static camera: the
// photon's momentum is a linear combination of the camera's orthonormal frame
// (ObserverFrame, built on the CPU), so it is null and has observed energy 1.
//
// Only rays with impact parameter b = L / E < 0.8 R_s are integrated, inside
// the radius R_s; outside it the weak field is handled analytically in these
// (areal radius) coordinates: the swept angle of a leg (weak_leg_angle), the
// exact entry point into R_s, and the deflection of rays that pass farther out
// (weak_pass_deflection). The CPU reference of the whole ray is
// trace_camera_ray(). The lensed region ends where the deflection drops below
// about half a pixel (b_max), with a soft edge blending into the normally
// drawn starfield.

#ifndef ASTRAXIS_BLACK_HOLE_HLSLI
#define ASTRAXIS_BLACK_HOLE_HLSLI

#include "common.hlsli"
#include "noise.hlsli"

FRAGMENT_TEXTURE(TextureCube, u_sky, u_sky_sampler, 0);     // Milky Way and stars in ICRF directions, mipmapped (where lensing squeezes the sky)
FRAGMENT_TEXTURE(Texture2D, u_blackbody, u_lut_sampler, 1); // log-temperature -> linear color

cbuffer Uniforms : register(b0, space3)
{
    float4 u_cam;         // xyz = camera position (hole frame, units of M), w = spin a
    float4 u_right;       // covariant right axis of the camera's frame * tan(fov/2) * aspect (w = t component)
    float4 u_up;          // covariant up axis * tan(fov/2)
    float4 u_forward;     // covariant forward axis
    float4 u_observer;    // covariant 4-velocity of the camera
    float4 u_right_dir;   // xyz = coordinate direction of the right axis * tan(fov/2) * aspect, w = tan(fov/2) * aspect
    float4 u_up_dir;      // xyz = ... of the up axis * tan(fov/2), w = tan(fov/2)
    float4 u_forward_dir; // xyz = ... of the forward axis, w = u^t of the camera
    float4x4 u_to_sky;    // hole frame -> ICRF (rotation)
    float4 u_radii;    // x = horizon r+, y = disk inner radius, z = disk outer radius, w = integration radius R_s
    float4 u_params;   // x = b_max, y = depth of the hole, z = disk brightness, w = disk time (M)
    float4 u_disk;     // x = peak temperature (K), y = 1 / max of the temperature profile, z = max steps, w = edge fade start (fraction of b_max)
    float4 u_targets;  // xy = size of the traced (reduced) targets, zw = scene size, in pixels
    float4 u_stars;    // x = sky grid cells per axis, y = star brightness, z = Milky Way brightness (0: none)
    float4 u_sky_params; // x = largest star dot half size (pixels), y = largest magnification of a star image,
                         // z = 1: uniform disk temperature, w = 1: disk frequency shifts on
    float4 u_flow0;      // hot flow (RIAF): x = 1: on (no thin disk), y = outer radius, z = ISCO radius
    float4 u_flow1;      // x = optical depth scale, y = nu_c scale (GHz), z = T_e0 (1e10 K), w = observed frequency (GHz)
    float4 u_flow2;      // x = display gain, y = far-view fade
    float4 u_spot0;      // flare hot spot: x = orbit radius (M), y = Gaussian sigma (M), z = phase now (rad), w = Omega (per M)
    float4 u_spot1;      // x = tau0 x strength (0: none), y = source temperature (1e10 K), z = |camera| (M)
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

// v rotated by angle about the unit axis (v perpendicular to it).
float3 rotate_in_plane(float3 v, float3 axis, float angle)
{
    return v * cos(angle) + cross(axis, v) * sin(angle);
}

// Angle swept between infinity and r = b / w along a monotonic leg (w <= ~0.8),
// to second order in M / r. See weak_leg_angle() in kerr_null.cpp.
float weak_leg_angle(float w, float b)
{
    float th = asin(w);
    float c = cos(th);
    float one_minus_c = w * w / (1.0 + c);
    float k1 = one_minus_c * one_minus_c / c;
    float t = w / c;
    float k2 = t * t * t / 3.0 - 2.0 * t + 2.5 * th - 0.5 * w * c;
    return th - k1 / b + 1.5 * k2 / (b * b);
}

// Deflection still to come from the tangent at radius r (b >~ 100); see
// weak_pass_deflection() in kerr_null.cpp.
float weak_pass_deflection(float r, float b, bool inward)
{
    float s_in = -4.0 / (b + sqrt(b * b + 8.0));
    float phi_in = asin(s_in);
    float phi_out = kPi - phi_in;
    float s = 2.0 * (b * b / r - 2.0) / (b + sqrt(max(b * b * (1.0 - 4.0 / r) + 8.0, 0.0)));
    s = clamp(s, -1.0, 1.0);
    float phi = inward ? asin(s) : kPi - asin(s);
    float u = sin(phi) / b + (1.0 + cos(phi) * cos(phi)) / (b * b);
    float du = cos(phi) * (1.0 / b - 2.0 * sin(phi) / (b * b));
    float first = phi_out - phi - atan2(u, -du);
    return first + (phi_out - phi) / (phi_out - phi_in) * 15.0 * kPi / (4.0 * b * b);
}

float3 blackbody(float t)
{
    float u = saturate(log(t / 1000.0) / log(40.0)); // 1000 K .. 40000 K
    return u_blackbody.SampleLevel(u_lut_sampler, float2(u, 0.5), 0).rgb;
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
    // u_forward_dir.w = u^t of the (static) camera: p.u_camera = p_t u^t.
    float g = pt * u_forward_dir.w / (ut * (pt - omega * c.y * p.x + omega * c.x * p.y));
    g = clamp(g, 0.0, 5.0);
    // Off: the disk as filmed for Interstellar, without Doppler and gravitational shifts.
    g = u_sky_params.w > 0.5 ? g : 1.0;

    // Thin-disk temperature profile T ~ r^-3/4 (1 - sqrt(r_in / r))^1/4, or uniform.
    float profile = u_sky_params.z > 0.5 ? 1.0
                                         : pow(rc, -0.75) * pow(max(1.0 - sqrt(r_in / rc), 0.0), 0.25) * u_disk.y;
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

// Hot accretion flow (RIAF): the formulas of src/ephem/riaf.cpp, in units that keep
// float in range - frequencies in GHz, temperatures in 1e10 K, optical depth through
// u_flow1.x = n_e0 sqrt(2) pi e^2 c / (6 k) M_cm 1e-19 (computed on the CPU).

// K_2 from Abramowitz & Stegun 9.8.1-9.8.8, as bessel_k2() in riaf.cpp.
float flow_bessel_k2(float x)
{
    float k0;
    float k1;
    if (x <= 2.0) {
        float y = 0.25 * x * x;
        float t = (x / 3.75) * (x / 3.75);
        float i0 = 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.360768e-1 + t * 0.45813e-2)))));
        float i1 = x * (0.5 + t * (0.87890594 + t * (0.51498869 + t * (0.15084934 + t * (0.2658733e-1 + t * (0.301532e-2 + t * 0.32411e-3))))));
        k0 = -log(0.5 * x) * i0 + (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.3488590e-1 + y * (0.262698e-2 + y * (0.10750e-3 + y * 0.74e-5))))));
        k1 = log(0.5 * x) * i1 + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.1919402e-1 + y * (-0.110404e-2 + y * (-0.4686e-4)))))));
    } else {
        float y = 2.0 / x;
        float e = exp(-x) * rsqrt(x);
        k0 = e * (1.25331414 + y * (-0.7832358e-1 + y * (0.2189568e-1 + y * (-0.1062446e-1 + y * (0.587872e-2 + y * (-0.251540e-2 + y * 0.53208e-3))))));
        k1 = e * (1.25331414 + y * (0.23498619 + y * (-0.3655620e-1 + y * (0.1504268e-1 + y * (-0.780353e-2 + y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
    return k0 + 2.0 / x * k1;
}

// Gas 4-velocity (contravariant, w = u^t) in Kerr-Schild coordinates: Keplerian
// shells outside the ISCO, a plunge at fixed theta inside (riaf_four_velocity()).
float4 flow_velocity(float3 x, float r, float a)
{
    float a2 = a * a;
    float cos_t = clamp(x.z / r, -1.0, 1.0);
    float cos2 = cos_t * cos_t;
    float sin2 = 1.0 - cos2;
    float rs = max(r, u_flow0.z);
    float omega = 1.0 / (rs * sqrt(rs) + a);
    float sigma_s = rs * rs + a2 * cos2;
    float delta_s = rs * rs - 2.0 * rs + a2;
    float big_a_s = (rs * rs + a2) * (rs * rs + a2) - a2 * delta_s * sin2;
    float g_tt = -(1.0 - 2.0 * rs / sigma_s);
    float g_tp_s = -2.0 * a * rs / sigma_s;
    float g_pp_s = big_a_s / sigma_s;
    float ut_s = rsqrt(max(-(g_tt + sin2 * (2.0 * omega * g_tp_s + omega * omega * g_pp_s)), 1e-6));
    if (r >= u_flow0.z) {
        return float4(-omega * ut_s * x.y, omega * ut_s * x.x, 0.0, ut_s);
    }
    float u_t = ut_s * (g_tt + sin2 * omega * g_tp_s);
    float l_hat = ut_s * (g_tp_s + omega * g_pp_s);
    float sigma = r * r + a2 * cos2;
    float delta = r * r - 2.0 * r + a2;
    float big_a = (r * r + a2) * (r * r + a2) - a2 * delta * sin2;
    float sd = sigma * delta;
    float gi_tt = -big_a / sd;
    float gi_tp = -2.0 * a * r / sd;
    float gi_pp_s = (delta - a2 * sin2) / sd;
    float up_bl = gi_tp * u_t + gi_pp_s * l_hat;
    float ut_bl = gi_tt * u_t + gi_tp * sin2 * l_hat;
    float rest = -1.0 - (gi_tt * u_t * u_t + 2.0 * gi_tp * sin2 * u_t * l_hat + gi_pp_s * sin2 * l_hat * l_hat);
    float ur = -sqrt(max(delta / sigma * rest, 0.0));
    float s = r * r + a2;
    float3 l = float3((r * x.x + a * x.y) / s, (r * x.y - a * x.x) / s, x.z / r);
    float up_ks = up_bl + a / delta * ur;
    return float4(l * ur + up_ks * float3(-x.y, x.x, 0.0), ut_bl + 2.0 * r / delta * ur);
}

// One step of the transfer at the end point x of a step of affine length h:
// adds the gas's brightness temperature (1e10 K) seen through flow_transmit.
void flow_step(float3 x, float3 p, float pt, float a, float h, inout float flow_tb, inout float flow_transmit)
{
    float r = ks_r(x, a);
    if (r > u_flow0.y || r <= u_radii.x) {
        return;
    }
    float half_r = 0.5 * r;
    float s2 = max(r * r - x.z * x.z, 1e-6);
    float n_rel = pow(half_r, -1.1) * exp(-x.z * x.z / (2.0 * s2));
    float t10 = u_flow1.z * pow(half_r, -0.84);
    float theta = 1.6863699 * t10; // k T / m_e c^2
    if (theta < 0.02 || n_rel < 1e-12) { // as riaf.cpp; below ~0.0116 K_2 and the emission underflow (0 / 0)
        return;
    }
    float4 u = flow_velocity(x, r, a);
    float pu = pt * u.w + dot(p, u.xyz); // nu_emit / nu_obs
    if (!(pu > 0.0)) {
        return;
    }
    float nu = u_flow1.w * pu;
    float nu_c = u_flow1.y * sqrt(n_rel / half_r);
    float nu_s = (2.0 / 9.0) * nu_c * theta * theta * 0.78539816;
    float xs = nu / nu_s;
    float shape = sqrt(xs) + 1.8877486 * pow(xs, 1.0 / 6.0); // 2^(11/12)
    float emission = shape * shape * exp(-pow(xs, 1.0 / 3.0));
    float dtau = u_flow1.x * n_rel * nu_s * emission / (flow_bessel_k2(1.0 / theta) * nu * nu * t10) * pu * h;
    float g = 1.0 / pu;
    flow_tb += flow_transmit * g * t10 * (1.0 - exp(-dtau));
    flow_transmit *= exp(-dtau);
}

// False color of a brightness temperature (1e10 K), linear HDR.
static const float kFlowTbLo = 0.3;  // 3e9 K: the bottom of the color ramp
static const float kFlowTbHi = 10.0; // 1e11 K: its top (Doppler-boosted gas seen edge-on)
float3 flow_color(float tb)
{
    static const float3 kRamp[5] = {
        float3(0.06, 0.0, 0.0), float3(0.45, 0.05, 0.0), float3(0.85, 0.35, 0.02), float3(1.0, 0.75, 0.2), float3(1.0, 0.97, 0.85)};
    float u = saturate(log(max(tb, 1e-6) / kFlowTbLo) / log(kFlowTbHi / kFlowTbLo)) * 4.0;
    int i = min((int)u, 3);
    float3 c = srgb_to_linear(lerp(kRamp[i], kRamp[i + 1], u - i));
    // Compressed: the beamed side outshines the rest ~100x. u_flow2.y fades a flow only a
    // few pixels across (riaf_far_view_fade()).
    return c * (u_flow2.x * u_flow2.y * 0.5 * sqrt(tb / kFlowTbHi));
}

// p.p - f (l.p)^2 over 2: grows away from 0 when the RK4 runs away (ray_runaway()).
float ray_hamiltonian(float3 x, float3 p, float pt, float a)
{
    float r = ks_r(x, a);
    float a2 = a * a;
    float f = 2.0 * r * r * r / (r * r * r * r + a2 * x.z * x.z);
    float s = r * r + a2;
    float3 l = float3((r * x.x + a * x.y) / s, (r * x.y - a * x.x) / s, x.z / r);
    float lp = -pt + dot(l, p);
    return 0.5 * (-pt * pt + dot(p, p) - f * lp * lp);
}

// dt/dl of Kerr-Schild coordinate time along the (backward) ray: -p_t + f (l.p), as
// geodesic_derivatives() in kerr_null.cpp.
float ray_time_rate(float3 x, float3 p, float pt, float a)
{
    float r = ks_r(x, a);
    float a2 = a * a;
    float f = 2.0 * r * r * r / (r * r * r * r + a2 * x.z * x.z);
    float s = r * r + a2;
    float3 l = float3((r * x.x + a * x.y) / s, (r * x.y - a * x.x) / s, x.z / r);
    return -pt + f * (-pt + dot(l, p));
}

// A flare's hot spot (HotSpot in riaf.cpp) where it was when this light left it:
// ray_time (<= 0) since the camera plus |camera| for "now" at the hole.
void spot_step(float3 x, float3 p, float pt, float a, float h, float ray_time, inout float flow_tb,
               inout float flow_transmit)
{
    float phi = u_spot0.z + u_spot0.w * (ray_time + u_spot1.z);
    float c = cos(phi);
    float s = sin(phi);
    float3 d = x - float3(u_spot0.x * c - a * s, u_spot0.x * s + a * c, 0.0);
    float sigma = u_spot0.y;
    float d2 = dot(d, d);
    if (d2 > 25.0 * sigma * sigma) {
        return;
    }
    float4 u = flow_velocity(x, ks_r(x, a), a);
    float pu = pt * u.w + dot(p, u.xyz); // nu_emit / nu_obs
    if (!(pu > 0.0)) {
        return;
    }
    float dtau = u_spot1.x * exp(-0.5 * d2 / (sigma * sigma)) / (sigma * 2.5066283) * pu * h;
    flow_tb += flow_transmit * u_spot1.y / pu * (1.0 - exp(-dtau));
    flow_transmit *= exp(-dtau);
}

// What the camera sees through the pixel at uv (0..1 over the scene, y down).
struct Trace
{
    bool lensed;    // inside the lensed region; elsewhere the normal starfield shows
    float3 disk;    // disk light reaching the camera
    float transmit; // fraction of the sky seen through the disk (0: hidden by the hole)
    float fade;     // coverage, falling to 0 at the edge of the lensed region
    float3 sky;     // ICRF direction the sky is seen in (zero where transmit = 0)
    bool opaque;    // hit the hole or the disk: drawn at the hole's depth
};

Trace trace_pixel(float2 uv)
{
    Trace result = (Trace)0;
    float a = u_cam.w;
    float2 ndc = float2(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0);
    float3 o = u_cam.xyz;
    float dist = length(o);

    // Photon from this pixel's direction in the camera's frame (traced backward):
    // its covariant momentum, coordinate direction and impact parameter L / E.
    float2 screen = ndc * float2(u_right_dir.w, u_up_dir.w);
    float4 pc = (u_forward + ndc.x * u_right + ndc.y * u_up) * rsqrt(1.0 + dot(screen, screen)) - u_observer;
    float3 d = normalize(u_forward_dir.xyz + ndc.x * u_right_dir.xyz + ndc.y * u_up_dir.xyz);
    float3 l_vec = cross(o, pc.xyz);
    float b = length(l_vec) / pc.w;
    float3 axis = l_vec / max(length(l_vec), 1e-30);
    float cos_in = -dot(o, d) / dist; // > 0: heading inward
    bool inward = cos_in > 0.0;

    float R = u_radii.w;
    bool integrate = b < 0.8 * R && (inward || dist <= R);
    // Lensed region: where the deflection still to come, ~(2 / b) (1 + cos_in), exceeds
    // about half a pixel (b_eff = 4 / deflection < b_max).
    float b_eff = integrate ? b : 2.0 * b / max(1.0 + cos_in, 1e-6);
    if (b_eff > u_params.x) {
        return result;
    }
    float fade = 1.0 - smoothstep(u_disk.w * u_params.x, u_params.x, b_eff);

    float3 color = 0.0;
    float transmit = 1.0;
    bool opaque = false;
    float3 sky_dir = d;

    if (!integrate) {
        // Weak field: passing beyond 0.8 R_s, or receding from a camera outside R_s.
        sky_dir = b >= 0.8 * R ? rotate_in_plane(d, axis, weak_pass_deflection(dist, b, inward))
                               : rotate_in_plane(o / dist, axis, weak_leg_angle(b / dist, b));
    } else {
        float3 x = o;
        float3 p = pc.xyz;
        float pt = pc.w;
        if (dist > R) {
            // Enter R_s where the ray crosses it: swept angle of the leg, tangential p = L / R,
            // radial p from H = 0 (Schwarzschild, f = 2 / R), then p_t made exactly null.
            float3 xr = rotate_in_plane(o / dist, axis, weak_leg_angle(b / R, b) - weak_leg_angle(b / dist, b));
            float f = 2.0 / R;
            float p_tan = b * pt / R;
            float p_rad = (-f * pt - sqrt(max(pt * pt - (1.0 - f) * p_tan * p_tan, 0.0))) / (1.0 - f);
            x = R * xr;
            p = p_rad * xr + p_tan * cross(axis, xr);
            pt = null_pt(x, p, a);
        }
        float r_plus = u_radii.x;
        bool escaped = false;
        int max_steps = (int)u_disk.z;
        bool flow = u_flow0.x > 0.5;
        float flow_tb = 0.0;       // brightness temperature of the hot flow (1e10 K)
        float flow_transmit = 1.0; // the flow's own transparency at its frequency (not the sky's)
        bool spot = flow && u_spot1.x > 0.0;
        // Coordinate time since the camera; the leg outside R_s as a straight line.
        float ray_time = dist > R ? -(dist - R) : 0.0;
        float time_rate = spot ? ray_time_rate(x, p, pt, a) : 0.0;

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

            // A ray spiralling into a fast-spinning hole can make the RK4 run away just
            // outside the horizon (ray_runaway() in kerr_null.cpp): it was falling in.
            if (ks_r(x, a) < 2.0 && abs(ray_hamiltonian(x, p, pt, a)) > 0.01 * pt * pt) {
                opaque = true;
                transmit = 0.0;
                break;
            }
            if (spot) {
                float next_rate = ray_time_rate(x, p, pt, a);
                ray_time += 0.5 * h * (time_rate + next_rate);
                time_rate = next_rate;
                spot_step(x, p, pt, a, h, ray_time, flow_tb, flow_transmit);
            }
            if (flow) {
                flow_step(x, p, pt, a, h, flow_tb, flow_transmit);
            }
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
                // The rest of the way out, analytically.
                float3 l_out = cross(x, p);
                float rx = length(x);
                float b_out = length(l_out) / pt;
                sky_dir = rotate_in_plane(x / rx, l_out / max(length(l_out), 1e-30),
                                          weak_leg_angle(min(b_out / rx, 0.95), b_out));
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
        if (flow) {
            color += flow_color(flow_tb); // transparent in visible light: the sky shows through
        }
    }

    result.lensed = true;
    result.disk = color;
    result.transmit = transmit;
    result.fade = fade;
    if (transmit > 0.0) {
        result.sky = mul((float3x3)u_to_sky, normalize(sky_dir));
    }
    result.opaque = opaque;
    return result;
}

#endif
