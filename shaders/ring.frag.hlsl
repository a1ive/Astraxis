// Planetary ring as a thin slab of normal optical depth tau, single scattering
// (single-scattering albedo w, phase function P normalized to a mean of 1):
//   lit face:   I/F = w P / 4 * mu0 / (mu + mu0) * (1 - exp(-tau (1/mu + 1/mu0)))
//   unlit face: I/F = w P / 4 * mu0 / (mu - mu0) * (exp(-tau/mu) - exp(-tau/mu0))
// with mu, mu0 the |cos| of the emission and incidence angles. For tau << 1 both
// become w P tau / (4 mu) (Jupiter's dusty rings); thick parts (Saturn's B ring)
// saturate on the lit face and go dark on the unlit one. mu is floored by the
// profile's thickness / width, so a thick halo seen edge-on stays finite.
// Alpha is the direct transmission exp(-tau/mu): whatever lies behind (planet,
// stars) shows through by that fraction. The planet's shadow is cast by its
// ellipsoid with a penumbra from the sun's disk. Everything is in body-fixed
// coordinates, where the ring plane is z = 0.
//
// Narrow eccentric and inclined rings (Uranus') are instead summed band by
// band: each edge is an ellipse with the planet at a focus, and the band lies
// at height z = cross(node, u) in direction u. Along the line of sight through
// this pixel's point on the z = 0 plane, the band's height is reached a little
// farther on (one fixed-point step). Its optical depth is filtered with a tent
// two pixels wide (like the profile's trilinear mipmaps), so a ring much
// narrower than a pixel keeps its equivalent width and its edges stay smooth.

#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_profile, u_sampler, 0); // r = normal optical depth, g = floor for mu

cbuffer Uniforms : register(b0, space3)
{
    float4 u_sun;    // xyz = unit direction to the sun (body frame), w = sun angular radius (rad)
    float4 u_camera; // xyz = camera position (body frame, km)
    float4 u_color;  // rgb = sRGB tint (albedo), a = gain
    float4 u_params; // x = phase asymmetry g, y = equatorial radius, z = polar radius (km)
    float4 u_radii;  // x = inner, y = outer radius of the profile (km), z = number of bands drawn one by one
    float4 u_band_edges[16];  // x, y = semi-major axis of the inner, outer edge (km); z, w = their eccentricity
    float4 u_band_shape[16];  // xy = unit vector to periapsis; zw = toward the ascending node, length a sin i (km)
    float4 u_band_optics[16]; // x = mean tau x mean width (km), y = floor for mu
};

// A view along the plane would shift a ring's crossing without bound.
static const float kMaxSlope = 1000.0;

struct PSInput
{
    float4 position : SV_Position;
    float3 local    : TEXCOORD0;
};

float henyey_greenstein(float cos_theta, float g)
{
    float d = 1.0 + g * g - 2.0 * g * cos_theta;
    return (1.0 - g * g) / (d * sqrt(d));
}

// Fraction of the sun's disk not hidden by the planet (scaled to a sphere).
float planet_shadow(float3 p, float3 L)
{
    float a = u_params.y;
    float stretch = a / u_params.z; // polar axis scaled up to the equatorial radius
    float3 q = float3(p.xy, p.z * stretch);
    float3 l = normalize(float3(L.xy, L.z * stretch));
    float along = dot(-q, l);
    if (along <= 0.0) {
        return 1.0; // the planet is behind, as seen toward the sun
    }
    float perp = length(-q - along * l);
    float penumbra = max(along * u_sun.w, 1.0);
    return smoothstep(a - penumbra, a + penumbra, perp);
}

// Integral of a unit-area tent of half-width h, centered on 0, up to x.
float tent_cdf(float x, float h)
{
    float t = clamp(x / h, -1.0, 1.0);
    return t < 0.0 ? 0.5 * (1.0 + t) * (1.0 + t) : 1.0 - 0.5 * (1.0 - t) * (1.0 - t);
}

float4 main(PSInput input) : SV_Target0
{
    float3 p = input.local;
    float r = length(p.xy);
    float u = (r - u_radii.x) / (u_radii.y - u_radii.x);
    // Line of sight: horizontal shift per km of height.
    float3 ray = p - u_camera.xyz;
    float run = max(abs(ray.z), length(ray.xy) / kMaxSlope + 1e-6);
    float2 slope = ray.xy / (ray.z < 0.0 ? -run : run);
    float2 radial = p.xy / max(r, 1e-6);
    // Before any branch: implicit derivatives (mip selection, footprints) need uniform control flow.
    float2 profile = u_profile.Sample(u_sampler, float2(u, 0.5)).rg;
    // Per pixel, not per 2x2 quad: coarse derivatives bead a ring seen at a low angle.
    float footprint = abs(ddx_fine(r)) + abs(ddy_fine(r));
    float shift = dot(slope, radial);
    float footprint_per_height = abs(ddx_fine(shift)) + abs(ddy_fine(shift));

    float tau = 0.0;
    float mu_floor = 0.0;
    float3 shade_at = p;
    int bands = (int)u_radii.z;
    if (bands == 0) {
        tau = (u >= 0.0 && u <= 1.0) ? profile.r : 0.0;
        mu_floor = profile.g;
    } else {
        float floor_weighted = 0.0;
        float3 weighted = 0.0;
        for (int k = 0; k < bands; ++k) {
            float4 edges = u_band_edges[k];
            float4 shape = u_band_shape[k];
            float2 q = p.xy;
            float z = 0.0;
            [unroll] for (int iteration = 0; iteration < 2; ++iteration) {
                float2 dir = normalize(q);
                z = shape.z * dir.y - shape.w * dir.x;
                q = p.xy + slope * z;
            }
            float rq = length(q);
            float c = dot(q / rq, shape.xy);
            float inner = edges.x * (1.0 - edges.z * edges.z) / (1.0 + edges.z * c);
            float outer = edges.y * (1.0 - edges.w * edges.w) / (1.0 + edges.w * c);
            float half_width = footprint + abs(z) * footprint_per_height + 1e-3;
            float covered = tent_cdf(outer - rq, half_width) - tent_cdf(inner - rq, half_width);
            // tau x width is kept where the ring narrows or widens.
            float t = u_band_optics[k].x / (outer - inner) * covered;
            tau += t;
            floor_weighted += t * u_band_optics[k].y;
            weighted += t * float3(q, z);
        }
        if (tau > 0.0) {
            mu_floor = floor_weighted / tau;
            shade_at = weighted / tau;
        }
    }

    float3 L = u_sun.xyz;
    float3 V = normalize(u_camera.xyz - p);
    float mu = max(abs(V.z), max(mu_floor, 1e-4));
    float mu0 = max(abs(L.z), 1e-4);

    float reflectance;
    if (V.z * L.z > 0.0) {
        reflectance = mu0 / (mu + mu0) * (1.0 - exp(-tau * (1.0 / mu + 1.0 / mu0)));
    } else {
        float d = mu - mu0;
        reflectance = abs(d) > 1e-4 ? mu0 / d * (exp(-tau / mu) - exp(-tau / mu0)) : tau / mu0 * exp(-tau / mu0);
    }

    float phase = henyey_greenstein(dot(-L, V), u_params.x);
    float3 radiance = srgb_to_linear(u_color.rgb) * u_color.a * 0.25 * phase * reflectance *
                      planet_shadow(shade_at, L);
    return float4(radiance, exp(-tau / mu));
}
