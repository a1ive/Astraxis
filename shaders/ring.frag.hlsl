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

#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_profile, u_sampler, 0); // r = normal optical depth, g = floor for mu

cbuffer Uniforms : register(b0, space3)
{
    float4 u_sun;    // xyz = unit direction to the sun (body frame), w = sun angular radius (rad)
    float4 u_camera; // xyz = camera position (body frame, km)
    float4 u_color;  // rgb = sRGB tint (albedo), a = gain
    float4 u_params; // x = phase asymmetry g, y = equatorial radius, z = polar radius (km)
    float4 u_radii;  // x = inner, y = outer radius of the profile (km)
};

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

float4 main(PSInput input) : SV_Target0
{
    float r = length(input.local.xy);
    float u = (r - u_radii.x) / (u_radii.y - u_radii.x);
    // Before any branch: implicit derivatives (mip selection) need uniform control flow.
    float2 profile = u_profile.Sample(u_sampler, float2(u, 0.5)).rg;
    float tau = (u >= 0.0 && u <= 1.0) ? profile.r : 0.0;

    float3 L = u_sun.xyz;
    float3 V = normalize(u_camera.xyz - input.local);
    float mu = max(abs(V.z), max(profile.g, 1e-4));
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
                      planet_shadow(input.local, L);
    return float4(radiance, exp(-tau / mu));
}
