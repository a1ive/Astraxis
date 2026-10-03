// Optically thin dusty ring, single scattering:
//   I/F = albedo * tau * P(theta) / (4 mu)
// (tau: normal optical depth, P: Henyey-Greenstein phase function normalized
// to a mean of 1, mu: |cos| of the emission angle, floored by the band's
// thickness / width so the edge-on view stays finite). The same holds for the
// lit and the unlit face when tau << 1. The planet's shadow is cast by its
// ellipsoid with a penumbra from the sun's disk. Everything is in body-fixed
// coordinates, where the ring plane is z = 0.

#include "common.hlsli"

#define MAX_BANDS 8

cbuffer Uniforms : register(b0, space3)
{
    float4 u_sun;              // xyz = unit direction to the sun (body frame), w = sun angular radius (rad)
    float4 u_camera;           // xyz = camera position (body frame, km)
    float4 u_color;            // rgb = sRGB tint (albedo), a = gain
    float4 u_params;           // x = band count, y = phase asymmetry g, z = equatorial radius, w = polar radius
    float4 u_bands[MAX_BANDS]; // inner, outer radius (km), optical depth, thickness (km)
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
    float a = u_params.z;
    float stretch = a / u_params.w; // polar axis scaled up to the equatorial radius
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
    float aa = fwidth(r); // before any branch: derivatives need uniform control flow

    float3 V = normalize(u_camera.xyz - input.local);
    float mu = abs(V.z);
    float sum = 0.0;
    int count = (int)u_params.x;
    for (int k = 0; k < count; ++k) {
        float4 b = u_bands[k];
        float width = b.y - b.x;
        float soft = max(0.02 * width, aa);
        float inside = smoothstep(b.x - soft, b.x + soft, r) * (1.0 - smoothstep(b.y - soft, b.y + soft, r));
        sum += inside * b.z / (4.0 * max(mu, b.w / width));
    }

    float3 L = u_sun.xyz;
    float phase = henyey_greenstein(dot(-L, V), u_params.y);
    float3 radiance = srgb_to_linear(u_color.rgb) * u_color.a * sum * phase * planet_shadow(input.local, L);
    return float4(radiance, 0.0);
}
