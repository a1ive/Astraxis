// Direct solar image refracted through the atmosphere. The CPU table holds
// curved-ray Beer-Lambert transmission and bending in a spherical atmosphere.
// Uniforms are declared by atmosphere.frag.hlsl.

float4 refracted_profile(float impact)
{
    float index = saturate((impact - u_refraction.x) / (u_refraction.y - u_refraction.x)) * 127.0;
    int lo = min((int)index, 126);
    return lerp(u_refracted[lo], u_refracted[lo + 1], index - lo);
}

float3 refracted_sun(float3 dir, float3 perpendicular, float footprint)
{
    if (u_refraction.z == 0.0 || dot(u_camera.xyz, dir) >= 0.0) return 0.0;
    float impact = length(perpendicular);
    float width = max(footprint, 1e-7);
    float lo = max(u_refraction.x, impact - 0.5 * width);
    float hi = min(u_refraction.y, impact + 0.5 * width);
    if (hi <= lo) return 0.0;
    float3 outward = perpendicular / max(impact, 1e-6);
    float distance = max(-dot(u_camera.xyz, dir), 1.0);
    float3 sum = 0.0;
    // Integrate the pixel's radial footprint: the red image is often much
    // narrower than a pixel. Sampling only its center breaks the ring.
    const int steps = 24;
    for (int i = 0; i < steps; ++i) {
        float b = lerp(lo, hi, (i + 0.5) / steps);
        float4 optical = refracted_profile(b);
        float3 view_ray = normalize(dir + outward * ((b - impact) / distance));
        float3 source = normalize(u_sun_unit.xyz - outward * (b * u_refraction.w));
        float3 bent = cos(optical.w) * view_ray - sin(optical.w) * outward;
        float separation = length(bent - source);
        // The ordinary Sun pass supplies the direct disk. Add only the image
        // outside that disk, so unbent rays do not duplicate its radiance.
        if (separation < u_sun.w && length(view_ray - source) > u_sun.w) {
            float mu = sqrt(saturate(1.0 - separation * separation / (u_sun.w * u_sun.w)));
            sum += optical.rgb * (0.4 + 0.6 * mu);
        }
    }
    return u_sun_color.rgb * sum * ((hi - lo) / (width * steps));
}
