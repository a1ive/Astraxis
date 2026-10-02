// Shared helpers. Scene rendering is in linear HDR; color constants authored
// in sRGB are converted with srgb_to_linear.

#ifndef ASTRAXIS_COMMON_HLSLI
#define ASTRAXIS_COMMON_HLSLI

// HLSL 2021 has no vector ternary; use select().
float3 srgb_to_linear(float3 c)
{
    return select(c <= 0.04045, c / 12.92, pow((max(c, 0.0) + 0.055) / 1.055, 2.4));
}

float3 linear_to_srgb(float3 c)
{
    c = saturate(c);
    return select(c <= 0.0031308, c * 12.92, 1.055 * pow(c, 1.0 / 2.4) - 0.055);
}

float luminance(float3 c)
{
    return dot(c, float3(0.2126, 0.7152, 0.0722));
}

#endif
