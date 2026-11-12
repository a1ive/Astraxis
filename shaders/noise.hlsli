// Shared hashes and value noise for procedural surfaces and gas.

#ifndef ASTRAXIS_NOISE_HLSLI
#define ASTRAXIS_NOISE_HLSLI

// A hash of a 2D cell in [0, 1).
float hash12(float2 q)
{
    float3 p3 = frac(float3(q.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.x + p3.y) * p3.z);
}

// A hash of a 3D cell in [0, 1).
float hash13(float3 p)
{
    p = frac(p * 0.1031);
    p += dot(p, p.yzx + 33.33);
    return frac((p.x + p.y) * p.z);
}

// Value noise in [0, 1], periodic in x with integer period `period`.
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

// Value noise in [-1, 1] around a circle without a seam: `turns` is the
// azimuth (periodic), with `period` cells per turn; `across` runs freely.
float turn_noise(float turns, float across, float period)
{
    return 2.0 * periodic_noise(float2(frac(turns) * period, across), period) - 1.0;
}

// Value noise in [-1, 1] over 3D space.
float value_noise3(float3 p)
{
    const float3 i = floor(p);
    const float3 t = smoothstep(0.0, 1.0, frac(p));
    const float a = lerp(hash13(i), hash13(i + float3(1, 0, 0)), t.x);
    const float b = lerp(hash13(i + float3(0, 1, 0)), hash13(i + float3(1, 1, 0)), t.x);
    const float c = lerp(hash13(i + float3(0, 0, 1)), hash13(i + float3(1, 0, 1)), t.x);
    const float d = lerp(hash13(i + float3(0, 1, 1)), hash13(i + float3(1, 1, 1)), t.x);
    return 2.0 * lerp(lerp(a, b, t.y), lerp(c, d, t.y), t.z) - 1.0;
}

#endif
