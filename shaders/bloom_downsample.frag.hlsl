// Bloom downsample: 13-tap filter from Jimenez, "Next Generation Post Processing
// in Call of Duty: Advanced Warfare" (SIGGRAPH 2014). The first pass uses a
// Karis average (weighting by 1 / (1 + luma)) to suppress fireflies.

#include "common.hlsli"

Texture2D u_source : register(t0, space2);
SamplerState u_sampler : register(s0, space2);

cbuffer Uniforms : register(b0, space3)
{
    float4 u_params; // xy = source texel size, z = 1 for Karis average, w = threshold (first pass only, 0 = off)
};

// Soft threshold with a quadratic knee: only light above ~threshold feeds the
// bloom, so lit planets stay crisp while the sun glares.
float3 apply_threshold(float3 c, float threshold)
{
    float knee = threshold * 0.5;
    float l = max(max(c.r, c.g), c.b);
    float soft = clamp(l - threshold + knee, 0.0, 2.0 * knee);
    soft = soft * soft / (4.0 * knee + 1e-5);
    float contribution = max(soft, l - threshold) / max(l, 1e-5);
    return c * contribution;
}

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

float3 tap(float2 uv, float2 offset)
{
    return u_source.SampleLevel(u_sampler, uv + offset * u_params.xy, 0).rgb;
}

float karis_weight(float3 c)
{
    return 1.0 / (1.0 + luminance(c));
}

float4 main(PSInput input) : SV_Target0
{
    float2 uv = input.uv;
    float3 a = tap(uv, float2(-2, 2));
    float3 b = tap(uv, float2(0, 2));
    float3 c = tap(uv, float2(2, 2));
    float3 d = tap(uv, float2(-2, 0));
    float3 e = tap(uv, float2(0, 0));
    float3 f = tap(uv, float2(2, 0));
    float3 g = tap(uv, float2(-2, -2));
    float3 h = tap(uv, float2(0, -2));
    float3 i = tap(uv, float2(2, -2));
    float3 j = tap(uv, float2(-1, 1));
    float3 k = tap(uv, float2(1, 1));
    float3 l = tap(uv, float2(-1, -1));
    float3 m = tap(uv, float2(1, -1));

    float3 result;
    if (u_params.z > 0.5) {
        // Five 2x2 box groups, each Karis-weighted.
        float3 g0 = (j + k + l + m) * 0.25;
        float3 g1 = (a + b + d + e) * 0.25;
        float3 g2 = (b + c + e + f) * 0.25;
        float3 g3 = (d + e + g + h) * 0.25;
        float3 g4 = (e + f + h + i) * 0.25;
        float w0 = karis_weight(g0) * 0.5;
        float w1 = karis_weight(g1) * 0.125;
        float w2 = karis_weight(g2) * 0.125;
        float w3 = karis_weight(g3) * 0.125;
        float w4 = karis_weight(g4) * 0.125;
        result = (g0 * w0 + g1 * w1 + g2 * w2 + g3 * w3 + g4 * w4) / (w0 + w1 + w2 + w3 + w4);
    } else {
        result = e * 0.125;
        result += (a + c + g + i) * 0.03125;
        result += (b + d + f + h) * 0.0625;
        result += (j + k + l + m) * 0.125;
    }
    if (u_params.w > 0.0) {
        result = apply_threshold(result, u_params.w);
    }
    return float4(max(result, 0.0), 1.0);
}
