// Final composite: HDR scene + bloom -> exposure -> ACES tonemap -> sRGB, with
// dithering to avoid banding in dark gradients.

#include "common.hlsli"

Texture2D u_scene : register(t0, space2);
Texture2D u_bloom : register(t1, space2);
SamplerState u_scene_sampler : register(s0, space2);
SamplerState u_bloom_sampler : register(s1, space2);

cbuffer Uniforms : register(b0, space3)
{
    float4 u_params; // x = exposure, y = bloom strength
};

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

// ACES filmic curve fit by Krzysztof Narkowicz (2015), "ACES Filmic Tone Mapping Curve".
float3 aces_fitted(float3 x)
{
    const float a = 2.51;
    const float b = 0.03;
    const float c = 2.43;
    const float d = 0.59;
    const float e = 0.14;
    return saturate((x * (a * x + b)) / (x * (c * x + d) + e));
}

// Interleaved gradient noise (Jimenez 2014).
float ign(float2 pixel)
{
    return frac(52.9829189 * frac(dot(pixel, float2(0.06711056, 0.00583715))));
}

float4 main(PSInput input) : SV_Target0
{
    float3 scene = u_scene.SampleLevel(u_scene_sampler, input.uv, 0).rgb;
    float3 bloom = u_bloom.SampleLevel(u_bloom_sampler, input.uv, 0).rgb;
    float3 hdr = scene + bloom * u_params.y;

    float3 color = linear_to_srgb(aces_fitted(hdr * u_params.x));
    color += (ign(input.position.xy) - 0.5) / 255.0;
    return float4(color, 1.0);
}
