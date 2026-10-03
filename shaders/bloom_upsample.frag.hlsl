// Bloom upsample: 3x3 tent filter, additively blended onto the next larger level.

#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_source, u_sampler, 0);

cbuffer Uniforms : register(b0, space3)
{
    float4 u_params; // xy = filter radius in uv units
};

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

float3 tap(float2 uv, float2 offset)
{
    return u_source.SampleLevel(u_sampler, uv + offset * u_params.xy, 0).rgb;
}

float4 main(PSInput input) : SV_Target0
{
    float2 uv = input.uv;
    float3 sum = tap(uv, float2(0, 0)) * 4.0;
    sum += (tap(uv, float2(-1, 0)) + tap(uv, float2(1, 0)) + tap(uv, float2(0, -1)) + tap(uv, float2(0, 1))) * 2.0;
    sum += tap(uv, float2(-1, -1)) + tap(uv, float2(1, -1)) + tap(uv, float2(-1, 1)) + tap(uv, float2(1, 1));
    return float4(sum / 16.0, 1.0);
}
