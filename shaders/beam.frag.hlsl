// Pulsar beam: soft across, fading along its length. Linear HDR, additive.

#include "common.hlsli"

cbuffer Uniforms : register(b0, space3)
{
    float4 u_color; // rgb = sRGB color, a = intensity
};

struct PSInput
{
    float4 position : SV_Position;
    float2 coord    : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    float across = exp(-5.0 * input.coord.x * input.coord.x);
    float t = saturate(input.coord.y);
    float along = pow(1.0 - t, 3.0) * smoothstep(0.0, 0.02, t);
    return float4(srgb_to_linear(u_color.rgb) * (u_color.a * across * along), 1.0);
}
