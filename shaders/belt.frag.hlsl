// Belt point: a soft round dot, added to the HDR target.

#include "common.hlsli"

cbuffer Uniforms : register(b0, space3)
{
    float4 u_color; // rgb = sRGB tint, a = brightness (linear radiance scale)
};

struct PSInput
{
    float4 position : SV_Position;
    float2 offset   : TEXCOORD0;
    float weight    : TEXCOORD1;
};

float4 main(PSInput input) : SV_Target0
{
    float falloff = saturate(1.0 - dot(input.offset, input.offset));
    return float4(srgb_to_linear(u_color.rgb) * (u_color.a * input.weight * falloff * falloff), 0.0);
}
