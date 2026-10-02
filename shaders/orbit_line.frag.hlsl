// Orbit line: fades from the body (head) toward the tail, anti-aliased edges.

#include "common.hlsli"

cbuffer Uniforms : register(b0, space3)
{
    float4 u_color;  // rgb = sRGB color, a = max opacity
    float4 u_params; // x = width in pixels, y = tail opacity factor
};

struct PSInput
{
    float4 position : SV_Position;
    float fade      : TEXCOORD0;
    float side      : TEXCOORD1;
};

float4 main(PSInput input) : SV_Target0
{
    float half_width = u_params.x * 0.5;
    float edge = saturate(half_width + 0.5 - abs(input.side));
    float trail = lerp(1.0, u_params.y, pow(saturate(input.fade), 0.7));
    return float4(srgb_to_linear(u_color.rgb), u_color.a * trail * edge);
}
