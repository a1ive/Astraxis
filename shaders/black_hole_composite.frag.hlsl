// Composites the (possibly reduced-resolution) black hole trace into the scene:
// color with its soft-edge coverage, and the depth the trace asked for.

#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_color, u_linear, 0);
FRAGMENT_TEXTURE(Texture2D, u_depth, u_point, 1);

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

struct PSOutput
{
    float4 color : SV_Target0;
    float depth  : SV_Depth;
};

PSOutput main(PSInput input)
{
    float4 c = u_color.SampleLevel(u_linear, input.uv, 0);
    if (c.a < 1.0 / 255.0) {
        discard;
    }
    PSOutput output;
    output.color = c; // premultiplied, so the bilinear blend with empty texels is exact
    output.depth = u_depth.SampleLevel(u_point, input.uv, 0).r;
    return output;
}
