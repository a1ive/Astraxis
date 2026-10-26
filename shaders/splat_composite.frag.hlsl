// A nebula's splats (a pinwheel's dust, the Homunculus), drawn off-screen at
// reduced resolution (render/splat_pass.hpp), added to the scene: the nebula
// is smooth, and thousands of overlapping splats blended into every MSAA
// sample cost far too much fill rate.

#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D<float4>, u_dust, u_dust_sampler, 0);

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    return float4(u_dust.Sample(u_dust_sampler, input.uv).rgb, 0.0);
}
