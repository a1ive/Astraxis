#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_star_depth, u_depth_sampler, 0);

#include "disk.hlsli"

struct PSInput {
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    return disk_integrate(disk_ray(input.uv));
}
