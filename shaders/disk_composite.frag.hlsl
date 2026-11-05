#include "common.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_disk, u_sampler, 0);
FRAGMENT_TEXTURE(Texture2D, u_star_depth, u_depth_sampler, 1);

#include "disk.hlsli"

struct PSInput {
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    float3 ray = disk_ray(input.uv);
    if (disk_edge(ray)) {
        return disk_integrate(ray);
    }
    return u_disk.SampleLevel(u_sampler, input.uv, 0);
}
