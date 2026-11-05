#include "common.hlsli"
#include "disk.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_disk, u_sampler, 0);

struct PSInput {
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    float3 ray = disk_ray(input.uv);
    float3 color;
    if (disk_edge(ray)) {
        color = disk_integrate(ray);
    } else {
        color = u_disk.SampleLevel(u_sampler, input.uv, 0).rgb;
    }
    return float4(color, 0.0);
}
