#include "common.hlsli"
#include "corona.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_corona, u_sampler, 0);

struct PSInput {
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    float3 ray = corona_ray(input.uv);
    float3 color;
    if (corona_edge(ray)) {
        color = u_color.rgb * corona_integrate(ray);
    } else {
        color = u_corona.SampleLevel(u_sampler, input.uv, 0).rgb;
    }
    return float4(color, 0.0);
}
