#include "corona.hlsli"

struct PSInput {
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    return float4(u_color.rgb * corona_integrate(corona_ray(input.uv)), 0.0);
}
