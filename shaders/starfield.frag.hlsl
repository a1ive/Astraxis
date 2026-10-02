// Background stars: soft gaussian dot, additively blended.

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
    float3 color    : TEXCOORD1;
};

float4 main(PSInput input) : SV_Target0
{
    float r2 = dot(input.uv, input.uv);
    float a = exp(-r2 * 3.5);
    return float4(input.color * a, 1.0);
}
