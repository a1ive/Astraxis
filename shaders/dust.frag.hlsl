// Dust splat: a Gaussian, added to what lies behind.

struct PSInput
{
    float4 position : SV_Position;
    float2 q        : TEXCOORD0;
    float3 radiance : TEXCOORD1;
};

float4 main(PSInput input) : SV_Target0
{
    return float4(input.radiance * exp(-0.5 * dot(input.q, input.q)), 1.0);
}
