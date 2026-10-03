// Planetary ring: a flat annulus in the planet's equatorial plane.
// Positions are camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // camera-relative translation * body-fixed axes (km)
    float4x4 u_view_proj;
    float4 u_radii;       // x = inner, y = outer radius of the mesh (km)
};

struct VSInput
{
    float2 direction : TEXCOORD0; // unit vector in the equatorial plane
    float edge       : TEXCOORD1; // 0 = inner edge, 1 = outer edge
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 local    : TEXCOORD0; // body-fixed position (km)
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float radius = lerp(u_radii.x, u_radii.y, input.edge);
    float3 local = float3(input.direction * radius, 0.0);
    output.position = mul(u_view_proj, mul(u_model, float4(local, 1.0)));
    output.local = local;
    return output;
}
