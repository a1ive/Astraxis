// Habitable zone: a flat annulus in its plane, in au around the zone's center.
// Positions are camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // camera-relative translation * plane axes (km per au)
    float4x4 u_view_proj;
    float4 u_radii;       // x = inner, y = outer radius of the annulus (au)
};

struct VSInput
{
    float2 direction : TEXCOORD0; // unit vector in the plane
    float edge       : TEXCOORD1; // 0 = inner edge, 1 = outer edge
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 local    : TEXCOORD0; // in the plane frame (au)
    float3 view     : TEXCOORD1; // camera-relative (km)
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float radius = lerp(u_radii.x, u_radii.y, input.edge);
    float3 local = float3(input.direction * radius, 0.0);
    float4 world = mul(u_model, float4(local, 1.0));
    output.position = mul(u_view_proj, world);
    output.local = local;
    output.view = world.xyz;
    return output;
}
