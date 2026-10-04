// Plume bounding box (a unit cube scaled to the box in the plume's local
// frame). Positions are camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // local frame (km) -> camera-relative world: translation * rotation
    float4x4 u_view_proj;
    float4 u_box_min;     // xyz: box corners in the local frame (km)
    float4 u_box_max;
};

struct VSInput
{
    float3 position : TEXCOORD0; // unit cube corner (0 or 1)
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0; // camera-relative world position (km)
    float3 local    : TEXCOORD1; // in the local frame (km)
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float3 local = lerp(u_box_min.xyz, u_box_max.xyz, input.position);
    float4 world = mul(u_model, float4(local, 1.0));
    output.position = mul(u_view_proj, world);
    output.rel_pos = world.xyz;
    output.local = local;
    return output;
}
