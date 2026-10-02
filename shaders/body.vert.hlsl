// Planet / moon ellipsoid. Positions are camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // translation (camera-relative) * rotation * scale
    float4x4 u_view_proj;
    float4x4 u_rotation;  // body-fixed axes in world space
    float4 u_inv_scale;   // xyz = 1 / radii, for ellipsoid normals
};

struct VSInput
{
    float3 position : TEXCOORD0; // unit sphere
    float2 uv       : TEXCOORD1; // equirectangular (u east from mesh longitude 0, v down from north)
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0; // camera-relative world position (km)
    float3 normal   : TEXCOORD1; // world-space normal
    float3 local    : TEXCOORD2; // body-fixed unit-sphere position
    float2 uv       : TEXCOORD3;
    float3 axis_x   : TEXCOORD4; // body-fixed axes in world space, for normals
    float3 axis_y   : TEXCOORD5; // perturbed in the body frame
    float3 axis_z   : TEXCOORD6;
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float4 world = mul(u_model, float4(input.position, 1.0));
    output.position = mul(u_view_proj, world);
    output.rel_pos = world.xyz;
    output.normal = mul(u_rotation, float4(input.position * u_inv_scale.xyz, 0.0)).xyz;
    output.local = input.position;
    output.uv = input.uv;
    output.axis_x = mul(u_rotation, float4(1.0, 0.0, 0.0, 0.0)).xyz;
    output.axis_y = mul(u_rotation, float4(0.0, 1.0, 0.0, 0.0)).xyz;
    output.axis_z = mul(u_rotation, float4(0.0, 0.0, 1.0, 0.0)).xyz;
    return output;
}
