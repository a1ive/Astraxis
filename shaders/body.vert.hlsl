// Planet / moon ellipsoid, or an irregular body's mesh. Positions are
// camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // translation (camera-relative) * rotation * scale
    float4x4 u_view_proj;
    float4x4 u_rotation;  // body-fixed axes in world space
    float4 u_inv_scale;   // xyz = 1 / scale, for normals; w = map u of longitude 0 (turns)
};

static const float kPi = 3.14159265;

struct VSInput
{
    float3 position : TEXCOORD0; // unit sphere (ellipsoids) or km (meshes), body-fixed
    float3 normal   : TEXCOORD1; // before scaling
    float2 uv       : TEXCOORD2; // x = east longitude / 360 deg (y: colatitude / 180 deg, unused)
    float  albedo   : TEXCOORD3; // relative to the body color
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0; // camera-relative world position (km)
    float3 normal   : TEXCOORD1; // world-space normal
    float3 local    : TEXCOORD2; // body-fixed position before scaling
    float2 uv       : TEXCOORD3;
    float3 axis_x   : TEXCOORD4; // body-fixed axes in world space, for normals
    float3 axis_y   : TEXCOORD5; // perturbed in the body frame
    float3 axis_z   : TEXCOORD6;
    float  albedo   : TEXCOORD7;
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float4 world = mul(u_model, float4(input.position, 1.0));
    output.position = mul(u_view_proj, world);
    output.rel_pos = world.xyz;
    output.normal = mul(u_rotation, float4(input.normal * u_inv_scale.xyz, 0.0)).xyz;
    output.local = input.position;

    // Maps are in planetocentric longitude and latitude: the direction of the
    // scaled point from the center, not the unit-sphere parameters (they differ
    // on a triaxial body by up to 7 deg for Vesta). The longitude correction is
    // small and wrapped, so u stays continuous across the seam column. A mesh is
    // unscaled, so its map u (that of its vertex direction) passes through.
    float3 q = input.position / u_inv_scale.xyz;
    float dlon = atan2(q.y, q.x) - atan2(input.position.y, input.position.x);
    dlon -= 2.0 * kPi * round(dlon / (2.0 * kPi));
    output.uv = float2(input.uv.x + dlon / (2.0 * kPi) + u_inv_scale.w,
                       0.5 - asin(clamp(q.z / max(length(q), 1e-20), -1.0, 1.0)) / kPi);

    output.axis_x = mul(u_rotation, float4(1.0, 0.0, 0.0, 0.0)).xyz;
    output.axis_y = mul(u_rotation, float4(0.0, 1.0, 0.0, 0.0)).xyz;
    output.axis_z = mul(u_rotation, float4(0.0, 0.0, 1.0, 0.0)).xyz;
    output.albedo = input.albedo;
    return output;
}
