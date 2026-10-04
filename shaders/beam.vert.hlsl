// Pulsar beam: a camera-facing quad along the beam axis whose width grows
// with distance from the star (a cone seen from the side). Positions are
// camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4 u_apex;   // xyz = camera-relative star position (km), w = beam length (km)
    float4 u_axis;   // xyz = unit beam direction, w = tan(half-angle)
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 coord    : TEXCOORD0; // x = across the beam (-1..1), y = along it (0..1)
};

static const float2 kCorners[6] = {
    float2(-1, 0), float2(1, 0), float2(1, 1),
    float2(-1, 0), float2(1, 1), float2(-1, 1),
};

VSOutput main(uint vertex_id : SV_VertexID)
{
    VSOutput output;
    float2 c = kCorners[vertex_id];
    float3 axis = u_axis.xyz;
    float along = c.y * u_apex.w;
    float3 center = u_apex.xyz + axis * along;
    // Across the beam, perpendicular to the axis and to the line of sight.
    float3 side = cross(axis, center);
    float len = length(side);
    side = len > 1e-6 ? side / len : float3(0, 0, 0);
    float3 p = center + side * (c.x * along * u_axis.w);
    output.position = mul(u_view_proj, float4(p, 1.0));
    output.coord = c;
    return output;
}
