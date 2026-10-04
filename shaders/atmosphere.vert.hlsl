// Atmosphere shell: a sphere around the planet, scaled to the top of the
// atmosphere. Positions are camera-relative (camera at origin).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_model;     // translation (camera-relative) * rotation * scale (planet radii)
    float4x4 u_view_proj;
    float4 u_top;         // x = radius of the shell top (planet radii)
};

struct VSInput
{
    float3 position : TEXCOORD0; // unit sphere
};

struct VSOutput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0; // camera-relative world position (km)
    float3 unit_pos : TEXCOORD1; // in the planet's unit-sphere frame
};

VSOutput main(VSInput input)
{
    VSOutput output;
    float3 unit_pos = input.position * u_top.x;
    float4 world = mul(u_model, float4(unit_pos, 1.0));
    output.position = mul(u_view_proj, world);
    output.rel_pos = world.xyz;
    output.unit_pos = unit_pos;
    return output;
}
