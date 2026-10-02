// The Sun as a camera-facing sprite at infinity (clip z = 0, reversed-Z far).
// Occluded by bodies through the depth test.

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj; // rotation-only view * projection
    float4 u_direction;   // xyz = camera-relative position (w = 1) or unit direction at infinity (w = 0)
    float4 u_sprite;      // x = sprite half-size in pixels, zw = viewport size
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 offset_px : TEXCOORD0; // pixel offset from the sun center
};

static const float2 kCorners[6] = {
    float2(-1, -1), float2(1, -1), float2(1, 1),
    float2(-1, -1), float2(1, 1), float2(-1, 1),
};

VSOutput main(uint vertex_id : SV_VertexID)
{
    VSOutput output;
    float2 corner = kCorners[vertex_id];

    float4 clip = mul(u_view_proj, float4(u_direction.xyz, u_direction.w));
    float2 offset_px = corner * u_sprite.x;
    clip.xy += offset_px * 2.0 / u_sprite.zw * clip.w;
    if (u_direction.w < 0.5) {
        clip.z = 0.0; // at infinity: the reversed-Z far plane
    }
    if (clip.w <= 0.0) {
        clip = float4(0, 0, -1, 1);
    }

    output.position = clip;
    output.offset_px = offset_px;
    return output;
}
