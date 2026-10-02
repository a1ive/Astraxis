// Background stars: one camera-facing quad per instance, sized in pixels and
// placed at infinity (clip z = 0, the reversed-Z far plane).

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj; // rotation-only view * projection
    float4 u_viewport;    // xy = size in pixels
    float4 u_params;      // x = brightness scale, y = size scale
};

struct VSInput
{
    float4 dir_size        : TEXCOORD0; // xyz = unit direction (ICRF), w = size in pixels
    float4 color_intensity : TEXCOORD1;
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
    float3 color    : TEXCOORD1;
};

static const float2 kCorners[6] = {
    float2(-1, -1), float2(1, -1), float2(1, 1),
    float2(-1, -1), float2(1, 1), float2(-1, 1),
};

VSOutput main(VSInput input, uint vertex_id : SV_VertexID)
{
    VSOutput output;
    float2 corner = kCorners[vertex_id];

    float4 clip = mul(u_view_proj, float4(input.dir_size.xyz, 0.0));
    float size = input.dir_size.w * u_params.y;
    clip.xy += corner * size / u_viewport.xy * clip.w;
    clip.z = 0.0;
    if (clip.w <= 0.0) {
        clip = float4(0, 0, -1, 1); // behind the camera: clipped
    }

    output.position = clip;
    output.uv = corner;
    output.color = input.color_intensity.rgb * input.color_intensity.a * u_params.x;
    return output;
}
