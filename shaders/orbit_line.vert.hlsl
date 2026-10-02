// Screen-space-width polyline. Points come from a storage buffer; each segment
// expands to a quad (6 vertices), clipped against the camera plane.

// SDL's D3D12 backend binds storage buffers as raw views.
// Each point is 16 bytes: xyz = camera-relative position (km), w = fade (0 = head, 1 = tail).
ByteAddressBuffer u_points : register(t0, space0);

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4 u_viewport; // xy = size in pixels
    float4 u_params;   // x = first point index, y = width in pixels
};

struct VSOutput
{
    float4 position : SV_Position;
    float fade      : TEXCOORD0;
    float side      : TEXCOORD1;
};

static const float kMinW = 1e-4;

VSOutput main(uint vertex_id : SV_VertexID)
{
    // Corner table: (end 0/1, side -1/+1).
    static const float2 kCorners[6] = {
        float2(0, -1), float2(0, 1), float2(1, 1),
        float2(0, -1), float2(1, 1), float2(1, -1),
    };

    uint segment = vertex_id / 6;
    float2 corner = kCorners[vertex_id % 6];
    uint first = (uint)u_params.x;

    float4 a = asfloat(u_points.Load4((first + segment) * 16u));
    float4 b = asfloat(u_points.Load4((first + segment + 1u) * 16u));

    float4 c0 = mul(u_view_proj, float4(a.xyz, 1.0));
    float4 c1 = mul(u_view_proj, float4(b.xyz, 1.0));
    float f0 = a.w;
    float f1 = b.w;

    VSOutput output;
    if (c0.w < kMinW && c1.w < kMinW) {
        output.position = float4(0, 0, -1, 1); // fully behind the camera
        output.fade = 1.0;
        output.side = 0.0;
        return output;
    }
    // Clip the part behind the camera.
    if (c0.w < kMinW) {
        float t = (kMinW - c0.w) / (c1.w - c0.w);
        c0 = lerp(c0, c1, t);
        f0 = lerp(f0, f1, t);
    } else if (c1.w < kMinW) {
        float t = (kMinW - c1.w) / (c0.w - c1.w);
        c1 = lerp(c1, c0, t);
        f1 = lerp(f1, f0, t);
    }

    float2 s0 = c0.xy / c0.w * u_viewport.xy;
    float2 s1 = c1.xy / c1.w * u_viewport.xy;
    float2 dir = s1 - s0;
    float len = length(dir);
    dir = len > 1e-6 ? dir / len : float2(1, 0);
    float2 normal = float2(-dir.y, dir.x);

    float4 c = corner.x < 0.5 ? c0 : c1;
    // Half width in pixels -> NDC (NDC spans 2 units per viewport).
    float half_width = u_params.y * 0.5 + 0.5;
    c.xy += normal * corner.y * half_width * 2.0 / u_viewport.xy * c.w;

    output.position = c;
    output.fade = corner.x < 0.5 ? f0 : f1;
    output.side = corner.y * half_width; // signed distance in pixels
    return output;
}
