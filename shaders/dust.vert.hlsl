// Dust grains of comet tails as Gaussian splats facing the camera. Each splat
// is a normalized 2D Gaussian of standard deviation sigma (km, at the grain;
// at least u_params.y pixels), so its total flux is its luminosity / distance^2,
// in the same units as the comet pass's columns.

#include "common.hlsli"

// SDL's D3D12 backend binds storage buffers as raw views. Each grain is 32 bytes:
// xyz = camera-relative position (km), w = sigma (km); rgb = color (sRGB), w = luminosity.
ByteAddressBuffer u_grains : register(t0, space0);

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4 u_params; // x = pixels per radian, y = minimum sigma (pixels), zw = viewport size (pixels)
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 q        : TEXCOORD0; // offset from the center in sigmas
    float3 radiance : TEXCOORD1; // at the center
};

static const float kReach = 3.0; // splat half-size in sigmas

VSOutput main(uint vertex_id : SV_VertexID)
{
    static const float2 kCorners[6] = {
        float2(-1, -1), float2(1, -1), float2(1, 1),
        float2(-1, -1), float2(1, 1), float2(-1, 1),
    };
    uint grain = vertex_id / 6;
    float2 corner = kCorners[vertex_id % 6];
    float4 a = asfloat(u_grains.Load4(grain * 32u));
    float4 b = asfloat(u_grains.Load4(grain * 32u + 16u));

    VSOutput output;
    float4 clip = mul(u_view_proj, float4(a.xyz, 1.0));
    float distance = length(a.xyz);
    if (clip.w <= 0.0 || distance <= 0.0) {
        output.position = float4(0, 0, -1, 1); // behind the camera
        output.q = 0.0;
        output.radiance = 0.0;
        return output;
    }
    float sigma = max(a.w, u_params.y * distance / u_params.x);
    float sigma_px = sigma / distance * u_params.x;
    clip.xy += corner * kReach * sigma_px * 2.0 / u_params.zw * clip.w;
    output.position = clip;
    output.q = corner * kReach;
    output.radiance = srgb_to_linear(b.rgb) * b.w / (6.2831853 * sigma * sigma);
    return output;
}
