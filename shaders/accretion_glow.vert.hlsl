// The glowing gas of an accretion system (the stream from L1, the rim's hot
// spot and the jets) as Gaussian splats facing the camera, drawn off-screen at
// reduced resolution (render/splat_pass.hpp) and added to the scene. The
// target has no depth, so the stars of the system hide the splats here: each
// is an ellipsoid in the disk's frame (the gainer, the donor's Roche lobe), and
// a splat behind one fades out over its own size near the limb. The fragment
// shader is dust.frag.hlsl.

#include "common.hlsli"

// SDL's D3D12 backend binds storage buffers as raw views. Each splat is 32 bytes
// (DustSplat): xyz = camera-relative position (km), w = sigma (km); rgb = color
// (sRGB), w = luminosity.
ByteAddressBuffer u_grains : register(t0, space0);

#define MAX_OCCLUDERS 2

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4 u_screen;      // xy = target size (px), z = pixels per radian, w = min sigma (px)
    float4 u_to_local[3]; // rows of the world -> disk rotation (xyz); u_to_local[0].w = first splat of this draw
    float4 u_center;      // xyz = the gainer relative to the camera (km), w = occluder count
    float4 u_occluder_centers[MAX_OCCLUDERS]; // xyz = center in the disk frame (km)
    float4 u_occluder_radii[MAX_OCCLUDERS];   // xyz = semi-axes along the disk's axes (km)
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 q        : TEXCOORD0; // offset from the center in sigmas
    float3 radiance : TEXCOORD1; // at the center
};

static const float kReach = 3.0; // splat half-size in sigmas

float3 to_local(float3 v)
{
    return float3(dot(u_to_local[0].xyz, v), dot(u_to_local[1].xyz, v), dot(u_to_local[2].xyz, v));
}

VSOutput culled()
{
    VSOutput output;
    output.position = float4(0, 0, -1, 1);
    output.q = 0.0;
    output.radiance = 0.0;
    return output;
}

// How much of a splat at `p` (disk frame) the ellipsoid leaves visible from `eye`.
float visible_past(float3 eye, float3 p, float sigma, float3 center, float3 radii)
{
    const float3 o = (eye - center) / radii;
    const float3 q = (p - center) / radii;
    if (dot(o, o) <= 1.0) {
        return 1.0; // the camera inside the star: nothing to hide behind
    }
    const float w = clamp(sigma / min(radii.x, min(radii.y, radii.z)), 0.02, 0.5);
    const float3 d = q - o;
    const float t = saturate(-dot(o, d) / dot(d, d));
    const float closest = length(o + t * d);
    const float behind = (1.0 - t) * length(d); // from the closest approach to the splat
    const float covered = (1.0 - smoothstep(1.0 - w, 1.0 + w, closest)) * smoothstep(0.0, w, behind);
    return smoothstep(1.0 - w, 1.0, length(q)) * (1.0 - covered);
}

VSOutput main(uint vertex_id : SV_VertexID)
{
    static const float2 kCorners[6] = {
        float2(-1, -1), float2(1, -1), float2(1, 1),
        float2(-1, -1), float2(1, 1), float2(-1, 1),
    };
    uint grain = uint(u_to_local[0].w) + vertex_id / 6;
    float2 corner = kCorners[vertex_id % 6];
    float4 a = asfloat(u_grains.Load4(grain * 32u));
    float4 b = asfloat(u_grains.Load4(grain * 32u + 16u));

    float4 clip = mul(u_view_proj, float4(a.xyz, 1.0));
    float distance = length(a.xyz);
    if (clip.w <= 0.0 || distance <= 0.0) {
        return culled();
    }
    float sigma = max(a.w, u_screen.w * distance / u_screen.z);

    const float3 eye = to_local(-u_center.xyz);
    const float3 p = to_local(a.xyz - u_center.xyz);
    float visible = 1.0;
    for (int k = 0; k < (int)u_center.w && k < MAX_OCCLUDERS; ++k) {
        visible *= visible_past(eye, p, sigma, u_occluder_centers[k].xyz, u_occluder_radii[k].xyz);
    }
    if (visible <= 1e-3) {
        return culled();
    }

    VSOutput output;
    float sigma_px = sigma / distance * u_screen.z;
    clip.xy += corner * kReach * sigma_px * 2.0 / u_screen.xy * clip.w;
    output.position = clip;
    output.q = corner * kReach;
    output.radiance = visible * srgb_to_linear(b.rgb) * b.w / (6.2831853 * sigma * sigma);
    return output;
}
