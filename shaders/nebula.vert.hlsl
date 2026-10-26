// Patches of an expanding nebula's dusty shell (scene/nebula.hpp, eta
// Carinae's Homunculus), splatted as flat 3D Gaussians like a pinwheel's dust
// (pinwheel.vert.hlsl): round along the shell, thin across it, projected to a
// 2D Gaussian on the screen (EWA splatting), so the shell seen edge-on gets
// its column density (limb brightening). The whole shell is scaled from the
// star by the expansion since the reference time. Each patch scatters the
// star's light: 1 / d^2 from the star, the column thinning as 1 / scale^2,
// and the Henyey-Greenstein phase function, compressed for display. The
// fragment shader is dust.frag.hlsl.

#include "common.hlsli"

// SDL's D3D12 backend binds storage buffers as raw views. Each patch is 32
// bytes (NebulaGpuGrain): xyz = position at the reference time relative to the
// star (ICRF, km), w = weight; xyz = the shell's outward normal (unit, ICRF).
ByteAddressBuffer u_grains : register(t0, space0);

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4x4 u_view;    // rotation only: display -> camera
    float4 u_center;    // xyz = the star relative to the camera (km, display axes), w = scale since the reference time
    float4 u_display_x; // rows of the display-from-ICRF rotation (xyz)
    float4 u_display_y;
    float4 u_display_z;
    float4 u_shape;     // x = patch sigma (km), y = wall sigma (km), both at the reference time; z = min sigma (px), w = g
    float4 u_light;     // x = face-on radiance at the reference distance, y = display gamma, z = reference distance (km), w = cull radiance
    float4 u_color;     // rgb = color of the scattered light (sRGB)
    float4 u_screen;    // xy = viewport (px), z = pixels per radian
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 q        : TEXCOORD0; // offset from the center in sigmas
    float3 radiance : TEXCOORD1; // at the center
};

static const float kReach = 2.5;      // splat half-size in sigmas
static const float kMaxOffView = 1.5; // patches centered beyond this (NDC) are skipped
static const float kTwoPi = 6.2831853;

float3 to_display(float3 icrf)
{
    return float3(dot(u_display_x.xyz, icrf), dot(u_display_y.xyz, icrf), dot(u_display_z.xyz, icrf));
}

float3 to_camera(float3 display)
{
    return mul(u_view, float4(display, 0.0)).xyz;
}

VSOutput culled()
{
    VSOutput output;
    output.position = float4(0, 0, -1, 1);
    output.q = 0.0;
    output.radiance = 0.0;
    return output;
}

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

    float scale = u_center.w;
    float3 offset = to_display(a.xyz) * scale; // from the star
    float d = length(offset);
    float3 world = u_center.xyz + offset;       // relative to the camera, display axes
    float3 p = to_camera(world);
    float depth = -p.z;
    if (scale <= 0.0 || d <= 0.0 || depth <= 0.0) {
        return culled();
    }

    // Light scattered toward the camera: starlight falling off from the star,
    // spread over a column thinning as the shell grows, and the phase function
    // at the scattering angle (cos = 1 when the light goes straight on to the
    // camera), all compressed by the display gamma.
    float q = u_light.z / d;
    q = q * q / (scale * scale);
    float g = u_shape.w;
    float cos_theta = dot(offset / d, -normalize(world));
    float h = 1.0 + g * g - 2.0 * g * cos_theta;
    float phase = (1.0 - g * g) / (h * sqrt(h));
    float surface = u_light.x * pow(q * phase, u_light.y) * a.w;

    float3 normal = to_camera(to_display(b.xyz)); // across the shell
    float3 helper = abs(normal.x) < 0.9 ? float3(1, 0, 0) : float3(0, 1, 0);
    float3 along = normalize(cross(normal, helper));
    float3 across = cross(normal, along);

    // Screen covariance (px^2) of the patch: J S J^T with J the projection's
    // Jacobian at p, plus the minimum size.
    float f = u_screen.z;
    float s_patch = u_shape.x * scale;
    float s_wall = u_shape.y * scale;
    float2 ja = f / depth * (along.xy + along.z * p.xy / depth) * s_patch;
    float2 jc = f / depth * (across.xy + across.z * p.xy / depth) * s_patch;
    float2 jn = f / depth * (normal.xy + normal.z * p.xy / depth) * s_wall;
    float min2 = u_shape.z * u_shape.z;
    float cxx = ja.x * ja.x + jc.x * jc.x + jn.x * jn.x + min2;
    float cyy = ja.y * ja.y + jc.y * jc.y + jn.y * jn.y + min2;
    float cxy = ja.x * ja.y + jc.x * jc.y + jn.x * jn.y;
    float mid = 0.5 * (cxx + cyy);
    float rad = sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);
    float l1 = mid + rad;
    float l2 = max(mid - rad, 1e-4);
    float2 e1 = abs(cxy) > 1e-6 * mid ? normalize(float2(cxy, l1 - cxx)) : (cxx >= cyy ? float2(1, 0) : float2(0, 1));
    float2 e2 = float2(-e1.y, e1.x);

    // The patch's flux (its area, 2 sigma by 2 sigma along the shell, seen
    // face-on at this distance) spread over the screen Gaussian.
    float area = 4.0 * s_patch * s_patch;
    float flux = surface * area / (depth * depth);
    float peak = flux * f * f / (kTwoPi * sqrt(l1 * l2));

    float3 radiance = srgb_to_linear(u_color.rgb) * peak;
    if (max(radiance.r, max(radiance.g, radiance.b)) < u_light.w) {
        return culled();
    }

    VSOutput output;
    float4 clip = mul(u_view_proj, float4(world, 1.0));
    // The projection is linearized at the patch's center: far outside the
    // view (a patch beside a camera inside the shell) it no longer holds, and
    // the ellipse could reach into the picture where the patch does not.
    if (any(abs(clip.xy) > kMaxOffView * clip.w)) {
        return culled();
    }
    float2 corner_px = kReach * (corner.x * sqrt(l1) * e1 + corner.y * sqrt(l2) * e2);
    clip.xy += corner_px * 2.0 / u_screen.xy * clip.w;
    output.position = clip;
    output.q = corner * kReach;
    output.radiance = radiance;
    return output;
}
