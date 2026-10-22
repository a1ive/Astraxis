// Dust of a colliding-wind binary's spiral (scene/pinwheel.hpp). Each grain
// stands for a patch of the expanding shell it lies on: a flat 3D Gaussian,
// spread along its ring and toward the neighboring rings (in proportion to its
// distance r from the center, so that the patches tile the shell) and thin
// across the shell (the cone's wall). The 3D Gaussian is projected to a 2D one
// on the screen (EWA splatting: the projection linearized at the grain), so
// a wall seen edge-on gets its column density (limb brightening). The age of
// the grain, its distance and its light follow pinwheel_grain_age and
// pinwheel_grain_brightness on the CPU; the fragment shader is dust.frag.hlsl.

#include "common.hlsli"

// SDL's D3D12 backend binds storage buffers as raw views. Each grain is 32 bytes:
// xyz = direction (unit, ICRF), w = phase; xyz = tangent along its ring (unit, ICRF), w = turn.
ByteAddressBuffer u_grains : register(t0, space0);

#define LUT_SIZE 64

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4x4 u_view;      // rotation only: display -> camera
    float4 u_center;      // xyz = the pair's barycenter relative to the camera (km)
    float4 u_display_x;   // rows of the display-from-ICRF rotation (xyz)
    float4 u_display_y;
    float4 u_display_z;
    float4 u_cycle;       // x = frac((t - t_ref) / period / turns), y = turns, z = step (km per period), w = r_nuc (km)
    float4 u_shape;       // x, y = patch sigma along the ring and across the rings / r, z = wall sigma (km), w = min sigma (px)
    float4 u_light;       // x = face-on radiance at r_nuc, y = display gamma, z = 1 / ln(r_end / r_nuc), w = cull radiance
    float4 u_screen;      // xy = viewport (px), z = pixels per radian
    float4 u_colors[LUT_SIZE]; // at ln(r / r_nuc) / ln(r_end / r_nuc) = i / (LUT_SIZE - 1): linear RGB, brightness factor
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 q        : TEXCOORD0; // offset from the center in sigmas
    float3 radiance : TEXCOORD1; // at the center
};

static const float kReach = 2.5; // splat half-size in sigmas
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

    // Age (in periods) and distance from the center, as pinwheel_grain_age.
    float turns = u_cycle.y;
    float age = turns * frac(u_cycle.x - (a.w + b.w) / turns);
    float r = age * u_cycle.z;
    float r_nuc = u_cycle.w;
    if (r <= r_nuc) {
        return culled();
    }
    // Faded in beyond r_nuc and out over the last turn (pinwheel_grain_brightness).
    float fade = min((r - r_nuc) / (0.25 * r_nuc), 1.0) * saturate(turns - age);
    // A face-on wall's surface brightness falls as the light taken up per
    // grain (r^-2) times the grains per area (r^-2), compressed for display.
    float q2 = r_nuc / r;
    q2 *= q2;
    float u = saturate(log(r / r_nuc) * u_light.z) * (LUT_SIZE - 1);
    uint i = min(uint(u), LUT_SIZE - 2);
    float4 color = lerp(u_colors[i], u_colors[i + 1], u - float(i));
    float surface = u_light.x * pow(q2 * q2, u_light.y) * color.a * fade;

    float3 world = u_center.xyz + to_display(a.xyz) * r; // relative to the camera, display axes
    float3 p = to_camera(world);
    float3 normal = to_camera(to_display(a.xyz)); // across the shell
    float3 along = to_camera(to_display(b.xyz));  // along the ring
    float3 across = cross(normal, along);
    float depth = -p.z;
    if (depth <= 0.0) {
        return culled();
    }

    // Screen covariance (px^2) of the patch: J S J^T with J the projection's
    // Jacobian at p, plus the minimum size.
    float f = u_screen.z;
    float s_along = u_shape.x * r;
    float s_across = u_shape.y * r;
    float s_wall = u_shape.z;
    float2 ja = f / depth * (along.xy + along.z * p.xy / depth) * s_along;
    float2 jc = f / depth * (across.xy + across.z * p.xy / depth) * s_across;
    float2 jn = f / depth * (normal.xy + normal.z * p.xy / depth) * s_wall;
    float min2 = u_shape.w * u_shape.w;
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
    float area = 4.0 * s_along * s_across;
    float flux = surface * area / (depth * depth);
    float peak = flux * f * f / (kTwoPi * sqrt(l1 * l2));

    float3 radiance = color.rgb * peak;
    if (max(radiance.r, max(radiance.g, radiance.b)) < u_light.w) {
        return culled();
    }

    VSOutput output;
    float4 clip = mul(u_view_proj, float4(world, 1.0));
    float2 offset = kReach * (corner.x * sqrt(l1) * e1 + corner.y * sqrt(l2) * e2); // px
    clip.xy += offset * 2.0 / u_screen.xy * clip.w;
    output.position = clip;
    output.q = corner * kReach;
    output.radiance = radiance;
    return output;
}
