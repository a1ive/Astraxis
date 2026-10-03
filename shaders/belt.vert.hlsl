// Point cloud of small bodies on two-body heliocentric orbits (the main asteroid
// belt, the Kuiper belt). Each object expands to a screen-space square (6
// vertices); its position is propagated here from J2000 ecliptic elements by
// solving Kepler's equation (see belt_position() in src/scene/scene.cpp).

// 7 floats (28 bytes) per object: a (au), e, i, node, arg_peri, mean anomaly at
// the epoch (rad), absolute magnitude H.
ByteAddressBuffer u_objects : register(t0, space0);

cbuffer Uniforms : register(b0, space1)
{
    float4x4 u_view_proj;
    float4 u_display_x; // rows of the display-from-ecliptic rotation (xyz)
    float4 u_display_y;
    float4 u_display_z;
    float4 u_sun;       // xyz = sun relative to the camera (km)
    float4 u_params;    // x = days since the epoch, y = point size (px), z = sqrt(GM sun) (au^1.5 / day), w = au (km)
    float4 u_viewport;  // xy = size in pixels, z = reference H (weight 1)
};

struct VSOutput
{
    float4 position : SV_Position;
    float2 offset   : TEXCOORD0; // -1..1 across the point
    float weight    : TEXCOORD1; // relative brightness from H
};

static const float kTwoPi = 6.28318530718;

VSOutput main(uint vertex_id : SV_VertexID)
{
    static const float2 kCorners[6] = {
        float2(-1, -1), float2(1, -1), float2(1, 1),
        float2(-1, -1), float2(1, 1), float2(-1, 1),
    };
    uint object = vertex_id / 6;
    float2 corner = kCorners[vertex_id % 6];

    uint base = object * 28u;
    float a = asfloat(u_objects.Load(base));
    float e = asfloat(u_objects.Load(base + 4u));
    float inc = asfloat(u_objects.Load(base + 8u));
    float node = asfloat(u_objects.Load(base + 12u));
    float peri = asfloat(u_objects.Load(base + 16u));
    float m0 = asfloat(u_objects.Load(base + 20u));
    float h = asfloat(u_objects.Load(base + 24u));

    // Mean anomaly now, wrapped to [0, 2 pi).
    float n = u_params.z / (a * sqrt(a)); // rad/day
    float m = m0 + n * u_params.x;
    m -= kTwoPi * floor(m / kTwoPi);

    // Kepler's equation by Newton's method (start at pi for very eccentric orbits).
    float big_e = e < 0.8 ? m + e * sin(m) : 3.14159265;
    [unroll] for (int k = 0; k < 8; ++k) {
        big_e -= (big_e - e * sin(big_e) - m) / (1.0 - e * cos(big_e));
    }
    float px = a * (cos(big_e) - e);
    float py = a * sqrt(1.0 - e * e) * sin(big_e);

    // Perifocal -> ecliptic: Rz(node) Rx(inc) Rz(peri).
    float cw = cos(peri), sw = sin(peri), ci = cos(inc), si = sin(inc), cn = cos(node), sn = sin(node);
    float x1 = cw * px - sw * py;
    float y1 = sw * px + cw * py;
    float3 ecl = float3(cn * x1 - sn * ci * y1, sn * x1 + cn * ci * y1, si * y1) * u_params.w;

    float3 rel = u_sun.xyz + float3(dot(u_display_x.xyz, ecl), dot(u_display_y.xyz, ecl), dot(u_display_z.xyz, ecl));
    float4 clip = mul(u_view_proj, float4(rel, 1.0));

    VSOutput output;
    output.offset = corner;
    // Brighter (larger) objects count more, compressed: x2.5 per 5 magnitudes.
    output.weight = clamp(pow(10.0, -0.08 * (h - u_viewport.z)), 0.3, 4.0);
    if (clip.w <= 1e-4) {
        output.position = float4(0, 0, -1, 1); // behind the camera
        return output;
    }
    clip.xy += corner * u_params.y / u_viewport.xy * clip.w;
    output.position = clip;
    return output;
}
