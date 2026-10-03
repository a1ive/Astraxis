// Planet / moon shading: Lambert lighting, soft eclipse shadows from other
// bodies (sun treated as a disk), the shadow of the planet's own rings,
// procedural gas-giant bands (Jupiter-like, Saturn-like) and a procedural
// battle station.

#include "common.hlsli"

#define MAX_OCCLUDERS 8
#define STYLE_GAS_GIANT 1
#define STYLE_STAR 2       // self-luminous: u_color is linear HDR radiance
#define STYLE_BLACK_HOLE 3 // event horizon: absorbs everything
#define STYLE_DEATH_STAR 4
#define STYLE_SATURN 5

static const float kPi = 3.14159265;

FRAGMENT_TEXTURE(Texture2D, u_albedo, u_sampler, 0); // sRGB texture: samples are linear
FRAGMENT_TEXTURE(Texture2D, u_ring_profile, u_ring_sampler, 1); // r = ring normal optical depth

cbuffer Uniforms : register(b0, space3)
{
    float4 u_sun;                      // xyz = unit direction to sun, w = sun angular radius (rad)
    float4 u_color;                    // rgb = base color (sRGB), a = style
    float4 u_params;                   // x = occluder count, y = ambient, z = 1 if textured, w = 1 to flip u
    float4 u_occluders[MAX_OCCLUDERS]; // xyz = camera-relative center, w = radius
    float4 u_ring_center;              // xyz = camera-relative center of the ring plane, w = 1 if it has rings
    float4 u_ring_normal;              // xyz = ring plane normal (unit)
    float4 u_ring_radii;               // x = inner, y = outer radius of the profile (km), z = profile samples
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 normal   : TEXCOORD1;
    float3 local    : TEXCOORD2;
    float2 uv       : TEXCOORD3;
    float3 axis_x   : TEXCOORD4;
    float3 axis_y   : TEXCOORD5;
    float3 axis_z   : TEXCOORD6;
    float  albedo   : TEXCOORD7; // relative (shape models; 1 on ellipsoids)
};

float eclipse_factor(float3 p, float3 L)
{
    float lit = 1.0;
    int count = (int)u_params.x;
    for (int k = 0; k < count; ++k) {
        float3 oc = u_occluders[k].xyz - p;
        float along = dot(oc, L);
        if (along <= 0.0) {
            continue;
        }
        float radius = u_occluders[k].w;
        float perp = length(oc - along * L);
        float penumbra = max(along * u_sun.w, 1e-3);
        // Fraction of the sun disk still visible, roughly.
        float s = smoothstep(radius - penumbra, radius + penumbra, perp);
        // Antumbra: an occluder smaller than the sun's disk never fully blocks it.
        float max_block = saturate((radius * radius) / (penumbra * penumbra));
        lit *= lerp(1.0 - max_block, 1.0, s);
    }
    return lit;
}

float3 gas_giant_color(float3 local)
{
    float lat = asin(clamp(local.z, -1.0, 1.0));
    float lon = atan2(local.y, local.x);

    float w = lat + 0.012 * sin(lon * 6.0 + lat * 20.0) + 0.006 * sin(lon * 13.0 - lat * 35.0);
    float b1 = 0.5 + 0.5 * sin(w * 22.0);
    float b2 = 0.5 + 0.5 * sin(w * 57.0 + 1.7);
    float t = saturate(b1 * 0.65 + b2 * 0.35);

    float3 light = float3(0.93, 0.88, 0.80);
    float3 dark = float3(0.70, 0.53, 0.40);
    float3 c = lerp(dark, light, t);

    // Bright equatorial zone, grayish polar regions.
    c = lerp(c, light, 0.5 * exp(-(lat / 0.12) * (lat / 0.12)));
    c = lerp(c, float3(0.62, 0.58, 0.55), smoothstep(0.9, 1.25, abs(lat)));

    // Cosmetic Great Red Spot (position is not tracked).
    float2 d = float2((lon - 0.9) / 0.20, (lat + 0.39) / 0.09);
    float spot = exp(-dot(d, d) * 1.5);
    c = lerp(c, float3(0.80, 0.46, 0.33), spot * 0.85);
    return c;
}

// Saturn: pale, low-contrast belts and zones around a brighter equatorial zone,
// and the hexagonal jet around the north pole (~78 deg N) with a darker,
// bluish-gray interior. Illustrative, tinted by `base` (linear).
float3 saturn_color(float3 local, float3 base)
{
    float lat = asin(clamp(local.z, -1.0, 1.0));
    float lon = atan2(local.y, local.x);

    float w = lat + 0.004 * sin(lon * 5.0 + lat * 30.0);
    float shade = 1.0 + 0.09 * sin(w * 26.0) + 0.05 * sin(w * 71.0 + 0.8) + 0.03 * sin(w * 150.0 + 2.1);
    float3 c = base * shade;
    c = lerp(c, base * 1.1, 0.6 * exp(-(lat / 0.14) * (lat / 0.14)));

    // Hexagon: the edge at 12 deg from the pole mid-side (13.9 deg at the corners).
    float colat = 1.5707963 - lat;
    float sector = fmod(lon + 6.2831853, 1.0471976) - 0.5235988;
    float edge = 0.2094395 / cos(sector);
    float inside = 1.0 - smoothstep(edge - 0.008, edge + 0.008, colat);
    c = lerp(c, base * float3(0.70, 0.76, 0.80), 0.8 * inside);
    c *= 1.0 - 0.15 * exp(-((colat - edge) / 0.01) * ((colat - edge) / 0.01));

    c = lerp(c, base * 0.85, smoothstep(1.2, 1.45, -lat)); // grayer south polar region
    return c;
}

// Where the ray from p toward the sun crosses the ring plane: the position in
// the ring profile (0..1 across [inner, outer]; < 0 if the plane is behind).
float ring_profile_u(float3 p, float3 L)
{
    float3 n = u_ring_normal.xyz;
    float3 d = p - u_ring_center.xyz;
    float denom = dot(L, n);
    float t = -dot(d, n) / (abs(denom) > 1e-6 ? denom : 1e-6);
    float r = length(d + t * L);
    return t > 0.0 ? (r - u_ring_radii.x) / (u_ring_radii.y - u_ring_radii.x) : -1.0;
}

float hash21(float2 q)
{
    float3 p3 = frac(float3(q.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.x + p3.y) * p3.z);
}

struct DeathStarSurface
{
    float3 albedo;
    float3 normal;    // body frame
    float3 glow;      // always emitted (superlaser focus)
    float3 lights;    // emitted, visible on the night side (windows)
};

// Hull panels in latitude bands at three scales, an equatorial trench and a
// concave superlaser dish at 28 deg N on the prime meridian. Details smaller
// than a pixel fade to their mean to avoid shimmer.
DeathStarSurface death_star(float3 p)
{
    DeathStarSurface s;
    s.normal = p;
    s.glow = 0.0;
    s.lights = 0.0;

    float lat = asin(clamp(p.z, -1.0, 1.0));
    float lon = atan2(p.y, p.x);

    float shade = 1.0;
    float seam = 0.0;
    float window = 0.0;
    const float kBands[3] = {14.0, 56.0, 224.0}; // bands per radian of latitude
    for (int k = 0; k < 3; ++k) {
        float y = lat * kBands[k];
        float w = fwidth(y); // pixel footprint in cells; also used for x (cells are ~square)
        float detail = saturate(1.5 - 2.0 * w);
        float row = floor(y);
        float cells = max(1.0, floor(cos((row + 0.5) / kBands[k]) * kBands[k] * 2.0 * kPi));
        float x = (lon / (2.0 * kPi) + 0.5) * cells;
        float2 cell = float2(floor(x), row + 997.0 * k);
        float h = hash21(cell);
        const float kContrast[3] = {0.05, 0.10, 0.16}; // subtle large sectors, busier fine plating
        shade *= lerp(1.0, 1.0 + kContrast[k] * (2.0 * h - 1.0), detail);
        float2 f = abs(frac(float2(x, y)) - 0.5);
        // The coarse level only draws latitude lines (the hull bands); finer levels draw panel edges.
        float edge = smoothstep(0.5 - max(w, 0.02), 0.5, k == 0 ? f.y : max(f.x, f.y));
        const float kSeam[3] = {0.6, 0.3, 0.2};
        seam = max(seam, edge * detail * kSeam[k]);
        if (k == 2) {
            // Lit windows on 3.5% of the fine panels; their mean glow when too small to resolve.
            window = lerp(0.035, step(0.965, hash21(cell + 31.7)), detail);
        }
    }
    s.albedo = 0.42 * shade * (1.0 - 0.45 * seam);
    s.lights = float3(1.0, 0.85, 0.6) * window * 0.25 * shade * shade * shade; // uneven by sector

    // Equatorial trench (dark, with lit rims).
    float trench_w = fwidth(lat);
    float trench = 1.0 - smoothstep(0.010, 0.010 + trench_w, abs(lat));
    float rims = exp(-pow((abs(lat) - 0.012) / max(0.002, trench_w), 2.0));
    s.albedo = lerp(s.albedo, 0.06, trench) + 0.15 * rims;
    s.lights *= 1.0 - trench;

    // Superlaser dish: a spherical-cap depression (angular radius 0.26 rad).
    const float kDishLat = 0.49;
    const float kDishRadius = 0.26;
    float3 c = float3(cos(kDishLat), 0.0, sin(kDishLat));
    float t = acos(clamp(dot(p, c), -1.0, 1.0)) / kDishRadius;
    float tw = fwidth(t); // outside the branch: derivatives need uniform control flow
    if (t < 1.05) {
        float inside = 1.0 - smoothstep(1.0 - tw, 1.0, t);
        float3 to_axis = c * dot(p, c) - p;
        float len = length(to_axis);
        if (len > 1e-6) {
            // Concave: normals tilt toward the dish axis, more toward the rim.
            float3 dish_n = normalize(p + to_axis / len * (0.9 * t));
            s.normal = normalize(lerp(p, dish_n, inside));
        }
        float rings = lerp(1.0, 0.85 + 0.15 * sin(t * 90.0), saturate(1.5 - 60.0 * tw));
        s.albedo = lerp(s.albedo, 0.30 * rings, inside);
        s.albedo += 0.2 * exp(-pow((t - 1.0) / max(0.02, tw), 2.0)); // raised rim
        s.lights *= 1.0 - inside;
        // Focus lens of the eight tributary beams.
        s.glow = float3(0.25, 1.0, 0.4) * 1.5 * exp(-pow(t / 0.05, 2.0)) * inside;
    }
    return s;
}

float4 main(PSInput input) : SV_Target0
{
    float3 n = normalize(input.normal);
    float3 L = u_sun.xyz;
    float3 V = normalize(-input.rel_pos);

    int style = (int)u_color.a;
    if (style == STYLE_BLACK_HOLE) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }
    if (style == STYLE_DEATH_STAR) {
        DeathStarSurface s = death_star(normalize(input.local));
        float3 ns = normalize(input.axis_x * s.normal.x + input.axis_y * s.normal.y + input.axis_z * s.normal.z);
        float ndl_s = dot(ns, L);
        float lit_s = saturate(ndl_s) * smoothstep(-0.05, 0.05, dot(n, L)) * eclipse_factor(input.rel_pos, L);
        float night = 1.0 - smoothstep(-0.15, 0.1, dot(n, L));
        float3 albedo_s = srgb_to_linear(u_color.rgb) / 0.42 * s.albedo; // u_color tints the hull
        return float4(albedo_s * (lit_s + u_params.y) + s.lights * night + s.glow, 1.0);
    }
    if (style == STYLE_STAR) {
        // Linear limb darkening, as for the sun sprite.
        float mu = saturate(dot(n, V));
        return float4(u_color.rgb * (1.0 - 0.6 * (1.0 - mu)), 1.0);
    }

    bool gas_giant = (int)u_color.a == STYLE_GAS_GIANT;
    bool saturn = (int)u_color.a == STYLE_SATURN;
    float3 albedo;
    if (u_params.z > 0.5) {
        float2 uv = input.uv;
        if (u_params.w > 0.5) {
            uv.x = 1.0 - uv.x;
        }
        albedo = u_albedo.Sample(u_sampler, uv).rgb;
    } else if (saturn) {
        albedo = saturn_color(normalize(input.local), srgb_to_linear(u_color.rgb));
    } else {
        albedo = srgb_to_linear(gas_giant ? gas_giant_color(normalize(input.local)) : u_color.rgb);
    }
    albedo *= input.albedo;
    float limb = gas_giant || saturn ? lerp(0.65, 1.0, pow(saturate(dot(n, V)), 0.4)) : 1.0;

    // Sunlight through the planet's own rings (mip level from the footprint on the profile).
    float ring_u = u_ring_center.w > 0.5 ? ring_profile_u(input.rel_pos, L) : -1.0;
    float ring_lod = log2(max(fwidth(ring_u) * u_ring_radii.z, 1.0));
    float ring_tau = u_ring_profile.SampleLevel(u_ring_sampler, float2(ring_u, 0.5), ring_lod).r;
    float ring_t = ring_u >= 0.0 && ring_u <= 1.0 ? exp(-ring_tau / max(abs(dot(L, u_ring_normal.xyz)), 1e-3)) : 1.0;

    float ndl = dot(n, L);
    float diffuse = saturate(ndl) * smoothstep(-0.05, 0.05, ndl);
    float lit = diffuse * eclipse_factor(input.rel_pos, L) * ring_t;

    float3 color = albedo * (lit * limb + u_params.y);
    return float4(color, 1.0);
}
