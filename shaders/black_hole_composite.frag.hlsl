// Composites the black hole into the scene at the scene's resolution, from the
// reduced-resolution prepass (black_hole.frag.hlsl):
// - where the four prepass rays around a pixel agree, their disk light,
//   transmission and sky direction are interpolated; where they disagree - the
//   hole's edge, or sky directions more than ~8 deg apart (close to the photon
//   ring) - the ray is traced again for this pixel, as in NPGS
//   (BlackHole_composite.frag.glsl);
// - the sky is then looked up per pixel along that direction (sky_light):
//   the stars as points where they fall (with the lens map linearized over the
//   pixel's neighbourhood), brightened by the magnification, and the Milky Way
//   from its own map; where lensing squeezes the sky so much that a dot would
//   reach farther than the star grid is searched (close to the photon ring),
//   from the mipmapped cube, in which the stars average out.
// Writes the hole's depth where the hole or the disk is hit; the lensed sky
// only replaces empty background.

#include "black_hole.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_trace_color, u_trace_color_sampler, 2); // point sampled
FRAGMENT_TEXTURE(Texture2D, u_trace_sky, u_trace_sky_sampler, 3);     // point sampled
FRAGMENT_TEXTURE(Texture2D, u_milky_way, u_milky_way_sampler, 4);     // equirectangular ICRF map (sRGB), mipmapped
// StarfieldPass::sky_stars() (32 bytes a star) and sky_cells(), after the samplers.
ByteAddressBuffer u_star_data : register(t5, space2);
ByteAddressBuffer u_star_cells : register(t6, space2);

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

struct PSOutput
{
    float4 color : SV_Target0;
    float depth  : SV_Depth;
};

static const float kEdgeCos = 0.99;      // neighbouring rays' sky directions this far apart: trace again
static const float kEdgeTransmit = 0.5;  // sky transmission differing this much (the hole's edge): trace again
static const float kSkyTransmit = 0.1;   // below this the sky barely shows (behind the disk): its direction does not matter

float4 trace_color(int2 p)
{
    p = clamp(p, int2(0, 0), int2(u_targets.xy) - 1);
    return u_trace_color.SampleLevel(u_trace_color_sampler, (float2(p) + 0.5) / u_targets.xy, 0);
}

float4 trace_sky(int2 p)
{
    p = clamp(p, int2(0, 0), int2(u_targets.xy) - 1);
    return u_trace_sky.SampleLevel(u_trace_sky_sampler, (float2(p) + 0.5) / u_targets.xy, 0);
}

bool sees_sky(float4 sky)
{
    return dot(sky.xyz, sky.xyz) > 0.25;
}

// Change of the sky direction per scene pixel along `step` (one prepass texel),
// from the neighbour on whichever side the rays agree with `base`.
float3 sky_gradient(int2 p, int2 step, float3 base, float texels_per_pixel)
{
    float4 ahead = trace_sky(p + step);
    float4 behind = trace_sky(p - step);
    float3 g = 0.0;
    float best = 1e30;
    if (sees_sky(ahead) && dot(ahead.xyz, base) > kEdgeCos) {
        g = ahead.xyz - base;
        best = dot(g, g);
    }
    if (sees_sky(behind) && dot(behind.xyz, base) > kEdgeCos) {
        float3 gb = base - behind.xyz;
        if (dot(gb, gb) < best) {
            g = gb;
        }
    }
    return g * texels_per_pixel;
}

// The Milky Way map along dir (gx, gy: change of dir per pixel), as the
// starfield draws it: RA 0h at the centre, increasing to the left.
float3 milky_way_light(float3 dir, float3 gx, float3 gy)
{
    if (u_stars.z <= 0.0) {
        return 0.0;
    }
    float rxy2 = max(dir.x * dir.x + dir.y * dir.y, 1e-12);
    float2 uv = float2(frac(0.5 - atan2(dir.y, dir.x) / (2.0 * kPi)), 0.5 - asin(clamp(dir.z, -1.0, 1.0)) / kPi);
    // d ra = (x dy - y dx) / (x^2 + y^2), d dec = dz / cos(dec): no seam at RA 0h.
    float cos_dec = sqrt(rxy2);
    float2 dx = float2(-(dir.x * gx.y - dir.y * gx.x) / rxy2 / (2.0 * kPi), -gx.z / cos_dec / kPi);
    float2 dy = float2(-(dir.x * gy.y - dir.y * gy.x) / rxy2 / (2.0 * kPi), -gy.z / cos_dec / kPi);
    return u_milky_way.SampleGrad(u_milky_way_sampler, uv, dx, dy).rgb * u_stars.z;
}

// The stars seen in this pixel, drawn like the starfield's dots (same size in
// pixels and profile) at the place each star's image falls: the lens map is
// linearized as sky offset = gx u + gy v (u, v in pixels). Their brightness is
// multiplied by the magnification, the pixel's solid angle without the hole
// over the sky's solid angle it now sees (capped: point images on a caustic).
float3 point_stars(float3 dir, float3 gx, float3 gy, float2 ndc)
{
    float a = dot(gx, gx);
    float b = dot(gx, gy);
    float c = dot(gy, gy);
    float det = a * c - b * b; // |gx x gy|^2
    if (det <= 0.0) {
        return 0.0;
    }
    float2 screen = ndc * float2(u_right_dir.w, u_up_dir.w);
    float len2 = 1.0 + dot(screen, screen);
    float pixel = 4.0 * u_right_dir.w * u_up_dir.w / (u_targets.z * u_targets.w) / (len2 * sqrt(len2));
    float magnification = min(pixel / sqrt(det), u_sky_params.y);

    // The 2 x 2 x 2 grid cells nearest to dir hold every star within half a cell
    // (1 / grid) of it; sky_light() hands over to the cube before a dot reaches farther.
    int grid = (int)u_stars.x;
    int3 c0 = int3(floor((dir * 0.5 + 0.5) * u_stars.x - 0.5));
    float3 sum = 0.0;
    [loop]
    for (int k = 0; k < 8; ++k) {
        int3 cell = c0 + int3(k & 1, (k >> 1) & 1, k >> 2);
        if (any(cell < 0) || any(cell >= grid)) {
            continue;
        }
        uint index = ((uint)cell.z * (uint)grid + (uint)cell.y) * (uint)grid + (uint)cell.x;
        uint first = u_star_cells.Load(index * 4);
        uint last = u_star_cells.Load(index * 4 + 4);
        [loop]
        for (uint i = first; i < last; ++i) {
            float4 dir_size = asfloat(u_star_data.Load4(i * 32));
            float4 color = asfloat(u_star_data.Load4(i * 32 + 16));
            float3 delta = dir_size.xyz - dir;
            float px = dot(gx, delta);
            float py = dot(gy, delta);
            float u = (c * px - b * py) / det;
            float v = (a * py - b * px) / det;
            // starfield.vert/.frag: a quad size px wide, exp(-3.5 r^2) with r = 1 at its edges.
            float half_size = 0.5 * dir_size.w;
            if (max(abs(u), abs(v)) < half_size) {
                sum += color.rgb * color.a * exp(-3.5 * (u * u + v * v) / (half_size * half_size));
            }
        }
    }
    return sum * (u_stars.y * magnification);
}

// The sky seen through a pixel along dir (gx, gy: change of dir per pixel).
float3 sky_light(float3 dir, float3 gx, float3 gy, float2 ndc)
{
    // A star's dot reaches as far as the largest stretch of the lens map times the
    // largest dot radius from its image; point_stars() finds stars within 1 / grid.
    float a = dot(gx, gx);
    float b = dot(gx, gy);
    float c = dot(gy, gy);
    float stretch = sqrt(0.5 * (a + c) + sqrt(0.25 * (a - c) * (a - c) + b * b));
    float cube = smoothstep(0.6, 1.0, stretch * u_sky_params.x * u_stars.x);
    float3 light = 0.0;
    if (cube > 0.0) {
        light += cube * u_sky.SampleGrad(u_sky_sampler, dir, gx, gy).rgb;
    }
    if (cube < 1.0) {
        light += (1.0 - cube) * (milky_way_light(dir, gx, gy) + point_stars(dir, gx, gy, ndc));
    }
    return light;
}

PSOutput main(PSInput input)
{
    // The four prepass rays around this pixel and their bilinear weights.
    float2 pos = input.uv * u_targets.xy - 0.5;
    float2 f = pos - floor(pos);
    int2 p0 = int2(floor(pos));
    const int2 offsets[4] = {int2(0, 0), int2(1, 0), int2(0, 1), int2(1, 1)};
    float weights[4] = {(1.0 - f.x) * (1.0 - f.y), f.x * (1.0 - f.y), (1.0 - f.x) * f.y, f.x * f.y};
    float4 c[4];
    float4 s[4];
    bool lensed = false;
    bool all_sky = true; // all four rays see the sky: the gradients come from them
    [unroll]
    for (int i = 0; i < 4; ++i) {
        c[i] = trace_color(p0 + offsets[i]);
        s[i] = trace_sky(p0 + offsets[i]);
        lensed = lensed || s[i].w != 0.0;
        all_sky = all_sky && sees_sky(s[i]);
    }
    if (!lensed) {
        discard;
    }

    bool edge = false;
    [unroll]
    for (int j = 0; j < 4; ++j) {
        [unroll]
        for (int k = j + 1; k < 4; ++k) {
            if (s[j].w != 0.0 && s[k].w != 0.0) {
                edge = edge || abs(c[j].a - c[k].a) > kEdgeTransmit;
                edge = edge || (min(c[j].a, c[k].a) > kSkyTransmit && sees_sky(s[j]) && sees_sky(s[k]) &&
                                dot(s[j].xyz, s[k].xyz) < kEdgeCos);
            }
        }
    }

    float4 color = 0.0; // rgb = disk light * fade, a = sky transmission * fade
    float fade = 0.0;
    float3 dir = 0.0;
    bool opaque = false;
    if (edge) {
        Trace t = trace_pixel(input.uv);
        if (!t.lensed) {
            discard;
        }
        color = float4(t.disk * t.fade, t.transmit * t.fade);
        fade = t.fade;
        dir = t.sky;
        opaque = t.opaque;
    } else {
        float nearest = -1.0;
        [unroll]
        for (int i = 0; i < 4; ++i) {
            color += weights[i] * c[i];
            fade += weights[i] * abs(s[i].w);
            dir += weights[i] * c[i].a * s[i].xyz; // the rays that see the sky, as much as they do
            if (weights[i] > nearest) {
                nearest = weights[i];
                opaque = s[i].w < 0.0;
            }
        }
    }

    float3 sky = 0.0;
    if (color.a > 0.0 && dot(dir, dir) > 1e-12) {
        dir = normalize(dir);
        float2 texels_per_pixel = u_targets.xy / u_targets.zw;
        float3 gx;
        float3 gy;
        if (!edge && all_sky) {
            gx = 0.5 * ((s[1].xyz - s[0].xyz) + (s[3].xyz - s[2].xyz)) * texels_per_pixel.x;
            gy = 0.5 * ((s[2].xyz - s[0].xyz) + (s[3].xyz - s[1].xyz)) * texels_per_pixel.y;
        } else {
            int2 pn = int2(floor(input.uv * u_targets.xy)); // prepass texel under this pixel
            float4 sn = trace_sky(pn);
            float3 base = dir;
            if (sees_sky(sn) && dot(sn.xyz, dir) > kEdgeCos) {
                base = sn.xyz;
            }
            gx = sky_gradient(pn, int2(1, 0), base, texels_per_pixel.x);
            gy = sky_gradient(pn, int2(0, 1), base, texels_per_pixel.y);
        }
        float2 ndc = float2(input.uv.x * 2.0 - 1.0, 1.0 - input.uv.y * 2.0);
        sky = sky_light(dir, gx, gy, ndc);
    }

    PSOutput output;
    output.color = float4(color.rgb + color.a * sky, fade); // premultiplied
    output.depth = opaque ? u_params.y : 0.0;
    return output;
}
