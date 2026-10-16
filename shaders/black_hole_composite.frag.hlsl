// Composites the black hole into the scene at the scene's resolution, from the
// reduced-resolution prepass (black_hole.frag.hlsl):
// - where the four prepass rays around a pixel agree, their disk light,
//   transmission and sky direction are interpolated and the sky is sampled per
//   pixel along that direction, with mip level and anisotropy from how fast the
//   direction changes across the screen (lensing squeezes the sky near the
//   photon ring); sampling the sky in the prepass instead blurred the stars
//   and lost those falling between its rays. (Splitting footprints longer than
//   the 8:1 of anisotropic filtering into several taps cost ~2 ms at 4.7 Mpx and
//   made no visible difference.)
// - where they disagree - the hole's edge, or sky directions more than ~8 deg
//   apart (close to the photon ring) - the ray is traced again for this pixel,
//   as in NPGS (BlackHole_composite.frag.glsl).
// Writes the hole's depth where the hole or the disk is hit; the lensed sky
// only replaces empty background.

#include "black_hole.hlsli"

FRAGMENT_TEXTURE(Texture2D, u_trace_color, u_trace_color_sampler, 2); // point sampled
FRAGMENT_TEXTURE(Texture2D, u_trace_sky, u_trace_sky_sampler, 3);     // point sampled

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
        sky = u_sky.SampleGrad(u_sky_sampler, dir, gx, gy).rgb;
    }

    PSOutput output;
    output.color = float4(color.rgb + color.a * sky, fade); // premultiplied
    output.depth = opaque ? u_params.y : 0.0;
    return output;
}
