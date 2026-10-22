// Earth's aurora: light emitted in a shell above the planet, integrated along
// the view ray in the planet's unit-sphere frame (lengths there in planet
// radii, converted to km for the emission profiles).
//
// Emission = horizontal pattern * vertical profile. Vertically, two Chapman
// layers exp(-z - e^-z) / H (z = (h - h_peak) / H, unit integral in km): the
// green, blue and violet lines low, the red line high. The red line follows
// the local 427.8 nm brightness by the observed ratio (red_over_blue in
// scene/aurora.cpp): faint aurora is redder. Horizontally, the
// auroral oval between its poleward and equatorward boundaries (looked up at
// the point's magnetic local time) and, dimmer, the diffuse aurora down to
// its own boundary, in both hemispheres. Within the oval, illustrative
// structure that keeps the mean: east-west arcs a few tens of km wide,
// field-aligned rays along them and brighter and fainter stretches, which fade
// to their mean where a pixel covers them. They are drawn from a tiling noise
// texture laid on the plane across the dipole axis (the aurora stays within
// 43 deg of the poles, so the projection has no seam). Night side only: the
// aurora fades out as the sun rises at the point below (it is there by day,
// but far fainter than the sunlit ground).

#include "common.hlsli"

static const int kSteps = 48;
static const int kOvalSamples = 48;
static const float kPi = 3.14159265;
static const float kNoiseCells = 32.0; // lattice cells across the noise tile

FRAGMENT_TEXTURE(Texture2D, u_noise, u_noise_sampler, 0); // four independent value noises (rgba), tiling

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_unit; // camera-relative world vector -> unit-sphere frame (3x3)
    float4 u_camera;    // xyz = camera in the unit-sphere frame, w = 1 to start the rays there
    float4 u_pole;      // xyz = north geomagnetic pole (unit), w = km per unit
    float4 u_noon;      // xyz = magnetic noon (unit, across the pole), w = shell bottom radius
    float4 u_sun;       // xyz = toward the sun (unit), w = shell top radius
    float4 u_low;       // rgb = column through the oval, low layer; w = animation time (s)
    float4 u_red;       // rgb = 1 kR of the red line; w = radians per pixel
    float4 u_blue;      // x = 427.8 nm of a column through the oval (kR)
    float4 u_layers;    // low peak (km), low scale height, high peak, high scale height
    float4 u_misc;      // x = diffuse aurora's flux relative to the oval's
    float4 u_oval[kOvalSamples]; // colatitudes (rad): x = diffuse, y = equatorward, z = poleward
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 unit_pos : TEXCOORD1;
};

float gradient_noise(float2 p)
{
    return frac(52.9829189 * frac(dot(p, float2(0.06711056, 0.00583715))));
}

// The noises (0..1, mean ~0.5) at a point of the polar plane (km), with
// lattice cells cell_km across.
float4 noise(float2 km, float cell_km)
{
    return u_noise.SampleLevel(u_noise_sampler, km / (cell_km * kNoiseCells), 0.0);
}

// I630.0 / I427.8 for a column of blue_kr: 2 at 0.1 kR, 0.4 at 3 kR.
float red_over_blue(float blue_kr)
{
    return 2.0 * pow(clamp(blue_kr, 0.1, 3.0) / 0.1, -0.4732);
}

float chapman(float h, float peak, float scale)
{
    float z = (h - peak) / scale;
    return exp(-z - exp(-z)) / scale;
}

// Brightness of the aurora at unit direction n (relative to the oval's mean
// flux), seen with a pixel footprint of footprint_km.
float aurora_pattern(float3 n, float footprint_km)
{
    float3 pole = u_pole.xyz;
    float m = dot(n, pole);
    float colat = acos(min(abs(m), 1.0));
    if (colat > 0.75) { // 43 deg: beyond any diffuse aurora of the model
        return 0.0;
    }
    float3 dusk = cross(pole, u_noon.xyz);
    float2 q = float2(dot(n, u_noon.xyz), dot(n, dusk)); // magnetic local time frame
    float mlt = 12.0 + atan2(q.y, q.x) * (12.0 / kPi);   // 0 .. 24 h
    float x = mlt * (kOvalSamples / 24.0);
    int k0 = (int)floor(x) % kOvalSamples;
    int k1 = (k0 + 1) % kOvalSamples;
    float3 oval = lerp(u_oval[k0].xyz, u_oval[k1].xyz, frac(x));

    const float edge = 0.01; // soft boundaries, ~0.6 deg
    float inside_pole = smoothstep(oval.z - edge, oval.z + edge, colat);
    float beyond_oval = smoothstep(oval.y - edge, oval.y + edge, colat);
    float band = inside_pole * (1.0 - beyond_oval);
    float diffuse = beyond_oval * (1.0 - smoothstep(oval.x - 2.0 * edge, oval.x + 2.0 * edge, colat));
    if (band + diffuse <= 0.0) {
        return 0.0;
    }

    // The polar plane (km); the southern cap takes another part of the noise.
    float2 p = q * u_pole.w + select(m < 0.0, float2(5100.0, 7300.0), float2(0.0, 0.0));
    float t = u_low.w;
    float4 coarse = noise(p + float2(t * 0.3, 0.0), 1500.0);
    float4 medium = noise(p + float2(0.0, t * 0.6), 400.0);
    // Brighter and fainter stretches along the oval (mean 1), slowly changing.
    float stretch = 0.3 + 1.4 * coarse.r;

    // Arcs: about three across the oval, folded on two scales, Gaussian in
    // cross-section; each brightens and fades along its length.
    float u = saturate((colat - oval.z) / max(oval.y - oval.z, 1e-3));
    float arc = u * 3.0 + 2.4 * (coarse.g - 0.5) + 0.6 * (medium.g - 0.5);
    float d = frac(arc) - 0.5;
    float width = 0.03 + 0.05 * medium.a; // of the spacing between arcs
    float strength = noise(p + floor(arc) * float2(370.0, 610.0) + float2(0.0, t * 0.5), 600.0).b;
    float arcs = exp(-d * d / (2.0 * width * width)) / (width * 2.5066283) * (strength * strength * strength / 0.19);
    arcs = lerp(1.0, arcs, saturate(1.5 - footprint_km / 30.0));
    // Rays along the arcs (field-aligned, so vertical columns in the shell).
    float rays = 1.0;
    float ray_detail = saturate(1.5 - footprint_km / 10.0);
    if (ray_detail > 0.0) {
        rays = lerp(1.0, 0.2 + 1.6 * noise(p + float2(t * 8.0, 0.0), 6.0).r, ray_detail);
    }
    float patches = 0.4 + 1.2 * medium.b;
    return stretch * (band * arcs * rays + u_misc.x * diffuse * patches);
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_unit, input.rel_pos));
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.unit_pos;

    float top = u_sun.w;
    float bottom = u_noon.w;
    float b = dot(o, dir);
    float p2 = dot(o, o) - b * b;
    if (p2 >= top * top) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }
    float q = sqrt(top * top - p2);
    float t0 = max(-b - q, 0.0);
    float t1 = -b + q;
    if (p2 < 1.0) {
        float ground = -b - sqrt(1.0 - p2);
        if (ground > 0.0) {
            t1 = min(t1, ground);
        }
    }
    // Below the shell's bottom nothing glows: march [t0, ta] and [tb, t1].
    float ta = t1;
    float tb = t1;
    if (p2 < bottom * bottom) {
        float qb = sqrt(bottom * bottom - p2);
        ta = clamp(-b - qb, t0, t1);
        tb = clamp(-b + qb, t0, t1);
    }
    float l1 = ta - t0;
    float length_total = l1 + (t1 - tb);
    if (length_total <= 0.0) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }

    float km = u_pole.w;
    float ds = length_total / kSteps;
    float jitter = gradient_noise(input.position.xy);
    float3 sum = 0.0;
    for (int i = 0; i < kSteps; ++i) {
        float s = (i + jitter) * ds;
        float t = s < l1 ? t0 + s : tb + (s - l1);
        float3 p = o + t * dir;
        float r = length(p);
        float3 n = p / r;
        // Fades in from sunset at the point below to the sun 12 deg below the horizon.
        float night = 1.0 - smoothstep(-0.2079, 0.0, dot(n, u_sun.xyz));
        if (night <= 0.0) {
            continue;
        }
        float distance_km = length(p - u_camera.xyz) * km;
        float pattern = aurora_pattern(n, distance_km * u_red.w);
        if (pattern <= 0.0) {
            continue;
        }
        float h = (r - 1.0) * km;
        float blue = u_blue.x * pattern;
        float3 emission = u_low.rgb * (pattern * chapman(h, u_layers.x, u_layers.y)) +
                          u_red.rgb * (blue * red_over_blue(blue) * chapman(h, u_layers.z, u_layers.w));
        sum += emission * (night * ds * km);
    }
    return float4(sum, 1.0);
}
