// Planet / moon shading: Lambert lighting by one or two stars (e.g. both suns
// of a circumbinary planet), soft eclipse shadows from other bodies (each star
// treated as a disk), the shadow of the planet's own rings,
// procedural gas-giant bands (Jupiter-like, Saturn-like) and a procedural
// battle station. Under an atmosphere, sunlight is reddened by the column it
// crosses (the twilight band at the terminator). City lights come on at dusk
// and shine up through the atmosphere. A relief map tilts the normal by the
// surface's slopes, and its water mirrors the sun; a cloud map lays a smooth white
// layer over all of it (scene/relief.hpp).

#include "atmosphere.hlsli"
#include "common.hlsli"
#include "noise.hlsli"

#define MAX_OCCLUDERS 8
#define STYLE_GAS_GIANT 1
#define STYLE_STAR 2       // self-luminous: u_color is linear HDR radiance
#define STYLE_BLACK_HOLE 3 // event horizon: absorbs everything
#define STYLE_DEATH_STAR 4
#define STYLE_SATURN 5
#define STYLE_ACCRETION_DISK 6
#define STYLE_ACCRETION_DONOR 7

static const float kPi = 3.14159265;

FRAGMENT_TEXTURE(Texture2D, u_albedo, u_sampler, 0); // sRGB texture: samples are linear
FRAGMENT_TEXTURE(Texture2D, u_ring_profile, u_ring_sampler, 1); // r = ring normal optical depth
FRAGMENT_TEXTURE(Texture2D, u_night_lights, u_night_sampler, 2); // sRGB gray map of emitted light
FRAGMENT_TEXTURE(Texture2D, u_field_map, u_field_sampler, 3); // r: radial magnetic field, 0.5 = none, compressed
FRAGMENT_TEXTURE(Texture2D, u_relief_map, u_relief_sampler, 4); // rg: slopes east and north, compressed; b: water
FRAGMENT_TEXTURE(Texture2D, u_cloud_map, u_cloud_sampler, 5); // r: cloud cover

cbuffer Uniforms : register(b0, space3)
{
    float4 u_sun;                      // xyz = unit direction to sun, w = sun angular radius (rad)
    float4 u_light2;                   // the same for a second star
    float4 u_light2_color;             // rgb = its light relative to the sun's (linear; 0 = none)
    float4 u_color;                    // rgb = base color (sRGB), a = style
    float4 u_params;                   // x = occluder count, y = ambient, z = 1 if textured, w = 1 to flip u
    float4 u_occluders[MAX_OCCLUDERS]; // xyz = camera-relative center, w = radius
    float4 u_ring_center;              // xyz = camera-relative center of the ring plane, w = 1 if it has rings
    float4 u_ring_normal;              // xyz = ring plane normal (unit)
    float4 u_ring_radii;               // x = inner, y = outer radius of the profile (km), z = profile samples
    float4 u_atmo_rayleigh;            // atmosphere: rgb = vertical optical depth, w = scale height (radii; 0 = none)
    float4 u_atmo_haze;                // rgb = vertical attenuation optical depth, w = scale height
    float4 u_night;                    // night lights: rgb = color (sRGB), w = I/F of a white pixel (0 = none)
    float4 u_surface;                  // x = scale of the reflected light, y = 1 to show the field map, w = its u offset,
                                       // z = 1 if the map is of the strength instead of Br
    float4 u_craft;                    // x = 1 if the vertex albedo is baked sunlight (a spacecraft's mesh),
                                       // y = illustrative fill light on the unlit side
    float4 u_spots[4];                 // a star's dark spots: xyz = direction in the body's axes, w = angular radius
    float4 u_spot_color;               // rgb = their linear radiance, a = how many
    float4 u_flare;                    // a star's flare: xyz = site in the body's axes, w = angular radius (0 = none)
    float4 u_flare_color;              // rgb = its linear radiance
    float4 u_spot_tints[4];            // a planet's dark spots (their places in u_spots, how many in u_spot_color.a):
                                       // rgb = relative change of the albedo at the center
    float4 u_relief;                   // x = slope of a full channel (0 = no map), y = the sea's mean square slope
                                       // (0 = no glint), z = its reflectance at normal incidence,
                                       // w = the albedo of the clouds (0 = none)
    float4 u_cloud_flow;               // the clouds' drift: x = place in the cycle (0..1), y = the cycle (days),
                                       // z, w = deg of longitude per day in the tropics and the middle latitudes
    float4 u_accretion;                // rgb = inner disk color, w = 1 / outer radius (km)
    float4 u_accretion_flow;           // x = turns of the disk texture
};

// A wobble of an edge on a star's surface that changes smoothly over the sphere.
float surface_wobble(float3 p)
{
    return sin(9.1 * p.x + 1.3) * sin(7.7 * p.y + 0.4) + 0.6 * sin(17.3 * p.z + 2.1) * sin(13.9 * p.y + 5.2);
}

// How much of a star's surface at unit vector p (body axes) is dark spot. The soft,
// uneven edge is only for the look: real spots are not circles. The light curve
// takes the circles (their area is about the same).
float star_spot_cover(float3 p)
{
    int count = (int)u_spot_color.a;
    float wobble = surface_wobble(p);
    float cover = 0.0;
    for (int i = 0; i < count; ++i) {
        float angle = acos(clamp(dot(p, u_spots[i].xyz), -1.0, 1.0));
        float x = angle / u_spots[i].w + 0.14 * wobble;
        cover = max(cover, 1.0 - smoothstep(0.78, 1.22, x));
    }
    return cover;
}

// A planet's dark spots: u_spots[i] = (east longitude, latitude, east-west width, north-south
// width), radians, widths at half maximum. dark_spot_weight in scene.cpp is the same.
float3 dark_spot_albedo(float3 p)
{
    float lat = asin(clamp(p.z, -1.0, 1.0));
    float lon = atan2(p.y, p.x);
    int count = (int)u_spot_color.a;
    float3 factor = 1.0;
    for (int i = 0; i < count; ++i) {
        float dlon = lon - u_spots[i].x;
        dlon -= 2.0 * kPi * round(dlon / (2.0 * kPi));
        float x = dlon * cos(u_spots[i].y) / u_spots[i].z;
        float y = (lat - u_spots[i].y) / u_spots[i].w;
        factor *= 1.0 + u_spot_tints[i].rgb * exp(-2.7725887 * (x * x + y * y));
    }
    return factor;
}

// A strength map's false color: viridis, dark where the field is weak (sRGB).
static const float3 kFieldWeak = float3(0.267, 0.005, 0.329);
static const float3 kFieldMiddle = float3(0.128, 0.567, 0.551);
static const float3 kFieldStrong = float3(0.993, 0.906, 0.144);
static const float kFieldStrengthOpacity = 0.8;

// The field map's false color, a map and not light, so it shows at night too: red
// where the field leaves the planet, blue where it enters.
static const float3 kFieldOut = float3(0.92, 0.26, 0.18);  // sRGB
static const float3 kFieldIn = float3(0.20, 0.47, 0.96);   // sRGB
static const float kFieldNight = 0.35;                     // brightness without sunlight
static const float kFieldDim = 0.6;                        // the surface under the map

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

float eclipse_factor(float3 p, float3 L, float light_radius)
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
        float penumbra = max(along * light_radius, 1e-3);
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

// How much of the light of a star in direction L (angular radius light_radius) reaches p:
// eclipses by other bodies and the planet's own rings.
float star_visibility(float3 p, float3 L, float light_radius)
{
    // Starlight through the rings (mip level from the footprint on the profile).
    float ring_u = u_ring_center.w > 0.5 ? ring_profile_u(p, L) : -1.0;
    float ring_lod = log2(max(fwidth(ring_u) * u_ring_radii.z, 1.0));
    float ring_tau = u_ring_profile.SampleLevel(u_ring_sampler, float2(ring_u, 0.5), ring_lod).r;
    float ring_t = ring_u >= 0.0 && ring_u <= 1.0 ? exp(-ring_tau / max(abs(dot(L, u_ring_normal.xyz)), 1e-3)) : 1.0;

    return eclipse_factor(p, L, light_radius) * ring_t;
}

// The cloud cover at a place of the map. The picture drifts along the parallels: two copies
// half a cycle apart, each sheared by the drift and faded out before it starts over;
// cloud_copy in scene/relief.hpp is the same.
float cloud_cover(float2 uv)
{
    float lat = (0.5 - uv.y) * kPi;
    float c = cos(lat);
    float s2 = sin(2.0 * lat);
    float c2 = c * c;
    float turns = (u_cloud_flow.w * s2 * s2 + u_cloud_flow.z * c2 * c2 * c2) * u_cloud_flow.y / 360.0;
    if (u_params.w > 0.5) {
        turns = -turns; // the map's east is to the left
    }
    float p0 = frac(u_cloud_flow.x);
    float p1 = frac(u_cloud_flow.x + 0.5);
    // The picture moves east: the cloud now here was to the west. One set of derivatives for
    // both lookups (frac jumps, and the two differ only by a smooth shift).
    float2 dx = ddx(uv);
    float2 dy = ddy(uv);
    float a = u_cloud_map.SampleGrad(u_cloud_sampler, float2(uv.x - turns * (p0 - 0.5), uv.y), dx, dy).r;
    float b = u_cloud_map.SampleGrad(u_cloud_sampler, float2(uv.x - turns * (p1 - 0.5), uv.y), dx, dy).r;
    return lerp(b, a, 1.0 - abs(2.0 * p0 - 1.0));
}

// Lambert light on a surface whose normal is n, tilted by the relief from the smooth
// surface's n0: the smooth surface's horizon still hides the star (there are no cast
// shadows, so a slope facing a star that has set would otherwise stay lit).
float lambert(float3 n, float3 n0, float3 L)
{
    return saturate(dot(n, L)) * smoothstep(-0.05, 0.05, dot(n0, L));
}

// Smith's shadowing of Gaussian slopes (the fit of Walter et al. 2007); sea_shadowing in
// scene/relief.hpp.
float sea_shadowing(float mu, float variance)
{
    float a = mu / max(sqrt(variance * max(1.0 - mu * mu, 0.0)), 1e-9);
    return a < 1.6 ? (3.535 * a + 2.181 * a * a) / (1.0 + 2.276 * a + 2.577 * a * a) : 1.0;
}

// The glint of a star off the sea, in the units of albedo x Lambert light; sea_glint in
// scene/relief.hpp is the same.
float sea_glint(float3 n, float3 L, float3 V)
{
    float mu_sun = dot(n, L);
    float mu_view = dot(n, V);
    if (mu_sun <= 0.0 || mu_view <= 0.0) {
        return 0.0;
    }
    float variance = u_relief.y;
    float3 h = normalize(L + V);
    float c = max(dot(n, h), 1e-4);
    float c2 = c * c;
    float facets = exp(-(1.0 - c2) / (c2 * variance)) / (variance * c2 * c2);
    float fresnel = u_relief.z + (1.0 - u_relief.z) * pow(1.0 - min(dot(V, h), 1.0), 5.0);
    return fresnel * facets * sea_shadowing(mu_sun, variance) * sea_shadowing(mu_view, variance) / (4.0 * mu_view);
}

// Fraction of sunlight (per channel) that reaches the surface through the
// atmosphere, with mu the cosine of the sun's zenith angle.
float3 atmosphere_transmission(float mu)
{
    if (u_atmo_rayleigh.w <= 0.0) {
        return 1.0;
    }
    return exp(-sun_optical_depth(1.0, mu, u_atmo_rayleigh, u_atmo_haze));
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
        float h = hash12(cell);
        const float kContrast[3] = {0.05, 0.10, 0.16}; // subtle large sectors, busier fine plating
        shade *= lerp(1.0, 1.0 + kContrast[k] * (2.0 * h - 1.0), detail);
        float2 f = abs(frac(float2(x, y)) - 0.5);
        // The coarse level only draws latitude lines (the hull bands); finer levels draw panel edges.
        float edge = smoothstep(0.5 - max(w, 0.02), 0.5, k == 0 ? f.y : max(f.x, f.y));
        const float kSeam[3] = {0.6, 0.3, 0.2};
        seam = max(seam, edge * detail * kSeam[k]);
        if (k == 2) {
            // Lit windows on 3.5% of the fine panels; their mean glow when too small to resolve.
            window = lerp(0.035, step(0.965, hash12(cell + 31.7)), detail);
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

float4 main(PSInput input, bool front_face : SV_IsFrontFace) : SV_Target0
{
    float3 n = normalize(input.normal);
    // Only a spacecraft's mesh is drawn from both sides: a sheet is lit on the side seen.
    if (!front_face) {
        n = -n;
    }
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
        float lit_s = saturate(ndl_s) * smoothstep(-0.05, 0.05, dot(n, L)) * eclipse_factor(input.rel_pos, L, u_sun.w);
        float night = 1.0 - smoothstep(-0.15, 0.1, dot(n, L));
        float3 albedo_s = srgb_to_linear(u_color.rgb) / 0.42 * s.albedo; // u_color tints the hull
        return float4(albedo_s * (lit_s + u_params.y) + s.lights * night + s.glow, 1.0);
    }
    if (style == STYLE_ACCRETION_DISK) {
        const float3 p = input.local * u_accretion.w;
        const float r = length(p.xy);
        const float turns = atan2(p.y, p.x) / (2.0 * kPi) - u_accretion_flow.x;
        const float face = pow(abs(dot(n, input.axis_z)), 2.0);
        const float shear = turns + 0.34 * log(r + 0.16);
        const float roll = turn_noise(turns, 4.0 * p.z + 6.0, 8.0);
        const float bands = turn_noise(shear, 29.0 * r + 0.9 * roll, 13.0);
        const float fine = turn_noise(shear, 92.0 * r + 3.0 * bands, 31.0);
        const float rim = turn_noise(turns + 0.08 * p.z, 33.0 * p.z + roll + 19.0, 12.0);
        const float grain = lerp(0.65 * rim + 0.35 * roll, (0.7 * bands + 0.3 * fine) * smoothstep(0.06, 0.25, r), face);
        const float hot = saturate(0.28 * log2(max(input.albedo, 1.0)));
        const float3 color = lerp(u_color.rgb, u_accretion.rgb, hot) * input.albedo;
        const float mu = saturate(dot(n, V));
        return float4(color * (1.0 - 0.6 * (1.0 - mu)) * (1.0 + 0.28 * grain), 1.0);
    }
    if (style == STYLE_STAR || style == STYLE_ACCRETION_DONOR) {
        // Linear limb darkening, as for the sun sprite; a Roche-lobe star's
        // mesh carries its gravity darkening and irradiation in the albedo.
        float mu = saturate(dot(n, V));
        float3 surface = u_color.rgb * input.albedo;
        if (style == STYLE_ACCRETION_DONOR) {
            const float3 p = normalize(input.local);
            const float mottling = 0.73 * value_noise3(18.0 * p) + 0.27 * value_noise3(41.0 * p);
            surface *= 1.0 + 0.15 * mottling;
        }
        if (u_spot_color.a > 0.5) {
            surface = lerp(surface, u_spot_color.rgb, star_spot_cover(normalize(input.local)));
        }
        if (u_flare.w > 0.0) {
            // Soft and uneven at the edge; the area is what sets the light.
            float3 p = normalize(input.local);
            float x = acos(clamp(dot(p, u_flare.xyz), -1.0, 1.0)) / u_flare.w + 0.12 * surface_wobble(p);
            surface = lerp(surface, u_flare_color.rgb, 1.0 - smoothstep(0.72, 1.28, x));
        }
        return float4(surface * (1.0 - 0.6 * (1.0 - mu)), 1.0);
    }

    bool gas_giant = (int)u_color.a == STYLE_GAS_GIANT;
    bool saturn = (int)u_color.a == STYLE_SATURN;
    float3 albedo;
    float2 uv = input.uv;
    if (u_params.w > 0.5) {
        uv.x = 1.0 - uv.x;
    }
    if (u_params.z > 0.5) {
        albedo = u_albedo.Sample(u_sampler, uv).rgb;
    } else if (saturn) {
        albedo = saturn_color(normalize(input.local), srgb_to_linear(u_color.rgb));
    } else {
        albedo = srgb_to_linear(gas_giant ? gas_giant_color(normalize(input.local)) : u_color.rgb);
    }
    if (u_spot_color.a > 0.5) {
        albedo *= dark_spot_albedo(normalize(input.local));
    }
    float baked_sunlight = 1.0;
    if (u_craft.x > 0.5) {
        baked_sunlight = input.albedo;
    } else {
        albedo *= input.albedo;
    }
    float limb = gas_giant || saturn ? lerp(0.65, 1.0, pow(saturate(dot(n, V)), 0.4)) : 1.0;

    // The relief tilts the normal that the sunlight falls on; horizons, the atmosphere and the
    // sea keep the smooth surface's.
    float3 relief = u_relief_map.Sample(u_relief_sampler, uv).rgb;
    float3 nr = n;
    if (u_relief.x > 0.0) {
        float2 s = 2.0 * relief.rg - 1.0;
        s *= abs(s) * u_relief.x;
        if (u_params.w > 0.5) {
            s.x = -s.x; // the map's east is to the left
        }
        float3 east = cross(input.axis_z, n);
        float len = length(east);
        if (len > 1e-4) {
            east /= len;
            nr = normalize(n - s.x * east - s.y * cross(n, east));
        }
    }

    // Clouds: a smooth white layer in place of the surface, as much as they cover it.
    float cover = cloud_cover(uv);
    if (u_relief.w > 0.0) {
        albedo = lerp(albedo, u_relief.w, cover);
        nr = normalize(lerp(nr, n, cover));
    } else {
        cover = 0.0;
    }

    float3 sun_reach = star_visibility(input.rel_pos, L, u_sun.w) * atmosphere_transmission(dot(n, L));
    float3 lit = lambert(nr, n, L) * sun_reach;
    float3 glint = 0.0;
    bool sea = u_relief.x > 0.0 && u_relief.y > 0.0;
    if (sea) {
        glint = sea_glint(n, L, V) * sun_reach;
    }
    if (any(u_light2_color.rgb > 0.0)) {
        float3 reach2 = u_light2_color.rgb * star_visibility(input.rel_pos, u_light2.xyz, u_light2.w) *
                        atmosphere_transmission(dot(n, u_light2.xyz));
        lit += lambert(nr, n, u_light2.xyz) * reach2;
        if (sea) {
            glint += sea_glint(n, u_light2.xyz, V) * reach2;
        }
    }

    lit *= baked_sunlight;
    // The fill comes from about the viewer's side, so the unlit side still shows its shape.
    float ambient = u_params.y + u_craft.y * (0.35 + 0.65 * saturate(dot(n, V)));
    float3 color = albedo * (lit * limb + ambient) * u_surface.x;
    // The mirrored sunlight keeps the sun's color, and leaves through the atmosphere.
    // It crosses the clouds twice.
    color += glint * (relief.b * (1.0 - cover) * (1.0 - cover) * baked_sunlight * u_surface.x) *
             atmosphere_transmission(saturate(dot(n, V)));
    if (u_night.w > 0.0) {
        // Lights come on between sunset and the sun 6 deg below the horizon,
        // and reach the viewer through the atmosphere.
        float dusk = 1.0 - smoothstep(-0.1045, 0.0, dot(n, L));
        float lights = u_night_lights.Sample(u_night_sampler, uv).r;
        color += srgb_to_linear(u_night.rgb) * (u_night.w * lights * dusk * (1.0 - cover)) *
                 atmosphere_transmission(saturate(dot(n, V)));
    }
    if (u_surface.y > 0.5) {
        // From the vertex's u, before any mirroring of a west-positive texture.
        float value = 2.0 * u_field_map.Sample(u_field_sampler, float2(input.uv.x + u_surface.w, input.uv.y)).r - 1.0;
        float3 tint = srgb_to_linear(value > 0.0 ? kFieldOut : kFieldIn);
        float shade = lerp(kFieldNight, 1.0, saturate(dot(lit, float3(0.2126, 0.7152, 0.0722))));
        // The map holds the square root of the strength. Mixed in linear light, that root would
        // look like a flat field of color wherever the field is more than a trace: its 1.5th
        // power (strength^0.75) keeps the weak fields faint. The surface is dimmed everywhere
        // alike: dimming with the strength would draw bright lines where the field changes sign.
        float mix_amount = pow(saturate(abs(value)), 1.5);
        if (u_surface.z > 0.5) {
            // Mixed in sRGB, as the ramp is uniform there.
            float t = saturate(0.5 * (value + 1.0));
            tint = srgb_to_linear(t < 0.5 ? lerp(kFieldWeak, kFieldMiddle, 2.0 * t) : lerp(kFieldMiddle, kFieldStrong, 2.0 * t - 1.0));
            mix_amount = kFieldStrengthOpacity;
        }
        color = lerp(color * kFieldDim, tint * shade, mix_amount);
    }
    return float4(color, 1.0);
}
