// Aurora: light emitted in a shell above a body, integrated along the view ray
// in the body's unit-sphere frame (lengths there in equatorial radii, converted
// to km for the emission profiles).
//
// Emission = horizontal pattern * vertical profile, shown as the lighting
// allows: three independent choices (u_mode), described with their parameters
// in scene/aurora_model.hpp.

#include "common.hlsli"

static const int kSteps = 48;
static const int kOvalSamples = 48;
static const int kJovianTails = 6;
static const int kTailPoints = 20; // JovianAurora::kTailPoints
static const float kPi = 3.14159265;
static const float kNoiseCells = 32.0; // lattice cells across the noise tile

// scene/aurora_model.hpp: AuroraPattern, AuroraVertical, AuroraLighting.
#define PATTERN_EARTH_OVAL 0
#define PATTERN_JOVIAN_OVAL 1
#define PATTERN_JOVIAN_FOOTPRINTS 2
#define PATTERN_SATURN_OVAL 3
#define PATTERN_BODY_MAP_SHARES 4
#define PATTERN_UNIFORM 5
#define PATTERN_EUROPA 6
#define PATTERN_IO_FLANKS 7
#define PATTERN_BODY_MAP_SEASONS 8
#define PATTERN_BODY_MAP_DUSK 9
#define VERTICAL_EARTH_CHAPMAN 0
#define VERTICAL_CHAPMAN_PAIR 1
#define VERTICAL_EXPONENTIAL 2
#define VERTICAL_GAUSSIAN 3
#define VERTICAL_UNIFORM_COLUMN 4
#define VERTICAL_TABLE 5
#define LIGHTING_ALWAYS 0
#define LIGHTING_NIGHT_SIDE 1
#define LIGHTING_NIGHT_OR_ECLIPSED 2

FRAGMENT_TEXTURE(Texture2D, u_noise, u_noise_sampler, 0); // four independent value noises (rgba), tiling
FRAGMENT_TEXTURE(Texture2D, u_ring, u_ring_sampler, 1);
FRAGMENT_TEXTURE(Texture2D, u_map, u_map_sampler, 2);

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_unit; // camera-relative world vector -> unit-sphere frame (3x3)
    float4 u_camera;    // xyz = camera in the unit-sphere frame, w = 1 to start the rays there
    float4 u_pole;      // xyz = the pattern's north axis (AuroraDrawItem::pole), w = km per unit
    float4 u_noon;      // xyz = its noon (AuroraDrawItem::noon)
    float4 u_sun;       // xyz = toward the sun (unit), w = its angular radius
    float4 u_low;       // rgb = low layer, w = animation time (s)
    float4 u_high;      // rgb = high layer
    float4 u_layers;    // see AuroraVertical
    float4 u_shell;     // bottom radius, top radius, polar ratio, radians per pixel
    float4 u_params;    // 427.8 nm of a column through the oval (kR), diffuse flux, map gain, sigma (km)
    int4 u_mode;        // pattern, vertical, lighting
    float4 u_parent;    // the body's parent: xyz = center (unit-sphere frame), w = radius
    float4 u_oval[kOvalSamples]; // AuroraProfile::oval
    float4 u_south_center;
    float4 u_south_east;
    float4 u_jovian_pattern[128]; // oval pairs or the bounds and samples of the six footprint tails
    float4 u_footprints[14]; // physical unit direction, illustrative display gain: 6 main spots, 6 beam spots, Io's 2 reflected
    float4 u_table[4]; // AuroraProfile::table, four per vector
    float4 u_ring_radii; // inner/outer km (zero disables ring attenuation)
    float4 u_map_sun; // AuroraProfile::map_sun
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

// Sample i of 256 of the contours: north, south.
float2 jovian_oval(int i)
{
    float4 pair = u_jovian_pattern[i / 2];
    return (i & 1) != 0 ? pair.zw : pair.xy;
}

// Head 2024 reference contours in System III. The two hemispheres have
// separate centers and shapes. The latitude is planetocentric, so undo
// the unit-sphere frame's polar stretching before locating the contour.
float jovian_pattern(float3 n, float footprint_km)
{
    n = normalize(n * float3(1.0, 1.0, u_shell.z));
    if (abs(n.z) < 0.65) return 0.0;
    bool south = n.z < 0.0;
    float3 center = south ? u_south_center.xyz : u_pole.xyz;
    float3 east = south ? u_south_east.xyz : u_noon.xyz;
    float3 north = cross(center, east);
    float2 q = float2(dot(n, east), dot(n, north));
    float rho = atan2(length(q), dot(n, center));
    float angle = atan2(q.y, q.x);
    float x = frac(angle / (2.0 * kPi) + 1.0) * 256.0;
    int i = (int)floor(x);
    float2 a = jovian_oval(i);
    float2 b = jovian_oval((i + 1) % 256);
    float radius = south ? lerp(a.y, b.y, frac(x)) : lerp(a.x, b.x, frac(x));
    float slope = (south ? b.y - a.y : b.x - a.x) * (256.0 / (2.0 * kPi));
    float distance = (rho - radius) * u_pole.w /
        sqrt(1.0 + pow(slope / max(sin(radius), 0.01), 2.0));
    // Gaussian pixel filtering keeps a subpixel oval's integrated brightness.
    float sigma = sqrt(u_params.w * u_params.w + footprint_km * footprint_km / 12.0);
    float band = exp(-0.5 * pow(distance / sigma, 2.0)) * u_params.w / sigma;
    if (band < 0.0001) return 0.0;
    // Smooth illustrative modulation. Mean dusk/dawn contrast is 2 (Head
    // 2024); the texture detail is not a reconstruction of an observation.
    float2 dusk = float2(-u_sun.y, u_sun.x);
    float local = dot(n.xy, dusk) / max(length(n.xy) * length(dusk), 0.001);
    float2 p = q * u_pole.w;
    float detail = 0.7 + 0.6 * noise(p + float2(u_low.w * 0.8, 0.0), 1500.0).r;
    return band * (1.0 + local / 3.0) * detail;
}

float footprint_tails(float3 n, float footprint_km)
{
    // A transverse Gaussian for all three wakes, based on Io's roughly 1000 km width.
    const float sigma = 400.0;
    float variance = sigma * sigma + footprint_km * footprint_km / 12.0;
    float spread = 4.0 * sqrt(variance) / u_pole.w;
    float tails = 0.0;
    int hemisphere = n.z < 0.0 ? 1 : 0;
    for (int k = hemisphere; k < kJovianTails; k += 2) {
        float4 cap = u_jovian_pattern[k];
        float3 delta = n - cap.xyz;
        float radius = cap.w + spread;
        if (cap.w <= 0.0 || dot(delta, delta) > radius * radius) continue;
        int start = kJovianTails + k * kTailPoints;
        float4 a = u_jovian_pattern[start];
        float best = spread * spread;
        float gain = 0.0;
        [loop]
        for (int i = 1; i < kTailPoints; ++i) {
            float4 b = u_jovian_pattern[start + i];
            float3 ab = b.xyz - a.xyz;
            float3 offset = n - a.xyz;
            float t = saturate(dot(offset, ab) / max(dot(ab, ab), 1e-12));
            float3 d = offset - t * ab;
            float distance2 = dot(d, d);
            if (distance2 < best) {
                best = distance2;
                gain = lerp(a.w, b.w, t);
            }
            a = b;
        }
        tails += gain * exp(-0.5 * best * u_pole.w * u_pole.w / variance) * sigma * rsqrt(variance);
    }
    return tails;
}

// Spots and wakes at the moons' footprints.
float footprints_pattern(float3 n, float footprint_km)
{
    n = normalize(n * float3(1.0, 1.0, u_shell.z));
    if (abs(n.z) < 0.65) return 0.0;
    float variance = u_params.w * u_params.w + footprint_km * footprint_km / 12.0;
    float spots = 0.0;
    for (int k = 0; k < 14; ++k) {
        float4 spot = u_footprints[k];
        if (spot.w <= 0.0 || dot(n, spot.xyz) < 0.5) continue;
        float3 d = cross(n, spot.xyz) * u_pole.w;
        // Two-dimensional filtering preserves unresolved spot flux.
        spots += spot.w * exp(-0.5 * dot(d, d) / variance) *
                 (u_params.w * u_params.w / variance);
    }
    return spots + footprint_tails(n, footprint_km);
}

// Carbary 2012's statistical oval, in solar local time.
float saturn_pattern(float3 n, float footprint_km)
{
    n = normalize(n * float3(1.0, 1.0, u_shell.z));
    float colat = acos(min(abs(n.z), 1.0));
    if (colat > 0.65) return 0.0;
    float3 dusk = cross(float3(0.0, 0.0, 1.0), u_noon.xyz);
    float hour = 12.0 + atan2(dot(n, dusk), dot(n, u_noon.xyz)) * (12.0 / kPi);
    float x = hour * (kOvalSamples / 24.0);
    int k = (int)floor(x) % kOvalSamples;
    float4 oval = lerp(u_oval[k], u_oval[(k + 1) % kOvalSamples], frac(x));
    float radius = n.z < 0.0 ? oval.y : oval.x;
    float sigma_km = oval.z * u_pole.w;
    float sigma = sqrt(sigma_km * sigma_km + footprint_km * footprint_km / 12.0);
    float band = exp(-0.5 * pow((colat - radius) * u_pole.w / sigma, 2.0)) * sigma_km / sigma;
    // Statistical peak columns ~40/60 kR north/south near 05 h; dusk ~10 kR.
    // No System III rotation, stochastic storm or terrestrial Kp modulation.
    return band * (10.0 + (n.z < 0.0 ? 50.0 : 30.0) * oval.w) / 40.0;
}

// AuroraMap: west longitude, south-to-north rows, cell centers, one layer per channel
// (seasons: r = the map, g = what is left of it without sunlight; shares: the long-term
// map, then above, inside and below the current sheet); n is the physical direction.
float4 map_layers(float3 n, float footprint_km)
{
    float2 uv = float2(-atan2(n.y, n.x) / (2.0 * kPi), asin(clamp(n.z, -1.0, 1.0)) / kPi + 0.5);
    float width, height;
    u_map.GetDimensions(width, height);
    float cell_km = u_pole.w * (kPi / height);
    float lod = max(0.0, log2(max(footprint_km / cell_km, 1.0)));
    return u_map.SampleLevel(u_map_sampler, uv, lod) * u_params.z;
}

// aurora_map_pattern in scene/aurora_model.cpp does the same.
float map_pattern(float3 n, float footprint_km)
{
    n = normalize(n * float3(1.0, 1.0, u_shell.z));
    float4 layers = map_layers(n, footprint_km);
    if (u_mode.x == PATTERN_BODY_MAP_SHARES) {
        return dot(layers, u_map_sun);
    }
    if (u_mode.x == PATTERN_BODY_MAP_SEASONS) {
        // The sun's zenith angle at this latitude's noon: u_sun is in the body's frame.
        float zenith = abs(asin(clamp(n.z, -1.0, 1.0)) - asin(clamp(u_sun.z, -1.0, 1.0)));
        float sunlit = zenith <= u_map_sun.y
            ? lerp(1.0, u_map_sun.w, saturate((zenith - u_map_sun.x) / (u_map_sun.y - u_map_sun.x)))
            : u_map_sun.w * saturate((u_map_sun.z - zenith) / (u_map_sun.z - u_map_sun.y));
        return lerp(layers.y, layers.x, sunlit);
    }
    if (u_mode.x == PATTERN_BODY_MAP_DUSK) {
        float across = length(u_sun.xy);
        return max(layers.x + u_map_sun.x * dot(n.xy, float2(-u_sun.y, u_sun.x)) / max(across, 1e-6), 0.0);
    }
    return layers.x;
}

// Earth's aurora, relative to the oval's mean flux: the auroral oval between
// its poleward and equatorward boundaries (looked up at the point's magnetic
// local time) and, dimmer, the diffuse aurora down to its own boundary, in both
// hemispheres. Within the oval, illustrative structure that keeps the mean:
// east-west arcs a few tens of km wide, field-aligned rays along them and
// brighter and fainter stretches, which fade to their mean where a pixel covers
// them. They are drawn from a tiling noise texture laid on the plane across the
// dipole axis (the aurora stays within 43 deg of the poles, so the projection
// has no seam).
float earth_pattern(float3 n, float footprint_km)
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
    return stretch * (band * arcs * rays + u_params.y * diffuse * patches);
}

// europa_pattern in scene/satellite_aurora.hpp does the same: u_pole.xyz is the field line.
float europa_pattern(float3 n)
{
    float3 afternoon = float3(-u_sun.y, u_sun.x, 0.0);
    float across = length(afternoon);
    float dusk = across > 1e-6 ? dot(n, afternoon) / across : 0.0;
    float along = dot(n, u_pole.xyz);
    return (1.0 + u_map_sun.w * dusk) *
           exp(n.z * (u_map_sun.x + u_map_sun.y * n.y) + u_map_sun.z * (along * along - n.z * n.z));
}

// Brightness across the surface at unit direction n, seen with a pixel
// footprint of footprint_km.
float aurora_pattern(float3 n, float footprint_km)
{
    switch (u_mode.x) {
    case PATTERN_EARTH_OVAL: return earth_pattern(n, footprint_km);
    case PATTERN_JOVIAN_OVAL: return jovian_pattern(n, footprint_km);
    case PATTERN_JOVIAN_FOOTPRINTS: return footprints_pattern(n, footprint_km);
    case PATTERN_SATURN_OVAL: return saturn_pattern(n, footprint_km);
    case PATTERN_BODY_MAP_SHARES:
    case PATTERN_BODY_MAP_SEASONS:
    case PATTERN_BODY_MAP_DUSK: return map_pattern(n, footprint_km);
    case PATTERN_EUROPA: return europa_pattern(n);
    case PATTERN_IO_FLANKS: return 1.0 + 0.75 * (n.x * n.x - n.z * n.z);
    default: return 1.0; // PATTERN_UNIFORM
    }
}

// AuroraProfile::table at height h (km).
float table_layer(float h)
{
    float x = (h - u_layers.x) / u_layers.y;
    if (x <= 0.0 || x >= 15.0) return 0.0;
    int k = (int)floor(x);
    float a = u_table[k / 4][k % 4];
    float b = u_table[(k + 1) / 4][(k + 1) % 4];
    return lerp(a, b, frac(x));
}

float exponential_layer(float h, float scale, float top)
{
    return exp(-h / scale) / (scale * (1.0 - exp(-top / scale)));
}

// Emission per km at height h (km) where the pattern is `pattern`.
float3 layer_emission(float h, float pattern)
{
    switch (u_mode.y) {
    case VERTICAL_EARTH_CHAPMAN: {
        // Two Chapman layers: the green, blue and violet lines low, the red
        // line high. The red line follows the local 427.8 nm brightness by the
        // observed ratio (red_over_blue in scene/aurora.cpp): faint aurora is redder.
        float blue = u_params.x * pattern;
        return u_low.rgb * (pattern * chapman(h, u_layers.x, u_layers.y)) +
               u_high.rgb * (blue * red_over_blue(blue) * chapman(h, u_layers.z, u_layers.w));
    }
    case VERTICAL_CHAPMAN_PAIR:
        // Independent columns; for Saturn's visible colors the exposure and
        // profiles are display assumptions, not a UV-to-visible conversion.
        return pattern * (u_low.rgb * chapman(h, u_layers.x, u_layers.y) +
                          u_high.rgb * chapman(h, u_layers.z, u_layers.w));
    case VERTICAL_EXPONENTIAL:
        // The high layer is a separate, smooth corona (Io's sodium, a fixed
        // eclipse-average spectrum).
        return u_low.rgb * (pattern * exponential_layer(h, u_layers.x, u_layers.y)) +
               u_high.rgb * exponential_layer(h, u_layers.z, u_layers.y);
    case VERTICAL_GAUSSIAN: {
        // The shell's own limits cut it off.
        float z = (h - u_layers.x) / (h < u_layers.x ? u_layers.y : u_layers.z);
        return u_low.rgb * (pattern * exp(-0.5 * z * z));
    }
    case VERTICAL_TABLE:
        return u_low.rgb * (pattern * table_layer(h));
    default: // VERTICAL_UNIFORM_COLUMN: Ganymede's published zenith correction uses a uniform column.
        return u_low.rgb * (pattern / u_layers.x);
    }
}

// How much of the emission shows at the point p, above the unit direction n.
float lighting(float3 n, float3 p)
{
    if (u_mode.z == LIGHTING_ALWAYS) {
        return 1.0;
    }
    // Fades in from sunset at the point below to the sun 12 deg below the horizon
    // (it is there by day, but far fainter than the sunlit ground). The normal is
    // the ellipsoid's, in physical body coordinates.
    float3 normal = normalize(n / float3(1.0, 1.0, u_shell.z));
    float night = 1.0 - smoothstep(-0.2079, 0.0, dot(normal, u_sun.xyz));
    if (u_mode.z == LIGHTING_NIGHT_OR_ECLIPSED) {
        // Same spherical parent/soft penumbra approximation as body shading.
        float3 oc = u_parent.xyz - p;
        float along = dot(oc, u_sun.xyz);
        if (along > 0.0 && u_parent.w > 0.0) {
            float penumbra = max(along * u_sun.w, 1e-6);
            float miss = length(oc - along * u_sun.xyz);
            float shadow = 1.0 - smoothstep(u_parent.w - penumbra, u_parent.w + penumbra, miss);
            shadow *= saturate(pow(u_parent.w / penumbra, 2.0));
            night = max(night, shadow);
        }
    }
    return night;
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_unit, input.rel_pos));
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.unit_pos;

    // Sample the ring before divergent shell clipping. Explicit gradients
    // select the same radial mip footprint as the ring's visible surface.
    float ring_t = 0.0;
    float transmission = 1.0;
    if (u_ring_radii.y > u_ring_radii.x) {
        float dz = abs(dir.z) > 1e-7 ? dir.z : 1e-7;
        ring_t = -o.z / dz;
        float3 crossing = o + ring_t * dir;
        float radius = length(crossing.xy) * u_pole.w;
        float u = (radius - u_ring_radii.x) / (u_ring_radii.y - u_ring_radii.x);
        float2 profile = u_ring.SampleGrad(u_ring_sampler, float2(u, 0.5),
                                          float2(ddx_fine(u), 0.0), float2(ddy_fine(u), 0.0)).rg;
        float mu = abs(normalize(dir * float3(1.0, 1.0, u_shell.z)).z);
        if (abs(dir.z) > 1e-7 && u >= 0.0 && u <= 1.0 && dot(crossing - u_camera.xyz, dir) > 0.0)
            transmission = exp(-profile.r / max(mu, max(profile.g, 1e-4)));
    }

    float bottom = u_shell.x;
    float top = u_shell.y;
    float b = dot(o, dir);
    float3 perpendicular = cross(o, dir);
    float p2 = dot(perpendicular, perpendicular);
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
    float ray_km = km * length(dir * float3(1.0, 1.0, u_shell.z));
    float ds = length_total / kSteps;
    float jitter = gradient_noise(input.position.xy);
    float3 sum = 0.0;
    for (int i = 0; i < kSteps; ++i) {
        float s = (i + jitter) * ds;
        float t = s < l1 ? t0 + s : tb + (s - l1);
        float3 p = o + t * dir;
        float r = length(p);
        float3 n = p / r;
        float shown = lighting(n, p);
        if (shown <= 0.0) {
            continue;
        }
        float distance_km = length(p - u_camera.xyz) * km;
        float pattern = aurora_pattern(n, distance_km * u_shell.w);
        if (pattern <= 0.0) {
            continue;
        }
        float h = (r - 1.0) * km * length(n * float3(1.0, 1.0, u_shell.z));
        float3 emission = layer_emission(h, pattern);
        float ring_visibility = t > ring_t ? transmission : 1.0;
        sum += emission * (shown * ring_visibility * ds * ray_km);
    }
    return float4(sum, 1.0);
}
