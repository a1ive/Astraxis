// What a comet's impacts left on a planet, multiplied over its surface: dark debris high in the
// atmosphere, seen through twice (sunlight in, light out), so darker toward the limb and the
// terminator [BEE96 Sec. 1]. The patterns are the same as src/scene/impact.cpp: change them together.
// Everything is in the body frame (km). The box is around the whole planet; each ray meets the ellipsoid
// and writes its depth, so what stands in front of the planet is untouched.

#define MAX_SITES 24

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_local;  // camera-relative world vector -> body frame (3x3 rotation)
    float4x4 u_view_proj;
    float4 u_camera;      // xyz = camera (body frame)
    float4 u_sun;         // xyz = unit direction to the sun (body frame)
    float4 u_radii;       // x = equatorial, y = polar, z = site count, w = km per radian near the sites
    float4 u_band;        // x = depth of the band, y, z = latitudes of its spread (rad, south, north)
    float4 u_tint;        // rgb = visible absorption per unit of optical depth at 890 nm
    float4 u_site[MAX_SITES];   // x = east longitude, y = latitude (rad), z = age (s), w = depth at the darkest
    float4 u_shape[MAX_SITES];  // x = ejecta km, y = core km, z = expansion km/s, w = seed + 16 * rings
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 local    : TEXCOORD1;
};

struct PSOutput
{
    float4 color : SV_Target0;
    float depth  : SV_Depth;
};

static const float kTwoPi = 6.2831853;
static const float kDeg = 0.017453293;

float smooth_step(float a, float b, float x)
{
    float s = saturate((x - a) / (b - a));
    return s * s * (3.0 - 2.0 * s);
}

float scar_shape(float east, float north, float age, float ejecta_km, float core_km, float expansion, int rings, float seed)
{
    const float ax = 0.57357644; // sin 145 deg
    const float ay = -0.81915204; // cos 145 deg
    float r = sqrt(east * east + north * north);
    float x = east * ax + north * ay;
    float y = -east * ay + north * ax;
    float ejecta = ejecta_km + expansion * min(age, 172800.0);
    float landed = smooth_step(300.0, 1100.0, age);
    float sum = 0.0;

    float wave = 0.025 * ejecta_km;
    float xw = wave > 0.0 ? x + wave * sin(y / (3.5 * wave) + seed) : x;
    float yw = wave > 0.0 ? y + wave * sin(x / (2.8 * wave) + 1.7 * seed) : y;

    float grown = ejecta_km > 0.0 ? ejecta / ejecta_km : 1.0;
    float core = core_km * (1.0 + (grown - 1.0) * saturate(age / 172800.0));
    sum += exp(-0.5 * (xw * xw + yw * yw) / (core * core));

    if (ejecta_km > 0.0) {
        float across = yw / (0.035 * ejecta_km);
        sum += 0.9 * smooth_step(-0.02 * ejecta_km, 0.0, xw) * (1.0 - smooth_step(0.22 * ejecta_km, 0.36 * ejecta_km, xw)) *
               exp(-0.5 * across * across);
        float x0 = xw - 0.115 * ejecta_km;
        float rd = sqrt(x0 * x0 + yw * yw);
        float az = atan2(yw, x0);
        float spread = az / 1.6;
        float rr = (rd / ejecta - 0.76) / (0.1 * max(1.0 - 0.14 * az * az, 0.3));
        float rays = 1.0 + 0.25 * sin(11.0 * az + seed + 2.0 * sin(3.0 * az + 1.3 * seed));
        sum += 0.7 * rays * exp(-spread * spread - 0.5 * rr * rr);
        float halo = r / ejecta;
        sum += 0.12 * exp(-3.0 * halo * halo);
    }
    sum *= landed;

    if (rings > 0 && age < 36000.0) {
        float fade = exp(-age / 10800.0) * smooth_step(120.0, 600.0, age);
        float outer = (r - (586.0 + 0.454 * age)) / 120.0;
        sum += 0.35 * fade * exp(-outer * outer);
        if (rings > 1) {
            float inner = (r - (586.0 + 0.189 * age)) / 90.0;
            sum += 0.2 * fade * exp(-inner * inner) * smooth_step(0.0, 0.5, -x / max(r, 1.0));
        }
    }
    return sum;
}

PSOutput main(PSInput input)
{
    float3 world_dir = normalize(input.rel_pos);
    float3 dir = mul((float3x3)u_to_local, world_dir);
    float3 o = u_camera.xyz;

    // The ray and the ellipsoid: a sphere once z is stretched.
    float stretch = u_radii.x / u_radii.y;
    float3 os = float3(o.xy, o.z * stretch);
    float3 ds = float3(dir.xy, dir.z * stretch);
    float a = dot(ds, ds);
    float b = dot(os, ds);
    float c = dot(os, os) - u_radii.x * u_radii.x;
    float disc = b * b - a * c;
    if (disc <= 0.0) {
        discard;
    }
    float t = (-b - sqrt(disc)) / a;
    if (t <= 0.0) {
        discard;
    }
    float3 q = o + t * dir;

    PSOutput output;
    float4 clip = mul(u_view_proj, float4(world_dir * t, 1.0));
    // A little nearer than the surface: at the planet mesh's vertices the two agree to rounding, and
    // pixels that lost to it would stay unmarked.
    output.depth = 1.0005 * clip.z / clip.w;

    float lat = asin(clamp(q.z / length(q), -1.0, 1.0));
    float lon = atan2(q.y, q.x);
    float cos_lat = cos(lat);
    float circle = kTwoPi * u_radii.w * cos_lat;

    float depth = 0.0;
    int count = (int)u_radii.z;
    for (int i = 0; i < count; ++i) {
        float north = (lat - u_site[i].y) * u_radii.w;
        // A pattern reaches no further than this north and south, and the sites share a latitude.
        if (abs(north) > 40000.0) {
            continue;
        }
        float d = lon - u_site[i].x;
        d -= kTwoPi * round(d / kTwoPi);
        float age = u_site[i].z;
        float wind = clamp(4.0e-6 * north, -0.010, 0.020);
        float east = d * u_radii.w * cos_lat - wind * age;
        east -= circle * round(east / circle);
        int rings = (int)(u_shape[i].w / 16.0);
        float seed = u_shape[i].w - 16.0 * rings;
        depth += u_site[i].w * scar_shape(east, north, age, u_shape[i].x, u_shape[i].y, u_shape[i].z, rings, seed);
    }

    // The band: the main part, and what has spread beyond it, at a fifth of the depth.
    if (u_band.x > 0.0) {
        float deg = lat / kDeg;
        float main_part = smooth_step(-48.0, -46.1, deg) * (1.0 - smooth_step(-39.4, -37.5, deg));
        float south = u_band.y / kDeg;
        float north_limit = u_band.z / kDeg;
        float spread = smooth_step(south - 2.0, south, deg) * (1.0 - smooth_step(north_limit, north_limit + 2.0, deg));
        depth += u_band.x * max(main_part, 0.2 * spread);
    }

    float3 n = normalize(float3(q.xy / (u_radii.x * u_radii.x), q.z / (u_radii.y * u_radii.y)));
    float mu = max(dot(-dir, n), 0.1);
    float mu0 = max(dot(u_sun.xyz, n), 0.1);
    output.color = float4(exp(-depth * u_tint.rgb * (1.0 / mu + 1.0 / mu0)), 1.0);
    return output;
}
