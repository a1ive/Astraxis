#ifndef ASTRAXIS_DISK_HLSLI
#define ASTRAXIS_DISK_HLSLI

// A disk's light along each ray: see scene/disk.hpp and DiskUniforms. Lengths in
// au, ICRF axes, relative to the star. The including shader declares the table of
// optical depth from the star, u_star_depth with u_depth_sampler.

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_ray_to_icrf;
    float4 u_camera;    // w: view (0 scattered, 1 infrared, 2 millimeter, 3 gas, 4 scattered in a thick disk)
    float4 u_params;    // radiance at the reference, gamma, reference, components
    float4 u_range;     // ln inner, ln outer, outer, occluders
    float4 u_phase_g;   // w: a source's least size (radians)
    float4 u_phase_w;   // w: a thick disk's gain on the first scattering
    float4 u_tint;      // w: sources
    float4 u_gas;       // sqrt(GM) in km/s sqrt(au), full color's speed, clumps
    float4 u_display;   // x: 1 / display_reference_au^2 in an r^2-scaled view (0: 1 / d^2)
    float4 u_clump[2];  // cos, sin of its azimuth, radius, radial sigma
    float4 u_clump_shape[2]; // sigma ahead, tail behind (radians), gain
    float4 u_axis_x[4]; // w: sigma
    float4 u_axis_y[4]; // w: second sigma
    float4 u_axis_z[4]; // w: second Gaussian's share
    float4 u_color[4];  // w: flare radius
    float4 u_temperature[4]; // K, radius, index, flare index
    float4 u_asymmetry[4];   // start, full, broad-layer warp/r, contrast
    float4 u_occluders[8];   // relative to the camera; w: radius
    float4 u_sources[4];     // relative to the camera; w: sigma
    float4 u_source_light[4];
    float4 u_thermal[16];    // false color, even in ln T; w: ln of the infrared emission
    float4 u_table[128];     // each component's column per radius, or its vertical optical depth
};

static const float kThermalLnLo = 3.4011974;  // ln 30 K (kDiskThermalLoK)
static const float kThermalLnHi = 6.9077553;  // ln 1000 K (kDiskThermalHiK)

float3 disk_ray(float2 uv)
{
    return normalize(mul(u_ray_to_icrf, float4(uv * float2(2, -2) + float2(-1, 1), 1, 1)).xyz);
}

float disk_sphere_front(float3 center, float radius, float3 ray)
{
    float b = dot(center, ray);
    float3 side = cross(center, ray);
    float h = radius * radius - dot(side, side);
    if (h <= 0 || b <= 0) return 1e20;
    return max(0.0, b - sqrt(h));
}

float disk_table(int c, float r)
{
    float f = (log(max(r, 1e-6)) - u_range.x) / (u_range.y - u_range.x) * 127.0;
    if (f <= 0.0 || f >= 127.0) return 0.0;
    int i = min((int)f, 126);
    return lerp(u_table[i][c], u_table[i + 1][c], f - i);
}

float disk_flare(int c, float r)
{
    float x = r / u_color[c].w;
    return pow(1.0 + x * x * x * x, 0.25 * u_temperature[c].w);
}

float disk_vertical(int c, float r, float z, float broad_shift)
{
    float flare = disk_flare(c, r);
    float s1 = u_axis_x[c].w * flare;
    float a1 = z / s1;
    float out_ = (1.0 - u_axis_z[c].w) * exp(-0.5 * a1 * a1) / s1;
    if (u_axis_z[c].w > 0.0) {
        float s2 = u_axis_y[c].w * flare;
        float a2 = (z - broad_shift) / s2;
        out_ += u_axis_z[c].w * exp(-0.5 * a2 * a2) / s2;
    }
    return 0.3989423 * out_;
}

float disk_hg(float g, float cos_theta)
{
    float d = 1.0 + g * g - 2.0 * g * cos_theta;
    return (1.0 - g * g) / (d * sqrt(d));
}

float disk_phase(float cos_theta)
{
    return u_phase_w.x * disk_hg(u_phase_g.x, cos_theta) + u_phase_w.y * disk_hg(u_phase_g.y, cos_theta) +
           u_phase_w.z * disk_hg(u_phase_g.z, cos_theta);
}

// ln of the grains' temperature d from the star.
float disk_ln_temperature(int c, float d)
{
    return log(u_temperature[c].x) - u_temperature[c].z * (log(max(d, 1e-3)) - log(u_temperature[c].y));
}

// False color; w: the black-body light in the infrared band.
float4 disk_thermal(float ln_t)
{
    float f = saturate((ln_t - kThermalLnLo) / (kThermalLnHi - kThermalLnLo)) * 15.0;
    int i = min((int)f, 14);
    float4 v = lerp(u_thermal[i], u_thermal[i + 1], f - i);
    return float4(v.rgb, exp(v.w));
}

// disk_clump_factor in scene/disk.cpp, at q in the component's plane.
float disk_clumps(float3 q, float r)
{
    float out_ = 1.0;
    for (int s = 0; s < (int)u_gas.z; ++s) {
        float a = (r - u_clump[s].z) / u_clump[s].w;
        float delta = atan2(u_clump[s].x * q.y - u_clump[s].y * q.x, u_clump[s].x * q.x + u_clump[s].y * q.y);
        float ahead = delta / u_clump_shape[s].x;
        float shape = delta >= 0.0 ? exp(-0.5 * ahead * ahead) : exp(delta / u_clump_shape[s].y);
        out_ += u_clump_shape[s].z * exp(-0.5 * a * a) * shape;
    }
    return out_;
}

// Blue toward us, red away, white at rest: the speed along the ray (positive away)
// in units of full color. disk_line_of_sight_km_s in scene/disk.cpp.
float3 disk_gas_color(float3 q, float3 d, float r)
{
    float speed = u_gas.x / sqrt(max(r, 1e-3));
    float along = speed * (-q.y * d.x + q.x * d.y) / max(r, 1e-3) / u_gas.y;
    float3 toward = along < 0.0 ? float3(0.10, 0.28, 1.0) : float3(1.0, 0.14, 0.06);
    return lerp(float3(1.0, 1.0, 1.0), toward, saturate(abs(along)));
}

// The steps of march() in scene/disk.cpp: rgb, and the light they sum to in w.
float4 disk_march(int c, float3 origin, float3 ray, float t0, float t1)
{
    float3 o = float3(dot(origin, u_axis_x[c].xyz), dot(origin, u_axis_y[c].xyz), dot(origin, u_axis_z[c].xyz));
    float3 d = float3(dot(ray, u_axis_x[c].xyz), dot(ray, u_axis_y[c].xyz), dot(ray, u_axis_z[c].xyz));
    float4 sum = 0;
    float t = t0;
    int view = (int)u_camera.w;
    for (int i = 0; i < 48 && t < t1; ++i) {
        float3 p = o + t * d;
        float sigma = u_axis_x[c].w * disk_flare(c, length(p.xy));
        float dt = min(0.5 * max(abs(p.z), sigma) / max(abs(d.z), 1e-4), 0.3 * max(length(p), 1.0));
        dt = min(dt, t1 - t);
        float3 q = o + (t + 0.5 * dt) * d;
        float r = length(q.xy);
        float light = disk_table(c, r);
        if (light > 0.0) {
            float4 asymmetry = u_asymmetry[c];
            float skew = view == 0 && any(asymmetry.zw != 0.0)
                ? smoothstep(asymmetry.x, asymmetry.y, r) * q.x : 0.0;
            float value = light * disk_vertical(c, r, q.z, asymmetry.z * skew) * dt;
            value *= 1.0 + asymmetry.w * skew / max(r, 1e-6);
            float distance = length(q);
            float3 rgb;
            if (view == 0) {
                value *= disk_phase(-dot(q, d) / distance) * (u_display.x > 0.0 ? u_display.x : 1.0 / (distance * distance));
                rgb = u_color[c].rgb;
            } else if (view == 1) {
                float4 thermal = disk_thermal(disk_ln_temperature(c, distance));
                value *= thermal.w;
                rgb = thermal.rgb;
            } else if (view == 2) {
                value *= exp(disk_ln_temperature(c, distance));
                rgb = u_tint.rgb;
            } else {
                value *= disk_clumps(q, r);
                rgb = disk_gas_color(q, d, r);
            }
            sum += float4(rgb * value, value);
        }
        t += dt;
    }
    return sum;
}

// The extinction of a thick disk's components, per au.
float disk_extinction(float r, float z)
{
    float out_ = 0.0;
    for (int c = 0; c < (int)u_params.w; ++c) {
        out_ += disk_table(c, r) * disk_vertical(c, r, z, 0.0);
    }
    return out_;
}

// What a sample of a thick disk sends on, without its 1 / d^2 and phase function:
// x the starlight that reaches it (disk_star_depth in scene/disk.cpp), y the source
// of the later scatterings (disk_diffuse).
float2 disk_star_light(float r, float z)
{
    float u = saturate((log(max(r, 1e-6)) - u_range.x) / (u_range.y - u_range.x));
    float v = sqrt(atan2(abs(z), max(r, 1e-6)) / 1.5707963);
    float2 uv = (float2(u, v) * float2(255.0, 127.0) + 0.5) / float2(256.0, 128.0);
    float2 table = u_star_depth.SampleLevel(u_depth_sampler, uv, 0).rg;
    return float2(exp(1.0 - exp(table.r)), exp(table.g));
}

// The steps of march_thick() in scene/disk.cpp: the light, and what gets through in w.
float2 disk_march_thick(float3 origin, float3 ray, float t0, float t1)
{
    float3 o = float3(dot(origin, u_axis_x[0].xyz), dot(origin, u_axis_y[0].xyz), dot(origin, u_axis_z[0].xyz));
    float3 d = float3(dot(ray, u_axis_x[0].xyz), dot(ray, u_axis_y[0].xyz), dot(ray, u_axis_z[0].xyz));
    float light = 0.0;
    float through = 1.0;
    float t = t0;
    for (int i = 0; i < 96 && t < t1 && through > 0.005; ++i) {
        float3 p = o + t * d;
        float r = length(p.xy);
        float sigma = u_axis_x[0].w * disk_flare(0, max(r, 0.5 * abs(p.z)));
        float height = max(abs(p.z), sigma);
        float climb = max(abs(d.z), 1e-3);
        float extinction = disk_extinction(r, p.z);
        float to_star = max(length(p), 1.0);
        float dt;
        if (extinction * sigma * 2.5 > 1e-4) {
            dt = min(min(0.25 * sigma * min(1.0, sigma / height) / climb, 0.15 / extinction), 0.06 * to_star);
        } else {
            dt = min(max(0.5 * sigma, abs(p.z) - 5.0 * sigma) / climb, 0.1 * to_star);
        }
        dt = min(dt, t1 - t);
        float3 q = o + (t + 0.5 * dt) * d;
        float rq = length(q.xy);
        float depth = disk_extinction(rq, q.z) * dt;
        if (depth > 0.0) {
            float distance = length(q);
            float pass_ = exp(-depth);
            float2 star = disk_star_light(rq, q.z);
            light += through * (1.0 - pass_) * (u_phase_w.w * disk_phase(-dot(q, d) / distance) * star.x + star.y) /
                     (distance * distance);
            through *= pass_;
        }
        t += dt;
    }
    return float2(light, through);
}

// The disk's light, and in w what it lets through.
float4 disk_integrate(float3 ray)
{
    float3 origin = u_camera.xyz;
    float along = -dot(origin, ray);
    float3 closest = origin + along * ray;
    float extent2 = u_range.z * u_range.z - dot(closest, closest);
    float3 rgb = 0;
    float through = 1.0;
    // The small things shining in the disk.
    for (int s = 0; s < (int)u_tint.w; ++s) {
        float along_source = dot(u_sources[s].xyz, ray);
        if (along_source > 0.0) {
            float3 side = cross(u_sources[s].xyz, ray);
            float size = max(u_sources[s].w, u_phase_g.w * along_source);
            rgb += u_source_light[s].rgb * exp(-0.5 * dot(side, side) / (size * size));
        }
    }
    if (extent2 <= 0) return float4(rgb, through);
    float extent = sqrt(extent2);
    float t0 = max(0.0, along - extent);
    float t1 = along + extent;
    for (int j = 0; j < (int)u_range.w; ++j) {
        t1 = min(t1, disk_sphere_front(u_occluders[j].xyz, u_occluders[j].w, ray));
    }
    if (t1 <= t0) return float4(rgb, through);

    float4 sum = 0;
    if ((int)u_camera.w == 4) {
        float2 thick = disk_march_thick(origin, ray, t0, t1);
        sum = float4(u_color[0].rgb * thick.x, thick.x);
        through = thick.y;
    } else {
        for (int c = 0; c < (int)u_params.w; ++c) {
            sum += disk_march(c, origin, ray, t0, t1);
        }
    }
    // Each component's color, with the range of brightness compressed.
    if (sum.w > 0) rgb += sum.rgb / sum.w * u_params.x * pow(sum.w / u_params.z, u_params.y);
    return float4(rgb, through);
}

// Integrate again near opaque silhouettes instead of filtering light across them.
bool disk_edge(float3 ray)
{
    float pixel = 3.0 * length(ddx(ray)) + 3.0 * length(ddy(ray));
    for (int j = 0; j < (int)u_range.w; ++j) {
        float4 o = u_occluders[j];
        float d = length(o.xyz);
        if (dot(o.xyz, ray) > 0 && abs(length(cross(o.xyz, ray)) - o.w) < pixel * d) return true;
    }
    return false;
}

#endif
