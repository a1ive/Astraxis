#ifndef ASTRAXIS_DISK_HLSLI
#define ASTRAXIS_DISK_HLSLI

// An optically thin disk's light along each ray: see scene/disk.hpp and
// DiskUniforms. Lengths in au, ICRF axes, relative to the star.

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_ray_to_icrf;
    float4 u_camera;    // w: view (0 scattered, 1 infrared, 2 millimeter)
    float4 u_params;    // radiance at the reference, gamma, reference, components
    float4 u_range;     // ln inner, ln outer, outer, occluders
    float4 u_phase_g;
    float4 u_phase_w;
    float4 u_tint;
    float4 u_axis_x[4]; // w: sigma
    float4 u_axis_y[4]; // w: second sigma
    float4 u_axis_z[4]; // w: second Gaussian's share
    float4 u_color[4];  // w: flare radius
    float4 u_temperature[4]; // K, radius, index, unused
    float4 u_occluders[8];   // relative to the camera; w: radius
    float4 u_thermal[16];    // false color, even in ln T; w: ln of the infrared emission
    float4 u_table[128];     // each component's column, per radius
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
    return sqrt(sqrt(1.0 + x * x * x * x));
}

float disk_vertical(int c, float r, float z)
{
    float flare = disk_flare(c, r);
    float s1 = u_axis_x[c].w * flare;
    float a1 = z / s1;
    float out_ = (1.0 - u_axis_z[c].w) * exp(-0.5 * a1 * a1) / s1;
    if (u_axis_z[c].w > 0.0) {
        float s2 = u_axis_y[c].w * flare;
        float a2 = z / s2;
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
            float value = light * disk_vertical(c, r, q.z) * dt;
            float distance = length(q);
            float3 rgb;
            if (view == 0) {
                value *= disk_phase(-dot(q, d) / distance) / (distance * distance);
                rgb = u_color[c].rgb;
            } else if (view == 1) {
                float4 thermal = disk_thermal(disk_ln_temperature(c, distance));
                value *= thermal.w;
                rgb = thermal.rgb;
            } else {
                value *= exp(disk_ln_temperature(c, distance));
                rgb = u_tint.rgb;
            }
            sum += float4(rgb * value, value);
        }
        t += dt;
    }
    return sum;
}

float3 disk_integrate(float3 ray)
{
    float3 origin = u_camera.xyz;
    float along = -dot(origin, ray);
    float3 closest = origin + along * ray;
    float extent2 = u_range.z * u_range.z - dot(closest, closest);
    if (extent2 <= 0) return 0;
    float extent = sqrt(extent2);
    float t0 = max(0.0, along - extent);
    float t1 = along + extent;
    for (int j = 0; j < (int)u_range.w; ++j) {
        t1 = min(t1, disk_sphere_front(u_occluders[j].xyz, u_occluders[j].w, ray));
    }
    if (t1 <= t0) return 0;

    float4 sum = 0;
    for (int c = 0; c < (int)u_params.w; ++c) {
        sum += disk_march(c, origin, ray, t0, t1);
    }
    if (sum.w <= 0) return 0;
    // Each component's color, with the range of brightness compressed.
    return sum.rgb / sum.w * u_params.x * pow(sum.w / u_params.z, u_params.y);
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
