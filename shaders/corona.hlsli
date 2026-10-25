#ifndef ASTRAXIS_CORONA_HLSLI
#define ASTRAXIS_CORONA_HLSLI

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_ray_to_solar;
    float4 u_camera;
    float4 u_params; // gain, activity, outer radius, occulter count
    float4 u_color;
    float4 u_screen;
    float4 u_occluders[8];
    float4 u_radial[128];
};

float3 corona_ray(float2 uv)
{
    return normalize(mul(u_ray_to_solar, float4(uv * float2(2, -2) + float2(-1, 1), 1, 1)).xyz);
}

// Cross products avoid cancellation for the very small Moon seen from Earth.
float sphere_front(float3 center, float radius, float3 ray)
{
    float b = dot(center, ray);
    float3 side = cross(center, ray);
    float h = radius * radius - dot(side, side);
    if (h <= 0 || b <= 0) return 1e20;
    return max(0.0, b - sqrt(h));
}

// Density contrast on wrinkled sheets and polar plumes. All dimensions and
// amplitudes here are illustrative, not a magnetic-field reconstruction.
float corona_structure(float3 p, float r)
{
    float3 n = p / r;
    float a = u_params.y;
    float phi = atan2(n.y, n.x);
    float cusp = smoothstep(1.0, 3.0, r);
    float warp = (0.08 + 0.38 * a) * sin(2.0 * phi + 0.7) + 0.13 * a * sin(5.0 * phi - 0.4);
    float width = lerp(0.12, 0.025 + 0.015 * a, cusp);
    float sheet = exp(-pow((n.z - warp) / width, 2.0));
    float fan = exp(-pow((dot(n, float3(0.78, 0.36, 0.51)) - 0.16 * sin(3.0 * phi)) /
                         (width * 1.2), 2.0));
    float fine = 0.55 + 0.25 * sin(37.0 * phi + n.z * 9.0) + 0.20 * sin(83.0 * phi - n.z * 16.0);
    float plumes = pow(abs(n.z), 8.0) * (0.35 + 0.65 * fine * fine);
    float striae = 0.18 + 3.0 * pow(saturate(fine), 4.0);
    float envelope = 0.05 + 4.5 * sheet * striae + a * 3.0 * fan * striae;
    return (envelope + (1.0 - a) * 0.7 * plumes) * (1.0 - smoothstep(24.0, u_params.z, r));
}

float2 radial_scattering(float r)
{
    float f = sqrt(saturate(log(max(r, 1.0)) / log(u_params.z))) * 127.0;
    int i = min((int)f, 126);
    return lerp(u_radial[i].xy, u_radial[i + 1].xy, f - i);
}

float corona_integrate(float3 ray)
{
    float3 o = u_camera.xyz;
    float z_camera = dot(o, ray);
    float3 closest = o - z_camera * ray;
    float b2 = dot(closest, closest);
    float extent2 = u_params.z * u_params.z - b2;
    if (extent2 <= 0) return 0;
    float extent = sqrt(extent2);
    float begin = max(z_camera, -extent);
    float end = extent;
    float hit = sphere_front(-o, 1.0, ray);
    for (int j = 0; j < (int)u_params.w; ++j) {
        hit = min(hit, sphere_front(u_occluders[j].xyz, u_occluders[j].w, ray));
    }
    end = min(end, z_camera + hit);
    if (end <= begin) return 0;

    // z=b*sinh(t) concentrates samples near the point of closest approach.
    float b = max(sqrt(b2), 0.25);
    float x0 = begin / b;
    float x1 = end / b;
    float t0 = sign(x0) * log(abs(x0) + sqrt(x0 * x0 + 1.0));
    float t1 = sign(x1) * log(abs(x1) + sqrt(x1 * x1 + 1.0));
    const int steps = 96;
    float dt = (t1 - t0) / steps;
    float sum = 0;
    for (int i = 0; i < steps; ++i) {
        float t = t0 + (i + 0.5) * dt;
        float z = b * sinh(t);
        float3 p = closest + z * ray;
        float r = length(p);
        float2 c = radial_scattering(r);
        float mu = z / r;
        sum += (c.x + c.y * mu * mu) * corona_structure(p, r) * b * cosh(t) * dt;
    }
    return pow(max(0.0, sum * u_params.x), u_color.w);
}

// Reintegrate near opaque silhouettes instead of filtering light across them.
bool corona_edge(float3 ray)
{
    float pixel = 3.0 * length(ddx(ray)) + 3.0 * length(ddy(ray));
    float3 sun = -u_camera.xyz;
    if (dot(sun, ray) > 0 && abs(length(cross(sun, ray)) - 1.0) < pixel * length(sun)) return true;
    for (int j = 0; j < (int)u_params.w; ++j) {
        float4 o = u_occluders[j];
        if (dot(o.xyz, ray) > 0 && abs(length(cross(o.xyz, ray)) - o.w) < pixel * length(o.xyz)) return true;
    }
    return false;
}

#endif
