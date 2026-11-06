// A comet fragment's impact: the flash and the ballistic plume, as src/scene/impact.cpp has them (change
// them together). In the local frame (km): origin at the site on the surface, x east, y north, z up.
//
// The plume is a swarm thrown from the origin into a cone, then falling freely, so what is at p at time t
// left with velocity v = (p + g t^2 / 2 z) / t and the density is f(v) / t^3: nothing is stored. Radiance is
// in I/F units, single scattering; it is lit above the planet's shadow and glows for a while in it.
// The flash is a Gaussian around the fireball, drawn over at least a few pixels, integrated in closed form.

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_local;  // camera-relative world vector -> local frame (3x3 rotation)
    float4 u_camera;      // xyz = camera (local), w = 1 to start the rays there
    float4 u_sun;         // xyz = unit direction to the sun (local), w = its angular radius
    float4 u_box_min;
    float4 u_box_max;
    float4 u_planet;      // xyz = planet center (local), w = radius (km)
    float4 u_params;      // x = age (s), y = extinction per unit of density (km^2)
    float4 u_flash;       // x = height of the fireball, y = sigma (km), z = radiance at its center, w = unhidden
    float4 u_flash_color; // rgb, by its temperature
    float4 u_glow;        // rgb = radiance per unit of optical depth in the shadow
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 local    : TEXCOORD1;
};

static const int kSteps = 64;
static const float kVMin = 4.5;
static const float kVMax = 11.81;
static const float kAlpha = 1.55;
static const float kGravity = 0.02325;
static const float kPi = 3.14159265;
static const float kHalfCone = 0.65449847; // 37.5 deg
static const float kSoft = 0.10471976;     // 6 deg
static const float3 kAxis = float3(0.28678822, -0.40957602, 0.8660254);
static const float kVanguardWidth = 0.08;
static const float kVanguardShare = 0.22425;
// kPlumeAlbedo in impact.hpp; slightly red, the color illustrative.
static const float3 kAlbedo = float3(0.74, 0.70, 0.64);

float jitter_noise(float2 p)
{
    return frac(52.9829189 * frac(dot(p, float2(0.06711056, 0.00583715))));
}

float density(float3 p, float age)
{
    float3 v = (p + float3(0.0, 0.0, 0.5 * kGravity * age * age)) / age;
    float speed = length(v);
    if (speed < 0.9 * kVMin || speed > 1.25 * kVMax) {
        return 0.0;
    }
    float cone = smoothstep(cos(kHalfCone + kSoft), cos(kHalfCone - kSoft), dot(v / speed, kAxis));
    if (cone <= 0.0) {
        return 0.0;
    }
    float solid = 2.0 * kPi * (1.0 - cos(kHalfCone));
    float f = 0.0;
    if (speed >= kVMin && speed <= kVMax) {
        f = pow(kVMin, kAlpha) * kAlpha * pow(speed, -kAlpha - 3.0) / solid;
    }
    float width = kVanguardWidth * kVMax;
    float d = (speed - kVMax) / width;
    f += kVanguardShare * exp(-0.5 * d * d) / (solid * kVMax * kVMax * width * 2.5066283);
    return cone * f / (age * age * age);
}

float sunlit(float3 p)
{
    float3 oc = u_planet.xyz - p;
    float along = dot(oc, u_sun.xyz);
    if (along <= 0.0) {
        return 1.0;
    }
    float perp = length(oc - along * u_sun.xyz);
    float penumbra = max(along * u_sun.w, 1e-3 * u_planet.w);
    return smoothstep(u_planet.w - penumbra, u_planet.w + penumbra, perp);
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_local, input.rel_pos));
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.local;

    float3 inv = 1.0 / select(abs(dir) > 1e-8, dir, 1e-8);
    float3 ta = (u_box_min.xyz - o) * inv;
    float3 tb = (u_box_max.xyz - o) * inv;
    float3 tn = min(ta, tb);
    float3 tf = max(ta, tb);
    float t0 = max(max(max(tn.x, tn.y), tn.z), 0.0);
    float t1 = min(min(tf.x, tf.y), tf.z);

    // The planet stops the ray.
    float3 oc = o - u_planet.xyz;
    float b = dot(oc, dir);
    float3 perp = oc - b * dir;
    float p2 = dot(perp, perp);
    float r2 = u_planet.w * u_planet.w;
    if (p2 < r2) {
        float ground = -b - sqrt(r2 - p2);
        if (ground > 0.0) {
            t1 = min(t1, ground);
        }
    }

    float3 light = float3(0.0, 0.0, 0.0);

    // The flash. Whether the planet hides it is decided for its center: it is usually drawn far wider than
    // it is, and cutting that by the limb would cut it in half.
    if (u_flash.z > 0.0) {
        float3 center = float3(0.0, 0.0, u_flash.x);
        float along = dot(center - o, dir);
        float3 miss = (o + along * dir) - center;
        float sigma2 = u_flash.y * u_flash.y;
        light += u_flash_color.rgb * u_flash.z * u_flash.w * exp(-0.5 * dot(miss, miss) / sigma2);
    }
    if (t1 <= t0) {
        return float4(light, 1.0);
    }

    float tau = 0.0;
    float age = u_params.x;
    if (u_params.y > 0.0 && age > 0.0) {
        float cos_angle = dot(dir, u_sun.xyz);
        // Henyey-Greenstein, a mean of 1. Small grains, nearly isotropic (illustrative).
        const float g = 0.2;
        float phase = (1.0 - g * g) / pow(1.0 + g * g - 2.0 * g * cos_angle, 1.5);
        float dt = (t1 - t0) / kSteps;
        float jitter = jitter_noise(input.position.xy);
        float scattered = 0.0;
        float glow = 0.0;
        for (int i = 0; i < kSteps; ++i) {
            float3 p = o + (t0 + (i + jitter) * dt) * dir;
            float altitude = p.z + dot(p.xy, p.xy) / (2.0 * u_planet.w);
            if (altitude < 0.0) {
                continue;
            }
            float sigma = u_params.y * density(float3(p.xy, altitude), age);
            if (sigma > 0.0) {
                float through = exp(-tau - 0.5 * sigma * dt) * sigma * dt;
                scattered += through * sunlit(p);
                glow += through;
                tau += sigma * dt;
            }
        }
        light += 0.25 * phase * kAlbedo * scattered + u_glow.rgb * glow;
    }
    return float4(light, exp(-tau));
}
