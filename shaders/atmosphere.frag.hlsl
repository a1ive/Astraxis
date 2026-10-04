// Atmosphere by single scattering along the view ray through the shell, in
// the planet's unit-sphere frame (an ellipsoidal planet scaled to the unit
// sphere; lengths there are in planet radii). Gas scatters with the Rayleigh
// phase function, haze with Henyey-Greenstein; both have exponential density
// profiles from the base (the surface, or the top of an opaque deck). Light
// that haze scatters into its narrow forward peak mostly keeps going, so haze
// attenuates with the scaled optical depth tau (1 - w g) (similarity relation),
// though it scatters with the full tau w. Sunlight reaching each sample is
// attenuated by the column toward the sun (Chapman approximation) and cut off
// by the planet's shadow, which gives the twilight arc and the bright ring of a
// backlit planet. Like the rings, radiance is in I/F units (phase functions
// with a mean of 1, hence the 1/4):
//   I/F = 1/4 * integral of (sigma_s P)(s) * exp(-tau_view(s) - tau_sun(s)) ds
// Alpha is the view transmission (mean over the channels), which dims what
// lies behind (ground, stars).

#include "atmosphere.hlsli"

static const int kSteps = 40;

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_to_unit;   // camera-relative world vector -> unit-sphere frame (3x3)
    float4 u_camera;      // xyz = camera in the unit-sphere frame, w = 1 to start the rays there
    float4 u_sun;         // xyz = unit direction to the sun (world), w = its angular radius (rad)
    float4 u_sun_unit;    // xyz = unit direction to the sun (unit-sphere frame), w = shell top radius
    float4 u_rayleigh;    // rgb = vertical optical depth, w = scale height (planet radii)
    float4 u_haze;        // rgb = vertical attenuation optical depth, tau (1 - w g), w = scale height
    float4 u_haze_scatter; // rgb = vertical scattering optical depth, tau w; w = asymmetry g
    float4 u_params;      // x = gain
};

struct PSInput
{
    float4 position : SV_Position;
    float3 rel_pos  : TEXCOORD0;
    float3 unit_pos : TEXCOORD1;
};

// Interleaved gradient noise (Jimenez 2014): per-pixel jitter of the samples.
float gradient_noise(float2 p)
{
    return frac(52.9829189 * frac(dot(p, float2(0.06711056, 0.00583715))));
}

float4 main(PSInput input) : SV_Target0
{
    float3 dir = normalize(mul((float3x3)u_to_unit, input.rel_pos));
    // Outside the shell the ray starts at the shell's surface (no precision
    // lost to a distant camera); inside, at the camera.
    float3 o = u_camera.w > 0.5 ? u_camera.xyz : input.unit_pos;

    float top = u_sun_unit.w;
    float b = dot(o, dir);
    float3 perp = o - b * dir;
    float p2 = dot(perp, perp);
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
    if (t1 <= t0) {
        return float4(0.0, 0.0, 0.0, 1.0);
    }

    float cos_theta = dot(normalize(input.rel_pos), u_sun.xyz);
    float phase_r = rayleigh_phase(cos_theta);
    float phase_h = henyey_greenstein_phase(cos_theta, u_haze_scatter.w);
    float3 sigma_r = u_rayleigh.rgb / u_rayleigh.w; // at the base, per planet radius
    float3 sigma_h = u_haze.rgb / u_haze.w;
    float3 sigma_hs = u_haze_scatter.rgb / u_haze.w;

    float dt = (t1 - t0) / kSteps;
    float jitter = gradient_noise(input.position.xy);
    float3 tau = 0.0;
    float3 sum = 0.0;
    for (int i = 0; i < kSteps; ++i) {
        float3 p = o + (t0 + (i + jitter) * dt) * dir;
        float r = length(p);
        float mu = dot(p, u_sun_unit.xyz) / r;
        float density_r = exp(-(r - 1.0) / u_rayleigh.w);
        float density_h = exp(-(r - 1.0) / u_haze.w);
        float3 extinction = sigma_r * density_r + sigma_h * density_h;
        float lit = ray_lit(r, mu, u_sun.w);
        if (lit > 0.0) {
            float3 t = exp(-(tau + 0.5 * dt * extinction) - sun_optical_depth(r, mu, u_rayleigh, u_haze));
            float3 scattering = sigma_r * density_r * phase_r + sigma_hs * density_h * phase_h;
            sum += scattering * t * (lit * dt);
        }
        tau += extinction * dt;
    }

    float3 transmission = exp(-tau);
    return float4(0.25 * u_params.x * sum, dot(transmission, 1.0 / 3.0));
}
