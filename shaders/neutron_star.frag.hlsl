// Backward ray map for a static observer outside a spherical opaque star.
// Only the two angles from the radial tables depend on the Schwarzschild
// metric. Each ray stays in the plane through the camera and star centre.
#include "common.hlsli"

FRAGMENT_TEXTURE(TextureCube, u_sky, u_sky_sampler, 0);
FRAGMENT_TEXTURE(Texture2D, u_milky_way, u_milky_way_sampler, 1);
ByteAddressBuffer u_star_data : register(t2, space2);
ByteAddressBuffer u_star_cells : register(t3, space2);

static const float kPi = 3.14159265;
static const int kSamples = 129;

cbuffer Uniforms : register(b0, space3)
{
    float4 u_right;
    float4 u_up;
    float4 u_forward;
    float4x4 u_to_sky;
    float4 u_observer;
    float4 u_geometry;
    float4 u_surface;
    float4 u_hot_spot;
    float4 u_stars;
    float4 u_sky_params;
    float4 u_bending[kSamples];
};

#include "lensed_sky.hlsli"

struct PSInput
{
    float4 position : SV_Position;
    float2 uv : TEXCOORD0;
};

struct PSOutput
{
    float4 color : SV_Target0;
    float depth : SV_Depth;
};

float2 bending(float fraction)
{
    float index = saturate(fraction) * (kSamples - 1);
    int lo = min((int)index, kSamples - 2);
    return lerp(u_bending[lo].xy, u_bending[lo + 1].xy, index - lo);
}

PSOutput main(PSInput input)
{
    float2 ndc = float2(input.uv.x * 2.0 - 1.0, 1.0 - input.uv.y * 2.0);
    float3 screen_ray = u_forward.xyz + ndc.x * u_right.xyz + ndc.y * u_up.xyz;
    float len2 = dot(screen_ray, screen_ray);
    float3 ray = screen_ray * rsqrt(len2);
    float3 radial = u_observer.xyz;
    float cos_angle = -dot(radial, ray);
    float3 across = cross(radial, ray);
    float sin_angle = length(across);
    float3 axis = across / max(sin_angle, 1e-10);
    float angle = atan2(sin_angle, cos_angle);
    float limb = u_geometry.x;
    float limb_width = max(fwidth(angle), 1e-6);
    float cover = saturate(0.5 + (limb - angle) / limb_width);

    // Evaluate the remapping in all lanes before derivatives. At the limb,
    // extrapolate the tangent ray over the small antialiasing strip.
    float sky_angle = sqrt(saturate(log(max(angle, limb) / limb) / log(kPi / limb)));
    float turn = bending(sky_angle).y;
    float3 sky_dir = mul((float3x3)u_to_sky, cos(turn) * ray + sin(turn) * cross(axis, ray));
    float3 gx = ddx_fine(sky_dir);
    float3 gy = ddy_fine(sky_dir);

    // Fade the lensed sky out once its displacement is subpixel. The limb
    // always stays opaque, even when the normal starfield is hidden.
    float sky_alpha = u_surface.w * smoothstep(0.35, 1.0, turn / u_sky_params.z);
    if (cover == 0.0 && sky_alpha == 0.0) {
        discard;
    }

    float3 sky = 0.0;
    if (cover < 1.0 && sky_alpha > 0.0 && u_stars.w > 0.0) {
        float pixel = u_geometry.w / (len2 * sqrt(len2));
        sky = u_stars.w * sky_alpha * sky_light(sky_dir, gx, gy, pixel);
    }
    float3 surface = 0.0;
    float depth = 0.0;
    if (cover > 0.0) {
        float b = saturate(sin_angle / u_geometry.z);
        float mu = sqrt(1.0 - b * b);
        float theta = bending(mu).x;
        float3 normal = cos(theta) * radial + sin(theta) * cross(axis, radial);
        float pole = abs(dot(normal, u_hot_spot.xyz));
        float hot = smoothstep(0.96126, 0.99144, pole) * u_hot_spot.w; // 16 to 7.5 degrees
        surface = u_surface.rgb * (0.4 + 0.6 * mu) * (1.0 + 2.8 * hot);
        depth = saturate(u_geometry.y / dot(u_forward.xyz, normal - u_observer.w * radial));
    }

    PSOutput output;
    output.color = float4(lerp(sky, surface, cover), lerp(sky_alpha, 1.0, cover));
    output.depth = depth;
    return output;
}
