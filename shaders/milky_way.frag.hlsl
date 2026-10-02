// Milky Way background: an equirectangular map in ICRF right ascension /
// declination (RA 0h at the centre, increasing to the left as seen from
// inside the sphere), looked up per pixel and added under the stars.

static const float kPi = 3.14159265;

Texture2D u_map : register(t0, space2); // sRGB texture: samples are linear
SamplerState u_sampler : register(s0, space2);

cbuffer Uniforms : register(b0, space3)
{
    float4x4 u_inv_view_proj; // inverse of the rotation-only (ICRF) view * projection
    float4 u_params;          // x = brightness, y = mip level
};

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    // A point on the near plane (reversed-Z: ndc z = 1) along this pixel's ray.
    float2 ndc = float2(input.uv.x * 2.0 - 1.0, 1.0 - input.uv.y * 2.0);
    float4 p = mul(u_inv_view_proj, float4(ndc, 1.0, 1.0));
    float3 d = normalize(p.xyz / p.w);

    float ra = atan2(d.y, d.x);
    float dec = asin(clamp(d.z, -1.0, 1.0));
    float2 uv = float2(frac(0.5 - ra / (2.0 * kPi)), 0.5 - dec / kPi);
    float3 c = u_map.SampleLevel(u_sampler, uv, u_params.y).rgb;
    return float4(c * u_params.x, 1.0);
}
