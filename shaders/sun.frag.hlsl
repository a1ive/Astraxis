// Sun disk with limb darkening plus a faint corona-like falloff. Values are
// HDR; the bloom pass turns the very bright disk into glare.

cbuffer Uniforms : register(b0, space3)
{
    float4 u_params; // x = disk radius in pixels, y = disk intensity, z = glow intensity
    float4 u_color;  // rgb = linear color
};

struct PSInput
{
    float4 position : SV_Position;
    float2 offset_px : TEXCOORD0;
};

float4 main(PSInput input) : SV_Target0
{
    float r = length(input.offset_px) / u_params.x; // in units of the disk radius
    float disk = 0.0;
    if (r < 1.0) {
        float mu = sqrt(1.0 - r * r);
        disk = 1.0 - 0.6 * (1.0 - mu); // linear limb darkening
    }
    // Anti-aliased edge for small disks.
    disk *= saturate((1.0 - r) * u_params.x + 0.5);

    float glow = u_params.z / (1.0 + r * r * 4.0) * saturate(1.0 - r / 12.0);
    return float4(u_color.rgb * (disk * u_params.y + glow), 1.0);
}
