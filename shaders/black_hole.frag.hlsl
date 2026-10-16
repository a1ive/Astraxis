// Black hole prepass at a reduced resolution (rays are expensive): what each
// pixel sees (black_hole.hlsli), stored for black_hole_composite.frag.hlsl,
// which samples the sky per scene pixel along the interpolated directions
// and traces again where neighbouring rays disagree.

#include "black_hole.hlsli"

struct PSInput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

struct PSOutput
{
    float4 color : SV_Target0; // rgb = disk light * fade, a = sky transmission * fade (premultiplied: interpolates exactly)
    float4 sky   : SV_Target1; // xyz = ICRF sky direction (zero: hidden), w = fade, negative where opaque
};

// Outside the lensed region nothing is written: both targets stay cleared to zero.
PSOutput main(PSInput input)
{
    Trace t = trace_pixel(input.uv);
    if (!t.lensed) {
        discard;
    }
    PSOutput output;
    output.color = float4(t.disk * t.fade, t.transmit * t.fade);
    output.sky = float4(t.sky, t.opaque ? -t.fade : t.fade);
    return output;
}
