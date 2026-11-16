// Habitable zone (scene/habitable_zone.hpp): at each point of the plane, every
// limit's insolation over S_eff summed over the stars, sum_i w_i / r_i^2 with
// w = L / S_eff (au^2). The optimistic zone has recent Venus < 1 < early Mars,
// the conservative one runaway greenhouse < 1 < maximum greenhouse. Edges are
// soft in the logarithm of that sum (at least a pixel and a half wide), and the
// light is added, so stars and orbits behind it stay as they are.

cbuffer Uniforms : register(b0, space3)
{
    float4 u_stars[4];     // xyz = star position in the plane frame (au)
    float4 u_weights[4];   // L / S_eff (au^2): x = recent Venus, y = runaway greenhouse,
                           // z = maximum greenhouse, w = early Mars
    float4 u_conservative; // rgb = linear radiance of the conservative zone
    float4 u_optimistic;   // rgb = of the optimistic one beyond it
    float4 u_params;       // x = number of stars, y, z = camera distance over which it fades in (km),
                           // w = edge softness (natural log of the sum)
};

struct PSInput
{
    float4 position : SV_Position;
    float3 local    : TEXCOORD0;
    float3 view     : TEXCOORD1;
};

float4 main(PSInput input) : SV_Target0
{
    float4 sum = 0.0;
    [unroll]
    for (int i = 0; i < 4; ++i) {
        if (i < (int)u_params.x) {
            float3 d = input.local - u_stars[i].xyz;
            sum += u_weights[i] / max(dot(d, d), 1e-12);
        }
    }
    float4 ln_s = log(max(sum, 1e-30));
    float4 soft = max(u_params.w, 1.5 * fwidth(ln_s));
    float optimistic = smoothstep(-soft.x, soft.x, -ln_s.x) * smoothstep(-soft.w, soft.w, ln_s.w);
    float conservative = smoothstep(-soft.y, soft.y, -ln_s.y) * smoothstep(-soft.z, soft.z, ln_s.z);
    // Close to the camera the plane would wash over the view (a planet in the zone).
    float near = smoothstep(u_params.y, u_params.z, length(input.view));
    float3 color = lerp(u_optimistic.rgb, u_conservative.rgb, conservative) * optimistic;
    return float4(color * near, 0.0);
}
