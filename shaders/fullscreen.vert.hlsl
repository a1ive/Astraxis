// Fullscreen triangle (3 vertices, no vertex buffer).

struct VSOutput
{
    float4 position : SV_Position;
    float2 uv       : TEXCOORD0;
};

VSOutput main(uint vertex_id : SV_VertexID)
{
    VSOutput output;
    float2 p = float2((vertex_id << 1) & 2, vertex_id & 2); // (0,0) (2,0) (0,2)
    output.position = float4(p * float2(2.0, -2.0) + float2(-1.0, 1.0), 0.0, 1.0);
    output.uv = p;
    return output;
}
