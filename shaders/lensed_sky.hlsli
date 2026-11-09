// Shared sky lookup. The lens shader supplies the sky textures, star grid,
// brightness and star-size limits in its own bindings and uniforms.

// The Milky Way map along dir (gx, gy: change of dir per pixel), as the
// starfield draws it: RA 0h at the centre, increasing to the left.
float3 milky_way_light(float3 dir, float3 gx, float3 gy)
{
    if (u_stars.z <= 0.0) {
        return 0.0;
    }
    float rxy2 = max(dir.x * dir.x + dir.y * dir.y, 1e-12);
    float2 uv = float2(frac(0.5 - atan2(dir.y, dir.x) / (2.0 * kPi)), 0.5 - asin(clamp(dir.z, -1.0, 1.0)) / kPi);
    // d ra = (x dy - y dx) / (x^2 + y^2), d dec = dz / cos(dec): no seam at RA 0h.
    float cos_dec = sqrt(rxy2);
    float2 dx = float2(-(dir.x * gx.y - dir.y * gx.x) / rxy2 / (2.0 * kPi), -gx.z / cos_dec / kPi);
    float2 dy = float2(-(dir.x * gy.y - dir.y * gy.x) / rxy2 / (2.0 * kPi), -gy.z / cos_dec / kPi);
    return u_milky_way.SampleGrad(u_milky_way_sampler, uv, dx, dy).rgb * u_stars.z;
}

// The stars seen in this pixel, drawn like the starfield's dots (same size in
// pixels and profile) at the place each star's image falls: the lens map is
// linearized as sky offset = gx u + gy v (u, v in pixels). Their brightness is
// multiplied by the magnification, the pixel's solid angle without the hole
// over the sky's solid angle it now sees (capped: point images on a caustic).
float3 point_stars(float3 dir, float3 gx, float3 gy, float pixel_solid_angle)
{
    float a = dot(gx, gx);
    float b = dot(gx, gy);
    float c = dot(gy, gy);
    float det = a * c - b * b; // |gx x gy|^2
    if (det <= 0.0) {
        return 0.0;
    }
    float magnification = min(pixel_solid_angle / sqrt(det), u_sky_params.y);

    // The 2 x 2 x 2 grid cells nearest to dir hold every star within half a cell
    // (1 / grid) of it; sky_light() hands over to the cube before a dot reaches farther.
    int grid = (int)u_stars.x;
    int3 c0 = int3(floor((dir * 0.5 + 0.5) * u_stars.x - 0.5));
    float3 sum = 0.0;
    [loop]
    for (int k = 0; k < 8; ++k) {
        int3 cell = c0 + int3(k & 1, (k >> 1) & 1, k >> 2);
        if (any(cell < 0) || any(cell >= grid)) {
            continue;
        }
        uint index = ((uint)cell.z * (uint)grid + (uint)cell.y) * (uint)grid + (uint)cell.x;
        uint first = u_star_cells.Load(index * 4);
        uint last = u_star_cells.Load(index * 4 + 4);
        [loop]
        for (uint i = first; i < last; ++i) {
            float4 dir_size = asfloat(u_star_data.Load4(i * 32));
            float4 color = asfloat(u_star_data.Load4(i * 32 + 16));
            float3 delta = dir_size.xyz - dir;
            float px = dot(gx, delta);
            float py = dot(gy, delta);
            float u = (c * px - b * py) / det;
            float v = (a * py - b * px) / det;
            // starfield.vert/.frag: a quad size px wide, exp(-3.5 r^2) with r = 1 at its edges.
            float half_size = 0.5 * dir_size.w;
            if (max(abs(u), abs(v)) < half_size) {
                sum += color.rgb * color.a * exp(-3.5 * (u * u + v * v) / (half_size * half_size));
            }
        }
    }
    return sum * (u_stars.y * magnification);
}

// The sky seen through a pixel along dir (gx, gy: change of dir per pixel).
float3 sky_light(float3 dir, float3 gx, float3 gy, float pixel_solid_angle)
{
    // A star's dot reaches as far as the largest stretch of the lens map times the
    // largest dot radius from its image; point_stars() finds stars within 1 / grid.
    float a = dot(gx, gx);
    float b = dot(gx, gy);
    float c = dot(gy, gy);
    float stretch = sqrt(0.5 * (a + c) + sqrt(0.25 * (a - c) * (a - c) + b * b));
    float cube = smoothstep(0.6, 1.0, stretch * u_sky_params.x * u_stars.x);
    float3 light = 0.0;
    if (cube > 0.0) {
        light += cube * u_sky.SampleGrad(u_sky_sampler, dir, gx, gy).rgb;
    }
    if (cube < 1.0) {
        light += (1.0 - cube) * (milky_way_light(dir, gx, gy) + point_stars(dir, gx, gy, pixel_solid_angle));
    }
    return light;
}
