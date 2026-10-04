// Exponential atmospheres in units of the planet radius (the planet is the
// unit sphere): column densities and phase functions shared by the atmosphere
// shell (atmosphere.frag.hlsl) and the sunlight reaching the ground
// (body.frag.hlsl).

#ifndef ASTRAXIS_ATMOSPHERE_HLSLI
#define ASTRAXIS_ATMOSPHERE_HLSLI

// Column density from radius r toward direction cosine mu (to the local
// vertical) to infinity, for a density exp(-(r - 1) / H), in units of H times
// the density at r = 1. Uses the Chapman grazing-incidence function
// Ch(x, mu) ~ c / ((c - 1) mu + 1), c = sqrt(pi x / 2) (exact at mu = 0 and 1,
// ~ 1/mu in between); below the horizon, twice the column from the periapsis
// minus the column on the far side (after C. Schüler, "An Approximation to
// the Chapman Grazing-Incidence Function for Atmospheric Scattering", GPU Pro
// 3, 2012). Rays that hit the planet are not handled here (see ray_lit).
float chapman_column(float r, float mu, float H)
{
    float c = sqrt(1.5707963 * r / H);
    if (mu >= 0.0) {
        return exp(-(r - 1.0) / H) * c / ((c - 1.0) * mu + 1.0);
    }
    float r0 = r * sqrt(saturate(1.0 - mu * mu)); // periapsis radius
    float c0 = sqrt(1.5707963 * r0 / H);
    return 2.0 * exp(-max(r0 - 1.0, 0.0) / H) * c0 - exp(-(r - 1.0) / H) * c / ((c - 1.0) * (-mu) + 1.0);
}

// Fraction of the sun's disk (angular radius sun_radius) seen above the
// planet from radius r at direction cosine mu.
float ray_lit(float r, float mu, float sun_radius)
{
    if (mu >= 0.0) {
        return 1.0;
    }
    float r0 = r * sqrt(saturate(1.0 - mu * mu));
    float penumbra = max(sun_radius * r * -mu, 1e-5); // the disk's width at the periapsis
    return smoothstep(1.0 - penumbra, 1.0 + penumbra, r0);
}

// Optical depth (per channel) toward the sun from radius r, direction cosine
// mu: rayleigh / haze hold the vertical (attenuation) optical depths (rgb) and
// scale heights (w).
float3 sun_optical_depth(float r, float mu, float4 rayleigh, float4 haze)
{
    return rayleigh.rgb * chapman_column(r, mu, rayleigh.w) + haze.rgb * chapman_column(r, mu, haze.w);
}

// Phase functions normalized to a mean of 1 over the sphere.
float rayleigh_phase(float cos_theta)
{
    return 0.75 * (1.0 + cos_theta * cos_theta);
}

float henyey_greenstein_phase(float cos_theta, float g)
{
    float d = 1.0 + g * g - 2.0 * g * cos_theta;
    return (1.0 - g * g) / (d * sqrt(d));
}

#endif
