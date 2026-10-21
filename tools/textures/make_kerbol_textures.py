"""Make the procedural surfaces of the Kerbol scene (assets/scenes/kerbol.toml).

    python -E -P make_kerbol_textures.py <output_dir> [name ...] [--width N]

Writes <name>.jpg (2:1 simple cylindrical, east longitude increasing to the right
from -180 deg at the left edge, as the scene's texture lookup expects) for every
world of the scene (see MAKERS), or only the names given; --width overrides the
map widths (2048, or 1024 for the smallest moons).

These are invented surfaces, not the game's: they keep a few traits that make
each world recognizable (Kerbin's blue seas, green land and ice caps; the Mun's
grey craters; Duna's rust plains and polar caps; Jool's green bands; Eve's purple
land and dark seas; Minmus' mint color and flat plains; Laythe's islands; Vall's
and Eeloo's cracked ice; Dres' canyons; Tylo's patches...) and nothing of the
game's actual geography or maps. Everything is fixed by the seeds below,
and the noise uses its own integer hash (no numpy random generator), so the
output does not depend on the numpy version.

Method: every pixel is a point on the unit sphere and all noise is 3-D gradient
noise sampled there, so the maps have no seam at 180 deg and no pinching at the
poles. Octaves stop where their wavelength reaches about two pixels. Colors are
mixed in linear light and written as sRGB. The scene shades without normal maps,
so craters are drawn in albedo only (darker floors, bright rims and ejecta).

Requires numpy and Pillow.
"""

import math
import os
import sys

import numpy as np
from PIL import Image

U32 = np.uint32


# ---------------------------------------------------------------------------
# Noise

def _hash(ix, iy, iz, seed):
    """32-bit integer hash of lattice coordinates (int64 arrays) and a seed."""
    h = (ix.astype(np.int64) * 0x8DA6B343 ^ iy.astype(np.int64) * 0xD8163841 ^
         iz.astype(np.int64) * 0xCB1AB31F ^ (seed * 0x9E3779B1)) & 0xFFFFFFFF
    h = h.astype(U32)
    h ^= h >> U32(16)
    h *= U32(0x7FEB352D)
    h ^= h >> U32(15)
    h *= U32(0x846CA68B)
    h ^= h >> U32(16)
    return h


_GRADIENTS = np.array([[1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
                       [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
                       [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1]], dtype=np.float64)


def gradient_noise(p, seed):
    """Perlin-style gradient noise at points p (..., 3), roughly in [-1, 1]."""
    cell = np.floor(p)
    f = p - cell
    i = cell.astype(np.int64)
    w = f * f * f * (f * (f * 6.0 - 15.0) + 10.0)
    total = 0.0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                h = _hash(i[..., 0] + dx, i[..., 1] + dy, i[..., 2] + dz, seed)
                g = _GRADIENTS[(h % U32(12)).astype(np.int64)]
                d = f - np.array([dx, dy, dz], dtype=np.float64)
                dot = np.einsum('...k,...k->...', g, d)
                wx = w[..., 0] if dx else 1.0 - w[..., 0]
                wy = w[..., 1] if dy else 1.0 - w[..., 1]
                wz = w[..., 2] if dz else 1.0 - w[..., 2]
                total = total + wx * wy * wz * dot
    return total


def fbm(p, seed, base, octaves, gain=0.5):
    """Fractal sum of gradient noise, normalized to about [-1, 1]."""
    total = np.zeros(p.shape[:-1])
    amplitude, norm, frequency = 1.0, 0.0, base
    for k in range(octaves):
        total += amplitude * gradient_noise(p * frequency, seed + 101 * k)
        norm += amplitude
        amplitude *= gain
        frequency *= 2.0
    return total / norm * 1.6


def octaves_for(base, width):
    """Octaves of fbm(base) whose wavelength stays above ~2 pixels."""
    finest = width / (2.0 * math.pi) / 2.0  # cycles per radian at 2 pixels per cycle
    return max(1, int(math.floor(math.log2(finest / base))) + 1)


# ---------------------------------------------------------------------------
# Map geometry and color helpers

def sphere_grid(width):
    height = width // 2
    lat = (0.5 - (np.arange(height) + 0.5) / height) * math.pi
    lon = ((np.arange(width) + 0.5) / width - 0.5) * 2.0 * math.pi
    lon2, lat2 = np.meshgrid(lon, lat)
    p = np.stack([np.cos(lat2) * np.cos(lon2), np.cos(lat2) * np.sin(lon2), np.sin(lat2)], axis=-1)
    return p, lat2, lon2


def weighted_quantile(values, lat, q):
    """Quantile of `values` weighted by pixel area (cos latitude)."""
    order = np.argsort(values, axis=None)
    weights = np.cos(lat).ravel()[order]
    cum = np.cumsum(weights)
    return values.ravel()[order][np.searchsorted(cum, q * cum[-1])]


def srgb_to_linear(c):
    c = np.asarray(c, dtype=np.float64) / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def linear_to_srgb8(c):
    c = np.clip(c, 0.0, 1.0)
    s = np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1.0 / 2.4) - 0.055)
    return np.round(s * 255.0).astype(np.uint8)


def color(hex_or_rgb):
    if isinstance(hex_or_rgb, str):
        v = int(hex_or_rgb.lstrip('#'), 16)
        hex_or_rgb = ((v >> 16) & 255, (v >> 8) & 255, v & 255)
    return srgb_to_linear(hex_or_rgb)


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def mix(a, b, t):
    t = np.asarray(t)[..., None]
    return a * (1.0 - t) + b * t


def ramp(stops, x):
    """Piecewise-linear color ramp: stops = [(x0, color), (x1, color), ...]."""
    out = np.zeros(x.shape + (3,))
    xs = [s[0] for s in stops]
    for k in range(3):
        out[..., k] = np.interp(x, xs, [s[1][k] for s in stops])
    return out


# ---------------------------------------------------------------------------
# Craters (albedo only)

class Rng:
    """Deterministic stream of uniforms in [0, 1) from the same hash."""

    def __init__(self, seed):
        self.seed = seed
        self.n = 0

    def uniform(self, count):
        idx = np.arange(self.n, self.n + count, dtype=np.int64)
        self.n += count
        return _hash(idx, idx * 0 + 7, idx * 0 + 13, self.seed).astype(np.float64) / 4294967296.0


def add_craters(albedo, p, lat, rng, count, r_min_deg, r_max_deg, slope, floor, rim, ejecta, rays=0):
    """Multiplies albedo by crater patterns. Radii follow a power law (cumulative
    N(>r) ~ r^-slope); each crater gets a random age that fades its contrast.
    The `rays` largest young craters also get bright rays."""
    height, width = lat.shape
    u = rng.uniform(count * 4)
    z = 2.0 * u[0::4] - 1.0
    phi = 2.0 * math.pi * u[1::4]
    centers = np.stack([np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z], axis=-1)
    a, b = r_min_deg ** -slope, r_max_deg ** -slope
    radii = np.radians((a + u[2::4] * (b - a)) ** (-1.0 / slope))
    ages = u[3::4]
    order = np.argsort(-radii)
    ray_set = set(order[np.argsort(ages[order[:max(rays * 4, 1)]])[:rays]]) if rays else set()
    pixel = math.pi / height
    for n in order:
        c, r, age = centers[n], radii[n], ages[n]
        reach = r * (6.0 if n in ray_set else 2.5)
        clat = math.asin(c[2])
        lat0, lat1 = clat - reach, clat + reach
        rows = slice(max(0, int((0.5 - lat1 / math.pi) * height)), min(height, int((0.5 - lat0 / math.pi) * height) + 2))
        if abs(clat) + reach >= math.pi / 2 or reach / max(math.cos(abs(clat) + reach), 1e-6) > math.pi:
            cols = [slice(0, width)]
        else:
            clon = math.atan2(c[1], c[0])
            dlon = reach / math.cos(abs(clat) + reach)
            c0 = int(((clon - dlon) / (2 * math.pi) + 0.5) * width)
            c1 = int(((clon + dlon) / (2 * math.pi) + 0.5) * width) + 2
            cols = [slice(max(c0, 0), min(c1, width))]
            if c0 < 0:
                cols.append(slice(width + c0, width))
            if c1 > width:
                cols.append(slice(0, c1 - width))
        fresh = (1.0 - age) ** 2
        for cs in cols:
            q = p[rows, cs]
            d = np.arccos(np.clip(q @ c, -1.0, 1.0))
            x = d / r
            # Antialias: soften features narrower than a pixel.
            soft = max(pixel / r, 0.04)
            f = np.ones(x.shape)
            f *= 1.0 - floor * (0.4 + 0.6 * fresh) * (1.0 - smoothstep(0.75 - soft, 0.85 + soft, x))
            f *= 1.0 + rim * (0.3 + 0.7 * fresh) * np.exp(-((x - 0.95) / (0.1 + soft)) ** 2)
            f *= 1.0 + ejecta * fresh * smoothstep(0.95, 1.1, x) * np.exp(-(x - 1.0) / 0.6) * smoothstep(2.5, 2.0, x)
            if n in ray_set:
                e1 = np.cross(c, [0.0, 0.0, 1.0] if abs(c[2]) < 0.9 else [1.0, 0.0, 0.0])
                e1 /= np.linalg.norm(e1)
                e2 = np.cross(c, e1)
                az = np.arctan2(q @ e2, q @ e1)
                k = 7.0 + 6.0 * age
                streak = np.maximum(0.0, np.cos(k * az + 10.0 * age) * np.cos(2.3 * k * az + 3.0)) ** 2
                f *= 1.0 + 0.25 * streak * smoothstep(1.0, 1.3, x) * np.exp(-(x - 1.0) / 1.8)
            albedo[rows, cs] *= f[..., None]


# ---------------------------------------------------------------------------
# Bodies

def continents(p, lat, seed, ocean_fraction, width):
    """Warped fractal terrain with ridged mountains, cut by a sea level that floods
    `ocean_fraction` of the surface: (land mask, elevation 0..1 above the sea,
    depth 0..1 below it)."""
    oct_fine = octaves_for(1.6, width)
    warp = np.stack([fbm(p, seed + 10 + k, 1.2, 4) for k in range(3)], axis=-1)
    q = p + 0.35 * warp
    h = fbm(q, seed, 1.6, oct_fine) + 0.35 * (1.0 - np.abs(fbm(q, seed + 1, 3.0, oct_fine - 1))) ** 3
    sea = weighted_quantile(h, lat, ocean_fraction)
    land = h > sea
    top = np.percentile(h[land], 99.5)
    elev = np.clip((h - sea) / (top - sea), 0.0, 1.0)
    depth = np.clip((sea - h) / (sea - np.percentile(h[~land], 1.0)), 0.0, 1.0)
    return land, elev, depth


def mottled(p, seed, dark, light, width, base=1.4, edges=(-0.15, 0.25), detail=0.1):
    """Two tones in large patches, with fine brightness noise."""
    tone = smoothstep(edges[0], edges[1], fbm(p, seed, base, 5))
    albedo = mix(color(dark), color(light), tone)
    return albedo * (1.0 + detail * fbm(p, seed + 1, 6.0, octaves_for(6.0, width)))[..., None]


def lines(p, lat, seed, base, octaves, thickness_deg, warp=0.3):
    """A network of curving lines (cracks, canyons): the zero contours of warped
    noise, about `thickness_deg` wide everywhere; 1 on a line, 0 away from it.
    The distance to the contour is |n| / |grad n|, with the gradient taken on the
    map (east-west steps shortened by cos(latitude))."""
    q = p + warp * np.stack([fbm(p, seed + 1 + k, base * 0.7, 4) for k in range(3)], axis=-1)
    n = fbm(q, seed, base, octaves)
    deg_per_px = 180.0 / lat.shape[0]
    gy = np.gradient(n, axis=0)
    gx = (np.roll(n, -1, axis=1) - np.roll(n, 1, axis=1)) * 0.5 / np.maximum(np.cos(lat), 0.05)
    distance_deg = np.abs(n) / np.maximum(np.hypot(gx, gy), 1e-9) * deg_per_px
    return 1.0 - smoothstep(0.5 * thickness_deg, thickness_deg + deg_per_px, distance_deg)


def make_kerbin(width):
    p, lat, _ = sphere_grid(width)
    oct_fine = octaves_for(1.6, width)
    land, elev, depth = continents(p, lat, 1, 0.64, width)   # two-thirds ocean

    moisture = fbm(p, 5, 2.0, 6) + 0.25 * np.cos(2.2 * lat)
    alat = np.abs(np.degrees(lat))
    desert = smoothstep(-0.05, -0.25, moisture) * smoothstep(8.0, 18.0, alat) * smoothstep(45.0, 32.0, alat)

    ocean = mix(color((24, 60, 104)), color((8, 22, 52)), smoothstep(0.0, 0.35, depth))
    ocean = mix(ocean, color((40, 96, 130)), smoothstep(0.06, 0.0, depth) * 0.7)
    ground = ramp([(0.0, color((70, 104, 48))), (0.18, color((52, 86, 40))), (0.45, color((88, 92, 58))),
                   (0.7, color((118, 104, 82))), (0.86, color((138, 134, 128))), (1.0, color((150, 148, 146)))],
                  elev)
    ground = mix(ground, color((198, 172, 116)), desert * (1.0 - smoothstep(0.5, 0.75, elev)))
    ground = mix(ground, color((176, 160, 118)), smoothstep(0.025, 0.0, elev) * 0.6)  # beaches
    ground *= (1.0 + 0.12 * fbm(p, 7, 24.0, max(1, oct_fine - 3)))[..., None]
    albedo = np.where(land[..., None], ground, ocean)

    snow_line = 0.94 - 0.45 * smoothstep(45.0, 72.0, alat)
    snow = smoothstep(snow_line, snow_line + 0.06, elev + 0.04 * fbm(p, 8, 10.0, 4)) * land
    cap_edge = 73.0 + 9.0 * fbm(p, 9, 2.0, 6)
    cap = smoothstep(cap_edge - 1.5, cap_edge + 1.5, alat)
    white = color((232, 238, 244))
    albedo = mix(albedo, white, np.maximum(snow, cap) * (1.0 - 0.08 * fbm(p, 10, 16.0, 3)))
    return albedo


def make_mun(width):
    p, lat, _ = sphere_grid(width)
    oct_fine = octaves_for(2.0, width)
    tone = smoothstep(-0.15, 0.25, fbm(p, 21, 1.4, 5))
    albedo = mix(color((136, 136, 135)), color((164, 164, 162)), tone)
    albedo *= (1.0 + 0.1 * fbm(p, 22, 6.0, oct_fine - 1))[..., None]
    add_craters(albedo, p, lat, Rng(23), count=4500, r_min_deg=0.3, r_max_deg=25.0, slope=1.5,
                floor=0.16, rim=0.28, ejecta=0.14, rays=2)
    return albedo


def make_duna(width):
    p, lat, _ = sphere_grid(width)
    oct_fine = octaves_for(1.4, width)
    warp = np.stack([fbm(p, 31 + k, 1.0, 4) for k in range(3)], axis=-1)
    h = fbm(p + 0.3 * warp, 34, 1.4, oct_fine)
    highland = smoothstep(-0.1, 0.25, h)
    albedo = mix(color((196, 108, 72)), color((128, 60, 42)), highland)
    albedo = mix(albedo, color((214, 136, 96)), smoothstep(0.1, 0.5, fbm(p, 35, 3.0, 5)) * (1.0 - highland) * 0.5)
    albedo *= (1.0 + 0.1 * fbm(p, 36, 12.0, max(1, oct_fine - 3)))[..., None]
    add_craters(albedo, p, lat, Rng(37), count=700, r_min_deg=0.4, r_max_deg=12.0, slope=1.6,
                floor=0.08, rim=0.08, ejecta=0.04)
    alat = np.abs(np.degrees(lat))
    cap_edge = 79.0 + 6.0 * fbm(p, 38, 2.5, 6)
    cap = smoothstep(cap_edge - 1.0, cap_edge + 1.0, alat)
    return mix(albedo, color((232, 224, 216)), cap)


def make_jool(width):
    p, lat, lon = sphere_grid(width)
    # Zonal flow: noise squashed along latitude makes features long in longitude.
    stretched = p * np.array([1.0, 1.0, 6.0])
    turb = fbm(stretched, 41, 1.5, octaves_for(9.0, width))
    swirl = fbm(p + 0.5 * np.stack([fbm(p, 42 + k, 2.0, 4) for k in range(3)], axis=-1), 45, 3.0, 6)
    latw = np.degrees(lat) + 3.5 * turb + 2.5 * swirl
    # Bands of irregular width: 1-D noise across the warped latitude.
    axis = np.stack([latw / 9.0, np.full(lat.shape, 0.37), np.full(lat.shape, 0.61)], axis=-1)
    bands = 0.5 + 0.55 * fbm(axis, 47, 1.0, 4, gain=0.45)
    filaments = fbm(p * np.array([1.0, 1.0, 10.0]), 48, 3.0, octaves_for(30.0, width))
    bands += 0.1 * filaments + 0.06 * swirl
    rng = Rng(46)
    palette = [(0.0, color((46, 92, 40))), (0.3, color((82, 134, 58))), (0.55, color((124, 172, 82))),
               (0.8, color((166, 204, 118))), (1.0, color((196, 222, 152)))]
    albedo = ramp(palette, np.clip(bands, 0.0, 1.0))
    # Storms: a few ovals, long in longitude.
    for n in range(5):
        u = rng.uniform(4)
        slat, slon = math.radians(-50 + 100 * u[0]), math.pi * (2 * u[1] - 1)
        a, b = math.radians(3 + 6 * u[2]), math.radians(1.5 + 2.5 * u[2])
        dlon = np.angle(np.exp(1j * (lon - slon))) * np.cos(lat)
        r2 = (dlon / a) ** 2 + ((lat - slat) / b) ** 2
        spot = np.exp(-r2 ** 2)
        albedo = mix(albedo, color((210, 228, 170)) if u[3] > 0.5 else color((52, 96, 50)), 0.7 * spot)
    polar = smoothstep(60.0, 85.0, np.abs(np.degrees(lat)))
    return albedo * (1.0 - 0.25 * polar)[..., None]


def make_moho(width):
    """Dark brown, heavily cratered."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 61, (88, 70, 58), (132, 108, 90), width)
    add_craters(albedo, p, lat, Rng(63), count=4000, r_min_deg=0.3, r_max_deg=20.0, slope=1.5,
                floor=0.14, rim=0.24, ejecta=0.12, rays=2)
    return albedo


def make_eve(width):
    """Purple land around dark seas."""
    p, lat, _ = sphere_grid(width)
    land, elev, depth = continents(p, lat, 71, 0.35, width)
    sea = mix(color((44, 28, 78)), color((18, 10, 42)), smoothstep(0.0, 0.4, depth))
    ground = ramp([(0.0, color((92, 58, 112))), (0.35, color((112, 76, 136))), (0.7, color((140, 106, 160))),
                   (1.0, color((178, 156, 192)))], elev)
    ground *= (1.0 + 0.12 * fbm(p, 75, 12.0, octaves_for(12.0, width)))[..., None]
    return np.where(land[..., None], ground, sea)


def make_gilly(width):
    """A small brown-grey rock."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 81, (98, 84, 72), (138, 120, 102), width, base=2.0, detail=0.15)
    add_craters(albedo, p, lat, Rng(83), count=250, r_min_deg=1.5, r_max_deg=20.0, slope=1.4,
                floor=0.1, rim=0.12, ejecta=0.05)
    return albedo


def make_minmus(width):
    """Mint green, with flat pale "salt" plains whose edges are crisp."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 91, (122, 190, 164), (156, 216, 192), width)
    f = fbm(p + 0.3 * np.stack([fbm(p, 93 + k, 1.5, 4) for k in range(3)], axis=-1), 92, 1.8, 6)
    flat_level = weighted_quantile(f, lat, 0.22)
    flats = smoothstep(flat_level + 0.01, flat_level - 0.01, f)
    albedo = mix(albedo, color((198, 234, 218)), flats)
    add_craters(albedo, p, lat, Rng(97), count=180, r_min_deg=0.8, r_max_deg=8.0, slope=1.5,
                floor=0.05, rim=0.06, ejecta=0.02)
    return albedo


def make_ike(width):
    """Dark grey, cratered."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 101, (80, 80, 84), (116, 116, 120), width)
    add_craters(albedo, p, lat, Rng(103), count=2500, r_min_deg=0.35, r_max_deg=18.0, slope=1.5,
                floor=0.14, rim=0.24, ejecta=0.1, rays=1)
    return albedo


def make_dres(width):
    """Grey, cratered, cut by a few long dark canyons."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 111, (110, 104, 98), (150, 144, 138), width)
    add_craters(albedo, p, lat, Rng(113), count=2200, r_min_deg=0.35, r_max_deg=15.0, slope=1.5,
                floor=0.12, rim=0.2, ejecta=0.08)
    canyons = lines(p, lat, 115, 0.9, 4, 0.6, warp=0.2)
    canyons = canyons * smoothstep(0.1, 0.4, fbm(p, 119, 1.2, 4))  # segments, not a global net
    return albedo * (1.0 - 0.55 * canyons)[..., None]


def make_laythe(width):
    """An ocean world with scattered sandy islands."""
    p, lat, _ = sphere_grid(width)
    land, elev, depth = continents(p, lat, 121, 0.86, width)
    ocean = mix(color((22, 56, 100)), color((8, 22, 54)), smoothstep(0.0, 0.3, depth))
    ocean = mix(ocean, color((44, 104, 136)), smoothstep(0.05, 0.0, depth) * 0.7)
    ground = ramp([(0.0, color((184, 168, 128))), (0.3, color((160, 146, 114))), (0.7, color((128, 120, 108))),
                   (1.0, color((150, 146, 140)))], elev)
    ground *= (1.0 + 0.1 * fbm(p, 125, 16.0, octaves_for(16.0, width)))[..., None]
    return np.where(land[..., None], ground, ocean)


def make_vall(width):
    """Pale blue ice crossed by darker cracks."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 131, (168, 188, 204), (206, 222, 232), width)
    cracks = np.maximum(lines(p, lat, 133, 1.5, 4, 0.35, warp=0.25), 0.6 * lines(p, lat, 137, 3.0, 3, 0.2, warp=0.2))
    albedo = mix(albedo, color((112, 134, 156)), 0.75 * cracks)
    add_craters(albedo, p, lat, Rng(139), count=400, r_min_deg=0.4, r_max_deg=10.0, slope=1.6,
                floor=0.06, rim=0.1, ejecta=0.06)
    return albedo


def make_tylo(width):
    """Light grey with large dark patches, cratered."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 141, (118, 114, 110), (206, 200, 194), width, base=0.9, edges=(-0.08, 0.08))
    add_craters(albedo, p, lat, Rng(143), count=2000, r_min_deg=0.3, r_max_deg=16.0, slope=1.5,
                floor=0.12, rim=0.2, ejecta=0.1, rays=2)
    return albedo


def make_bop(width):
    """Rough dark brown."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 151, (88, 62, 48), (128, 94, 72), width, base=2.0, detail=0.2)
    add_craters(albedo, p, lat, Rng(153), count=500, r_min_deg=1.0, r_max_deg=18.0, slope=1.4,
                floor=0.12, rim=0.16, ejecta=0.06)
    return albedo


def make_pol(width):
    """Lumpy yellow-tan."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 161, (186, 160, 92), (222, 204, 140), width, base=2.4, detail=0.18)
    albedo *= (1.0 - 0.18 * smoothstep(0.2, 0.6, fbm(p, 163, 5.0, 5)))[..., None]
    add_craters(albedo, p, lat, Rng(165), count=150, r_min_deg=1.0, r_max_deg=12.0, slope=1.4,
                floor=0.06, rim=0.08, ejecta=0.03)
    return albedo


def make_eeloo(width):
    """White ice with dark brownish cracks and a few dark patches."""
    p, lat, _ = sphere_grid(width)
    albedo = mottled(p, 171, (204, 208, 210), (232, 234, 236), width)
    albedo = mix(albedo, color((150, 140, 130)), 0.35 * smoothstep(0.35, 0.7, fbm(p, 173, 1.6, 5)))
    cracks = np.maximum(lines(p, lat, 175, 1.2, 4, 0.4), 0.7 * lines(p, lat, 179, 2.6, 3, 0.25))
    albedo = mix(albedo, color((118, 98, 84)), 0.8 * cracks)
    add_craters(albedo, p, lat, Rng(181), count=300, r_min_deg=0.4, r_max_deg=8.0, slope=1.6,
                floor=0.05, rim=0.08, ejecta=0.05)
    return albedo


# Name -> (maker, default width). The smallest moons get 1024-pixel maps.
MAKERS = {"kerbin": (make_kerbin, 2048), "mun": (make_mun, 2048), "duna": (make_duna, 2048),
          "jool": (make_jool, 2048), "moho": (make_moho, 2048), "eve": (make_eve, 2048),
          "gilly": (make_gilly, 1024), "minmus": (make_minmus, 2048), "ike": (make_ike, 2048),
          "dres": (make_dres, 2048), "laythe": (make_laythe, 2048), "vall": (make_vall, 2048),
          "tylo": (make_tylo, 2048), "bop": (make_bop, 1024), "pol": (make_pol, 1024),
          "eeloo": (make_eeloo, 2048)}


def main():
    args = sys.argv[1:]
    width = None
    if "--width" in args:
        i = args.index("--width")
        width = int(args[i + 1])
        del args[i:i + 2]
    if not args:
        print(__doc__)
        sys.exit(1)
    out_dir, names = args[0], args[1:] or list(MAKERS)
    for name in names:
        maker, default_width = MAKERS[name]
        albedo = maker(width or default_width)
        path = os.path.join(out_dir, name + ".jpg")
        Image.fromarray(linear_to_srgb8(albedo), "RGB").save(path, "JPEG", quality=90, optimize=True)
        print(path)


if __name__ == "__main__":
    main()
