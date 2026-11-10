"""Relief maps: a body's surface slopes and its water, baked from a height map.

The scene's `[bodies.relief]` texture, in the projection of the body's albedo map
(2:1 simple cylindrical, rows from north to south, east longitude increasing to
the right). The body shader tilts the surface normal by the slopes, which shows
the relief near the terminator, and reflects the sun off the water.

    red, green  slope toward the east and toward the north (rise over run), as
                0.5 + 0.5 sign(s) sqrt(|s| / MAX_SLOPE): the square root keeps the
                gentle slopes, which are nearly all of a planet seen at kilometers
                per pixel, clear of the JPEG's errors
    blue        share of the pixel that is open water

relief_decode in src/scene/relief.hpp reads the same encoding.

Requires numpy and Pillow.
"""

import math

import numpy as np
from PIL import Image

MAX_SLOPE = 1.0  # kReliefMaxSlope in src/scene/relief.hpp


def smooth_rows(height):
    """Low-pass along each row, more toward the poles, where the columns of the
    map are closer together than its rows: without it the east-west slopes there
    come from differences over a fraction of the map's resolution."""
    rows, width = height.shape
    lat = (0.5 - (np.arange(rows) + 0.5) / rows) * math.pi
    k = np.fft.rfftfreq(width, d=1.0 / width)  # cycles around the row
    excess = 1.0 / np.maximum(np.cos(lat), 1e-3) ** 2 - 1.0
    window = np.exp(-np.outer(excess, (k / (0.5 * width)) ** 2))
    return np.fft.irfft(np.fft.rfft(height, axis=1) * window, n=width, axis=1)


def slopes(height_km, radius_km, coverage=None):
    """East and north slopes on a sphere, with optional valid area per cell."""
    rows, width = height_km.shape
    if coverage is None:
        h = smooth_rows(height_km)
    else:
        # Missing cells cannot pull the polar low-pass toward a fictitious height.
        weights = smooth_rows(coverage)
        weighted = smooth_rows(height_km * coverage)
        h = np.divide(weighted, weights, out=np.zeros_like(weighted), where=weights > 1e-6)
    lat = (0.5 - (np.arange(rows) + 0.5) / rows) * math.pi
    step = math.pi / rows  # radians per pixel, both ways
    east = (np.roll(h, -1, axis=1) - np.roll(h, 1, axis=1)) / (2.0 * step * radius_km)
    east /= np.maximum(np.cos(lat), 1e-3)[:, None]
    north = np.empty_like(h)
    north[1:-1] = (h[:-2] - h[2:]) / (2.0 * step * radius_km)
    if coverage is not None:
        # Fade the difference by the least observed cell in its stencil.
        east *= np.minimum(coverage, np.minimum(np.roll(coverage, 1, axis=1), np.roll(coverage, -1, axis=1)))
        north[1:-1] *= np.minimum(coverage[1:-1], np.minimum(coverage[:-2], coverage[2:]))
    north[0] = north[1]
    north[-1] = north[-2]
    if coverage is not None:
        north[0] *= coverage[0]
        north[-1] *= coverage[-1]
    return east, north


def encode(slope):
    s = np.clip(slope / MAX_SLOPE, -1.0, 1.0)
    return np.round(255.0 * (0.5 + 0.5 * np.sign(s) * np.sqrt(np.abs(s)))).astype(np.uint8)


def write(path, height_km, radius_km, water=None, coverage=None):
    """Writes the relief map of a height map; `water` is the share of each pixel
    that is open water (its surface is level whatever the heights say)."""
    east, north = slopes(np.asarray(height_km, dtype=np.float64), radius_km, coverage)
    if water is None:
        water = np.zeros(east.shape)
    water = np.clip(np.asarray(water, dtype=np.float64), 0.0, 1.0)
    east *= 1.0 - water
    north *= 1.0 - water
    rgb = np.stack([encode(east), encode(north), np.round(255.0 * water).astype(np.uint8)], axis=-1)
    # No chroma subsampling: the channels are unrelated quantities.
    Image.fromarray(rgb, "RGB").save(path, "JPEG", quality=92, subsampling=0, optimize=True)
    steep = np.degrees(np.arctan(np.hypot(east, north)))
    weight = np.cos((0.5 - (np.arange(east.shape[0]) + 0.5) / east.shape[0]) * math.pi)[:, None] * np.ones(east.shape)
    order = np.argsort(steep, axis=None)
    cum = np.cumsum(weight.ravel()[order])
    q = [steep.ravel()[order][np.searchsorted(cum, f * cum[-1])] for f in (0.5, 0.9, 0.99)]
    print("%s: slopes %.2f / %.2f / %.2f deg (median, 90%%, 99%%), steepest %.1f deg, water %.1f%%" %
          (path, q[0], q[1], q[2], steep.max(), 100.0 * (water * weight).sum() / weight.sum()))
