"""Downscale global surface maps to equirectangular textures for assets/textures/.

Sources and conventions are listed in assets/textures/SOURCES.md. Download the
originals there first, then:

    python prepare_textures.py <source_dir> <output_dir> [width]

Expected file names in <source_dir> (the default is 2048 px wide):
    PIA07782.jpg  -> jupiter.jpg
    world.200407.3x5400x2700.jpg -> earth.jpg
    lroc_color_poles_2k.tif      -> moon.jpg
    io.tif        -> io.jpg
    europa.tif    -> europa.jpg
    ganymede.tif  -> ganymede.jpg
    callisto.tif  -> callisto.jpg
    MI_170630_DLR_basemap_degrees.tif (from Cassini_DLR_Mimas.zip) -> mimas.jpg
    Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif -> enceladus.jpg
    Tethys_Cassini_mosaic_global_293m.tif            -> tethys.jpg
    Dione_Cassini_Voyager_mosaic_global_154m.tif     -> dione.jpg
    Rhea_Cassini_Voyager_mosaic_global_417m.tif      -> rhea.jpg
    Iapetus_Cassini_Voyager_mosaic_global_783m.tif   -> iapetus.jpg

Requires Pillow. The USGS mosaics are 100-200 MB GeoTIFFs (up to ~190 Mpx), so
the decompression-bomb guard is disabled and images are pre-reduced by an
integer factor before the final Lanczos resize.
"""

import os
import sys

from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# (source, output, fill polar no-data gaps). Gap filling is only for the USGS
# mosaics: on complete maps it would mistake genuinely dark pixels (e.g. polar
# ocean on Earth) for missing data.
MAPS = [
    ('PIA07782.jpg', 'jupiter.jpg', False),
    ('io.tif', 'io.jpg', True),
    ('europa.tif', 'europa.jpg', True),
    ('ganymede.tif', 'ganymede.jpg', True),
    ('callisto.tif', 'callisto.jpg', True),
    ('world.200407.3x5400x2700.jpg', 'earth.jpg', False),
    ('lroc_color_poles_2k.tif', 'moon.jpg', False),
    ('MI_170630_DLR_basemap_degrees.tif', 'mimas.jpg', True),
    ('Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif', 'enceladus.jpg', True),
    ('Tethys_Cassini_mosaic_global_293m.tif', 'tethys.jpg', True),
    ('Dione_Cassini_Voyager_mosaic_global_154m.tif', 'dione.jpg', True),
    ('Rhea_Cassini_Voyager_mosaic_global_417m.tif', 'rhea.jpg', True),
    ('Iapetus_Cassini_Voyager_mosaic_global_783m.tif', 'iapetus.jpg', True),
]


def fill_polar_gaps(img, max_value=2, polar_fraction=0.3):
    """Replaces no-data pixels (pure black) near the poles with the mean of the
    valid pixels in the same row, or of the nearest row that has valid pixels.
    Some mosaics have no coverage close to the poles, which would otherwise
    render as black caps."""
    w, h = img.size
    px = img.load()
    polar_rows = int(h * polar_fraction)
    rows = list(range(polar_rows)) + list(range(h - polar_rows, h))

    def is_gap(p):
        return max(p) <= max_value

    row_mean = {}
    for y in rows:
        valid = [px[x, y] for x in range(w) if not is_gap(px[x, y])]
        if len(valid) > w // 20:
            row_mean[y] = tuple(sum(c[i] for c in valid) // len(valid) for i in range(3))

    filled = 0
    for y in rows:
        mean = row_mean.get(y)
        if mean is None:
            nearest = min(row_mean, key=lambda r: abs(r - y), default=None)
            if nearest is None:
                continue
            mean = row_mean[nearest]
        for x in range(w):
            if is_gap(px[x, y]):
                px[x, y] = mean
                filled += 1
    return filled


def prepare(src, dst, width, fill_gaps):
    height = width // 2
    with Image.open(src) as img:
        print(f'{os.path.basename(src)}: {img.size[0]}x{img.size[1]} {img.mode}')
        img = img.convert('RGB') if img.mode not in ('RGB', 'L') else img
        factor = max(1, min(img.size[0] // (width * 2), img.size[1] // (height * 2)))
        if factor > 1:
            img = img.reduce(factor)
        img = img.resize((width, height), Image.Resampling.LANCZOS).convert('RGB')
        filled = fill_polar_gaps(img) if fill_gaps else 0
        img.save(dst, 'JPEG', quality=90, optimize=True)
    print(f'  -> {dst} ({os.path.getsize(dst) // 1024} KB, {filled} gap pixels filled)')


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    src_dir, out_dir = sys.argv[1], sys.argv[2]
    width = int(sys.argv[3]) if len(sys.argv) == 4 else 2048
    os.makedirs(out_dir, exist_ok=True)
    for src_name, dst_name, fill_gaps in MAPS:
        src = os.path.join(src_dir, src_name)
        if not os.path.exists(src):
            print(f'skip {src_name} (not found)')
            continue
        prepare(src, os.path.join(out_dir, dst_name), width, fill_gaps)


if __name__ == '__main__':
    main()
