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
    ura1vuu2.tif .. ura5vuu2.tif (JPL Solar System Simulator) -> ariel, umbriel,
                                                    titania, oberon, miranda.jpg
    Triton_Voyager2_ClrMosaic_GlobalFill_600m.tif    -> triton.jpg
    Pluto_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif  -> pluto.jpg
    Charon_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif -> charon.jpg
    Mercury_MESSENGER_mosaic_global_250m_2013.tif    -> mercury.jpg
    Mars_Viking_ClrMosaic_global_925m.tif            -> mars.jpg
    Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif        -> phobos.jpg
    Vesta_Dawn_FC_HAMO_Mosaic_Global_74ppd.tif       -> vesta.jpg
    Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif       -> ceres.jpg

Requires Pillow. The USGS mosaics are 100-200 MB GeoTIFFs (up to ~190 Mpx), so
the decompression-bomb guard is disabled and images are pre-reduced by an
integer factor before the final Lanczos resize.
"""

import os
import sys

from PIL import Image, ImageChops, ImageFilter, ImageMath, ImageStat

Image.MAX_IMAGE_PIXELS = None

# (source, output, gap filling): False for complete maps (filling would mistake
# genuinely dark pixels, e.g. polar ocean on Earth, for missing data), 'polar' for
# USGS mosaics with small polar gaps, 'unimaged' for maps with large unimaged
# areas (the Uranian moons' and Triton's northern hemispheres, Pluto and Charon
# south of ~30 S). ('unimaged', {...}) passes options to fill_unimaged: the
# Uranian maps have dark shading spikes along the limit and no genuinely black
# terrain, so everything darker than 40 counts as unimaged there.
MAPS = [
    ('PIA07782.jpg', 'jupiter.jpg', False),
    ('io.tif', 'io.jpg', 'polar'),
    ('europa.tif', 'europa.jpg', 'polar'),
    ('ganymede.tif', 'ganymede.jpg', 'polar'),
    ('callisto.tif', 'callisto.jpg', 'polar'),
    ('world.200407.3x5400x2700.jpg', 'earth.jpg', False),
    ('lroc_color_poles_2k.tif', 'moon.jpg', False),
    ('MI_170630_DLR_basemap_degrees.tif', 'mimas.jpg', 'polar'),
    ('Enceladus_Cassini_ISS_Global_Mosaic_100m_HPF.tif', 'enceladus.jpg', 'polar'),
    ('Tethys_Cassini_mosaic_global_293m.tif', 'tethys.jpg', 'polar'),
    ('Dione_Cassini_Voyager_mosaic_global_154m.tif', 'dione.jpg', 'polar'),
    ('Rhea_Cassini_Voyager_mosaic_global_417m.tif', 'rhea.jpg', 'polar'),
    ('Iapetus_Cassini_Voyager_mosaic_global_783m.tif', 'iapetus.jpg', 'polar'),
    ('ura1vuu2.tif', 'ariel.jpg', ('unimaged', {'max_value': 40})),
    ('ura2vuu2.tif', 'umbriel.jpg', ('unimaged', {'max_value': 40})),
    ('ura3vuu2.tif', 'titania.jpg', ('unimaged', {'max_value': 40})),
    ('ura4vuu2.tif', 'oberon.jpg', ('unimaged', {'max_value': 40})),
    ('ura5vuu2.tif', 'miranda.jpg', ('unimaged', {'max_value': 40})),
    ('Triton_Voyager2_ClrMosaic_GlobalFill_600m.tif', 'triton.jpg', 'unimaged'),
    ('Pluto_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif', 'pluto.jpg', 'unimaged'),
    ('Charon_NewHorizons_Global_Mosaic_300m_Jul2017_8bit.tif', 'charon.jpg', 'unimaged'),
    ('Mercury_MESSENGER_mosaic_global_250m_2013.tif', 'mercury.jpg', 'polar'),
    ('Mars_Viking_ClrMosaic_global_925m.tif', 'mars.jpg', False),
    ('Phobos_Viking_Mosaic_40ppd_DLRcontrol.tif', 'phobos.jpg', 'polar'),
    ('Vesta_Dawn_FC_HAMO_Mosaic_Global_74ppd.tif', 'vesta.jpg', False),
    ('Ceres_Dawn_FC_DLR_global_20ppd_Oct2015.tif', 'ceres.jpg', 'polar'),
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


def fill_unimaged(img, max_value=2, rim=15, radius_fraction=0.02, blend=0.3):
    """Fills large no-data areas (pure black) with what is known nearby, fading
    into the map's mean brightness farther away, so that an unimaged hemisphere
    renders as a plain surface instead of a black cap. Normalized convolution:
    blur(image * mask) / blur(mask) is the local mean of the valid pixels."""
    w, h = img.size
    mask = img.convert('L').point(lambda v: 255 if v > max_value else 0)
    # Close isolated black specks inside the imaged area (e.g. Pluto's dark
    # Cthulhu), then drop the dark, ragged rim along the limit of the imaging.
    mask = mask.filter(ImageFilter.MaxFilter(7)).filter(ImageFilter.MinFilter(7))
    mask = mask.filter(ImageFilter.MinFilter(rim))
    radius = max(2, int(w * radius_fraction))
    blurred_mask = mask.filter(ImageFilter.GaussianBlur(radius))
    mean = [int(round(v)) for v in ImageStat.Stat(img, mask).mean]
    bands = []
    for band, m in zip(img.split(), mean):
        known = ImageChops.multiply(band, mask).filter(ImageFilter.GaussianBlur(radius))
        local = ImageMath.lambda_eval(
            lambda a: a['convert'](a['min'](a['k'] * 255 / a['max'](a['bm'], 1), 255), 'L'),
            k=known, bm=blurred_mask)
        # Weight of the local estimate: 1 where enough valid pixels are near, 0 far away.
        weight = blurred_mask.point(lambda v: min(255, int(v / blend)))
        bands.append(Image.composite(local, Image.new('L', (w, h), m), weight))
    fill = Image.merge(img.mode, bands)
    seam = mask.filter(ImageFilter.GaussianBlur(2))
    filled = mask.histogram()[0]
    return Image.composite(img, fill, seam), filled


def prepare(src, dst, width, fill_gaps):
    height = width // 2
    with Image.open(src) as img:
        print(f'{os.path.basename(src)}: {img.size[0]}x{img.size[1]} {img.mode}')
        img = img.convert('RGB') if img.mode not in ('RGB', 'L') else img
        factor = max(1, min(img.size[0] // (width * 2), img.size[1] // (height * 2)))
        if factor > 1:
            img = img.reduce(factor)
        img = img.resize((width, height), Image.Resampling.LANCZOS).convert('RGB')
        filled = 0
        if fill_gaps == 'polar':
            filled = fill_polar_gaps(img)
        elif fill_gaps:
            mode, options = fill_gaps if isinstance(fill_gaps, tuple) else (fill_gaps, {})
            assert mode == 'unimaged', mode
            img, filled = fill_unimaged(img, **options)
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
