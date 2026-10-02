"""Convert the NASA SVS Deep Star Maps 2020 Milky Way layer to assets/textures/milky_way.jpg.

Source and license: assets/textures/SOURCES.md. Download milkyway_2020_4k.exr
there first, then:

    python convert_milkyway.py <milkyway_2020_4k.exr> <output.jpg> [--stats]

The EXR is plate carree in ICRF (J2000) right ascension / declination, linear
half-float RGB. The SVS page does not state its units; --stats shows the
values already span 0..1 (maximum exactly 1.0), so they are divided by WHITE
= 1.0 and sRGB encoded. The renderer multiplies by a per-scene brightness, so
WHITE only sets the quantization range, not the look.

Only the standard library and Pillow are needed: the reader handles the subset
of OpenEXR these files use (scanline, ZIP compression, HALF channels).
"""

import struct
import sys
import zlib

from PIL import Image

# Linear value mapped to 1.0 (see --stats and SOURCES.md).
WHITE = 1.0


def read_attributes(data):
    if data[:4] != b'\x76\x2f\x31\x01':
        raise ValueError('not an OpenEXR file')
    i = 8
    attrs = {}
    while True:
        j = data.index(b'\0', i)
        name = data[i:j].decode()
        if not name:
            return attrs, j + 1
        k = data.index(b'\0', j + 1)
        typ = data[j + 1:k].decode()
        size = struct.unpack_from('<i', data, k + 1)[0]
        attrs[name] = (typ, data[k + 5:k + 5 + size])
        i = k + 5 + size


def read_exr(path):
    """Returns (width, height, channels) with channels = {name: [float] * w * h}, top row first."""
    data = open(path, 'rb').read()
    attrs, pos = read_attributes(data)

    compression = attrs['compression'][1][0]
    if compression != 3:
        raise ValueError(f'unsupported compression {compression} (only ZIP = 3)')
    xmin, ymin, xmax, ymax = struct.unpack('<4i', attrs['dataWindow'][1])
    width = xmax - xmin + 1
    height = ymax - ymin + 1

    names = []
    chlist = attrs['channels'][1]
    c = 0
    while chlist[c] != 0:
        e = chlist.index(b'\0', c)
        name = chlist[c:e].decode()
        pixel_type = struct.unpack_from('<i', chlist, e + 1)[0]
        if pixel_type != 1:
            raise ValueError(f'channel {name}: only HALF is supported')
        names.append(name) # stored in this (alphabetical) order
        c = e + 1 + 16

    lines_per_block = 16
    blocks = (height + lines_per_block - 1) // lines_per_block
    offsets = struct.unpack_from(f'<{blocks}Q', data, pos)
    channels = {n: [0.0] * (width * height) for n in names}

    for offset in offsets:
        y, size = struct.unpack_from('<ii', data, offset)
        raw = zlib.decompress(data[offset + 8:offset + 8 + size])
        # Undo the predictor, then the even/odd byte split.
        t = bytearray(raw)
        for i in range(1, len(t)):
            t[i] = (t[i - 1] + t[i] - 128) & 0xFF
        half = (len(t) + 1) // 2
        out = bytearray(len(t))
        out[0::2] = t[:half]
        out[1::2] = t[half:]

        rows = min(lines_per_block, ymax + 1 - y)
        values = struct.unpack(f'<{len(out) // 2}e', out)
        k = 0
        for row in range(rows):
            base = (y - ymin + row) * width
            for n in names:
                channels[n][base:base + width] = values[k:k + width]
                k += width
    return width, height, channels


def srgb_encode(x):
    x = min(max(x, 0.0), 1.0)
    return 12.92 * x if x <= 0.0031308 else 1.055 * x ** (1.0 / 2.4) - 0.055


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    width, height, ch = read_exr(sys.argv[1])
    r, g, b = ch['R'], ch['G'], ch['B']
    print(f'{width}x{height}')

    if '--stats' in sys.argv:
        lum = sorted(0.2126 * r[i] + 0.7152 * g[i] + 0.0722 * b[i] for i in range(0, width * height, 7))
        for p in (0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 1.0):
            print(f'luminance p{p * 100:g}: {lum[min(int(p * len(lum)), len(lum) - 1)]:.6g}')

    lut = [round(255.0 * srgb_encode(i / 65535.0)) for i in range(65536)]
    pixels = bytearray(width * height * 3)
    scale = 65535.0 / WHITE
    for i in range(width * height):
        pixels[3 * i] = lut[min(max(int(r[i] * scale + 0.5), 0), 65535)]
        pixels[3 * i + 1] = lut[min(max(int(g[i] * scale + 0.5), 0), 65535)]
        pixels[3 * i + 2] = lut[min(max(int(b[i] * scale + 0.5), 0), 65535)]
    Image.frombytes('RGB', (width, height), bytes(pixels)).save(sys.argv[2], quality=92, subsampling=0)
    print(f'wrote {sys.argv[2]}')


if __name__ == '__main__':
    main()
