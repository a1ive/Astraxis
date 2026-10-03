"""Render the Astraxis logo SVGs to a multi-resolution Windows .ico.

Sources:
    res/astraxis.svg        full-detail mark, used for sizes >= 48 px
    res/astraxis-small.svg  simplified bold mark, used for sizes <= 40 px

Usage:
    python make_icon.py [--inkscape <path/to/inkscape>] [--out res/astraxis.ico]

Requires Inkscape (CLI) for rasterization; no other dependencies.
Every frame is stored PNG-compressed (supported by Windows Vista and later).
"""

import argparse
import os
import shutil
import struct
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SIZES = [16, 20, 24, 32, 40, 48, 64, 128, 256]
SMALL_MAX = 40


def find_inkscape(explicit):
    if explicit:
        return explicit
    found = shutil.which('inkscape')
    if found:
        return found
    default = r'C:\Program Files\Inkscape\bin\inkscape.exe'
    if os.path.exists(default):
        return default
    sys.exit('Inkscape not found; pass --inkscape')


def render(inkscape, svg, size, png):
    subprocess.run([inkscape, svg, '--export-type=png', f'--export-filename={png}',
                    '-w', str(size), '-h', str(size)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def write_ico(frames, dst):
    """frames: list of (size, png_bytes), written as an ICONDIR with PNG payloads."""
    header = struct.pack('<HHH', 0, 1, len(frames))
    offset = len(header) + 16 * len(frames)
    entries, payload = b'', b''
    for size, data in frames:
        dim = 0 if size >= 256 else size  # 0 means 256 in ICONDIRENTRY
        entries += struct.pack('<BBBBHHII', dim, dim, 0, 0, 1, 32, len(data), offset)
        payload += data
        offset += len(data)
    with open(dst, 'wb') as f:
        f.write(header + entries + payload)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inkscape')
    ap.add_argument('--out', default=os.path.join(ROOT, 'res', 'astraxis.ico'))
    args = ap.parse_args()

    inkscape = find_inkscape(args.inkscape)
    full = os.path.join(ROOT, 'res', 'astraxis.svg')
    small = os.path.join(ROOT, 'res', 'astraxis-small.svg')

    frames = []
    with tempfile.TemporaryDirectory() as tmp:
        for size in SIZES:
            png = os.path.join(tmp, f'{size}.png')
            render(inkscape, small if size <= SMALL_MAX else full, size, png)
            with open(png, 'rb') as f:
                frames.append((size, f.read()))
    write_ico(frames, args.out)
    print(f'wrote {args.out} ({", ".join(str(s) for s in SIZES)})')


if __name__ == '__main__':
    main()
