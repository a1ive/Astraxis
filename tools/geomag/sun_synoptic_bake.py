"""Write the Sun's photospheric field, one harmonic set per Carrington rotation.

    python -I tools/geomag/sun_synoptic_bake.py [--merge] <dir with the maps> assets/geomag/sun_hmi_synoptic.txt

With --merge, the rotations already in the output are kept and those of the maps given
replace them or are added (tools/refresh/monthly_refresh.py bakes only the new ones).

The directory holds hmi.Synoptic_Mr_polfil.<CR>.fits from
http://jsoc.stanford.edu/data/hmi/synoptic/ (SDO/HMI radial-field synoptic maps with
the poles filled, CC0 1.0); see assets/geomag/SOURCES.md. Each map (3600 Carrington
longitudes by 1440 steps in sine latitude, Mx/cm^2 = G) is averaged to 1 degree by
1/180 in sine latitude and projected on Schmidt quasi-normalized harmonics:

    g'_nm, h'_nm = (2n + 1) / (X Y) sum Br P_n^m(sin lat) (cos, sin)(m lon)

so that Br at the photosphere is sum P_n^m (g' cos m lon + h' sin m lon), the
convention of the Stanford PFSS (Sun, WSO pfss.pdf, eq. 10 and 15). The monopole
(the map's net flux) is left out. Only the standard library is used.
"""

import math
import multiprocessing
import sys
from pathlib import Path

DEGREE = 20
LON_BINS = 360
SINLAT_BINS = 180
PIXEL_DEG = 0.1  # CRPIX1 = 1800 at 180 deg, longitude rising 0.1 deg per pixel: pixel k (0-based) at 0.1 (k + 1)


def cards(data, start):
    """Header cards from `start` to END; returns (dict, offset of the data)."""
    out = {}
    i = start
    while True:
        card = data[i:i + 80].decode("ascii")
        i += 80
        if card.startswith("END"):
            break
        if card[8:10] == "= ":
            value = card[10:].strip()
            # A quoted string may hold a slash; otherwise one starts the comment.
            value = value[1:value.index("'", 1)].strip() if value.startswith("'") else value.split("/")[0].strip()
            out[card[:8].strip()] = value
    return out, (i + 2879) // 2880 * 2880


def rice_decode(tile, count, block=32):
    """FITS RICE_1 for 32-bit pixels (cfitsio fits_rdecomp)."""
    last = int.from_bytes(tile[:4], "big", signed=True)
    nbits = (len(tile) - 4) * 8
    bits = bin(int.from_bytes(tile[4:], "big") | (1 << nbits))[3:] if nbits else ""
    p = 0
    out = []
    while len(out) < count:
        fs = int(bits[p:p + 5], 2) - 1
        p += 5
        n = min(block, count - len(out))
        if fs < 0:
            out.extend([last] * n)
            continue
        for _ in range(n):
            if fs == 25:
                diff = int(bits[p:p + 32], 2)
                p += 32
            else:
                one = bits.find("1", p)
                diff = (one - p) << fs
                p = one + 1
                if fs:
                    diff |= int(bits[p:p + fs], 2)
                    p += fs
            diff = diff >> 1 if diff & 1 == 0 else ~(diff >> 1)
            last = (last + diff + 0x80000000) % 0x100000000 - 0x80000000
            out.append(last)
    return out


def read_map(path):
    """The map averaged to LON_BINS x SINLAT_BINS (G), south row first; and the blank fraction."""
    data = path.read_bytes()
    _, offset = cards(data, 0)
    head, offset = cards(data, offset)
    assert head["XTENSION"] == "BINTABLE" and head["ZCMPTYPE"] == "RICE_1" and head["ZBITPIX"] == "32"
    assert head["CTYPE1"] == "CRLN-CEA" and head["CTYPE2"] == "CRLT-CEA" and head["BUNIT"] == "Mx/cm^2"
    nx, ny = int(head["ZNAXIS1"]), int(head["ZNAXIS2"])
    assert (nx, ny) == (3600, 1440) and int(head["ZTILE2"]) == 1
    assert float(head["CRPIX1"]) == 1800.0 and abs(float(head["CDELT1"]) + 0.1) < 1e-9
    assert float(head["CRPIX2"]) == 720.5
    scale, zero = float(head.get("BSCALE", "1")), float(head.get("BZERO", "0"))
    blank = int(head.get("BLANK", "-2147483648"))
    row_bytes = int(head["NAXIS1"])
    heap = offset + int(head.get("THEAP", row_bytes * ny))
    sx, sy = nx // LON_BINS, ny // SINLAT_BINS
    sums = [[0.0] * LON_BINS for _ in range(SINLAT_BINS)]
    counts = [[0] * LON_BINS for _ in range(SINLAT_BINS)]
    blanks = 0
    for row in range(ny):
        at = offset + row * row_bytes
        size = int.from_bytes(data[at:at + 4], "big")
        start = int.from_bytes(data[at + 4:at + 8], "big")
        pixels = rice_decode(data[heap + start:heap + start + size], nx, int(head["ZVAL1"]))
        s, c = sums[row // sy], counts[row // sy]
        for k, v in enumerate(pixels):
            if v == blank:
                blanks += 1
                continue
            s[k // sx] += v * scale + zero
            c[k // sx] += 1
    grid = [[s / c if c else 0.0 for s, c in zip(srow, crow)] for srow, crow in zip(sums, counts)]
    return int(head["CAR_ROT"]), grid, blanks / (nx * ny)


def schmidt(degree, x):
    """P[n][m](x), Schmidt quasi-normalized, without the Condon-Shortley phase."""
    s = math.sqrt(max(1.0 - x * x, 0.0))
    p = [[0.0] * (degree + 1) for _ in range(degree + 1)]
    p[0][0] = 1.0
    for n in range(1, degree + 1):
        for m in range(n):
            above = math.sqrt((n - 1) ** 2 - m * m) * p[n - 2][m] if n - 1 > m else 0.0
            p[n][m] = ((2 * n - 1) * x * p[n - 1][m] - above) / math.sqrt(n * n - m * m)
        p[n][n] = s * p[n - 1][n - 1] * (1.0 if n == 1 else math.sqrt((2 * n - 1) / (2 * n)))
    return p


def coefficients(grid):
    # The rows' Fourier sums first: a[j][m], b[j][m] = sum over longitude of Br (cos, sin)(m lon).
    # Bin i averages pixels k = s i .. s i + s - 1, at 0.1 (k + 1) degrees.
    s = 3600 // LON_BINS
    lons = [math.radians(PIXEL_DEG * (s * i + (s + 1) / 2.0)) for i in range(LON_BINS)]
    trig = [[(math.cos(m * l), math.sin(m * l)) for l in lons] for m in range(DEGREE + 1)]
    g = {}
    h = {}
    for j, row in enumerate(grid):
        x = -1.0 + (j + 0.5) * 2.0 / SINLAT_BINS
        p = schmidt(DEGREE, x)
        for m in range(DEGREE + 1):
            a = sum(v * c for v, (c, _) in zip(row, trig[m]))
            b = sum(v * s for v, (_, s) in zip(row, trig[m]))
            for n in range(max(m, 1), DEGREE + 1):
                g[n, m] = g.get((n, m), 0.0) + a * p[n][m]
                h[n, m] = h.get((n, m), 0.0) + b * p[n][m]
    norm = 1.0 / (LON_BINS * SINLAT_BINS)
    return {k: (2 * k[0] + 1) * norm * v for k, v in g.items()}, {k: (2 * k[0] + 1) * norm * v for k, v in h.items()}


def bake(path):
    rotation, grid, blank = read_map(path)
    net = sum(sum(row) for row in grid) / (LON_BINS * SINLAT_BINS)
    g, h = coefficients(grid)
    values = []
    for n in range(1, DEGREE + 1):
        values.append(g[n, 0])
        for m in range(1, n + 1):
            values += [g[n, m], h[n, m]]
    return rotation, values, blank, net


HEADER = ("# The Sun's photospheric radial field from SDO/HMI synoptic maps with the poles filled\n"
          "# (hmi.Synoptic_Mr_polfil_720s, Sun et al. 2011; data CC0 1.0, courtesy of NASA/SDO and\n"
          "# the HMI science team). A value-added product of this project, not an HMI data product.\n"
          "# Per Carrington rotation: g'_nm, h'_nm (G) with Br(R_sun) = sum P_n^m(sin lat)\n"
          "# (g' cos m lon + h' sin m lon), Schmidt quasi-normalized, Carrington longitude east;\n"
          "# monopole left out. Order: n = 1.. degree; g'_n0, then g'_nm h'_nm for m = 1..n.\n"
          "# Written by tools/geomag/sun_synoptic_bake.py; see SOURCES.md.\n")


def read_rows(path):
    """{rotation: the line's values as written} of an existing output."""
    rows = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith("degree"):
                if int(line.split()[1]) != DEGREE:
                    sys.exit(f"{path}: degree {line.split()[1]}, the tool writes {DEGREE}")
                continue
            rotation, values = line.split(maxsplit=1)
            rows[int(rotation)] = values.strip()
    return rows


def main():
    args = sys.argv[1:]
    merge = "--merge" in args
    if merge:
        args.remove("--merge")
    if len(args) != 2:
        sys.exit(__doc__)
    maps = sorted(Path(args[0]).glob("hmi.Synoptic_Mr_polfil.*.fits"))
    if not maps:
        sys.exit(f"no maps in {args[0]}")
    with multiprocessing.Pool() as pool:
        baked = sorted(pool.map(bake, maps))
    rows = read_rows(args[1]) if merge else {}
    for rotation, values, _, _ in baked:
        rows[rotation] = " ".join(f"{v:.5g}" for v in values)
    rotations = sorted(rows)
    if rotations != list(range(rotations[0], rotations[-1] + 1)):
        sys.exit(f"rotations missing between {rotations[0]} and {rotations[-1]}")
    with open(args[1], "w", newline="\n") as out:
        out.write(HEADER)
        out.write(f"degree {DEGREE}\n")
        for rotation in rotations:
            out.write(f"{rotation} {rows[rotation]}\n")
    for rotation, values, blank, net in baked:
        if blank > 0.0 or abs(net) > 0.2:
            print(f"CR{rotation}: blank {blank:.4f}, net flux {net:+.3f} G")
    print(f"baked CR{baked[0][0]}..{baked[-1][0]}; the file holds {len(rotations)} rotations, "
          f"CR{rotations[0]}..{rotations[-1]}; g10 of the last baked {baked[-1][1][0]:+.3f} G")


if __name__ == "__main__":
    main()
