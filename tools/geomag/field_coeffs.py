"""Write the planetary field models used by [bodies.magnetic_field].

    python -I tools/geomag/field_coeffs.py <dir with the downloads> assets/geomag

The directory holds ah5.dat, gsfco8.dat, jrm33.dat, cassini11.dat, anderson2012.dat,
weber2022dip.dat from libinternalfield and igrf14coeffs.txt from IAGA; see
assets/geomag/SOURCES.md.
Only the standard library is used.
"""

import sys
from pathlib import Path

LIB = "libinternalfield, commit 9a6b6e5 (MIT, Copyright (c) 2022 Dr Matt James)"

# file stem, output, degree kept, source note.
MODELS = (
    ("ah5", "field_uranus_ah5.txt", 4,
     f"Uranus AH5, Herbert 2009, JGR 114, A11206, via {LIB}. Radius 25,600 km; ULS."),
    ("gsfco8", "field_neptune_o8.txt", 3,
     f"Neptune O8, Connerney et al. 1991, JGR 96, 19023, via {LIB}. Radius 24,765 km; NLS."),
    ("jrm33", "field_jupiter_jrm33.txt", 13,
     f"Jupiter JRM33, Connerney et al. 2022, JGR Planets 127, e2021JE007055, via {LIB}. "
     "Radius 71,492 km; System III."),
    ("cassini11", "field_saturn_cassini11.txt", 11,
     f"Saturn Cassini 11, Dougherty et al. 2018, Science 362, eaat5434, via {LIB}. "
     "Radius 60,268 km; axisymmetric by construction."),
    ("anderson2012", "field_mercury_anderson2012.txt", 4,
     f"Mercury, Anderson et al. 2012, JGR Planets 117, E00L12, via {LIB}. Radius 2,440 km; axisymmetric."),
    ("weber2022dip", "field_ganymede_weber2022.txt", 1,
     f"Ganymede dipole, Weber et al. 2022, GRL 49, e2022GL098633, via {LIB}. Radius 2,631 km; "
     "longitude 0 faces Jupiter, positive toward the trailing side."),
)


def header(note):
    return (f"# {note}\n"
            "# Schmidt quasi-normalized internal Gauss coefficients in nT, model frame right-handed\n"
            "# (longitude positive east of the model's prime meridian). Written by\n"
            "# tools/geomag/field_coeffs.py; see SOURCES.md.\n")


def library_model(path, degree):
    rows = []
    for line in path.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in ("g", "h") and int(parts[1]) <= degree:
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    if not any(r[:3] == ("g", 1, 0) for r in rows):
        sys.exit(f"{path}: no dipole term")
    return rows


def igrf_series(path):
    """The epochs, the year the secular variation runs to, and rows (g|h, n, m, values, rate)."""
    rows = []
    epochs = None
    for line in path.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if parts[:3] == ["g/h", "n", "m"]:
            # "1900.0 .. 2025.0", then the secular variation as "2025-30".
            epochs = parts[3:-1]
            first, last = parts[-1].split("-")
            until = float(first[:2] + last)
            if float(epochs[-1]) != float(first) or not until > float(first):
                sys.exit(f"{path}: unexpected secular variation column {parts[-1]}")
        elif epochs is not None and parts and parts[0] in ("g", "h"):
            if len(parts) != len(epochs) + 4:
                sys.exit(f"{path}: bad row {line}")
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3:-1], parts[-1]))
    if epochs is None or not rows:
        sys.exit(f"{path}: no coefficients")
    return epochs, until, rows


def write(path, note, rows):
    text = header(note) + "".join(f"{c} {n} {m} {v}\n" for c, n, m, v in rows)
    path.write_text(text, encoding="ascii", newline="\n")
    print(f"{path}: {len(rows)} coefficients")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source, out = Path(sys.argv[1]), Path(sys.argv[2])
    for stem, name, degree, note in MODELS:
        write(out / name, note, library_model(source / f"{stem}.dat", degree))
    epochs, until, rows = igrf_series(source / "igrf14coeffs.txt")
    path = out / "field_earth_igrf14.txt"
    text = ("# Earth IGRF-14, IAGA Working Group V-MOD "
            "(https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt). Radius 6,371.2 km; geocentric.\n"
            "# Schmidt quasi-normalized internal Gauss coefficients in nT at each epoch (linear between\n"
            "# them), then their rate in nT per year from the last epoch until rate_until; model frame\n"
            "# right-handed (longitude positive east). Written by tools/geomag/field_coeffs.py; see SOURCES.md.\n"
            f"epochs {' '.join(epochs)}\nrate_until {until:.1f}\n"
            + "".join(f"{c} {n} {m} {' '.join(v)} {rate}\n" for c, n, m, v, rate in rows))
    path.write_text(text, encoding="ascii", newline="\n")
    print(f"{path}: {len(rows)} coefficients at {len(epochs)} epochs, {epochs[0]} .. {epochs[-1]}, rates until {until:.1f}")


if __name__ == "__main__":
    main()
