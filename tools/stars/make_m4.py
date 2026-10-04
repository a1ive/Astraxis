"""Build the star field of the globular cluster M4 (NGC 6121) for the PSR B1620-26 scene.

Members: Gaia DR3 (CDS catalogue I/355/gaiadr3, queried through VizieR; ESA/Gaia/DPAC,
licence CC BY-NC 3.0 IGO, so the output is for non-commercial use only) within
20' of the cluster centre, G < 19, proper motion within 1.5 mas/yr of the
cluster's and parallax within 3 sigma of it.

Sources
  [BV21] Baumgardt & Vasiliev 2021, MNRAS 505, 5957, via the online table
         https://people.smp.uq.edu.au/HolgerBaumgardt/globular/orbits_table.txt:
         centre RA 245.896744, Dec -26.525749 deg; distance 1.85 kpc;
         proper motion (-12.511, -19.017) mas/yr.
  [HA10] Harris 1996 (2010 edition), https://physics.mcmaster.ca/~harris/mwgc.dat:
         King concentration c = 1.65, core radius r_c = 1.16'.
  [SI03] Sigurdsson et al. 2003, Science 301, 193, Table 1 note: A(F555W) = 1.31,
         taken as A_V (M4's extinction law is unusual, R_V ~ 3.7).
  [WC19] Wang & Chen 2019, ApJ 877, 116, Table 3: A_G / A_V = 0.789,
         A_BP / A_V = 1.002, A_RP / A_V = 0.589.
  [PM13] Pecaut & Mamajek 2013, ApJS 208, 9, online table
         https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
         (version 2022.04.16): G - V and B - V as functions of Bp - Rp (dwarfs;
         used for giants too, which is approximate).

Positions on the sky are Gaia's; the depth of each star along the line of sight
is not measured and is drawn from a King (1962) density profile at its
projected radius (fixed seed): the 3D structure is statistical, not real.
Stars fainter than G = 19 (M_V ~ 6.5, below ~0.65 M_sun) and members lost to
crowding in the core are missing.

Usage:
    python make_m4.py <output.csv>

Output columns: east_pc, north_pc, away_pc (relative to the cluster centre, on
the tangent-plane axes at the centre; away = along the line of sight), abs_vmag
(absolute V, extinction removed), bv (intrinsic B - V).
"""

import math
import random
import sys
import urllib.request

CENTER_RA = 245.896744
CENTER_DEC = -26.525749
DISTANCE_PC = 1850.0
PM_RA, PM_DEC = -12.511, -19.017
PM_RADIUS = 1.5            # mas/yr
RADIUS_ARCMIN = 20.0
G_LIMIT = 19.0
KING_C = 1.65
RC_ARCMIN = 1.16
A_V = 1.31
A_G = 0.789 * A_V
E_BP_RP = (1.002 - 0.589) * A_V
SEED = 20261004

VIZIER = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/355/gaiadr3'
          f'&-c={CENTER_RA}%20{CENTER_DEC}&-c.rm={RADIUS_ARCMIN}&-out.max=unlimited'
          '&-out=RA_ICRS,DE_ICRS,Gmag,BP-RP,Plx,e_Plx,pmRA,pmDE'
          f'&pmRA={PM_RA - PM_RADIUS}..{PM_RA + PM_RADIUS}'
          f'&pmDE={PM_DEC - PM_RADIUS}..{PM_DEC + PM_RADIUS}&Gmag=%3C{G_LIMIT}')
PM13 = 'https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt'


def fetch(url):
    with urllib.request.urlopen(url, timeout=600) as r:
        return r.read().decode('utf-8', errors='replace')


def gaia_members():
    rows = []
    header = None
    for line in fetch(VIZIER).splitlines():
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if header is None:
            header = [c.strip() for c in cols]
            continue
        if cols[0].startswith('-') or cols[0].strip() in ('deg',):
            continue
        rec = dict(zip(header, (c.strip() for c in cols)))
        if not rec.get('BP-RP'):
            continue
        pm_ra, pm_dec = float(rec['pmRA']), float(rec['pmDE'])
        if math.hypot(pm_ra - PM_RA, pm_dec - PM_DEC) > PM_RADIUS:
            continue
        if rec.get('Plx') and rec.get('e_Plx'):
            if abs(float(rec['Plx']) - 1000.0 / DISTANCE_PC) > 3.0 * float(rec['e_Plx']):
                continue
        rows.append((float(rec['RA_ICRS']), float(rec['DE_ICRS']), float(rec['Gmag']), float(rec['BP-RP'])))
    return rows


def pm13_colors():
    """(Bp-Rp, G-V, B-V) rows of the dwarf table where all three are given."""
    out = []
    lines = fetch(PM13).splitlines()
    header = None
    for line in lines:
        if line.startswith('#SpT'):
            header = line[1:].split()
            continue
        if header is None or not line.strip() or line.startswith('#'):
            if header is not None and out and not line.strip():
                break  # end of the first table
            continue
        cols = line.split()
        if len(cols) < len(header):
            continue
        rec = dict(zip(header, cols))
        try:
            out.append((float(rec['Bp-Rp']), float(rec['G-V']), float(rec['B-V'])))
        except (KeyError, ValueError):
            continue
    out.sort()
    return out


def interp(table, x, k):
    if x <= table[0][0]:
        return table[0][k]
    if x >= table[-1][0]:
        return table[-1][k]
    for a, b in zip(table, table[1:]):
        if a[0] <= x <= b[0]:
            f = (x - a[0]) / (b[0] - a[0]) if b[0] > a[0] else 0.0
            return a[k] + f * (b[k] - a[k])
    return table[-1][k]


def king_density(r, rc, rt):
    """King (1962) volume density, unnormalised."""
    if r >= rt:
        return 0.0
    z = math.sqrt((1.0 + (r / rc) ** 2) / (1.0 + (rt / rc) ** 2))
    return (math.acos(z) / z - math.sqrt(1.0 - z * z)) / (z * z)


def sample_depth(rng, R, rc, rt):
    """Line-of-sight offset at projected radius R, drawn from the King profile."""
    if R >= rt:
        return 0.0
    half = math.sqrt(rt * rt - R * R)
    peak = king_density(R, rc, rt)
    while True:
        l = rng.uniform(-half, half)
        if rng.random() * peak <= king_density(math.hypot(R, l), rc, rt):
            return l


def main(dst):
    members = gaia_members()
    colors = pm13_colors()
    pc_per_rad = DISTANCE_PC
    rc = RC_ARCMIN / 60.0 * math.pi / 180.0 * pc_per_rad
    rt = rc * 10.0 ** KING_C
    ra0 = math.radians(CENTER_RA)
    dec0 = math.radians(CENTER_DEC)
    rng = random.Random(SEED)
    dist_mod = 5.0 * math.log10(DISTANCE_PC / 10.0)

    rows = []
    for ra, dec, g, bp_rp in members:
        # Gnomonic projection onto the tangent plane at the centre.
        a, d = math.radians(ra), math.radians(dec)
        cosc = math.sin(dec0) * math.sin(d) + math.cos(dec0) * math.cos(d) * math.cos(a - ra0)
        xi = math.cos(d) * math.sin(a - ra0) / cosc
        eta = (math.cos(dec0) * math.sin(d) - math.sin(dec0) * math.cos(d) * math.cos(a - ra0)) / cosc
        east, north = xi * pc_per_rad, eta * pc_per_rad
        away = sample_depth(rng, math.hypot(east, north), rc, rt)
        g0 = g - A_G
        c0 = bp_rp - E_BP_RP
        vmag = g0 - interp(colors, c0, 1) - dist_mod
        bv = interp(colors, c0, 2)
        rows.append((east, north, away, vmag, bv))

    with open(dst, 'w', encoding='ascii', newline='\n') as out:
        out.write('# M4 (NGC 6121) members from Gaia DR3 (CDS I/355), generated by tools/stars/make_m4.py.\n')
        out.write('# Gaia data: ESA/Gaia/DPAC, licence CC BY-NC 3.0 IGO (non-commercial use only).\n')
        out.write(f'# Centre RA {CENTER_RA}, Dec {CENTER_DEC} deg, distance {DISTANCE_PC:.0f} pc (Baumgardt & Vasiliev 2021).\n')
        out.write(f'# Depths drawn from a King profile (c = {KING_C}, r_c = {rc:.3f} pc, r_t = {rt:.2f} pc; Harris 2010).\n')
        out.write(f'# Extinction removed: A_V = {A_V}, A_G = {A_G:.3f}, E(BP-RP) = {E_BP_RP:.3f}; V and B-V via Pecaut & Mamajek 2013.\n')
        out.write('east_pc,north_pc,away_pc,abs_vmag,bv\n')
        for e, n, a, v, bv in rows:
            out.write(f'{e:.4f},{n:.4f},{a:.4f},{v:.3f},{bv:.3f}\n')
    print(f'{len(rows)} members written to {dst} (r_c = {rc:.3f} pc, r_t = {rt:.2f} pc)')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1])
