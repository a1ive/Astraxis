"""Linear ephemeris of Algol's primary eclipses from TESS light curves.

Usage:
    python algol_ephemeris.py <TESS SPOC *_lc.fits> ...

Needs numpy. Reads the SPOC 2-minute light curves of TIC 346783960 (sectors 18,
58 and 85; download addresses in assets/ephem/SOURCES.md), finds each primary
minimum as the time of mirror symmetry of the SAP flux within +-0.25 d of the
deepest point (minima with gaps in that window are skipped), removes the
light-time effect of the eclipsing pair's orbit around the barycenter with
Algol C, and fits a linear ephemeris to the result. The scene
(assets/scenes/algol.toml) puts the eclipses on that ephemeris; its outer orbit
then puts the light-time effect back for an observer (tests/stellar_system_tests.cpp
compares the scene with the minima printed here).

Outer orbit: Baron et al. 2012, ApJ 752, 20, Table 5 (C relative to AB), with
the masses of their Table 4 (m_C / (m_A + m_B + m_C) = 1.76 / 5.63). The depth
convention is the scene's (VisualOrbit): the companion recedes at its
ascending node.
"""

import sys

import numpy as np

AU_KM = 149597870.7
C_KM_S = 299792.458
PARALLAX_MAS = 34.7          # Baron et al. 2012 (Zavala et al. 2010, Peterson et al. 2011)
OUTER = dict(t_peri_jd=2446927.22, period_d=680.168, a_mas=93.43, i_deg=83.66, e=0.227, w_deg=310.02)
C_FRACTION = 1.76 / (3.17 + 0.70 + 1.76)
TIDAK = (2445641.5135, 2.86730431)  # HJD, d (Jetsu 2021, ApJ 920, 137, eq. 1): for E numbers only

FITS_TYPES = {'D': '>f8', 'E': '>f4', 'J': '>i4', 'I': '>i2', 'K': '>i8', 'B': 'u1', 'L': 'u1', 'A': 'S1'}


def read_header(f):
    cards = {}
    while True:
        block = f.read(2880).decode('ascii')
        for i in range(0, 2880, 80):
            card = block[i:i + 80]
            if card.startswith('END'):
                return cards
            if card[8:10] == '= ':
                cards[card[:8].strip()] = card[10:].split('/')[0].strip().strip("'").strip()


def read_light_curve(path):
    """TIME (BTJD), SAP_FLUX and QUALITY of a SPOC light curve file."""
    with open(path, 'rb') as f:
        read_header(f)
        h = read_header(f)
        rows, row_bytes = int(h['NAXIS2']), int(h['NAXIS1'])
        columns, offset = {}, 0
        for k in range(1, int(h['TFIELDS']) + 1):
            form = h['TFORM%d' % k]
            count = int(form[:-1]) if len(form) > 1 else 1
            dtype = np.dtype(FITS_TYPES[form[-1]])
            columns[h['TTYPE%d' % k]] = (offset, dtype)
            offset += dtype.itemsize * count
        raw = np.frombuffer(f.read(rows * row_bytes), dtype='u1').reshape(rows, row_bytes)
    if h.get('TIMESYS') != 'TDB' or int(h.get('BJDREFI', 0)) != 2457000:
        sys.exit(f'{path}: unexpected time system')
    out = {}
    for name in ('TIME', 'SAP_FLUX', 'QUALITY'):
        offset, dtype = columns[name]
        out[name] = raw[:, offset:offset + dtype.itemsize].copy().view(dtype).ravel().astype(float)
    return out


def symmetry_time(t, f, guess, half=0.25):
    near = np.abs(t - guess) < half + 0.02
    t, f = t[near], f[near]
    if len(t) < 2 or np.max(np.diff(t)) > 0.01 or t[0] > guess - half or t[-1] < guess + half:
        return None  # a gap or an edge inside the window
    dt = np.linspace(0.01, half, 120)
    best = None
    for t0 in np.arange(guess - 0.03, guess + 0.03, 0.0001):
        s = np.mean((np.interp(t0 + dt, t, f) - np.interp(t0 - dt, t, f)) ** 2)
        if best is None or s < best[0]:
            best = (s, t0)
    return best[1]


def kepler_e(m, e):
    ea = m.copy()
    for _ in range(50):
        ea -= (ea - e * np.sin(ea) - m) / (1.0 - e * np.cos(ea))
    return ea


def light_time_days(jd):
    """Delay of the eclipsing pair's light by its depth around the barycenter."""
    o = OUTER
    a_km = o['a_mas'] / PARALLAX_MAS * AU_KM
    m = 2.0 * np.pi * (jd - o['t_peri_jd']) / o['period_d']
    ea = kepler_e(np.mod(m, 2.0 * np.pi), o['e'])
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + o['e']) * np.sin(ea / 2.0), np.sqrt(1.0 - o['e']) * np.cos(ea / 2.0))
    r = a_km * (1.0 - o['e'] * np.cos(ea))
    depth_c = r * np.sin(nu + np.radians(o['w_deg'])) * np.sin(np.radians(o['i_deg']))  # C relative to AB
    return -C_FRACTION * depth_c / C_KM_S / 86400.0


def main(paths):
    minima = []
    for path in paths:
        lc = read_light_curve(path)
        t, f, q = lc['TIME'], lc['SAP_FLUX'], lc['QUALITY']
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        t, f = t[ok], f[ok] / np.median(f[ok])
        picked = []
        for i in np.argsort(f):
            if f[i] > 0.75:
                break
            if all(abs(t[i] - p) > 1.0 for p in picked):
                picked.append(t[i])
        for guess in sorted(picked):
            t0 = symmetry_time(t, f, guess)
            if t0 is not None:
                minima.append(t0 + 2457000.0)
    jd = np.array(minima)
    epoch = np.round((jd - TIDAK[0]) / TIDAK[1])
    geometric = jd - light_time_days(jd)
    a = np.vstack([np.ones_like(epoch), epoch - np.round(epoch.mean())]).T
    (t_mid, period), *_ = np.linalg.lstsq(a, geometric, rcond=None)
    residual = geometric - a @ np.array([t_mid, period])
    e_mid = int(np.round(epoch.mean()))
    print('# observed BJD_TDB    light time (d)  residual (d)')
    for x, lt, res in zip(jd, light_time_days(jd), residual):
        print(f'{x:.5f}  {lt:+.5f}  {res:+.5f}')
    print(f'# {len(jd)} minima, rms {np.sqrt(np.mean(residual ** 2)) * 1440:.1f} min')
    print(f'# geometric primary minimum BJD_TDB {t_mid:.5f} (TIDAK E = {e_mid}), period {period:.7f} d')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
