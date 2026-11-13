"""Point clouds of the trans-Neptunian objects by dynamical class, for
assets/scenes/trans_neptunian.toml.

Usage (from the repository root):
    python tools/belts/make_tno_classes.py [output_dir]

Resonances: the Deep Ecliptic Survey (DES) class list of M. W. Buie
(https://www.boulder.swri.edu/~buie/kbo/desclass.html; Elliot et al. 2005, AJ 129,
1117), from 10 Myr integrations of each orbit and two clones 3 sigma away. The page
lists every object DES finds in a mean-motion resonance with Neptune and marks the
insecure classes (the three integrations disagree) in italics; only secure ones are
used. One request; the page is cached in tools/belts/cache/des/. (The non-resonant
classes are only given per object, on pages the server does not want crawled.)

The other objects are sorted by their osculating elements (no integration), after
Gladman et al. 2008 (The Solar System Beyond Neptune, 43): e < 0.24 classical; else
perihelion q < 40 au scattered disk, q >= 40 au detached. Orbits with SBDB condition
code 6 or worse are left unclassified.

Orbits: the JPL Small-Body Database (SBDB) Query API, classes TNO and CEN (Neptune's
Trojans can be either), as in make_belts.py: heliocentric osculating J2000 ecliptic
elements moved to one epoch with the two-body mean motion. Objects are matched by
designation (DES writes 2014 RM74 as 14RM74). Centaurs that DES does not place in a
resonance are left out.

Groups written (AXBELT1 files, the format of make_belts.py):
    tno_trojans.bin       1:1
    tno_plutinos.bin      3:2
    tno_twotinos.bin      2:1
    tno_resonant.bin      every other resonance (and objects in two at once)
    tno_classical.bin     not resonant, e < 0.24
    tno_scattered.bin     not resonant, e >= 0.24, q < 40 au
    tno_detached.bin      not resonant, e >= 0.24, q >= 40 au
    tno_unclassified.bin  insecure DES resonances and poorly known orbits
"""

import json
import math
import os
import re
import struct
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(__file__))
import make_belts  # noqa: E402  (SBDB API, epoch and mean motion)

DES_CLASSES = 'https://www.boulder.swri.edu/~buie/kbo/desclass.html'
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache', 'des')
USER_AGENT = 'astraxis-belt-tool (python urllib)'

GROUPS = ['trojans', 'plutinos', 'twotinos', 'resonant', 'classical', 'scattered', 'detached', 'unclassified']
CLASSICAL_MAX_E = 0.24
DETACHED_MIN_Q_AU = 40.0
WORST_CONDITION_CODE = 5


def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                return resp.read().decode('latin-1')
        except urllib.error.HTTPError as exc:
            sys.exit(f'{url}: {exc}')  # refused: do not insist
        except OSError as exc:
            print(f'  retry {attempt + 1} {url}: {exc}')
            time.sleep(5 * (attempt + 1))
    sys.exit(f'fetch failed: {url}')


def des_resonances():
    """{designation: (code, secure)} for the objects DES lists by name."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, 'desclass.html')
    if not os.path.exists(path):
        with open(path, 'w', encoding='latin-1') as f:
            f.write(get(DES_CLASSES))
    with open(path, encoding='latin-1') as f:
        page = f.read()
    generated = re.search(r'Table generated ([^<]*)', page)
    print(f'DES class list: {generated.group(1) if generated else "(no date)"}')
    out = {}
    rows = re.findall(r'<tr>\s*<td>(.*?)</td>\s*<td>\s*(\d+)</td>\s*<td>(.*?)</td>\s*</tr>', page, flags=re.S)
    for code, count, cell in rows:
        for name, italic in re.findall(r'<a href="astrom/([^"]+)\.html">(<i>)?', cell):
            d = designation(name)
            if d:
                out[d] = (code.strip(), not italic)
    return out


def designation(des_name):
    """DES name -> SBDB primary designation: '90377' stays, '14RM74' -> '2014 RM74'."""
    if des_name.isdigit():
        return des_name
    m = re.fullmatch(r'(\d\d)([A-Z]{2})(\d*)', des_name)
    if not m:
        return None
    yy = int(m.group(1))
    return f'{1900 + yy if yy >= 90 else 2000 + yy} {m.group(2)}{m.group(3)}'


def resonance_group(code):
    if ':' not in code:
        return None
    if '+' in code:
        return 'resonant'
    p, q = (int(re.match(r'\d+', s).group(0)) for s in code.split(':'))
    return {(1, 1): 'trojans', (3, 2): 'plutinos', (2, 1): 'twotinos'}.get((p, q), 'resonant')


def fetch(sb_class):
    rows, offset = [], 0
    fields = make_belts.FIELDS + ',pdes,class,condition_code'
    while True:
        params = {'fields': fields, 'sb-class': sb_class, 'limit': str(make_belts.PAGE), 'limit-from': str(offset)}
        d = json.loads(get(make_belts.API + '?' + urllib.parse.urlencode(params)))
        page = [dict(zip(d['fields'], r)) for r in d.get('data', [])]
        rows.extend(page)
        print(f'  SBDB {sb_class}: {len(rows)} / {d["count"]}')
        if len(page) < make_belts.PAGE:
            return rows
        offset += make_belts.PAGE


def elements(r):
    """(H, a, e, i, node, peri, M at the common epoch) or None, as make_belts.main does."""
    try:
        a, e = float(r['a']), float(r['e'])
        i, om, w, ma = (math.radians(float(r[k])) for k in ('i', 'om', 'w', 'ma'))
        epoch, h = float(r['epoch']), float(r['H'])
    except (TypeError, ValueError):
        return None
    if not (a > 0.0 and 0.0 <= e < 1.0):
        return None
    m = math.remainder(ma + make_belts.mean_motion_rad_per_day(a) * (make_belts.EPOCH_JD - epoch),
                       2.0 * math.pi) % (2.0 * math.pi)
    return (h, a, e, i, om, w, m)


def classify(r, el, resonances):
    code, secure = resonances.get(r['pdes'], ('', False))
    group = resonance_group(code) if code else None
    if group:
        return group if secure else 'unclassified'
    if r['class'] != 'TNO':
        return None  # a Centaur outside the resonances
    try:
        cc = int(r['condition_code'])
    except (TypeError, ValueError):
        cc = 9
    if cc > WORST_CONDITION_CODE:
        return 'unclassified'
    _, a, e = el[:3]
    if e < CLASSICAL_MAX_E:
        return 'classical'
    return 'detached' if a * (1.0 - e) >= DETACHED_MIN_Q_AU else 'scattered'


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', '..', 'assets', 'belts')
    os.makedirs(out_dir, exist_ok=True)
    resonances = des_resonances()
    print(f'DES resonant objects: {len(resonances)} ({sum(s for _, s in resonances.values())} secure)')
    groups = {g: [] for g in GROUPS}
    for sb_class in ('TNO', 'CEN'):
        for r in fetch(sb_class):
            el = elements(r)
            if el is None:
                continue
            g = classify(r, el, resonances)
            if g:
                groups[g].append(el)
    for g in GROUPS:
        objects = sorted(groups[g])
        path = os.path.normpath(os.path.join(out_dir, f'tno_{g}.bin'))
        with open(path, 'wb') as f:
            f.write(b'AXBELT1\0')
            f.write(struct.pack('<IId', len(objects), 0, make_belts.EPOCH_JD))
            for h, a, e, i, om, w, m in objects:
                f.write(struct.pack('<7f', a, e, i, om, w, m, h))
        print(f'tno_{g}.bin: {len(objects)} objects ({os.path.getsize(path) // 1024} KB)')


if __name__ == '__main__':
    main()
