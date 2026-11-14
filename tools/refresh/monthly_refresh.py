"""Monthly refresh of the data that ends in a forecast or at the latest measurement
(run by .github/workflows/refresh-data.yml; works the same by hand).

    python -I tools/refresh/monthly_refresh.py [--skip-omni] [--skip-ephem] [--body <file>]

1. OMNI2: downloads the yearly hourly files into build/omni (years already there are kept,
   the last three are always fetched again: recent values get revised), rebakes
   assets/geomag/solar_wind.bin with tools/geomag/solar_wind_bake.py and updates its
   numbers in assets/geomag/SOURCES.md.
2. Horizons: for each telescope in TARGETS, reads the end of its coverage from Horizons,
   moves `stop` in tools/bake/ephem.toml to the last whole hour at least a minute inside
   it, rebakes with tools/bake/horizons_bake.py and updates the table row in
   assets/ephem/SOURCES.md.
3. With --body, writes a pull request description: what changed, the trajectory notes
   Horizons gives for each target, and what is left to check by hand.

Every documented number it rewrites must be found exactly once; if a page was reworded,
the script stops instead of leaving the data and the page apart. Standard library only.
"""

import datetime
import hashlib
import os
import re
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
OMNI_URL = 'https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_{}.dat'
OMNI_DIR = os.path.join(ROOT, 'build', 'omni')
OMNI_FIRST_YEAR = 1963
SOLAR_WIND = os.path.join(ROOT, 'assets', 'geomag', 'solar_wind.bin')
GEOMAG_SOURCES = os.path.join(ROOT, 'assets', 'geomag', 'SOURCES.md')
EPHEM_CONFIG = os.path.join(ROOT, 'tools', 'bake', 'ephem.toml')
EPHEM_SOURCES = os.path.join(ROOT, 'assets', 'ephem', 'SOURCES.md')
HORIZONS_API = 'https://ssd.jpl.nasa.gov/api/horizons.api'

# The telescopes whose tables run into a forecast (assets/scenes/space_telescopes.toml).
TARGETS = ['jwst_geo', 'euclid_geo', 'roman_geo', 'tess_geo', 'chandra_geo']

MONTHS = {m: i + 1 for i, m in enumerate(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                                          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])}


# --- Text files -------------------------------------------------------------

def read_text(path):
    """Returns (text with LF line ends, the file's own line end)."""
    with open(path, encoding='utf-8', newline='') as f:
        text = f.read()
    return text.replace('\r\n', '\n'), '\r\n' if '\r\n' in text else '\n'


def write_text(path, text, newline):
    with open(path, 'w', encoding='utf-8', newline='') as f:
        f.write(text.replace('\n', newline))


def replace_once(text, pattern, repl, where):
    new, n = re.subn(pattern, repl, text)
    if n != 1:
        sys.exit(f'{where}: expected one match for {pattern!r}, found {n}; update the pattern or the page')
    return new


def run(*args):
    print('$', ' '.join(args), flush=True)
    subprocess.run([sys.executable, '-I', *args], cwd=ROOT, check=True)


def fetch(url, timeout=120):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code == 404 or attempt == 4:
                raise
        except Exception:
            if attempt == 4:
                raise
        time.sleep(10 * (attempt + 1))


# --- OMNI2 ------------------------------------------------------------------

def solar_wind_summary(path):
    with open(path, 'rb') as f:
        data = f.read()
    magic, (version, first_year, count) = data[:4], struct.unpack_from('<HHI', data, 4)
    if magic != b'AXSW' or version != 1:
        sys.exit(f'{path}: not an AXSW v1 file')
    hours = list(struct.iter_unpack('<HbB', data[12:12 + 4 * count]))
    start = datetime.datetime(first_year, 1, 1)

    def median(values):
        values = sorted(values)
        return values[len(values) // 2]
    pressures = [h[0] / 100 for h in hours if h[0] != 65535]
    speeds = [150 + 8 * h[2] for h in hours if h[2] != 255]
    return {
        'hours': count,
        'start': start,
        'end': start + datetime.timedelta(hours=count),
        'sha256': hashlib.sha256(data).hexdigest(),
        'median_pressure': median(pressures),
        'median_speed': median(speeds),
        'max_pressure': max(pressures),
    }


def refresh_omni(today):
    os.makedirs(OMNI_DIR, exist_ok=True)
    last_year = None
    for year in range(OMNI_FIRST_YEAR, today.year + 1):
        path = os.path.join(OMNI_DIR, f'omni2_{year}.dat')
        if os.path.exists(path) and year < today.year - 2:
            last_year = year
            continue
        try:
            data = fetch(OMNI_URL.format(year))
        except urllib.error.HTTPError as e:
            if e.code == 404 and year == today.year:  # early January: the new year's file is not out yet
                break
            raise
        with open(path + '.part', 'wb') as f:
            f.write(data)
        os.replace(path + '.part', path)
        print(f'downloaded omni2_{year}.dat ({len(data) // 1024} KB)')
        last_year = year

    old = solar_wind_summary(SOLAR_WIND)
    run(os.path.join('tools', 'geomag', 'solar_wind_bake.py'), OMNI_DIR, SOLAR_WIND)
    new = solar_wind_summary(SOLAR_WIND)

    text, nl = read_text(GEOMAG_SOURCES)
    where = 'assets/geomag/SOURCES.md'
    text = replace_once(text, r'`omni2_1963\.dat` 至 `omni2_\d{4}\.dat`', f'`omni2_1963.dat` 至 `omni2_{last_year}.dat`', where)
    text = replace_once(text, r'(OMNI 数据），)\d{4}-\d{2}-\d{2}( 下载)', rf'\g<1>{today:%Y-%m-%d}\g<2>', where)
    text = replace_once(text, r'seq 1963 \d{4}', f'seq 1963 {last_year}', where)
    text = replace_once(
        text,
        r'当前文件 [\d,]+ 小时，1963-01-01 至 \d{4}-\d{2}-\d{2} \d{2}:\d{2} UT，SHA-256 `[0-9a-f]{64}`；'
        r'流压中位数 [\d.]+ nPa、速度 \d+ km/s，最大的小时流压 [\d.]+ nPa',
        f'当前文件 {new["hours"]:,} 小时，1963-01-01 至 {new["end"]:%Y-%m-%d %H:%M} UT，SHA-256 `{new["sha256"]}`；'
        f'流压中位数 {new["median_pressure"]:.2f} nPa、速度 {new["median_speed"]} km/s，'
        f'最大的小时流压 {new["max_pressure"]:.2f} nPa',
        where)
    write_text(GEOMAG_SOURCES, text, nl)
    return old, new


# --- Horizons ---------------------------------------------------------------

def target_blocks(text):
    """{name: (start, end)} spans of the [[target]] blocks in ephem.toml."""
    starts = [m.start() for m in re.finditer(r'^\[\[target\]\]', text, re.M)] + [len(text)]
    blocks = {}
    for a, b in zip(starts, starts[1:]):
        m = re.search(r'^name = "([^"]+)"', text[a:b], re.M)
        if m:
            blocks[m.group(1)] = (a, b)
    return blocks


def block_value(text, span, key):
    m = re.search(rf'^{key} = "([^"]+)"', text[span[0]:span[1]], re.M)
    if not m:
        sys.exit(f'tools/bake/ephem.toml: no {key} in the block at offset {span[0]}')
    return m.group(1)


def horizons_coverage(command, center, start):
    """Asks for daily vectors to 2200 and reads where Horizons says its data ends.
    Returns (coverage end as a datetime in TDB, the trajectory notes of the object data)."""
    params = {
        'format': 'text', 'COMMAND': f"'{command}'", 'OBJ_DATA': "'YES'", 'MAKE_EPHEM': "'YES'",
        'EPHEM_TYPE': "'VECTORS'", 'CENTER': f"'{center}'", 'START_TIME': f"'{start}'",
        'STOP_TIME': "'2200-01-01'", 'STEP_SIZE': "'30 d'", 'TIME_TYPE': "'TDB'",
    }
    text = fetch(HORIZONS_API + '?' + urllib.parse.urlencode(params)).decode('utf-8', errors='replace')
    m = re.search(r'No ephemeris for target .* after A\.D\. (\d{4})-([A-Z]{3})-(\d{2}) (\d{2}):(\d{2}):(\d{2})', text)
    if not m:
        sys.exit(f'Horizons {command}: no end of coverage in the reply:\n{text[-1500:]}')
    y, mon, d, hh, mm, ss = m.groups()
    end = datetime.datetime(int(y), MONTHS[mon], int(d), int(hh), int(mm), int(ss))
    lines = text.splitlines()
    revised = next((ln.strip() for ln in lines if ln.strip().startswith('Revised:')), '')
    notes, keep = [revised] if revised else [], False
    for ln in lines:
        if 'TRAJECTORY' in ln.upper() and not keep:
            keep = True
        elif keep and ln.startswith('*****'):
            break
        if keep:
            notes.append(ln.rstrip())
    return end, '\n'.join(notes).strip()


def eph_knots(path):
    with open(path, 'rb') as f:
        head = f.read(20)
    return struct.unpack_from('<I', head, 16)[0]


def refresh_ephem():
    config, config_nl = read_text(EPHEM_CONFIG)
    blocks = target_blocks(config)
    missing = [t for t in TARGETS if t not in blocks]
    if missing:
        sys.exit(f'tools/bake/ephem.toml has no target {missing}')

    results = []
    for name in TARGETS:
        span = target_blocks(config)[name]
        command, center = block_value(config, span, 'command'), block_value(config, span, 'center')
        start, old_stop = block_value(config, span, 'start'), block_value(config, span, 'stop')
        end, notes = horizons_coverage(command, center, start)
        stop_time = (end - datetime.timedelta(minutes=1)).replace(minute=0, second=0)
        stop = f'{stop_time:%Y-%m-%d %H:%M}'
        print(f'{name}: Horizons {command} ends {end:%Y-%m-%d %H:%M:%S} TDB; stop {old_stop} -> {stop}')
        block = config[span[0]:span[1]]
        block = re.sub(r'^stop = "[^"]+"', f'stop = "{stop}"', block, count=1, flags=re.M)
        config = config[:span[0]] + block + config[span[1]:]
        path = os.path.join(ROOT, 'assets', 'ephem', name + '.eph')
        results.append({'name': name, 'command': command, 'start': start, 'old_stop': old_stop, 'stop': stop,
                        'end': end, 'notes': notes, 'old_knots': eph_knots(path), 'path': path})
    write_text(EPHEM_CONFIG, config, config_nl)

    run(os.path.join('tools', 'bake', 'horizons_bake.py'), EPHEM_CONFIG, *TARGETS)

    text, nl = read_text(EPHEM_SOURCES)
    for r in results:
        r['knots'] = eph_knots(r['path'])
        row = re.compile(rf'^\| `{re.escape(r["name"])}\.eph` \|.*$', re.M)
        found = row.findall(text)
        if len(found) != 1:
            sys.exit(f'assets/ephem/SOURCES.md: expected one table row for {r["name"]}.eph, found {len(found)}')
        cells = found[0].split('|')
        cells[4] = f' {r["start"]} .. {r["stop"]} '
        cells[-2] = f' {r["knots"]} '
        text = row.sub(lambda _: '|'.join(cells), text)
    write_text(EPHEM_SOURCES, text, nl)
    return results


# --- Pull request -----------------------------------------------------------

def pr_body(today, omni, ephem):
    out = [f'每月自动刷新（{today:%Y-%m-%d}，`tools/refresh/monthly_refresh.py`）。数据由脚本重新下载和烘焙，'
           '资料页里的数字已同步；文字说明需要人工核对后再合并。', '']
    if omni:
        old, new = omni
        out += ['## OMNI2 逐小时太阳风', '',
                '| | 小时数 | 记录末端（UT） | 流压中位数 | 速度中位数 | 最大流压 |', '|---|---:|---|---:|---:|---:|']
        for label, s in (('原有', old), ('刷新后', new)):
            out.append(f'| {label} | {s["hours"]:,} | {s["end"]:%Y-%m-%d %H:%M} | {s["median_pressure"]:.2f} nPa | '
                       f'{s["median_speed"]} km/s | {s["max_pressure"]:.2f} nPa |')
        out.append('')
    if ephem:
        out += ['## 空间望远镜星历', '',
                '| 文件 | 原结束 | 新结束 | Horizons 覆盖到（TDB） | 节点数 |', '|---|---|---|---|---|']
        for r in ephem:
            out.append(f'| `{r["name"]}.eph` | {r["old_stop"]} | {r["stop"]} | {r["end"]:%Y-%m-%d %H:%M:%S} | '
                       f'{r["old_knots"]} → {r["knots"]} |')
        out += ['', 'Horizons 的轨道说明（用于更新资料页里的跟踪截止日期和轨道文件名）：', '']
        for r in ephem:
            out += [f'<details><summary><code>{r["name"]}</code>（{r["command"]}）</summary>', '', '```',
                    r['notes'] or '(no trajectory notes)', '```', '', '</details>', '']
    out += ['## 合并前核对', '',
            '- [ ] Windows 和 Linux 的构建与测试通过（工作流已在本分支上触发 `build.yml`）。',
            '- [ ] `assets/ephem/SOURCES.md` 的说明文字：各望远镜的烘焙日期、跟踪数据截止日、预报来源。',
            '- [ ] `tools/bake/ephem.toml` 中望远镜一节的覆盖注释（“checked …”、轨道文件名）。',
            '- [ ] 结束时间变化较大时，打开 Space Telescopes 场景看最后一段轨迹，并核对测试里依赖预报段的跨度。',
            '- [ ] OMNI2 的年代覆盖率（`solar_wind_bake.py` 的输出）与 `assets/geomag/SOURCES.md` 一致。']
    return '\n'.join(out) + '\n'


def main():
    args = sys.argv[1:]
    body_path = None
    if '--body' in args:
        i = args.index('--body')
        if i + 1 >= len(args):
            sys.exit(__doc__)
        body_path = args[i + 1]
        del args[i:i + 2]
    unknown = [a for a in args if a not in ('--skip-omni', '--skip-ephem')]
    if unknown:
        sys.exit(__doc__)
    today = datetime.datetime.now(datetime.timezone.utc).date()
    omni = None if '--skip-omni' in args else refresh_omni(today)
    ephem = None if '--skip-ephem' in args else refresh_ephem()
    if body_path:
        with open(body_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(pr_body(today, omni, ephem))
        print(f'wrote {body_path}')


if __name__ == '__main__':
    main()
