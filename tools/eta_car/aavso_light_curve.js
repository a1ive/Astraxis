// Eta Carinae's light curve from the AAVSO International Database, reduced to
// assets/stars/eta_car_aavso.csv (see assets/stars/SOURCES.md).
//
// The AAVSO search pages sit behind a browser challenge, so this runs in a
// browser: open https://apps.aavso.org/v2/data/search/photometry/ (no login
// needed), paste this into the developer console and save what it prints.
//
// It fetches every visual and Johnson V observation of "eta Car" (1000 per
// page), drops duplicates (the search lists some twice) and limits ("<", ">"),
// and bins them in kBinYears: the median of the visual estimates before
// kCutoverYear and of the instrumental V measures (CCD, PEP, DSLR) from then
// on, each falling back on the other kind in a bin without it. Visual
// estimates of the nebulous star run ~0.3 mag fainter than V (1970s-1990s);
// taking whichever kind a bin has would zigzag between them. Around 2004 CCD
// monitoring takes over. The historical estimates in the database (1595 on)
// are Frew 2004's, reduced to the V scale, with Smith & Frew 2011's additions
// for 1838-1845.

const kBinYears = 0.25;
const kCutoverYear = 2004.0;
const kBase = 'https://apps.aavso.org/v2/data/search/photometry/?target=eta+Car&band=2&band=0&submit=Search' +
              '&per_page=1000&page=';

async function fetchAll() {
    const rows = [];
    for (let page = 1;; ++page) {
        const response = await fetch(kBase + page);
        if (!response.ok) {
            break;
        }
        const doc = new DOMParser().parseFromString(await response.text(), 'text/html');
        const trs = doc.querySelectorAll('table tbody tr');
        if (trs.length === 0) {
            break;
        }
        // #, target, observer, JD, UTC date, magnitude, uncertainty, band, type
        trs.forEach(tr => rows.push([...tr.querySelectorAll('td')].map(td => td.innerText.trim())));
        if (trs.length < 1000) {
            break;
        }
        await new Promise(resolve => setTimeout(resolve, 500)); // gently
    }
    return rows;
}

function median(values) {
    const v = values.slice().sort((a, b) => a - b);
    const n = v.length;
    return n % 2 ? v[(n - 1) / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

async function reduce() {
    const rows = await fetchAll();
    const seen = new Set();
    const bins = new Map();
    let kept = 0;
    for (const [, , observer, jdText, , magText, , band, type] of rows) {
        const key = observer + ' ' + jdText + ' ' + magText + ' ' + band;
        if (seen.has(key) || !/^-?\d+(\.\d+)?$/.test(magText)) {
            continue; // a duplicate, or a limit
        }
        seen.add(key);
        const jd = parseFloat(jdText);
        const year = 2000.0 + (jd - 2451545.0) / 365.25;
        const instrumental = band === 'Johnson V' && type !== 'Visual';
        if (band !== 'Visual' && !instrumental) {
            continue;
        }
        const k = Math.floor(year / kBinYears);
        if (!bins.has(k)) {
            bins.set(k, {v: [], visual: [], years: []});
        }
        const b = bins.get(k);
        (instrumental ? b.v : b.visual).push(parseFloat(magText));
        b.years.push({year, instrumental});
        ++kept;
    }
    const lines = ['# Eta Carinae, AAVSO International Database (visual and Johnson V), retrieved ' +
                       new Date().toISOString().slice(0, 10) + ' by tools/eta_car/aavso_light_curve.js.',
                   '# ' + kept + ' observations in ' + kBinYears + '-year bins: year (mean of those used), ' +
                       'median V, count, source (V: instrumental, vis: visual).',
                   'year,v,count,source'];
    for (const k of [...bins.keys()].sort((a, b) => a - b)) {
        const b = bins.get(k);
        const prefer = (k + 0.5) * kBinYears >= kCutoverYear;
        const instrumental = prefer ? b.v.length > 0 : b.visual.length === 0;
        const mags = instrumental ? b.v : b.visual;
        const years = b.years.filter(y => y.instrumental === instrumental).map(y => y.year);
        const mean = years.reduce((s, y) => s + y, 0) / years.length;
        lines.push(mean.toFixed(3) + ',' + median(mags).toFixed(2) + ',' + mags.length + ',' +
                   (instrumental ? 'V' : 'vis'));
    }
    return lines.join('\n');
}

reduce().then(text => console.log(text));
