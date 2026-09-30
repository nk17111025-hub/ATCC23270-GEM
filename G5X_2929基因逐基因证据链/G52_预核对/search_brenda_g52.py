"""Read BRENDA reference tables for every numeric G52 EC label.

This is a literature lead search. A same-species entry still requires
strain, sequence and measured reaction validation before gene-level credit.
"""
import concurrent.futures
import json
import re
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from lxml import html

DIR = Path(__file__).resolve().parent
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
ec = sorted({s.strip() for r in records for field in ('UniProt EC', 'KEGG live EC', 'D3 KEGG EC') for s in r[field].replace(',', ';').split(';') if re.fullmatch(r'\d+\.\d+\.\d+\.\d+', s.strip())})
path = DIR / 'g52_brenda_ec_references.json'
data = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {'by_ec': {}}

def get_one(ecno):
    url = 'https://www.brenda-enzymes.org/enzyme.php?ecno=' + ecno + '&onlyTable=Reference&showtm=0'
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'G52 evidence review/1.0'})
            with urllib.request.urlopen(req, timeout=40) as response:
                body = response.read()
            tree = html.fromstring(body)
            matches = []
            for tr in tree.iter('div'):
                if tr.get('class') != 'row' or not (tr.get('id') or '').startswith('tab30r'):
                    continue
                txt = ' '.join(' '.join(tr.itertext()).split())
                if 'Acidithiobacillus ferrooxidans' in txt or 'Ferrobacillus ferrooxidans' in txt:
                    matches.append(txt[:1500])
            return ecno, {'url': url, 'matching_reference_rows': matches, 'count': len(matches), 'html_bytes': len(body), 'parser_version': 4, 'boundary': 'BRENDA reference row is a species/EC literature lead, not a current-WP assay.'}
        except (urllib.error.URLError, OSError, ValueError) as exc:
            if attempt == 2:
                return ecno, {'url': url, 'error': str(exc)}
            time.sleep(attempt + 1)

missing = [s for s in ec if s not in data['by_ec'] or 'error' in data['by_ec'][s] or data['by_ec'][s].get('parser_version') != 4]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    futures = [pool.submit(get_one, e) for e in missing]
    for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
        e, result = future.result()
        data['by_ec'][e] = result
        if i % 10 == 0 or i == len(missing):
            data['accessed_utc'] = datetime.now(timezone.utc).isoformat()
            data['source'] = 'https://www.brenda-enzymes.org/'
            path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
            print(f'checked {i}/{len(missing)}; species-reference ECs {sum(bool(v.get("count")) for v in data["by_ec"].values())}', flush=True)

assert len(data['by_ec']) == len(ec)
print(json.dumps({'ec': len(ec), 'species_reference_ec': sum(bool(v.get('count')) for v in data['by_ec'].values()), 'errors': sum('error' in v for v in data['by_ec'].values())}))
