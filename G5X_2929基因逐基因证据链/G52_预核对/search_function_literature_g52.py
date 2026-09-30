"""Discovery-only Europe PMC search for G52 current product names.

Every current G52 locus is queried; a negative index result is not evidence
that no original experiment exists, and positive hits require paper review.
"""
import concurrent.futures
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
path = DIR / 'g52_function_literature_search.json'
data = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {'by_locus': {}}

def normalize(name):
    name = urllib.parse.unquote(name or '')
    name = re.sub(r'\s+', ' ', name).strip()
    return name

def search(record):
    loc, name = record['当前locus'], normalize(record['当前product'])
    query = f'"Acidithiobacillus ferrooxidans" AND "{name}"'
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': query, 'format': 'json', 'pageSize': 20})
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'G52-evidence-review/1.0'})
            with urllib.request.urlopen(req, timeout=35) as response:
                body = json.load(response)
            hits = [{'title': x.get('title', ''), 'doi': x.get('doi', ''), 'pmid': x.get('pmid', ''), 'year': x.get('pubYear', '')} for x in body.get('resultList', {}).get('result', [])]
            return loc, {'product_query': name, 'hit_count': body.get('hitCount', 0), 'hits': hits, 'url': url, 'boundary': 'Index hit only; exact WP identity, strain and experimental type unverified.'}
        except (urllib.error.URLError, OSError, ValueError) as exc:
            if attempt == 2:
                return loc, {'product_query': name, 'error': str(exc), 'url': url}
            time.sleep(attempt + 1)

missing = [r for r in records if r['当前locus'] not in data['by_locus'] or 'error' in data['by_locus'][r['当前locus']]]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for loc, result in pool.map(search, missing):
        data['by_locus'][loc] = result

data['accessed_utc'] = datetime.now(timezone.utc).isoformat()
data['source'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'loci': len(data['by_locus']), 'searched_now': len(missing), 'hit_loci': sum(bool(v.get('hit_count')) for v in data['by_locus'].values()), 'errors': sum('error' in v for v in data['by_locus'].values())}, ensure_ascii=False))
