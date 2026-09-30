import concurrent.futures
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
afe = sorted({r['旧AFE'] for r in records if r['旧AFE']})
path = DIR / 'g52_literature_search.json'
existing = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {'by_afe': {}}

def search(a):
    params = urllib.parse.urlencode({'query': '"' + a + '"', 'format': 'json', 'pageSize': 50})
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + params
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'G52-evidence-review/1.0 contact=research'})
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.load(response)
            hits = [{'title': x.get('title', ''), 'doi': x.get('doi', ''), 'pmid': x.get('pmid', ''), 'year': x.get('pubYear', ''), 'source': x.get('source', '')} for x in data.get('resultList', {}).get('result', [])]
            return a, {'hitCount': data.get('hitCount', 0), 'hits': hits, 'url': url}
        except (urllib.error.URLError, OSError, ValueError) as exc:
            if attempt == 2:
                return a, {'error': str(exc), 'url': url}
            time.sleep(attempt + 1)

missing = [a for a in afe if a not in existing['by_afe'] or 'error' in existing['by_afe'][a]]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for a, result in pool.map(search, missing):
        existing['by_afe'][a] = result

existing['accessed_utc'] = datetime.now(timezone.utc).isoformat()
existing['source'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
path.write_text(json.dumps(existing, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'afe': len(afe), 'searched': len(missing), 'hit_afe': sum(bool(v.get('hitCount')) for v in existing['by_afe'].values()), 'errors': sum('error' in v for v in existing['by_afe'].values())}))
