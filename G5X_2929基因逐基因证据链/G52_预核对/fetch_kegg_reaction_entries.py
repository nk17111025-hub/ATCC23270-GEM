import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
data = json.loads((DIR / 'g52_kegg.json').read_text(encoding='utf-8'))
ids = sorted({rid.removeprefix('rn:') for links in data['ko_reactions'].values() for rid in links})
out = {'accessed_utc': datetime.now(timezone.utc).isoformat(), 'source': 'https://rest.kegg.jp/get/', 'reactions': {}, 'failed': []}

def parse_entry(text):
    fields = {}
    key = ''
    for line in text.splitlines():
        if line == '///':
            break
        if line[:12].strip():
            key = line[:12].strip()
            fields.setdefault(key, []).append(line[12:].strip())
        elif key:
            fields[key].append(line[12:].strip())
    return {k: '; '.join(v) for k, v in fields.items() if k in {'ENTRY', 'NAME', 'DEFINITION', 'EQUATION', 'COMMENT', 'ENZYME', 'ORTHOLOGY', 'DBLINKS'}}

for start in range(0, len(ids), 10):
    chunk = ids[start:start + 10]
    url = 'https://rest.kegg.jp/get/' + '+'.join(chunk)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'G52-evidence-review/1.0'})
            with urllib.request.urlopen(req, timeout=45) as response:
                body = response.read().decode('utf-8')
            for entry in body.split('///'):
                parsed = parse_entry(entry)
                rid = (parsed.get('ENTRY', '').split() or [''])[0]
                if rid in chunk:
                    out['reactions'][rid] = parsed
            break
        except (urllib.error.URLError, OSError) as exc:
            if attempt == 2:
                out['failed'].append({'ids': chunk, 'error': str(exc)})
            else:
                time.sleep(attempt + 1)

(DIR / 'g52_kegg_reaction_entries.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'requested': len(ids), 'received': len(out['reactions']), 'failed_batches': len(out['failed'])}))
