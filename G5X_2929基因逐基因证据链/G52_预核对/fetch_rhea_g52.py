"""Fetch current official Rhea equations for all G52 cross-reference IDs."""
import csv
import io
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
ids = {s for r in records for field in ('UniProt reaction Rhea', 'D3 Rhea交叉引用') for s in r[field].split(';') if s.startswith('RHEA:')}
params = {'query': '', 'columns': 'rhea-id,equation,chebi-id,ec,pubmed,reaction-xref(KEGG),reaction-xref(MetaCyc)', 'format': 'tsv', 'limit': '100000'}
url = 'https://www.rhea-db.org/rhea/?' + urllib.parse.urlencode(params)
request = urllib.request.Request(url, headers={'User-Agent': 'G52-evidence-review/1.0'})
with urllib.request.urlopen(request, timeout=120) as response:
    body = response.read().decode('utf-8-sig')
reader = csv.DictReader(io.StringIO(body), delimiter='\t')
entries = {}
total = 0
for row in reader:
    total += 1
    key = row.get('Reaction identifier', '')
    if key in ids:
        entries[key] = row
out = {'source': url, 'accessed_utc': datetime.now(timezone.utc).isoformat(), 'requested_rhea_ids': len(ids), 'server_rows': total, 'matched': len(entries), 'missing': sorted(ids - set(entries)), 'reactions': entries, 'boundary': 'Rhea is curated chemical reaction definition; gene-specific catalysis requires independent evidence.'}
(DIR / 'g52_rhea_reaction_entries.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'requested': len(ids), 'matched': len(entries), 'server_rows': total, 'missing': len(ids - set(entries))}))
