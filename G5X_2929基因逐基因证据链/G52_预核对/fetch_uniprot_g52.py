import csv
import json
import re
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

DIR = Path(__file__).resolve().parent
with (DIR / 'g52_precheck.json').open(encoding='utf-8') as f:
    records = json.load(f)['records']
targets = {r['WP'] for r in records}
fields = 'accession,id,reviewed,protein_name,ec,xref_refseq,xref_interpro,xref_pfam,ft_act_site,ft_binding'
url = f'https://rest.uniprot.org/uniprotkb/search?query={quote("organism_id:243159")}&fields={fields}&format=tsv&size=500'
all_rows = []
pages = 0
while url:
    req = urllib.request.Request(url, headers={'User-Agent':'G52-evidence-review/1.0','Accept':'text/tab-separated-values'})
    with urllib.request.urlopen(req, timeout=60) as response:
        body = response.read().decode('utf-8-sig')
        header = response.headers.get('Link', '')
    rows = list(csv.DictReader(body.splitlines(), delimiter='\t'))
    all_rows.extend(rows)
    pages += 1
    link = re.search(r'<([^>]+)>;\s*rel="next"', header)
    url = link.group(1) if link else ''
    print('page', pages, 'rows', len(rows), 'total', len(all_rows), flush=True)
    if url:
        time.sleep(0.4)

by_wp = defaultdict(list)
for row in all_rows:
    for wp in re.findall(r'WP_\d+\.\d+', row.get('RefSeq','')):
        if wp in targets:
            by_wp[wp].append(row)

result = {}
for wp in sorted(targets):
    entries = by_wp.get(wp, [])
    result[wp] = [{
        'accession': r['Entry'],
        'entry_name': r['Entry Name'],
        'reviewed': r['Reviewed'],
        'protein_names': r['Protein names'],
        'ec': r['EC number'],
        'interpro': r['InterPro'],
        'pfam': r['Pfam'],
        'active_site': r['Active site'],
        'binding_site': r['Binding site'],
        'url': f"https://www.uniprot.org/uniprotkb/{r['Entry']}/entry",
    } for r in entries]

out = {'accessed_utc':datetime.now(timezone.utc).isoformat(),'query':'organism_id:243159','source':'https://rest.uniprot.org/uniprotkb/search','pages':pages,'total_taxon_entries':len(all_rows),'by_wp':result}
with (DIR/'g52_uniprot.json').open('w',encoding='utf-8') as f:
    json.dump(out,f,ensure_ascii=False,indent=2)
print(json.dumps({'targets':len(targets),'matched_wp':sum(bool(v) for v in result.values()),'unmatched_wp':sum(not v for v in result.values()),'reviewed_wp':sum(any(e['reviewed']=='reviewed' for e in v) for v in result.values())}))
