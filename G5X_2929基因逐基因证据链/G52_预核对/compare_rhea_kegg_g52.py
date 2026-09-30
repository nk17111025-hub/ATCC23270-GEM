"""Cross-check G52 Rhea and KEGG reaction xrefs without equating xref with assay."""
import json
from pathlib import Path

DIR = Path(__file__).resolve().parent
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
rhea = json.loads((DIR / 'g52_rhea_reaction_entries.json').read_text(encoding='utf-8'))['reactions']
out = []
for r in records:
    ids = sorted({s for f in ('UniProt reaction Rhea', 'D3 Rhea交叉引用') for s in r[f].split(';') if s.startswith('RHEA:')})
    ke = {s for s in r['KEGG live reaction'].split(';') if s}
    mapped = set()
    for rid in ids:
        mapped |= {s.removeprefix('KEGG:') for s in rhea[rid]['Cross-reference (KEGG)'].split(',') if s.startswith('KEGG:')}
    overlap = sorted(ke & mapped)
    if not ids or not ke:
        status = 'one-side-or-none'
    elif overlap:
        status = 'xref-overlap'
    else:
        status = 'different-granularity-or-chemistry-to-review'
    out.append({'locus': r['当前locus'], 'wp': r['WP'], 'rhea_ids': ids, 'kegg_ids': sorted(ke), 'rhea_kegg_xrefs': sorted(mapped), 'overlap': overlap, 'status': status, 'boundary': 'Shared xref is one annotation chain, not independent experimental confirmation.'})
summary = {'both': sum(bool(x['rhea_ids'] and x['kegg_ids']) for x in out), 'overlap': sum(x['status'] == 'xref-overlap' for x in out), 'different': sum(x['status'] == 'different-granularity-or-chemistry-to-review' for x in out)}
assert summary == {'both': 74, 'overlap': 70, 'different': 4}
(DIR / 'g52_rhea_kegg_xref_comparison.json').write_text(json.dumps({'summary': summary, 'rows': out}, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary))
