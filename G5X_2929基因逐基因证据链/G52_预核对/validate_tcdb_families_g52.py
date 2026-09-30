"""Compare selected current G52 WPs to TCDB within their annotated families.

Family selection comes from the current RefSeq product, not the all-TCDB top hit.
These are homology leads only: neither substrate nor transport direction is inferred.
"""
import json
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

DIR = Path(__file__).resolve().parent
ROOT = Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1')
records = json.loads((DIR / 'g52_precheck.json').read_text(encoding='utf-8'))['records']
by_locus = {r['当前locus']: r for r in records}

# The expected families are independently anchored to the current product name.
expected = {
    'RU820_RS01515': ('3.A.5.', 'SecE'),
    'RU820_RS01680': ('3.A.5.', 'SecY'),
    'RU820_RS01835': ('1.B.9.', 'FadL'),
    'RU820_RS02250': ('2.A.4.', 'CDF'),
    'RU820_RS02330': ('3.A.1.125.', 'LolC'),
    'RU820_RS02335': ('3.A.1.125.', 'LolD'),
    'RU820_RS02350': ('1.A.30.', 'ExbB/TolQ'),
    'RU820_RS02640': ('5.A.1.', 'DsbD'),
    'RU820_RS02705': ('2.A.66.', 'MurJ'),
    'RU820_RS02895': ('1.C.126.', 'CorC'),
}

wp_to_locus = {by_locus[loc]['WP']: loc for loc in expected}
queries = {}
for p in SeqIO.parse(ROOT / 'protein.faa', 'fasta'):
    if p.id in wp_to_locus:
        queries[wp_to_locus[p.id]] = str(p.seq).replace('*', '').upper()
assert set(queries) == set(expected)

refs = {}
for p in SeqIO.parse(DIR / 'tcdb_2026_09_25.faa', 'fasta'):
    tcid = p.id.rsplit('|', 1)[-1]
    seq = str(p.seq).replace('*', '').upper()
    if len(seq) >= 70:
        refs.setdefault(tcid, (p.description, seq))

aligner = PairwiseAligner(mode='local')
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def metrics(q, t):
    aln = aligner.align(q, t)[0]
    aligned = matches = 0
    for (q0, q1), (t0, t1) in zip(*aln.aligned):
        q_piece, t_piece = q[q0:q1], t[t0:t1]
        aligned += min(len(q_piece), len(t_piece))
        matches += sum(a == b for a, b in zip(q_piece, t_piece))
    return {
        'identity': round(matches / aligned, 4) if aligned else 0,
        'query_coverage': round(aligned / len(q), 4),
        'target_coverage': round(aligned / len(t), 4),
        'aligned_residues': aligned,
        'score': round(aln.score, 1),
    }

out = {}
for loc, (prefix, family) in expected.items():
    hits = []
    for tcid, (description, seq) in refs.items():
        if tcid.startswith(prefix):
            hits.append({'tcid': tcid, 'description': description, **metrics(queries[loc], seq)})
    hits.sort(key=lambda x: (x['identity'] * x['query_coverage'], x['score']), reverse=True)
    best = hits[0] if hits else None
    out[loc] = {
        'wp': by_locus[loc]['WP'],
        'current_product': by_locus[loc]['当前product'],
        'expected_family': family,
        'expected_tcid_prefix': prefix,
        'family_reference_count': len(hits),
        'best': best,
        'family_homology_supported': bool(best and best['identity'] >= .30 and best['query_coverage'] >= .60 and not (family == 'SecE' and 'secy' in best['description'].lower())),
        'evidence_boundary': 'Sequence-family support only; transport substrate, direction and compartment unresolved.',
    }

result = {
    'source': 'https://www.tcdb.org/public/tcdb',
    'retrieved_utc': datetime.now(timezone.utc).isoformat(),
    'selection_rule': 'Current RefSeq product family, followed by within-family local alignment',
    'by_locus': out,
}
(DIR / 'g52_tcdb_family_validation.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
for loc, v in out.items():
    b = v['best']
    print(loc, v['expected_family'], b['tcid'] if b else '-', b['identity'] if b else '-', b['query_coverage'] if b else '-', v['family_homology_supported'])
