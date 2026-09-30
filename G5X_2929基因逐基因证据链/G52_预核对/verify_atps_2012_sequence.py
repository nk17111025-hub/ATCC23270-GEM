"""Compare the 2012 ATCC 23270 assayed atpS sequence with current WP."""
import json
import urllib.request
from io import StringIO
from pathlib import Path

from Bio import SeqIO

DIR = Path(__file__).resolve().parent
ROOT = Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1')
fasta = DIR / 'FM177944.1_atpS_ATCC23270.fasta'
if not fasta.exists():
    with urllib.request.urlopen('https://www.ebi.ac.uk/ena/browser/api/fasta/FM177944?download=false', timeout=35) as r:
        fasta.write_bytes(r.read())
old = next(SeqIO.parse(fasta, 'fasta'))
current = next(x for x in SeqIO.parse(ROOT / 'cds_from_genomic.fna', 'fasta') if 'RU820_RS02615' in x.description)
a, b = str(old.seq), str(current.seq)
assert len(a) == len(b) == 1674
nt_differences = [{'position_1_based': i+1, 'FM177944': x, 'current': y} for i, (x, y) in enumerate(zip(a,b)) if x != y]
old_aa, new_aa = str(old.seq.translate(to_stop=True)), str(current.seq.translate(to_stop=True))
aa_differences = [{'position_1_based': i+1, 'FM177944': x, 'current': y} for i, (x, y) in enumerate(zip(old_aa,new_aa)) if x != y]
assert len(nt_differences)==7 and len(aa_differences)==4
same_length = []
for p in SeqIO.parse(ROOT / 'protein.faa', 'fasta'):
    seq = str(p.seq).replace('*','')
    if len(seq)==len(old_aa):
        same_length.append({'wp': p.id, 'mismatches': sum(x!=y for x,y in zip(seq,old_aa))})
same_length.sort(key=lambda x:x['mismatches'])
assert same_length[0]['wp']=='WP_012536192.1' and same_length[0]['mismatches']==4
out = {
    'paper_doi': '10.6026/97320630008695',
    'paper_assayed_gene_accession': 'FM177944.1',
    'current_locus': 'RU820_RS02615',
    'current_wp': 'WP_012536192.1',
    'legacy_afe': 'AFE_0539',
    'old_and_current_cds_length_nt': 1674,
    'old_and_current_protein_length_aa': len(old_aa),
    'nt_differences': nt_differences,
    'aa_differences': aa_differences,
    'protein_identity': round(1-len(aa_differences)/len(old_aa),6),
    'same_length_current_proteins_ranked': same_length[:5],
    'experimental_result': 'Crude E. coli extracts expressing cloned ATCC 23270 atpS formed ATP from APS and pyrophosphate; induced extract versus controls.',
    'evidence_boundary': 'Direct assay is of a near-identical, same-strain historical allele, not the exact current WP sequence. APS kinase domain activity and SFAT2 were not assayed.',
    'sequence_source': 'https://www.ebi.ac.uk/ena/browser/view/FM177944',
}
(DIR / 'g52_atps_2012_sequence_comparison.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'nt_differences':len(nt_differences),'aa_differences':len(aa_differences),'protein_identity':out['protein_identity'],'nearest':same_length[:3]},ensure_ascii=False))
