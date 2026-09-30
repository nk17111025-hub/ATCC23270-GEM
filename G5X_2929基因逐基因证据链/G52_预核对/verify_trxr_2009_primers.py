"""Check published ATCC 23270 TrxR primer ends against the current CDS.

Primer sequences: Wang et al. 2009, doi:10.1007/s00284-009-9390-2,
Materials and Methods, pp. 36–37. Only the gene-specific 3' tails are used.
"""
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

DIR = Path(__file__).resolve().parent
CDS = Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/cds_from_genomic.fna')
target = next(r for r in SeqIO.parse(CDS, 'fasta') if 'RU820_RS01820' in r.description)
seq = str(target.seq).upper()
forward_gene_tail = 'GATCCGAAAAGTGTCATCGATGAAAAAC'
reverse_gene_tail = 'GCCCTCTTGTTGCTGCTCCAGCCAACGCTCC'
reverse_rc = str(Seq(reverse_gene_tail).reverse_complement())
fpos, rpos = seq.find(forward_gene_tail), seq.find(reverse_rc)
assert fpos == 6 and rpos == len(seq) - len(reverse_rc) - 3
result = {
    'locus': 'RU820_RS01820',
    'wp': 'WP_009566072.1',
    'legacy_afe': 'AFE_0375',
    'paper_doi': '10.1007/s00284-009-9390-2',
    'paper_strain': 'ATCC 23270',
    'current_cds_length_nt': len(seq),
    'forward_gene_specific_tail': forward_gene_tail,
    'forward_match_0_based': fpos,
    'forward_mismatches': 0,
    'reverse_gene_specific_tail': reverse_gene_tail,
    'reverse_complement_match_0_based': rpos,
    'reverse_mismatches': 0,
    'both_primer_ends_match_current_cds': True,
    'paper_experiment': 'Recombinant TrxR directly reduced oxidized thioredoxin with NADPH; Cys142/Cys145 mutants lost most activity.',
    'evidence_boundary': 'Exact primer ends link the assayed ATCC 23270 gene to the current CDS; published assay does not repair the 2016 TRDR charge imbalance.',
}
(DIR / 'g52_trxr_2009_exact_identity.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False))
