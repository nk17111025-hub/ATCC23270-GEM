import csv
import pathlib
import re

root = pathlib.Path(__file__).parent
papers = list(csv.DictReader((root / 'G54_论文来源表.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
index = list(csv.DictReader((root / 'G54_基因索引_初核.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
by_locus = {r['当前 RU820 locus']: f"G54-{int(r['总序号']):04d}" for r in index}
fields = ['Source ID', 'Full citation', 'Year', 'DOI', 'Source type', 'Strain',
          'Experimental type', 'Genes / proteins studied',
          'Related Reaction IDs / candidate IDs', 'What was directly measured',
          'Applicable confidence level', 'Evidence limitation']
sources = [
    ('G54-R001', 'NCBI RefSeq GCF_049532655.1 / NZ_CP136162.1; local frozen GFF/CDS/protein', '2026', 'official genome', 'ATCC 23270', 'current locus/WP/product/coordinates', 'https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/', 'Annotation and sequence identity; no enzyme assay'),
    ('G54-R002', 'B1 exact RefSeq protein_id old/new locus mapping; local project table', '2026', 'derived mapping', 'ATCC 23270', 'old AFE tracing only', '', '233 exact, 56 unmapped, 4 multicopy; old AFE is not primary identity'),
    ('G54-R003', '2016 iMC507 original paper and mmc1.xls Table 1 (615 reaction rows)', '2016', 'model baseline', 'ATCC 23270', 'old reaction/GPR/formula/confidence/PMID/media bounds', 'https://doi.org/10.1016/j.meteno.2016.03.003', 'Model assignments and 1.8 proton coefficient are not direct single-enzyme measurements'),
    ('G54-R004', 'UniProtKB RefSeq_Protein ID mapping; InterPro/Pfam cross-references', '2026', 'protein annotation', 'strain tagged when available', 'WP protein and domain leads', 'https://rest.uniprot.org/idmapping/run', 'Unreviewed TrEMBL and homolog annotations remain inferred'),
    ('G54-R005', 'Europe PMC REST search for 293 RU820/WP/AFE identities', '2026', 'literature search', 'mixed', 'original paper discovery', 'https://www.ebi.ac.uk/europepmc/webservices/rest/', 'Abstract API hits are leads; strain and experimental fact require full-text review'),
    ('G54-R006', 'KEGG/Rhea primary reaction entries, with BioCyc links where accessible', '2026', 'reaction reference', 'cross species', '29 candidate reaction-ID loci', 'https://www.kegg.jp/;https://www.rhea-db.org/', 'Database cross-references share annotation inheritance and do not prove target-strain catalysis'),
    ('G54-R007', 'Rhea RHEA:26333 molybdopterin synthase', '2026', 'reaction reference', 'cross species', 'RU820_RS04675;RU820_RS04680', 'https://www.rhea-db.org/rhea/26333', 'Full equation has protein-bound sulfur carrier; ATCC 23270 enzyme assay absent'),
    ('G54-R008', 'Rhea RHEA:27758 glycine cleavage system', '2026', 'reaction reference', 'cross species', 'RU820_RS04695–RU820_RS04710;GHMT3 comparison', 'https://www.rhea-db.org/rhea/27758', 'Existing GHMT3 chemical matches but old GPR/EC differs; GcvL identity unresolved'),
    ('G54-R009', 'Rhea RHEA:15341 cobaltochelatase', '2026', 'reaction reference', 'Pseudomonas denitrificans source experiment', 'RU820_RS06225 candidate', 'https://www.rhea-db.org/rhea/15341', 'CobN/S/T complex; CobS candidate alone cannot carry whole reaction'),
    ('G54-R010', 'PDB 3KPK SQR chain compared with current/old RefSeq protein.faa', '2026', 'sequence reconciliation', 'ATCC 23270', 'RU820_RS06210 conflict; RU820_RS08260 match', 'https://www.rcsb.org/structure/3KPK', '3KPK C160A mutant matches WP_012536761.1 433/434; not WP_201763903.1'),
]
rows = []
for sid, title, year, stype, strain, genes, url, boundary in sources:
    rows.append({
        'Source ID': sid, 'Full citation': title + (f' [{url}]' if url else ''), 'Year': year,
        'DOI': '10.1016/j.meteno.2016.03.003' if sid == 'G54-R003' else '',
        'Source type': stype, 'Strain': strain, 'Experimental type': 'database or model audit',
        'Genes / proteins studied': genes, 'Related Reaction IDs / candidate IDs': genes,
        'What was directly measured': 'No target-strain single-enzyme assay' if stype != 'sequence reconciliation' else 'Sequence identity comparison',
        'Applicable confidence level': 'not scored alone', 'Evidence limitation': boundary,
    })
for p in papers:
    related = [by_locus[x] for x in dict.fromkeys(re.findall(r'RU820_RS\d+', p['Related RU820 loci'])) if x in by_locus]
    if p['Source ID'] in {'G54-S004', 'G54-S008'}:
        related.append('G54-1169 (historical-ID conflict only); RU820_RS08260 (out of group)')
    rows.append({
        'Source ID': p['Source ID'],
        'Full citation': p['Full citation'] + ' [' + p['Article URL'] + ']',
        'Year': p['Year'], 'DOI': p['DOI'], 'Source type': p['Source type'],
        'Strain': p['Strain'], 'Experimental type': p['Experimental type'],
        'Genes / proteins studied': p['Genes / proteins studied'],
        'Related Reaction IDs / candidate IDs': ';'.join(related),
        'What was directly measured': p['Directly measured/reported facts'],
        'Applicable confidence level': 'supports expression/context only' if p['Source ID'] in {'G54-S005','G54-S006','G54-S007'} else 'not scored alone',
        'Evidence limitation': p['Evidence limitation'],
    })
assert len({r['Source ID'] for r in rows}) == len(rows)
with (root / 'G54_证据来源表.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print(len(rows), 'sources')
