import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).parent

def read(name):
    with (ROOT / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

index = read('G54_基因索引_初核.tsv')
ids = {r['总序号']: (r['当前 RU820 locus'], r['当前 WP protein_id']) for r in index}
assert len(index) == len(ids) == 293
assert set(ids) == {str(i) for i in range(880, 1173)}

specs = [
    ('G54_293数据库核查.tsv', '总序号', '当前 RU820 locus', '当前 WP protein_id'),
    ('G54_293论文检索.tsv', '总序号', 'RU820', 'WP'),
    ('G54_2016逐基因对照.tsv', '总序号', 'RU820 locus', 'WP protein_id'),
]
for name, no, locus, wp in specs:
    rows = read(name)
    assert len(rows) == 293, (name, len(rows))
    got = {r[no]: (r[locus], r[wp]) for r in rows}
    assert got == ids, name
    print(f'{name}: 293 identities matched')

reactions = read('G54_反应候选_独立原始记录裁决.tsv')
assert len(reactions) == 29
assert len({r['关联行ID'] for r in reactions}) == 29
for row in reactions:
    number = row['关联行ID'].split('-')[-1]
    assert ids[str(int(number))][0] == row['RU820_locus']
print('Independent reaction candidates: 29 identities matched')

model = read('G54_2016逐基因对照.tsv')
direct = [r for r in model if r['Direct Reaction ID']]
assert len(direct) == 2
assert {(r['RU820 locus'], r['Direct Reaction ID']) for r in direct} == {
    ('RU820_RS04575', 'CYTBD'), ('RU820_RS04580', 'CYTBD')}
print('Direct 2016 GPR: 2 CYTBD genes')

alpha = read('G54_Alpha逐基因审查.tsv')
assert len(alpha) == 293
assert 'Fe²⁺ lb' in alpha[0] and 'Fe²⁺ ub' in alpha[0] and 'Fe2+ lb' not in alpha[0]
assert {r['Candidate ID'] for r in alpha} == {f'G54-{n:04d}' for n in range(880, 1173)}
assert all(r['Current locus'] == ids[str(int(r['Candidate ID'].split('-')[-1]))][0] for r in alpha)
assert all(r['2026 reaction score'] in {'不评分', '1', '2', '3', '4'} for r in alpha)
assert all(r['2026 reaction score'] or r['Score scope'] for r in alpha)
assert next(r for r in alpha if r['Candidate ID'] == 'G54-1169')['2026 reaction score'] == '不评分'
for n in (1028, 1058, 1092, 1143):
    a = next(r for r in alpha if r['Candidate ID'] == f'G54-{n:04d}')
    assert 'AFE_' in a['Legacy AFE locus'] and ';' in a['Legacy AFE locus']
sources = {r['Source ID'] for r in read('G54_证据来源表.tsv')}
assert list(read('G54_证据来源表.tsv')[0]) == [
    'Source ID', 'Full citation', 'Year', 'DOI', 'Source type', 'Strain',
    'Experimental type', 'Genes / proteins studied',
    'Related Reaction IDs / candidate IDs', 'What was directly measured',
    'Applicable confidence level', 'Evidence limitation']
assert all(set(r['Source ID'].split(';')) <= sources for r in alpha)
for n, sid in ((887, 'G54-S001'), (989, 'G54-S002'), (1061, 'G54-S003'),
               (923, 'G54-S009'), (1016, 'G54-S010')):
    row = next(r for r in alpha if r['Candidate ID'] == f'G54-{n:04d}')
    assert sid in row['Source ID'] and row['Primary evidence DOI']
for r in direct:
    a = next(a for a in alpha if a['Current locus'] == r['RU820 locus'])
    assert a['Reaction ID'] == r['Direct Reaction ID'] == 'CYTBD'
    assert a['Reaction Formula'] == r['Reaction formula']
    assert a['Gene-Reaction Association'] == r['GRA']
    for key in ('Fe2+ lb', 'Fe2+ ub', 'tetrathionate lb', 'tetrathionate ub', 'sulfur lb', 'sulfur ub'):
        alpha_key = key.replace('Fe2+ ', 'Fe²⁺ ')
        assert a[alpha_key] == r[key], (r['RU820 locus'], key)
print('Alpha: 293 distinct genes, source IDs resolved, direct 2016 A-P aligned')

comparison = read('G54_Alpha既有反应逐项对照.tsv')
assert len(comparison) == 13
assert all(r['Reaction ID'] and r['Reaction Formula'] for r in comparison)
gcv = [r for r in comparison if r['Reaction ID'] == 'GHMT3' and r['Current locus'] in {
    'RU820_RS04695', 'RU820_RS04700', 'RU820_RS04705', 'RU820_RS04710'}]
assert len(gcv) == 4
assert all(r['Gene-Reaction Association'] == 'AFE_0295' for r in gcv)
assert all(r['2026 reaction score'] == '不评分' for r in comparison)
print('Alpha 2016 comparison: 13 non-direct rows, including four GHMT3 original A-P records')
