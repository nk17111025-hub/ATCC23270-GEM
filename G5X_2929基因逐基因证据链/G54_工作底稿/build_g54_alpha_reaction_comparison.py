import csv
import pathlib

root = pathlib.Path(__file__).parent
def read(name):
    with (root / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

fields = list(read('G54_Alpha逐基因审查.tsv')[0])
chemical = read('G54_2016反应化学错配.tsv')
old = ['2016 comparator Reaction ID', 'Reaction name', 'Reaction formula', 'Confidence Level',
       'EC Number', 'PMID', 'Subsystem', 'GRA', 'GPRA', 'PRA', 'Fe2+ lb', 'Fe2+ ub',
       'tetrathionate lb', 'tetrathionate ub', 'sulfur lb', 'sulfur ub']
alpha = fields[:16]
rows = []
for r in chemical:
    if not r['2016 comparator Reaction ID']:
        continue
    parts = {c: r[c].split(' || ') for c in old}
    n = len(parts[old[0]])
    assert all(len(v) == n for v in parts.values()), r['RU820 locus']
    for j in range(n):
        if parts[old[0]][j] == 'CYTBD' and r['RU820 locus'] in {'RU820_RS04575', 'RU820_RS04580'}:
            continue  # direct GPR and original A-P are already in the main Alpha rows
        row = {k: '' for k in fields}
        for a, c in zip(alpha, old):
            row[a] = parts[c][j]
        row.update({
            'Record class': '2016 既有反应原行化学比较；非当前基因直接 GPR',
            'Candidate ID': f"G54-{int(r['总序号']):04d}-{parts[old[0]][j]}",
            'Candidate group': 'G54',
            '2016 linked Reaction ID': parts[old[0]][j],
            '2026 reaction object': r['Reaction-ID clue(s)'],
            'Legacy AFE locus': r['旧 AFE locus'],
            'Current locus': r['RU820 locus'],
            'Evidence type': '2016 原表 A-P；KEGG 化合物集合/反应线索比较',
            'Evidence boundary': r['same-chemistry decision'] + ' ' + r['scope/remaining verification'],
            '2026 action': '反应比较记录；不得将原 GPR 自动移给当前基因',
            '2026 reaction score': '不评分',
            'Source ID': 'G54-R003;G54-R006',
            'Score scope': '仅比较旧反应；本行不构成当前基因 reaction/GPR 赋值',
        })
        rows.append(row)
with (root / 'G54_Alpha既有反应逐项对照.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(len(rows), 'per-gene comparison rows')
