import csv
import json
import pathlib

root = pathlib.Path(__file__).parent
names = [
    ('Alpha审查', 'G54_Alpha逐基因审查.tsv'),
    ('基因索引', 'G54_逐基因索引与状态.tsv'),
    ('Alpha反应比较', 'G54_Alpha既有反应逐项对照.tsv'),
    ('2016对照', 'G54_2016逐基因对照.tsv'),
    ('反应比较', 'G54_2016反应化学错配.tsv'),
    ('反应原始记录', 'G54_反应候选_待原始记录核查.tsv'),
    ('原始反应裁决', 'G54_反应候选_独立原始记录裁决.tsv'),
    ('数据库核查', 'G54_293数据库核查.tsv'),
    ('论文检索', 'G54_293论文检索.tsv'),
    ('证据来源', 'G54_证据来源表.tsv'),
    ('冲突待验证', 'G54_冲突缺证与待验证.tsv'),
]
data = []
for tab, name in names:
    with (root / name).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        rows = list(reader)
    assert all(len(r) == len(rows[0]) for r in rows), name
    data.append({'name': tab, 'rows': rows})
folder = root / '_artifact_build'
folder.mkdir(exist_ok=True)
path = folder / 'g54_workbook_data.json'
path.write_text(json.dumps(data, ensure_ascii=False), encoding='utf-8')
print(path, [(x['name'], len(x['rows']) - 1, len(x['rows'][0])) for x in data])
