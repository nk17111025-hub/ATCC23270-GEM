import csv
import pathlib

root = pathlib.Path(__file__).parent
def read(name):
    with (root / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

index = read('G54_基因索引_初核.tsv')
audit = {r['Candidate ID']: r for r in read('G54_Alpha逐基因审查.tsv')}
rows = []
for gene in index:
    number = int(gene['总序号'])
    aid = f'G54-{number:04d}'
    a = audit[aid]
    issues = []
    if gene['映射状态'] != '完全匹配':
        issues.append(('历史编号', gene['映射状态'] + '：' + gene['AFE 候选'], '按 WP/序列与旧基因组再核对；不得选一个多拷贝 AFE'))
    if gene['是否有反应ID线索'] == '是' and a['2026 reaction score'] == '不评分':
        issues.append(('反应线索', gene['KEGG reaction'] + ';' + gene['BioCyc reaction ID'] + ';' + gene['Rhea 交叉引用'], a['2026 action']))
    if number in (881, 882):
        issues.append(('反应定义', 'CYTBD 1.8 H+/2e− 为模型拟合，非本株 CydAB 单酶测量', '保持现行系数待模型质量复核；更新证据层级'))
    if number in (886, 905, 906, 907, 908):
        issues.append(('GPR/EC', 'GCV 复合体 Lpd 身份未决；GHMT3 原 GPR/EC 指 GlyA', '确定 Lpd 并审 GHMT3 GPR/EC，不重复添加同化学反应'))
    if number in (901, 902):
        issues.append(('新反应定义', 'MoaD/MoaE 候选反应含蛋白载体，2016 无同式', '定义 MoaD 硫化/非硫化载体及再活化步骤'))
    if number == 1169:
        issues.append(('文献旧编号冲突', 'AFE_1293 SQR 论文结构序列与现 RU820_RS08260/AFE_1792 对应', '本基因不赋 SQR；将论文转交组外位点'))
    if number == 1172:
        issues.append(('亚基不完整', 'CobS 候选缺 CobN/T 同株伙伴证据', '核实蛋白域与完整钴螯合酶复合体'))
    for kind, fact, next_step in issues:
        rows.append({
            'Candidate ID': aid,
            'RU820': gene['当前 RU820 locus'],
            'WP': gene['当前 WP protein_id'],
            '冲突/缺证类型': kind,
            '已核事实或边界': fact,
            '下一步验证': next_step,
            '当前裁决': a['2026 action'],
        })
fields = list(rows[0])
with (root / 'G54_冲突缺证与待验证.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fields, delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(len(rows), 'issue rows', len({r['Candidate ID'] for r in rows}), 'genes')
