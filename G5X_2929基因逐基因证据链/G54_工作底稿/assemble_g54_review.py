import csv
import pathlib

ROOT = pathlib.Path(__file__).parent


def read(name):
    with (ROOT / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(name, rows, fields):
    with (ROOT / name).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


index = read('G54_基因索引_初核.tsv')
model = {r['总序号']: r for r in read('G54_2016逐基因对照.tsv')}
db = {r['总序号']: r for r in read('G54_293数据库核查.tsv')}
lit = {r['总序号']: r for r in read('G54_293论文检索.tsv')}
reaction = {str(int(r['关联行ID'].split('-')[-1])): r for r in read('G54_反应候选_独立原始记录裁决.tsv')}
assert len(index) == len(model) == len(db) == len(lit) == 293

base_fields = [c.replace('Fe2+ ', 'Fe²⁺ ') for c in read('G54_Alpha结构审查表_初核.tsv')[0]]
out = []
status = []
for r in index:
    n = r['总序号']
    mr, dr, lr = model[n], db[n], lit[n]
    row = {k: '' for k in base_fields}
    row['Record class'] = '逐基因审查'
    row['Candidate ID'] = f'G54-{int(n):04d}'
    row['Candidate group'] = 'G54'
    row['2016 linked Reaction ID'] = mr['Direct Reaction ID']
    row['Legacy AFE locus'] = r['AFE 候选'] if r['映射状态'] != '完全匹配' and r['AFE 候选'] else r['旧 AFE locus']
    row['Current locus'] = r['当前 RU820 locus']
    row['Evidence type'] = '当前 RefSeq GFF/CDS；UniProt/InterPro/Pfam；KEGG/Rhea/BioCyc/BRENDA/TCDB 适用性；Europe PMC；2016 原模型'
    row['Source ID'] = 'G54-R001;G54-R002;G54-R003;G54-R004;G54-R005'
    row['2026 reaction score'] = '不评分'
    row['Score scope'] = '无可归属的精确小分子反应；不评分'
    row['Evidence boundary'] = (
        f"当前注释：{r['当前 NCBI product']}；WP={r['当前 WP protein_id']}；"
        f"AFE映射={r['映射状态']}；2016直接GPR={mr['Direct 2016 GPR status']}；"
        f"论文={lr['证据状态']}。"
    )
    row['2026 action'] = '保留逐基因记录；当前未证实可加入或修改的精确反应/GPR'
    if n in reaction:
        a = reaction[n]
        row['2026 reaction object'] = a['候选反应ID']
        row['Evidence boundary'] += ' 反应原始记录裁决：' + a['反应与基因适配性/建议']
        row['2026 action'] = '保留候选并按原始反应记录裁决；未达到正式改模条件'
        row['Source ID'] += ';G54-R006'
    if mr['Direct Reaction ID']:
        remap = {
            'Reaction ID': 'Direct Reaction ID',
            'Reaction Name': 'Reaction name',
            'Reaction Formula': 'Reaction formula',
            'Confidence Level': '2016 confidence',
            'EC Number': 'EC',
            'PMID': 'PMID',
            'Subsystem': 'Subsystem',
            'Gene-Reaction Association': 'GRA',
            'Gene-Protein-Reaction Association': 'GPRA',
            'Protein-Reaction-Association': 'PRA',
            'Fe²⁺ lb': 'Fe2+ lb',
            'Fe²⁺ ub': 'Fe2+ ub',
            'tetrathionate lb': 'tetrathionate lb',
            'tetrathionate ub': 'tetrathionate ub',
            'sulfur lb': 'sulfur lb',
            'sulfur ub': 'sulfur ub',
        }
        for alpha_col, model_col in remap.items():
            row[alpha_col] = mr[model_col]
    if n in ('881', '882'):
        row['2026 reaction object'] = 'CYTBD; RHEA:40527 / KEGG R11325；原模型 1.8 H+ / 2e− 待单独校准'
        row['Primary evidence DOI'] = '10.1186/1471-2164-10-394;10.3389/fmicb.2019.00592;10.1016/j.meteno.2016.03.003'
        row['Evidence type'] = '同株转录/蛋白丰度；异株氧化酶研究；数据库同源；2016 全网络拟合'
        row['Evidence boundary'] = '本株 CydA/B 表达已报告；没有单酶反应与 1.8 H+/2e− 直接测量。Brasseur 2004 实验株 ATCC19859。'
        row['2026 action'] = '保留 CYTBD 双基因 GPR；更新证据边界；质子系数进入模型质量复核'
        row['2026 reaction score'] = '2'
        row['Score scope'] = '精确 2016 反应定义，含 1.8 质子系数；序列和表达支持核心功能，系数为模型拟合'
        row['Source ID'] += ';G54-S005;G54-S006;G54-S007'
    elif n in ('901', '902'):
        row['2026 reaction object'] = '候选 R09395 / RHEA:26333；2 硫化 MoaD 载体 + cyclic pyranopterin phosphate + H2O → molybdopterin + 2 非硫化 MoaD + 2 H+'
        row['Evidence type'] = '同株未审校 UniProt/InterPro/Pfam 双亚基注释；Rhea/KEGG 反应；异株机制论文'
        row['Evidence boundary'] = 'MoaD/MoaE 双亚基候选，Rhea 式涉及蛋白载体；本株无直接酶活，旧模型无同式。'
        row['2026 action'] = '候选新增钼蝶呤合酶反应；先定义蛋白载体形态与再活化步骤，再审双基因 GPR'
        row['2026 reaction score'] = '2'
        row['Score scope'] = '候选 Rhea 完整反应；跨菌种生化加本株同源，非本株直接测量'
        row['Source ID'] += ';G54-R007'
    elif n in ('905', '906', '907', '908'):
        row['2016 linked Reaction ID'] = ''
        row['2026 reaction object'] = 'R01221 / RHEA:27758；与旧模型 GHMT3 核心物质集合相同'
        row['Evidence boundary'] = 'GcvT/H/PA/PB 为系统成员；旧 GHMT3 关联 GlyA，EC 2.1.2.1 与 Rhea 1.4.1.27 不符；GcvL 身份未决。'
        row['2026 action'] = '核查 GHMT3 原 GPR 和 EC；暂不新增同化学反应，不写未证实的 GcvL/AND GPR'
        row['Score scope'] = 'THF 立体、区室及完整本株 GPR 未决；不评分'
        row['Source ID'] += ';G54-R008'
    elif n == '886':
        row['2016 linked Reaction ID'] = ''
        row['2026 reaction object'] = 'Lpd 候选；K00382 多反应，GcvL 与 PDH E3 均未唯一确定'
        row['2026 action'] = '保留 Lpd 候选；核实组外 RU820_RS08740 同源和复合体归属后定 GPR'
    elif n == '1172':
        row['2026 reaction object'] = 'CobS 候选 K09882 / R05227 / RHEA:15341'
        row['Evidence boundary'] = '现产品 AAA ATPase，旧注释 CobS；CobN/T 伙伴及本株酶活未确认，异株 Rhea 实验不能赋予单基因完整反应。'
        row['2026 action'] = '保留 CobS 候选；先核蛋白域及 CobN/T 伙伴，再评估反应和 GPR'
        row['Source ID'] += ';G54-R009'
    elif n == '1169':
        row['Evidence boundary'] = '旧 AFE_1293 的 SQR 论文编号与结构序列冲突：PDB 3KPK 匹配组外 RU820_RS08260/AFE_1792，不匹配当前 115 aa WP_201763903.1。'
        row['2026 action'] = '纠正旧文献错误归属；本基因保持 hypothetical，不赋 SQR 反应或分数'
        row['Source ID'] += ';G54-R010;G54-S004;G54-S008'
    elif n == '913':
        row['2016 linked Reaction ID'] = ''
        row['2026 action'] = '保留 R00188/R11188 功能冲突；PHP 域不足以替换 BPNT 原 GPR'
    elif n == '1046':
        row['2026 action'] = '核苷脱氧核糖基转移酶候选；待定义底物特异性，PUNP1 为磷解不等同'
    elif n in ('1002', '1005'):
        row['2026 action'] = '记录蛋白二硫键反应；不据此更改小分子代谢模型 TRDR'
    elif n in ('937', '957', '976', '1080', '1083', '1099', '1101', '1113'):
        row['2026 action'] = '记录核酸/移动元件加工功能；本轮无适用的小分子代谢模型反应'
    elif n in ('989', '1061'):
        row['Evidence boundary'] += ' 精确旧 AFE 文献实验为异株，不能作本株转运测量。'
        row['2026 action'] = '保留异株转运功能候选；底物、方向及本株证据不足以定反应'
        row['Primary evidence DOI'] = '10.1007/s12011-012-9589-0' if n == '989' else '10.1007/s12223-013-0244-8'
        row['Source ID'] += ';G54-S002' if n == '989' else ';G54-S003'
    elif n == '887':
        row['Primary evidence DOI'] = '10.3390/genes9070347'
        row['Source ID'] += ';G54-S001'
        row['2026 action'] = '保留 glucan biosynthesis 候选；论文实验株 ATCC 53993，缺本株反应实测'
    elif n == '923':
        row['Primary evidence DOI'] = '10.3390/ijms140816901'
        row['Source ID'] += ';G54-S009'
        row['2026 action'] = '保留转录调控候选；同株基因组论文仅作计算预测，无代谢反应'
    elif n == '1016':
        row['Primary evidence DOI'] = '10.1371/journal.pone.0112226'
        row['Source ID'] += ';G54-S010'
        row['2026 action'] = '保留 MazF 家族候选；同株论文为家族分类，无已测代谢反应'
    elif n in ('883', '885', '890', '910', '911', '940', '946', '951', '971', '975', '988', '1034', '1050', '1095', '1096', '1102', '1110', '1112'):
        row['2026 action'] = '保留当前候选酶或辅因子功能；底物/产物或催化身份未充分确认，不新增精确反应'
    row['Evidence boundary'] += ' 数据库覆盖：' + dr['核查覆盖结论']
    out.append(row)
    status.append({
        '总序号': n,
        'RU820': r['当前 RU820 locus'],
        'WP': r['当前 WP protein_id'],
        '当前产品': r['当前 NCBI product'],
        '旧 AFE': row['Legacy AFE locus'],
        '旧 AFE 映射状态': r['映射状态'],
        '2016 直接 GPR': mr['Direct Reaction ID'],
        '反应线索': row['2026 reaction object'],
        '2026 动作': row['2026 action'],
        '2026 分数': row['2026 reaction score'],
        '数据库核查行': f'G54-{int(n):04d}',
        '论文检索行': f'G54-{int(n):04d}',
        '未决点': r['未决点'],
    })

write('G54_Alpha逐基因审查.tsv', out, base_fields)
write('G54_逐基因索引与状态.tsv', status, list(status[0]))
print('assembled', len(out), 'Alpha rows; numeric scores', sum(r['2026 reaction score'].isdigit() for r in out))
