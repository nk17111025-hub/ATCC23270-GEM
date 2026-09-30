import csv
import json
import re
from urllib.parse import unquote
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / 'G56_结果'
OUT.mkdir(exist_ok=True)

def rows(rel):
    with (ROOT / rel).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def toks(s):
    return [v.strip() for v in re.split(r'[;|]', str(s or '')) if v.strip()]

genes = []
for line in (ROOT / 'G56.txt').read_text(encoding='utf-8').splitlines():
    if re.match(r'^\d{4}｜RU820_RS', line):
        n, locus, wp, coord, product = line.split('｜', 4)
        genes.append(dict(number=int(n), locus=locus, wp=wp, coord=coord, product=product))
assert len(genes) == 293 and [g['number'] for g in genes] == list(range(1466, 1759))
assert len({g['locus'] for g in genes}) == 293

gff_path = ROOT / 'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff'
gff_gene, gff_cds = {}, {}
for line in gff_path.open(encoding='utf-8'):
    if line.startswith('#'):
        continue
    a = line.rstrip('\n').split('\t')
    if len(a) < 9 or a[2] not in ('gene', 'CDS'):
        continue
    attrs = {k: unquote(v) for k, v in (part.split('=', 1) for part in a[8].split(';') if '=' in part)}
    locus = attrs.get('locus_tag')
    if not locus:
        continue
    dest = gff_gene if a[2] == 'gene' else gff_cds
    dest[locus] = dict(seq=a[0], start=int(a[3]), end=int(a[4]), strand=a[6], **attrs)

b1 = defaultdict(list)
for x in rows('B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv'):
    b1[x['当前GCF049位点']].append(x)
cross = {x['当前locus']: x for x in rows('D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv')}
evidence = defaultdict(list)
for x in rows('D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_注释证据表.tsv'):
    evidence[x['当前locus']].append(x)

uniprot = defaultdict(list)
uniprot_afe = defaultdict(list)
for x in rows('G5X_2929基因逐基因证据链/G56_uniprot_taxon243159_20260926.tsv'):
    for wp in re.findall(r'WP_\d+\.\d+', x['RefSeq']):
        uniprot[wp].append(x)
    for afe_locus in re.findall(r'AFE_\d+', x['Gene Names']):
        uniprot_afe[afe_locus].append(x)

rhea_equations = {}
for x in rows('D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/Rhea_全基因EC查询汇总.tsv'):
    rhea_equations[x['Reaction identifier']] = x.get('Equation', '')

alpha = ROOT / '00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx'
sheet = openpyxl.load_workbook(alpha, read_only=True, data_only=True).worksheets[0]
it = sheet.iter_rows(values_only=True)
next(it)
alpha_headers = list(next(it)[:30])
model_by_id = {}
model_by_afe = defaultdict(list)
for r in it:
    if r[16] != '2016 original':
        continue
    values = list(r[:30]) + [None] * max(0, 30-len(r))
    rid = str(values[0] or '').strip()
    if not rid:
        continue
    model_by_id[rid] = values
    for afe in set(re.findall(r'AFE_\d+', ' '.join(str(v or '') for v in values[7:10]))):
        model_by_afe[afe].append(rid)

paper_map = {
    'AFE_1677': ['G56-P01'], 'AFE_1678': ['G56-P01'],
    'AFE_1683': ['G56-P01'], 'AFE_1684': ['G56-P01'], 'AFE_1685': ['G56-P01'],
    'AFE_1686': ['G56-P01'], 'AFE_1687': ['G56-P01'],
    'AFE_1688': ['G56-P01','G56-P05'], 'AFE_1689': ['G56-P01','G56-P05'],
    'AFE_1690': ['G56-P01','G56-P05'], 'AFE_1691': ['G56-P01','G56-P05'], 'AFE_1692': ['G56-P01'],
    'AFE_1997': ['G56-P02'], 'AFE_1998': ['G56-P02'], 'AFE_1999': ['G56-P02'],
    'AFE_1799': ['G56-P04','G56-P05'],
    'AFE_1847': ['G56-P03'], 'AFE_1883': ['G56-P03'],
    'AFE_1792': ['G56-P08'],
}

index = []
review = []
for g in genes:
    locus, wp = g['locus'], g['wp']
    gene, cds = gff_gene.get(locus), gff_cds.get(locus)
    identity = bool(gene and cds and cds.get('protein_id') == wp and cds.get('product') == g['product'] and
                    g['coord'] == f"{cds['seq']}:{cds['start']}-{cds['end']}({cds['strand']})" and
                    gene['start'] == cds['start'] and gene['end'] == cds['end'] and gene['strand'] == cds['strand'])
    bx = b1[locus]
    exact_b1 = [x for x in bx if x['映射状态'] == '完全匹配' and x['当前蛋白ID'] == wp]
    afe = exact_b1[0]['旧AFE位点'] if len(exact_b1) == 1 else ''
    mapped = '完全匹配' if afe else ('映射歧义' if bx else '未确认')
    cx = cross.get(locus, {})
    ux = uniprot[wp]
    old_ux = [u for u in uniprot_afe.get(afe,[]) if u not in ux] if afe else []
    reviewed = [x for x in ux if x['Reviewed'] == 'reviewed']
    db = {
        'NCBI': f'https://www.ncbi.nlm.nih.gov/protein/{wp}',
        'UniProt': ';'.join(f"https://www.uniprot.org/uniprotkb/{u['Entry']}" for u in ux),
        'InterPro': ';'.join(f"https://www.ebi.ac.uk/interpro/entry/InterPro/{v}" for u in ux for v in toks(u['InterPro'])),
        'Pfam': ';'.join(f"https://www.ebi.ac.uk/interpro/entry/pfam/{v}" for u in ux for v in toks(u['Pfam'])),
        'KEGG': f"https://www.kegg.jp/entry/afr:{cx.get('KEGG afr gene')}" if cx.get('KEGG afr gene') else '',
        'BioCyc': f"https://biocyc.org/AFER243159/keyword-search?keyword={cx.get('BioCyc gene ID')}" if cx.get('BioCyc gene ID') else '',
        'Rhea': ';'.join(f'https://www.rhea-db.org/rhea/{v.replace("RHEA:", "")}' for v in toks(cx.get('Rhea精确交叉引用'))),
    }
    rids = sorted(set(model_by_afe.get(afe, []))) if afe else []
    id_source = 'G56-S01;G56-S02;G56-S03;G56-S04;G56-S05'
    if ux or old_ux: id_source += ';G56-S06'
    if cx.get('KEGG afr gene'): id_source += ';G56-S07'
    if cx.get('BioCyc gene ID'): id_source += ';G56-S08'
    if cx.get('Rhea精确交叉引用'): id_source += ';G56-S09'
    if cx.get('BRENDA EC记录'): id_source += ';G56-S10'
    pids = paper_map.get(afe, [])
    if pids: id_source += ';' + ';'.join(pids)
    exact_rxn = bool(rids or cx.get('Rhea精确交叉引用'))
    if not identity:
        action = '证据不足暂缓'
        issue = '官方GFF身份与任务清单不一致，须先解决版本/坐标/WP。'
    elif not afe and bx:
        action = '证据不足暂缓'
        issue = '旧AFE映射有多拷贝或歧义；不能归入2016 GPR。'
    elif rids:
        action = '无需修改'
        issue = '旧模型已关联反应；仅凭数据库注释不能改写原反应或GPR。'
    elif cx.get('KEGG reaction') or cx.get('BioCyc reaction ID') or cx.get('Rhea精确交叉引用'):
        action = '证据不足暂缓'
        issue = '数据库候选尚未逐反应证明精确化学式、区室及2016等价状态。'
    elif re.search(r'transporter|permease|efflux|antiporter|symporter|exporter|TonB|porin', g['product'], re.I):
        action = '证据不足暂缓'
        issue = '转运底物、方向、膜侧或能量耦联未由该基因专属证据确定。'
    elif re.search(r'oxidoreductase|dehydrogenase|oxygenase|synthase|synthetase|transferase|hydrolase|phosphatase|kinase|decarboxylase|aldolase|epimerase|isomerase|desaturase|glyoxalase|phospholipase|lyase|ligase', g['product'], re.I) and not re.search(r'DNA|RNA|tRNA|rRNA|transposase|protease|peptidase|sigma factor|replication|partition', g['product'], re.I):
        action = '证据不足暂缓'
        issue = '仅有酶名/家族线索，未取得该基因专属精确底物—产物及2016等价反应核对。'
    else:
        action = '无需修改'
        issue = '当前证据未给出可加入代谢模型的精确反应。'
    if afe == 'AFE_1677':
        issue = '当前NorD名称与旧文献CbbO1指称冲突；KEGG K02448/R00294不能直接作为本蛋白催化反应。'
        action = '证据不足暂缓'
    if afe == 'AFE_1688':
        issue = 'CsoS3/Can1为羧酶体壳碳酸酐酶候选；旧HCO3E已有AFE_0287，羧酶体区室与可合并性未证。'
        action = '证据不足暂缓'
    if afe == 'AFE_1883':
        issue = 'PPC已存在但2016 GPR为AFE_1810；AFE_1883同工酶能力需独立证据，暂不写OR。'
        action = '证据不足暂缓'
    if afe == 'AFE_1799':
        issue = '45项历史裁决GLCP为仅更新证据：转录/蛋白检出不能核实糖原链长及旧反应计量。'
        action = '仅更新证据或编号'
    if afe == 'AFE_1847':
        issue = '45项历史裁决B-005暂缓：C4二羧酸/亚碲酸盐底物、膜侧、方向与耦联均未证。'
        action = '证据不足暂缓'
    if afe == 'AFE_1999':
        issue = '同株及异源表达证实AHL产物范围；各acyl-ACP前体、逐物种反应式和区室仍未定。'
        action = '证据不足暂缓'
    if afe == 'AFE_1792':
        issue = 'SQRED1旧PMID为通路预测；同株2025过表达支持硫氧化生理作用，但旧式q8膜侧及Sre OR分支未证，三培养条件上下界全0另列QC。'
        action = '仅更新证据或编号'
    if afe in {'AFE_1690','AFE_1691','AFE_1906','AFE_1907','AFE_1908','AFE_1910','AFE_1987'}:
        issue = '2016高分的所引证据未直接测本蛋白的逐反应催化；反应级行将2026分数按间接证据评2，旧行保留。'
        action = '仅更新证据或编号'
    paper_boundary = ''
    if afe in ['AFE_1688','AFE_1690','AFE_1691']:
        paper_boundary = '同株测得CO2条件下表达/蛋白丰度；未直接测定该蛋白催化反应。'
    elif afe == 'AFE_1999':
        paper_boundary = '同株检出AHL及转录，异源表达产AHL；未逐项测定底物到产物的反应式。'
    elif afe == 'AFE_1799':
        paper_boundary = '同株转录/蛋白组提示表达；未直接测定GLCP精确反应。'
    elif afe in ['AFE_1847','AFE_1883']:
        paper_boundary = '同株基因组序列及功能预测；未直接测定转运底物或PPC催化。'
    elif afe == 'AFE_1792':
        paper_boundary = '同株过表达改变硫氧化与浸矿表型；未证明旧SQRED1的q8底物、膜侧及Sre复合体OR分支。'
    elif pids:
        paper_boundary = '同株研究涉及该操纵子；未直接证实此蛋白的独立代谢反应。'
    else:
        paper_boundary = '当前检索未定位该位点的同株直接反应实验；数据库/同源注释不等于直接验证。'
    score = '不评分'
    if rids:
        score = '逐反应沿用2016原分数；本轮不迁移到新GPR'
    rowid = f"G56-G{g['number']:04d}"
    record = [None]*30
    record[16:30] = ['2026 gene review', rowid, 'G56', ';'.join(rids), '', afe or '未确认', locus,
                     '10.1128/AEM.71.11.7033-7040.2005' if afe == 'AFE_1999' else ('10.1128/aem.00170-25' if afe == 'AFE_1792' else ''),
                     '官方注释；数据库预测；同源/结构域；原始论文另列' if pids else '官方注释；数据库预测；同源/结构域',
                     paper_boundary, action, '不评分', id_source, '仅对精确反应评分；gene行不评分']
    review.append(record)
    historical = {'AFE_1799':'GLCP','AFE_1847':'B-005','AFE_1883':'B-003'}.get(afe,'')
    index.append({
        '总序号': g['number'], '当前locus': locus, 'WP': wp, '坐标链向': g['coord'], '当前NCBI product': g['product'],
        'GFF身份核对': '通过' if identity else '不通过', '旧AFE': afe or '未确认', '映射状态': mapped,
        'B1映射方法': exact_b1[0]['映射方法'] if afe else '',
        'UniProt accession': ';'.join(u['Entry'] for u in ux),
        'UniProt旧AFE关联但WP未匹配': ';'.join(u['Entry'] for u in old_ux),
        'UniProt人工审校': ';'.join(u['Entry'] for u in reviewed),
        'UniProt功能名': ';'.join(u['Protein names'] for u in ux),
        'UniProt最近修改日期': ';'.join(u['Date of last modification'] for u in ux),
        'UniProt条目版本': ';'.join(u['Entry version'] for u in ux),
        'InterPro': ';'.join(sorted({v for u in ux for v in toks(u['InterPro'])})),
        'Pfam': ';'.join(sorted({v for u in ux for v in toks(u['Pfam'])})),
        'KEGG gene': cx.get('KEGG afr gene',''), 'KO': cx.get('KEGG KO',''), 'KEGG EC': cx.get('KEGG EC',''),
        'KEGG reaction': cx.get('KEGG reaction',''), 'BioCyc gene': cx.get('BioCyc gene ID',''),
        'BioCyc reaction': cx.get('BioCyc reaction ID',''), 'Rhea精确交叉引用': cx.get('Rhea精确交叉引用',''),
        'Rhea仅EC候选': cx.get('Rhea仅EC候选',''), 'BRENDA EC记录': cx.get('BRENDA EC记录',''),
        'Rhea交叉引用反应式': ';'.join(f'{v}: {rhea_equations.get(v, "")}' for v in toks(cx.get('Rhea精确交叉引用'))),
        'TCDB基因级编号': '未取得' if re.search(r'transporter|permease|efflux|antiporter|symporter|exporter|TonB|porin', g['product'], re.I) else '不适用',
        '数据库分歧': cx.get('数据库分歧',''), '2016关联反应': ';'.join(rids),
        '是否找到精确reaction': '2016已关联' if rids else ('数据库候选待核' if exact_rxn or cx.get('KEGG reaction') or cx.get('BioCyc reaction ID') else '否'),
        '2026动作': action, '关联行ID': rowid, '未决点': issue, '论文Source ID': ';'.join(pids),
        '证据边界': paper_boundary, '来源ID': id_source, 'NCBI链接': db['NCBI'],
        '45项历史裁决ID': historical,
        'UniProt链接': db['UniProt'], 'InterPro链接': db['InterPro'], 'Pfam链接': db['Pfam'],
        'KEGG链接': db['KEGG'], 'BioCyc链接': db['BioCyc'], 'Rhea链接': db['Rhea'],
        '数据库预测证据': f"KEGG KO={cx.get('KEGG KO','')}; EC={cx.get('KEGG EC','')}; reaction={cx.get('KEGG reaction','')}; BioCyc reaction={cx.get('BioCyc reaction ID','')}" if any(cx.get(k) for k in ('KEGG KO','KEGG EC','KEGG reaction','BioCyc reaction ID')) else '无该位点专属候选',
        '同源及结构域证据': f"精确WP匹配UniProt={';'.join(u['Entry'] for u in ux)}; InterPro={';'.join(sorted({v for u in ux for v in toks(u['InterPro'])}))}; Pfam={';'.join(sorted({v for u in ux for v in toks(u['Pfam'])}))}" if ux else ('仅旧AFE关联，当前WP未匹配' if old_ux else '未取得当前WP直接条目'),
        '实验论文实测及边界': paper_boundary,
    })

for ix in index:
    locus, afe = ix['当前locus'], ix['旧AFE']
    if afe == '未确认':
        continue
    for rid in toks(ix['2016关联反应']):
        base = model_by_id.get(rid)
        if not base:
            continue
        rec = list(base)
        score = 2 if base[3] in (2,3) else '不评分'
        score_boundary = '同株序列/同源与功能注释支持反应类别；精确底物、方向、区室及该蛋白直接周转未被引用PMID逐项证实。'
        score_source = 'G56-S01;G56-S02;G56-S04;G56-S05;G56-S06'
        if str(base[5]) == '238956':
            score_boundary = '2016所引PMID 238956测脂多糖脂肪酸组成，未测本Fab酶的该条链长反应；本轮按序列/功能注释评2。'
            score_source += ';G56-P06'
        elif rid == 'UHGA':
            score_boundary = '2016所引PMID 15044494直接研究GnnA/GnnB，不能证明LpxH/UHGA式；本轮仅按注释评2。'
            score_source += ';G56-P07'
        elif rid == 'SQRED1':
            score_boundary = '2009所引PMID为通路预测；2025同株AFE_1792过表达改变硫氧化，未测旧q8反应膜侧及Sre OR分支；界限全0另列QC。'
            score_source += ';G56-P08'
        elif rid == 'RUBISCO':
            score_boundary = '同株cbb1表达与蛋白组支持存在；G56未定位该CbbL1/S1复合体对旧式的直接酶学测定。'
            score_source += ';G56-P01;G56-P05'
        rec[16:30] = ['2026 reaction review', f"G56-R{ix['总序号']}-{rid}", 'G56', rid, rid,
                      afe, locus, '', '数据库预测；同源证据；原始论文另列', score_boundary,
                      '仅更新证据或编号' if rid in ('GLCP','SQRED1') or str(base[5]) in ('238956','15044494') else '无需修改',
                      score, score_source,
                      f'具体2016反应；2016原Confidence={base[3]}；2026按间接证据评{score}；不转给新GPR/区室/方向']
        review.append(rec)

# 关键未决反应保留反应级行。2016核心列仅在已有同式反应时复制；
# 未经接受的候选不写入原行，也不把数据库反应号当作模型Reaction ID。
for locus, linked, object_name, doi, boundary, src in [
    ('RU820_RS07800','HCO3E','羧酶体壳CsoS3催化的CO2/HCO3-互变；模型无羧酶体区室','10.3389/fmicb.2019.00603','同株表达证据；无CsoS3纯化酶反应测定，也不能把细胞质HCO3E直接赋给它。','G56-P01'),
    ('RU820_RS08690','PPC','AFE_1883拟补PPC GPR，历史B-003继续暂缓','10.1186/1471-2164-9-597','旧PPC的2016分数仅属既有反应；AFE_1883独立催化能力未证。','G56-P03'),
    ('RU820_RS09235','','AfeI多种长链AHL生成；尚不能定义单一精确反应','10.1128/AEM.71.11.7033-7040.2005','同株及异源表达测得AHL产物；各acyl-ACP底物与逐产物计量未直接测定。','G56-P02'),
]:
    ix = next(x for x in index if x['当前locus']==locus)
    base = model_by_id.get(linked)
    rec = list(base) if base else [None]*30
    rec[16:30] = ['2026 reaction candidate—hold', f'G56-H{ix["总序号"]}', 'G56', linked,
                  object_name, ix['旧AFE'], locus, doi, '数据库预测；同源证据；实验论文分列',
                  boundary, '证据不足暂缓', '不评分', f'G56-S01;G56-S02;G56-S05;{src}',
                  '未形成可接受的精确新反应；既有2016分数只见核心列']
    review.append(rec)

assert len(index) == 293 and len(review) >= 293
source_headers = ['Source ID','Full citation','Year','DOI','Source type','Strain','Experimental type',
                  'Genes / proteins studied','Related Reaction IDs / candidate IDs',
                  'What was directly measured','Applicable confidence level','Evidence limitation']
sources = [
['G56-S01','NCBI RefSeq GCF_049532655.1 genomic.gff; https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/',2026,'','Official annotation','ATCC 23270','Genome sequence/annotation','G56全部293位点','G56-G1466—G56-G1758','坐标、CDS、WP与product来自冻结官方文件','不能独立评分','product多为计算注释；文件SHA-256见报告。'],
['G56-S02','B1 旧新基因映射.tsv; 本项目已核验的WP一致/邻域映射',2026,'','Mapping evidence','ATCC 23270','Sequence mapping','G56全组','G56全组','旧AFE与现行WP的映射证据','不评分','歧义映射标未确认。'],
['G56-S03','G56任务文档; https://docs.google.com/document/d/1CCJTvgr4-2G5r7YnB3GRZc6W3tCWzreoGFzXru6FMHs',2026,'','Frozen task list','ATCC 23270','List specification','G56全部293位点','G56全组','责任清单及当前product','不评分','任务清单不是酶学验证。'],
['G56-S04','D3-1-R全基因跨数据库交叉引用.tsv及注释证据表.tsv; 项目快照',2026,'','Database cross-reference snapshot','ATCC 23270/旧AFE','Computational integration','G56全组','G56全组','数据库ID、KO/EC及候选反应交叉引用','2至多；须按反应核定','来源存在转注；MetaCyc未独立获取，数据库交叉引用不等于湿实验。'],
['G56-S05','Campodonico et al. (2016), iMC507 original supplement Table 1; https://doi.org/10.1016/j.meteno.2016.03.003',2016,'10.1016/j.meteno.2016.03.003','Published model baseline','ATCC 23270','Model reconstruction','2016旧AFE GPR','2016各已关联反应','原始Reaction ID、式、分数、GPR及三种培养条件上下界','保留原Confidence列','模型收录及原分数不是本轮新增GPR或化学反应的实验验证。'],
['G56-S06','UniProtKB taxon 243159 2026-09-26 TSV快照; https://rest.uniprot.org/uniprotkb/stream?query=organism_id%3A243159&format=tsv',2026,'','Protein annotation/database','ATCC 23270 / DSM 14882','Sequence annotation','精确WP匹配的G56蛋白','G56全组','reviewed状态、名称、InterPro/Pfam/EC交叉引用、条目修改日期与版本','2至多，按反应核定','多数为自动或同源注释；未匹配WP者不借旧AFE冒充精确蛋白。'],
['G56-S07','KEGG afr基因/KO/EC/Reaction本地导出; https://www.kegg.jp/kegg-bin/show_organism?org=afr',2026,'','Pathway database','旧ATCC 23270 AFE','Computational annotation','映射到afr的G56基因','G56全组','KO、EC与候选反应','2至多，按反应核定','KO多对多，不等于该基因催化全部关联反应。'],
['G56-S08','BioCyc AFER243159 2026-09项目原始导出; https://biocyc.org/AFER243159/organism-summary',2026,'','Tier-3 pathway database','旧ATCC 23270 AFE','Computational pathway reconstruction','映射到BioCyc的G56基因','G56全组','gene→protein→reaction ID链','2至多，按反应核定','Tier-3与MetaCyc可能共享继承链；需独立实验佐证。'],
['G56-S09','Rhea EC/KEGG/MetaCyc cross-reference 2026-09项目导出; https://www.rhea-db.org/',2026,'','Reaction chemistry database','通用反应','Curated reaction/EC mapping','G56有Rhea ID的位点','G56全组','反应式、ChEBI与交叉引用','2至多，基因归属另证','精确反应式不证明本株该基因、区室或方向。'],
['G56-S10','BRENDA organism/EC本地导出; https://www.brenda-enzymes.org/',2026,'','Enzyme database','混合菌株','EC/酶学记录','G56有EC记录的位点','G56全组','EC关联及文献线索','按原实验菌株与精确反应核定','非同株或未追到原实验者仅为候选。'],
['G56-P01','Esparza et al. (2019), Effect of CO2 Concentration on Uptake and Assimilation of Inorganic Carbon in the Extreme Acidophile Acidithiobacillus ferrooxidans; https://doi.org/10.3389/fmicb.2019.00603',2019,'10.3389/fmicb.2019.00603','Original experimental paper','ATCC 23270','CO2培养、RT-qPCR、蛋白分析','cbb1/carboxysome locus, including AFE_1688 context','G56-H1476及相邻基因','CO2条件下生长及部分cbb基因转录/蛋白丰度','2，不能直接给酶反应4','CsoS3的碳酸酐酶催化未直接测得；表中其他cso基因多为预测功能。'],
['G56-P02','Farah et al. (2005), Evidence for a Functional Quorum-Sensing Type AI-1 System in the Extremophilic Bacterium Acidithiobacillus ferrooxidans; https://doi.org/10.1128/AEM.71.11.7033-7040.2005',2005,'10.1128/AEM.71.11.7033-7040.2005','Original experimental paper','ATCC 23270; AfeI heterologous in E. coli','AHL LC-MS/MS, RT-PCR, heterologous expression','AFE_1997/1998/1999; AfeI','G56-H1750','同株AHL产生与转录；表达AfeI的E. coli产生多种AHL','特定产物支持功能；精确反应暂不评分','各acyl-ACP底物、逐产物计量未直接测定；异源表达须另标。'],
['G56-P03','Valdés et al. (2008), Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications; https://doi.org/10.1186/1471-2164-9-597',2008,'10.1186/1471-2164-9-597','Genome paper/computational annotation','ATCC 23270','Genome sequencing/pathway inference','AFE_1847、AFE_1883等','B-003;B-005','同株基因组序列；功能主要由同源预测','2，限于序列/功能推测','没有AFE_1883催化PPC或AFE_1847转运底物实测。'],
['G56-P04','Mamani et al. (2016), Insights into the Quorum Sensing Regulon of the Acidophilic Acidithiobacillus ferrooxidans; https://doi.org/10.3389/fmicb.2016.01365',2016,'10.3389/fmicb.2016.01365','Original experimental paper','ATCC 23270','QS类似物处理转录组','AFE_1799等','GLCP','RNA丰度变化','2，间接','不测GLCP底物链长、产物计量或活性。'],
['G56-P05','Bellenberg et al. (2019), Proteomics Reveal Enhanced Oxidative Stress Responses and Metabolic Adaptation in Acidithiobacillus ferrooxidans Biofilm Cells on Pyrite; https://doi.org/10.3389/fmicb.2019.00592',2019,'10.3389/fmicb.2019.00592','Original experimental paper','DSM 14882T = ATCC 23270','Fe(II)/黄铁矿生物膜蛋白组','AFE_1688/1689/1690/1691/1799等','GLCP;RUBISCO;G56-H1476','蛋白检出与相对丰度','2，间接','蛋白丰度不证明GLCP、Rubisco或CsoS3的精确化学式。'],
['G56-P06','Hirt & Vestal (1975), Physical and chemical studies of Thiobacillus ferrooxidans lipopolysaccharides; https://doi.org/10.1128/jb.123.2.642-650.1975',1975,'10.1128/jb.123.2.642-650.1975','Original experimental paper','T. ferrooxidans; 菌株与ATCC23270同一性未核实','LPS成分分析','脂多糖脂肪酸；非G56 Fab酶','3OAR1-17;3OAS1-15;3OAS27-28;MCOATA','LPS脂肪酸组成','不能用于Fab精确反应3分','不能证明AFE_1906/1907/1908/1910的逐链长催化；旧PMID误归因。'],
['G56-P07','Sweet, Ribeiro & Raetz (2004), Oxidation and Transamination of the 3-position of UDP-N-Acetylglucosamine by Enzymes from Acidithiobacillus ferrooxidans; https://doi.org/10.1074/jbc.M400596200',2004,'10.1074/jbc.M400596200','Original biochemical paper','A. ferrooxidans来源；原文需核菌株','重组GnnA/GnnB及产物NMR','GnnA/GnnB；非LpxH AFE_1987','UHGA','GnnA/GnnB协同形成UDP-GlcNAc3N','不能用于UHGA直接催化3/4分','旧UHGA PMID对象错配；不得迁移酶学强度。'],
['G56-P08','Jung, Inaba & Banta (2025), Overexpression of sulfide:quinone reductase (SQR) in Acidithiobacillus ferrooxidans enhances sulfur, pyrite, and pyrrhotite oxidation; https://doi.org/10.1128/aem.00170-25',2025,'10.1128/aem.00170-25','Original genetic/physiological paper','ATCC 23270及AFE_1792过表达株','过表达、硫氧化和浸矿表型','AFE_1792','SQRED1','过表达相关硫氧化/浸矿变化与产物测量','3限于生理功能；旧精确q8式评2','未测旧模型Sre OR分支、q8膜侧和旧反应全部计量。'],
]
assert len({r[0] for r in sources}) == len(sources)
with (OUT/'G56_证据来源.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(source_headers);w.writerows(sources)
with (OUT/'G56_审查数据.json').open('w', encoding='utf-8') as f:
    json.dump({'alpha_headers': alpha_headers, 'index': index, 'review': review,
               'sources': sources, 'source_headers': source_headers,
               'model_by_id': model_by_id, 'counts': dict(Counter(x['2026动作'] for x in index))}, f, ensure_ascii=False)
with (OUT/'G56_基因索引.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(index[0]),delimiter='\t');w.writeheader();w.writerows(index)
print('genes',len(index),'identity pass',sum(x['GFF身份核对']=='通过' for x in index),
      'B1 exact',sum(x['旧AFE']!='未确认' for x in index),'UniProt exact WP',sum(bool(x['UniProt accession']) for x in index),
      'model covered',sum(bool(x['2016关联反应']) for x in index),'actions',dict(Counter(x['2026动作'] for x in index)))
