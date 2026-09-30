import csv
import json
import re
import urllib.parse
from collections import Counter
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
GFF = ROOT / '01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/genomic.gff'
D3 = ROOT / 'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv'
B1 = ROOT / 'B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv'
G5 = ROOT / 'G5_全新基因识别与编号/G5_02_全新基因身份核验.tsv'
MODEL = ROOT / 'A0_复现2016与2024模型/A0_主任务/脚本与环境/中间/mmc1.xlsx'
UNIPROT = OUT / 'g52_uniprot.json'
UNIPROT_DETAILS = OUT / 'g52_uniprot_details.json'
KEGG = OUT / 'g52_kegg.json'
BELLENBERG = OUT / 'g52_bellenberg2019_hits.json'


def attrs(s):
    return {k: urllib.parse.unquote(v) for part in s.split(';') if '=' in part for k, v in [part.split('=', 1)]}


def tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


genes = []
cds = {}
with GFF.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) != 9:
            continue
        a = attrs(c[8])
        locus = a.get('locus_tag', '')
        if c[2] == 'gene' and a.get('gene_biotype') == 'protein_coding':
            genes.append({'locus': locus, 'seqid': c[0], 'start': int(c[3]), 'end': int(c[4]), 'strand': c[6], 'gene_name': a.get('gene', '')})
        elif c[2] == 'CDS' and locus:
            cds.setdefault(locus, []).append({'wp': a.get('protein_id', ''), 'product': a.get('product', ''), 'start': int(c[3]), 'end': int(c[4]), 'strand': c[6], 'partial': a.get('partial', ''), 'pseudo': a.get('pseudo', '')})

assert len(genes) == 2929, len(genes)
group = genes[293:586]
assert len(group) == 293
assert group[0]['locus'] == 'RU820_RS01515'
assert group[-1]['locus'] == 'RU820_RS03035'

d3 = {r['当前locus']: r for r in tsv(D3)}
b1 = {}
for r in tsv(B1):
    b1.setdefault(r['当前GCF049位点'], []).append(r)
g5 = {r['当前locus']: r for r in tsv(G5)}
uniprot = json.load(UNIPROT.open(encoding='utf-8'))['by_wp'] if UNIPROT.exists() else {}
uniprot_details = json.load(UNIPROT_DETAILS.open(encoding='utf-8'))['by_wp'] if UNIPROT_DETAILS.exists() else {}
kegg = json.load(KEGG.open(encoding='utf-8')) if KEGG.exists() else {'by_afe':{},'ko_reactions':{}}
bellenberg_hits = json.load(BELLENBERG.open(encoding='utf-8')) if BELLENBERG.exists() else []
bellenberg_by_locus = {}
for hit in bellenberg_hits:
    bellenberg_by_locus.setdefault(hit['locus'], []).append(hit)

w = load_workbook(MODEL, read_only=True, data_only=True)
old = {}
old_by_reaction = {}
for row in list(w['Table 1'].values)[2:]:
    if not row[0]:
        continue
    fields = ['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe2 lb','Fe2 ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub']
    rec = dict(zip(fields, row))
    old_by_reaction[str(row[0])] = rec
    for afe in set(re.findall(r'AFE_\d{4}', str(row[7]))):
        old.setdefault(afe, []).append(rec)

records = []
reaction_rows = []
for ordinal, gene in enumerate(group, 294):
    locus = gene['locus']
    cur_cds = cds.get(locus, [])
    d = d3.get(locus, {})
    mappings = b1.get(locus, [])
    accepted_b1 = [m for m in mappings if m['映射状态'] == '完全匹配']
    g = g5.get(locus, {})
    afe = accepted_b1[0]['旧AFE位点'] if len(accepted_b1) == 1 else ''
    mapping_basis = accepted_b1[0]['映射方法'] if afe else ''
    if not afe and g.get('是否进入母表') == '否' and g.get('本轮核验结论') == '重新识别为原始GCA旧AFE基因对应' and re.fullmatch(r'AFE_\d{4}', g.get('旧AFE候选', '')):
        afe = g['旧AFE候选']
        mapping_basis = 'G5原始AFE辅助核验：' + g.get('本轮核验结论', '')
    status = '已确认' if afe else '未确认'
    linked = old.get(afe, []) if afe else []
    wp = cur_cds[0]['wp'] if cur_cds and len({c['wp'] for c in cur_cds}) == 1 else ''
    product = cur_cds[0]['product'] if cur_cds and len({c['product'] for c in cur_cds}) == 1 else ''
    u = uniprot.get(wp, [])
    assert len(u) <= 1
    u = u[0] if u else {}
    ud = uniprot_details.get(wp, {})
    u_rxn = ud.get('reactions', [])
    k = kegg['by_afe'].get(afe, {}) if afe else {}
    ko = [x.removeprefix('ko:') for x in k.get('ko', [])]
    k_rxn = sorted({x.removeprefix('rn:') for ko_id in k.get('ko', []) for x in kegg['ko_reactions'].get(ko_id, [])})
    bh = bellenberg_by_locus.get(locus, [])
    issue = []
    if len(cur_cds) != 1:
        issue.append('分段CDS需复核：' + str(len(cur_cds)) + '段；同一WP=' + str(len({c['wp'] for c in cur_cds}) == 1))
    if cur_cds and (min(c['start'] for c in cur_cds),max(c['end'] for c in cur_cds),cur_cds[0]['strand']) != (gene['start'],gene['end'],gene['strand']):
        issue.append('gene/CDS边界不一致')
    if not afe:
        issue.append('AFE历史映射未确认')
        if g.get('旧AFE候选'):
            issue.append('G5候选=' + g['旧AFE候选'] + '；状态=' + g.get('本轮核验结论',''))
    if d and (d.get('当前protein ID') != wp or urllib.parse.unquote(d.get('NCBI当前功能', '')) != product):
        issue.append('D3汇总与官方CDS不一致')
    if bh and not any(x['wp_exact'] for x in bh):
        if afe == 'AFE_0545':
            issue.append('Bellenberg2019旧WP比当前多N端19 aa；其后597 aa与当前完全一致；蛋白检出可作序列追溯，边界需保留')
        else:
            issue.append('Bellenberg2019旧WP与当前WP不一致；蛋白检出不能直接投射')
    if bh and afe in {'AFE_0363','AFE_0434'}:
        issue.append('Bellenberg2019与当前NCBI/UniProt功能判定冲突；须核验反应归属')
    issue.extend(['UniProt/InterPro/Pfam原始条目及证据代码未逐条复核','反应化学及2016适用性未逐条复核','同株原始实验检索未完成'])
    record = {
        '总序号': ordinal, '当前locus': locus, 'WP': wp, '当前product': product,
        '染色体': gene['seqid'], '起点': gene['start'], '终点': gene['end'], '链向': gene['strand'],
        '旧AFE': afe, '映射状态': status, '映射依据': mapping_basis,
        'D3 B1状态': d.get('B1映射状态',''),
        'D3 KEGG KO': d.get('KEGG KO',''), 'D3 KEGG EC': d.get('KEGG EC',''), 'D3 KEGG reaction': d.get('KEGG reaction',''),
        'D3 BioCyc reaction ID': d.get('BioCyc reaction ID',''), 'D3 Rhea交叉引用': d.get('Rhea精确交叉引用',''),
        'KEGG live KO': ';'.join(ko), 'KEGG live EC': ';'.join(x.removeprefix('ec:') for x in k.get('enzyme',[])),
        'KEGG live reaction': ';'.join(k_rxn), 'KEGG live pathways': ';'.join(x.removeprefix('path:') for x in k.get('pathway',[])),
        'KEGG AFE URL': f'https://www.kegg.jp/entry/afr:{afe}' if afe else '',
        'UniProt accession': u.get('accession',''), 'UniProt reviewed': u.get('reviewed',''),
        'UniProt protein name': u.get('protein_names',''), 'UniProt EC': u.get('ec',''),
        'InterPro ID': u.get('interpro',''), 'Pfam ID': u.get('pfam',''),
        'UniProt active site': u.get('active_site',''), 'UniProt binding site': u.get('binding_site',''),
        'UniProt URL': u.get('url',''),
        'UniProt catalytic activity': '；'.join(x['name'] for x in u_rxn),
        'UniProt reaction Rhea': ';'.join(sorted({c['id'] for x in u_rxn for c in x['cross_refs'] if c['database']=='Rhea'})),
        'UniProt reaction evidence': ';'.join(sorted({e['evidenceCode'] + '|' + e.get('source','') + ':' + e.get('id','') for x in u_rxn for e in x['evidence']})),
        '2019蛋白组命中': '是' if bh else '否',
        '2019同WP命中': '是' if any(x['wp_exact'] for x in bh) else '否' if bh else '',
        '2019序列追溯': '旧蛋白N端多19 aa，其后597 aa与当前完全一致' if afe == 'AFE_0545' and bh else '',
        '2019补表原文': '；'.join(f"Table page {x['page']}: {x['snippet']}" for x in bh),
        '2019来源DOI': '10.3389/fmicb.2019.00592' if bh else '',
        '2016反应ID': ';'.join(x['Reaction ID'] for x in linked),
        '2016反应数': len(linked), '2016对照状态': '有AFE-GPR文本关联；逐反应需复核' if linked else '未找到已确认AFE的GPR文本关联；不等于功能缺失',
        '精确reaction确认': '未完成', '2026动作': '待裁决', '审查状态': '预核对；未完成完整证据链',
        '关联行ID': ';'.join('G52-2016-'+str(x['Reaction ID'])+'-'+locus for x in linked),
        '未决点': '；'.join(issue),
        '官方来源': f'https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/',
        'WP来源': f'https://www.ncbi.nlm.nih.gov/protein/{wp}' if wp else '',
        '任务来源': 'https://docs.google.com/document/d/1tUSXPtaTYA8u1SQthnR4YBR5GwnlToUPudalLySccd8/edit',
    }
    records.append(record)
    for oldrec in linked:
        reaction_rows.append({'Candidate ID': 'G52-2016-'+str(oldrec['Reaction ID'])+'-'+locus,'Current locus':locus,'Legacy AFE locus':afe, **oldrec})

assert len(records) == 293 and len({r['当前locus'] for r in records}) == 293
assert [r['总序号'] for r in records] == list(range(294,587))
assert all(r['WP'].startswith('WP_') for r in records)
OUT.mkdir(exist_ok=True)
with (OUT/'g52_precheck.json').open('w', encoding='utf-8') as f:
    json.dump({'records': records, 'reaction_rows': reaction_rows}, f, ensure_ascii=False, indent=2)
print(json.dumps({'genes':len(records),'linked_reaction_rows':len(reaction_rows),'loci_with_2016_links':sum(bool(r['2016反应数']) for r in records),'afe_confirmed':sum(r['映射状态']=='已确认' for r in records),'afe_unconfirmed':sum(r['映射状态']!='已确认' for r in records),'gff_cds_issues':sum('CDS数量需复核' in r['未决点'] or 'gene/CDS边界不一致' in r['未决点'] for r in records)}, ensure_ascii=False))
