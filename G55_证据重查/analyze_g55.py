import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
GFF = ROOT / '01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/genomic.gff'
D3 = ROOT / 'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正'
ALPHA = ROOT / '00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx'

def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

genes = []
cds = {}
with GFF.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) != 9:
            continue
        attrs = dict(part.split('=', 1) for part in cols[8].split(';') if '=' in part)
        if cols[2] == 'gene' and attrs.get('gene_biotype') == 'protein_coding':
            genes.append(dict(locus=attrs['locus_tag'], start=int(cols[3]), end=int(cols[4]), strand=cols[6]))
        elif cols[2] == 'CDS' and 'locus_tag' in attrs:
            cds[attrs['locus_tag']] = dict(wp=attrs.get('protein_id',''), product=attrs.get('product',''), start=int(cols[3]), end=int(cols[4]), strand=cols[6])

assert len(genes) == 2929, len(genes)
group = genes[1172:1465]
assert len(group) == 293 and group[0]['locus'] == 'RU820_RS06230' and group[-1]['locus'] == 'RU820_RS07745'
cross = {r['当前locus']:r for r in rows(D3/'D3-1-R_全基因跨数据库交叉引用.tsv')}
ev = defaultdict(list)
for r in rows(D3/'D3-1-R_注释证据表.tsv'):
    ev[r['当前locus']].append(r)
wb = openpyxl.load_workbook(ALPHA, read_only=True, data_only=True)
s = wb.worksheets[0]
alpha = []
for row in list(s.values)[2:]:
    if row[16] == '2016 original' and row[0]:
        alpha.append(row)
out = []
for num,g in enumerate(group,1173):
    loc=g['locus']; c=cds[loc]; x=cross[loc]
    afe=x['旧AFE locus']
    old=[{'id':r[0],'name':r[1],'formula':r[2],'score':r[3],'ec':r[4],'pmid':r[5],'subsystem':r[6],'gra':r[7],'gpra':r[8],'pra':r[9],'bounds':list(r[10:16])} for r in alpha if afe and re.search(r'(?<![A-Za-z0-9_])'+re.escape(afe)+r'(?![A-Za-z0-9_])',str(r[7])+' '+str(r[8]))]
    out.append({'num':num,**g,**c,'afe':afe,'cross':x,'evidence':ev[loc],'old_reactions':old})
print('group',len(out),'unique',len(set(z['locus'] for z in out)),'cds_match',sum(z['start']==z['cross'].get('坐标','') for z in []))
print('mapped',sum(bool(z['afe']) for z in out),'old_linked',sum(bool(z['old_reactions']) for z in out),'old_rxns',sum(len(z['old_reactions']) for z in out))
for k in ['旧模型未覆盖功能候选','转运候选','重点分类冲突','功能判定']:
    print(k,Counter(z['cross'].get(k,'') for z in out))
print('candidates')
for z in out:
    x=z['cross']
    if x['旧模型未覆盖功能候选']=='是' or x['重点分类冲突'] not in ('','否'):
        print(z['num'],z['locus'],z['product'],'old',len(z['old_reactions']),x['候选依据'][:130],x['重点分类冲突'])
(Path(__file__).parent/'g55_data.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
