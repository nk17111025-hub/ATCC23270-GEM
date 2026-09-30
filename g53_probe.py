import csv, re, collections
from pathlib import Path
from openpyxl import load_workbook

root=Path('.')
genes=[]
for line in (root/'G53.txt').read_text(encoding='utf-8-sig').splitlines():
    m=re.match(r'^(\d{4})｜(RU820_RS\d+)｜(WP_[\d.]+)｜([^｜]+)｜(.+)$',line)
    if m: genes.append(m.groups())
ids={x[1] for x in genes}
def readtsv(path,key):
    with open(path,encoding='utf-8-sig',newline='') as f:
        return {r[key]:r for r in csv.DictReader(f,delimiter='\t') if r.get(key) in ids}
cross=readtsv(root/'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv','当前locus')
maprows={}
with open(root/'B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv',encoding='utf-8-sig',newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['当前GCF049位点'] in ids: maprows.setdefault(r['当前GCF049位点'],[]).append(r)
model={}
w=load_workbook(root/'00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx',read_only=True,data_only=True)
for row in list(w.worksheets[0].values)[2:]:
    v=' '.join(str(x or '') for x in row[7:10])
    for afe in set(re.findall(r'AFE_\d+',v)):model.setdefault(afe,[]).append(row)
print('genes',len(genes),'cross',len(cross),'map',len(maprows),'unique',len(ids))
print('mapping status',collections.Counter(tuple(x['映射状态'] for x in maprows.get(g[1],[])) for g in genes))
print('2016 GPR',sum(bool(cross.get(g[1],{}).get('2016 GPR反应')) for g in genes))
print('KEGG',sum(bool(cross.get(g[1],{}).get('KEGG reaction')) for g in genes),'BioCyc',sum(bool(cross.get(g[1],{}).get('BioCyc reaction ID')) for g in genes))
for seq,locus,wp,coord,prod in genes:
    c=cross.get(locus,{})
    if c.get('2016 GPR反应') or c.get('KEGG reaction') or c.get('BioCyc reaction ID') or c.get('旧模型未覆盖功能候选')=='是':
        print(seq,locus,prod,'AFE',c.get('旧AFE locus'),'2016',c.get('2016 GPR反应'),'KEGG',c.get('KEGG reaction'),'BioCyc',c.get('BioCyc reaction ID'))
