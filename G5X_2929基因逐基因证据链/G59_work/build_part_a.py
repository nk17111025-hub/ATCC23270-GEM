from pathlib import Path
import csv,re,openpyxl,json
from collections import defaultdict
root=Path('.')
out=root/'G5X_2929基因逐基因证据链/G59_work/part_a.xlsx'
gff=root/'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff'
mapf=root/'B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv'
d3=root/'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正'
xlsx=root/'A0_复现2016与2024模型/A0_主任务/脚本与环境/中间/mmc1.xlsx'
def attrs(s):
 d={}
 for x in s.split(';'):
  if '=' in x:
   k,v=x.split('=',1);d[k]=v
 return d
def read_tsv(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
# Freeze ordered protein-coding genes by source GFF and join CDS on locus.
genes=[]; cds={}
for line in gff.open(encoding='utf-8'):
 if not line or line[0]=='#': continue
 f=line.rstrip('\n').split('\t'); a=attrs(f[8]); typ=f[2]
 if typ=='CDS' and a.get('locus_tag'): cds[a['locus_tag']]=a
 if typ=='gene' and a.get('gene_biotype')=='protein_coding':
  genes.append({'seqid':f[0],'start':int(f[3]),'end':int(f[4]),'strand':f[6],'locus':a.get('locus_tag',''),'gene_product':a.get('product','')})
assert len(genes)==2929, len(genes)
# Use the G59 responsibility document to avoid pseudogenes and noncoding loci in the genomic interval.
source=(root/'G5X_2929基因逐基因证据链/G59_work/G59_source.txt').read_text(encoding='utf-8-sig')
entries=[]
for line in source.splitlines():
 m=re.match(r'^(\d+)｜(RU820_RS\d+)｜(WP_[^｜]+)｜([^｜]+)｜(.+)$',line.strip())
 if m and 2345<=int(m.group(1))<=2442: entries.append((int(m.group(1)),m.group(2),m.group(3),m.group(4),m.group(5)))
assert len(entries)==98 and entries[0][1]=='RU820_RS12360' and entries[-1][0]==2442, (len(entries), entries[0], entries[-1])
gffby={g['locus']:g for g in genes}
seg=[]
for idx,loc,wp,coord,product in entries:
 g=gffby.get(loc)
 assert g is not None,(idx,loc)
 assert cds.get(loc,{}).get('protein_id')==wp,(idx,loc,wp,cds.get(loc,{}).get('protein_id'))
 g['wp']=wp;g['product']=product;g['gff_coord_doc']=coord;g['source_index']=idx;seg.append(g)
# authoritative GFF protein mapping.
for g in seg:
 a=cds.get(g['locus'],{}); g['wp']=a.get('protein_id','');g['product']=a.get('product',g['gene_product']);g['cds_join']='exact locus_tag CDS match' if a else 'CDS not found'
# B1 history maps only one-to-one exact current locus AND WP.
mp=read_tsv(mapf); mapby=defaultdict(list)
for r in mp:
 for loc in r.get('当前GCF049位点','').split('|'):
  if loc: mapby[loc].append(r)
# D3 all-db records + annotation per gene.
files={n:d3/n for n in ['D3-1-R_全基因跨数据库交叉引用.tsv','D3-1-R_注释证据表.tsv','D3-1-R_模型覆盖状态.tsv','D3-1-R_旧模型未覆盖功能候选.tsv','D3-1-R_数据库分歧.tsv']}
d3data={k:read_tsv(v) for k,v in files.items()}
# 2016 Table 1 reactions
wb0=openpyxl.load_workbook(xlsx,read_only=True,data_only=True); s=wb0['Table 1']; rows=list(s.iter_rows(values_only=True)); h=[str(x).strip() if x is not None else '' for x in rows[1]]
ridx={k:i for i,k in enumerate(h)}
reactions=[]
for row in rows[2:]:
 if not row[0]:continue
 reactions.append({k:(row[i] if i<len(row) else None) for k,i in ridx.items()})
# detect 2016 rows by mapped legacy AFE token.
for i,g in enumerate(seg,2345):
 g['idx']=g['source_index']
 candidates=mapby.get(g['locus'],[])
 exact=[r for r in candidates if r.get('当前蛋白ID','')==g['wp'] and r.get('映射状态','')=='完全匹配']
 g['legacy']='; '.join(dict.fromkeys(r.get('旧AFE位点','') for r in exact if r.get('旧AFE位点'))) if len(exact)==1 else ('未确认' if not exact else '映射歧义')
 g['map_note']=('; '.join(r.get('映射方法','')+'; '+r.get('映射状态','') for r in exact) if exact else ('non-unique candidate in B1, unresolved' if candidates else 'no exact B1 mapping row'))
 g['cross']= [r for r in d3data['D3-1-R_全基因跨数据库交叉引用.tsv'] if r.get('当前locus')==g['locus']]
 g['ann']=[r for r in d3data['D3-1-R_注释证据表.tsv'] if r.get('当前locus')==g['locus']]
 g['cover']=[r for r in d3data['D3-1-R_模型覆盖状态.tsv'] if r.get('当前locus')==g['locus']]
 g['cand']=[r for r in d3data['D3-1-R_旧模型未覆盖功能候选.tsv'] if r.get('当前locus')==g['locus']]
 g['div']=[r for r in d3data['D3-1-R_数据库分歧.tsv'] if r.get('当前locus')==g['locus']]
 g['rx16']=[r for r in reactions if g['legacy'] not in ('未确认','映射歧义') and g['legacy'] and g['legacy'] in str(r.get('Gene-Reaction Association',''))]
 # text digests with original-source identifiers and evidence grade as extracted by D3.
 g['dbsum']='; '.join(f"{r.get('数据库','')}:{r.get('记录ID','')} ({r.get('证据字段','')}; {r.get('证据等级/类型','')})" for r in g['ann'])
 g['crosssum']='; '.join(f"{r.get('数据库','')}={r.get('记录ID','')}" for r in g['cross'])
# make workbook
w=openpyxl.Workbook(); inv=w.active; inv.title='Gene index'
headers=['Total no.','Current locus','WP protein_id','GFF coordinates','Current NCBI product','Legacy AFE locus','B1 mapping basis/status','Exact reaction found?','2016 linked Reaction IDs','Action','Evidence summary','Unresolved points','Associated review row ID','Identity evidence source']
inv.append(headers)
for g in seg:
 # D3 is a historical cross-db aggregate; presence does not qualify as exact experimental reaction.
 exact=False
 rxids='; '.join(str(r.get('Reaction ID')) for r in g['rx16'])
 evidence=f"GFF/CDS exact identity; B1: {g['map_note']}; D3 annotations: {g['dbsum'] or 'none in local D3 table'}; xrefs: {g['crosssum'] or 'none in local D3 table'}"
 action='仅更新证据或编号' if g['legacy'] not in ('未确认','映射歧义') else '证据不足暂缓'
 unresolved='当前功能的精确底物/产物、辅因子、方向/区室及原始实验未逐项确证；不据NCBI product定反应。' if not exact else '未核实是否有本株直接实验；需补UniProt/InterPro/Pfam原始条目。'
 inv.append([g['idx'],g['locus'],g['wp'],g['gff_coord_doc'],g['product'],g['legacy'],g['map_note'],'未确认（D3注释/xref不是精确反应的直接验证）',rxids,action,evidence,unresolved,f"A-{g['idx']}",'GFF GCF_049532655.1 + same-locus CDS; B1 mapping TSV'])
# Alpha structural audit rows: one per gene review; fields match A-P then Q-AD. No invented reaction chemistry.
a=w.create_sheet('Alpha review rows')
alpha=['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe2 lb','Fe2 ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub','Record class','Candidate ID','Candidate group','2016 linked Reaction ID','2026 reaction object','Legacy AFE locus','Current locus','Primary evidence DOI','Evidence type','Evidence boundary','2026 action','2026 reaction score','Source ID','Score scope']
a.append(alpha)
for g in seg:
 if g['rx16']:
  for r in g['rx16']:
   a.append([r.get('Reaction ID'),r.get('Reaction Name'),r.get('Reaction Formula '),r.get('Confidence Level'),r.get('EC Number'),r.get('PMID'),r.get('Subsystem'),r.get('Gene-Reaction Association'),r.get('Gene-Protein-Reaction Association'),r.get('Protein-Reaction-Association'),r.get('lb'),r.get('ub'),r.get('lb'),r.get('ub'),r.get('lb'),r.get('ub'),'2016 linked reaction audit',f"G59-{g['idx']}-{r.get('Reaction ID')}",'G59',r.get('Reaction ID'),'Reaction existing in iMC507; chemical equation retained verbatim',g['legacy'],g['locus'],'','2016 model record; not direct gene-specific experimental revalidation','2016 reaction/GPR stated in iMC507; this review has not independently verified exact substrate/cofactor/localization claims.','仅更新证据或编号',r.get('Confidence Level'),'G59-SRC-2016','2016 original reaction score; no transfer to new chemistry/GPR'])
 else:
  a.append(['','','','','','','','','','','','','','','','','gene_review',f"G59-{g['idx']}",'G59','','No exact reaction assigned from evidence reviewed',g['legacy'],g['locus'],'','database annotation / homology / genome context only' if g['ann'] else 'identity and local annotation records only','Direct evidence not located in reviewed local records; exact reaction chemistry, direction and localization remain unconfirmed.','证据不足暂缓','不评分','G59-SRC-GFF; G59-SRC-B1; G59-SRC-D3','No reaction defined; no score'])
# Evidence rows from D3 per locus, preserving accession/source identity and original tier labels.
ev=w.create_sheet('Evidence ledger'); ev.append(['Source ID','Total no.','Current locus','Database/source','Record accession','Evidence field','Description / value','Evidence grade/source type','Source file','Source limit'])
ev.append(['G59-SRC-GFF','','','NCBI RefSeq GFF/CDS','GCF_049532655.1 / NZ_CP136162.1','assembly, gene and CDS identity','Official fixed assembly coordinate and CDS protein_id from project local RefSeq files','primary sequence annotation','B1 standardization data/GCF_049532655.1/{genomic.gff,protein.faa}','Product is annotation; no enzyme activity implied'])
ev.append(['G59-SRC-B1','','','B1 legacy mapping','旧新基因映射.tsv','protein_id exact match','Legacy AFE accepted only for one-to-one current locus and exact WP protein ID','sequence identity mapping','B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv','Multiple candidate mapping rows retained as unresolved'])
ev.append(['G59-SRC-D3','','','D3 local multi-database scan','D3-1-R aggregate tables','cross-reference leads','Per-gene evidence/index rows copied as leads for traceability','mixed database annotations; not independent verification','D3_全基因组候选功能发现/D3-1-R_全基因跨数据库交叉引用.tsv; 注释证据表.tsv','Aggregated database cross-references may share an annotation lineage; source records need primary-entry recheck'])
ev.append(['G59-SRC-2016','','','iMC507 original supplementary Table 1','mmc1.xlsx / Table 1','reaction, score, GPR, media bounds','Original rows are copied only where an exact B1 AFE mapping could be linked to a gene reaction association','original model evidence, not necessarily gene-specific experimental evidence','A0_复现2016与2024模型/A0_主任务/脚本与环境/中间/mmc1.xlsx','Confidence score applies to reaction as reported in 2016; this task did not reassign it to new chemistry'])
for g in seg:
 for r in g['ann']:
  ev.append([f"G59-D3-ANN-{g['idx']}-{r.get('数据库','')}-{r.get('记录ID','')}",g['idx'],g['locus'],r.get('数据库',''),r.get('记录ID',''),r.get('证据字段',''),r.get('记录内容',''),r.get('证据等级/类型',''),r.get('来源文件',''),'D3 aggregate; accession and original database page need reinspection; automated annotation is not a wet experiment'])
 for r in g['cross']:
  ev.append([f"G59-D3-XREF-{g['idx']}-{r.get('数据库','')}-{r.get('记录ID','')}",g['idx'],g['locus'],r.get('数据库',''),r.get('记录ID',''),r.get('关系类型',''),r.get('对象/交叉引用',''),r.get('证据类型',''),r.get('来源文件',''),'D3 cross-reference only; repeated provenance is one source chain'])
# unresolved list
p=w.create_sheet('Open issues'); p.append(['Total no.','Current locus','Topic','Status / next verification'])
for g in seg:
 p.append([g['idx'],g['locus'],'protein family, enzyme/substrate, cofactors, direction, compartment, complex/isoenzyme, strain-specific paper, 2016 coverage','未确认：需按总说明对UniProt/InterPro/Pfam、KEGG/KO/EC、Rhea/ChEBI、BioCyc/MetaCyc、BRENDA/TCDB原始记录与同株论文逐项回查；D3汇总仅作检索入口。'])
# 2016 Table 1 linked rows exact raw evidence for review
rsh=w.create_sheet('2016 linked rows'); rsh.append(h)
for g in seg:
 for r in g['rx16']: rsh.append([r.get(k) for k in h])
# QA/limits
q=w.create_sheet('Readme'); q.append(['Item','Result'])
for row in [
 ('Scope','G59 subassignment total sequence numbers 2345–2442 inclusive'),('Count',len(seg)),('First / last',f"{seg[0]['locus']} / {seg[-1]['locus']}"),('Genome baseline','GCF_049532655.1, NZ_CP136162.1; local GFF gene order and CDS join'),('Mapping rule','B1 exact current locus + exact protein_id + mapping status 完全匹配; ambiguous paralogs excluded'),('Evidence boundary','D3 tables are pre-existing aggregated scans, used as lookup trails; no claim they were independently verified during this pass'),('Reaction score','不评分 for unconfirmed reactions; original score reproduced only for linked iMC507 rows'),('Completion status','该区段 98 行身份索引已生成；多数独立数据库原条目/论文证据链仍未完成，行动按证据不足暂缓，不得据此合并为已验收完整审查'),('2016 full fields','Where linked rows exist, raw Table 1 fields copied; formulas, bounds and confidence score remain original')]: q.append(row)
for sh in w:
 sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions
 for col in sh.columns:
  letter=col[0].column_letter
  maxlen=min(55,max((len(str(c.value)) if c.value is not None else 0 for c in col),default=10))
  sh.column_dimensions[letter].width=max(12,maxlen+2)
w.save(out)
print(json.dumps({'path':str(out.resolve()),'n':len(seg),'first':seg[0]['locus'],'last':seg[-1]['locus'],'workbook_sheets':w.sheetnames},ensure_ascii=True))
