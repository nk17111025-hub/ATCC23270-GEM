import csv, re, os
from pathlib import Path
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
root=Path.cwd()
out=root/'G5X_2929基因逐基因证据链'/'G59_work'/'part_b.xlsx'
out.parent.mkdir(parents=True,exist_ok=True)
# canonical list from local coordinate-sorted GFF and exact G59 sequence positions
GFF=root/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'GCF_049532655.1'/'genomic.gff'
features=[]
with GFF.open(encoding='utf-8-sig') as f:
 for line in f:
  if line.startswith('#'): continue
  a=line.rstrip('\n').split('\t')
  if len(a)<9 or a[2]!='gene' or 'gene_biotype=protein_coding' not in a[8]: continue
  attrs=dict(re.findall(r'([^;=]+)=([^;]*)',a[8]))
  features.append((a[0],int(a[3]),int(a[4]),a[6],attrs))
features.sort(key=lambda x:(x[0],x[1],x[2]))
# Number is G59 position from 2345; requested 2443-2540
subset=features[2442:2540]
assert len(subset)==98
# extract cds products / proteins by locus from same GFF
cds={}
with GFF.open(encoding='utf-8-sig') as f:
 for line in f:
  if line.startswith('#'): continue
  a=line.rstrip('\n').split('\t')
  if len(a)<9 or a[2]!='CDS': continue
  at=dict(re.findall(r'([^;=]+)=([^;]*)',a[8]))
  locus=at.get('locus_tag')
  if locus: cds[locus]=(at.get('protein_id',''),at.get('product',''))
# TSV reader accommodating literal field quoting
def read_tsv(path):
 with open(path,encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
mp=read_tsv(root/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'旧新基因映射.tsv')
map_by={r.get('当前GCF049位点',''):r for r in mp}
d3dir=root/'D3_全基因组候选功能发现'/'D3-1-R_全基因组重扫与修正'
xref=read_tsv(d3dir/'D3-1-R_全基因跨数据库交叉引用.tsv')
xby={r.get('当前locus',''):r for r in xref}
ann=read_tsv(d3dir/'D3-1-R_注释证据表.tsv')
annby={}
for r in ann: annby.setdefault(r.get('当前locus',''),[]).append(r)
# model base
model_path=root/'A0_复现2016与2024模型'/'A0_主任务'/'脚本与环境'/'中间'/'mmc1.xlsx'
model_byAFE={}
if model_path.exists():
 wb0=load_workbook(model_path,read_only=True,data_only=True)
 for sh in wb0.worksheets:
  for row in sh.iter_rows(values_only=True):
   vals=[str(v) if v is not None else '' for v in row]
   for afe in re.findall(r'AFE_\d+', ' '.join(vals)):
    model_byAFE.setdefault(afe,[]).append(' | '.join(v for v in vals[:16] if v))
wb=Workbook(); idx=wb.active; idx.title='Gene index'
headers=['总序号','当前 RU820 locus','当前 WP protein_id','染色体坐标','当前 NCBI product','旧 AFE locus','AFE 映射依据','映射状态','KEGG KO','KEGG EC','KEGG reaction','BioCyc reaction ID','Rhea 精确交叉引用','2016 GPR反应线索','是否找到精确 reaction','动作','关联行 ID','审查状态','未决点','身份来源','数据库交叉引用来源']
idx.append(headers)
alpha=wb.create_sheet('Alpha structure review')
alpha.append(['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe2+ lb','Fe2+ ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub','Record class','Candidate ID','Candidate group','2016 linked Reaction ID','2026 reaction object','Legacy AFE locus','Current locus','Primary evidence DOI','Evidence type','Evidence boundary','2026 action','2026 reaction score','Source ID','Score scope'])
sources=wb.create_sheet('Sources'); sources.append(['Source ID','Full citation / record','Year','DOI','Source type','Strain','Experimental type','Genes / proteins studied','Related Reaction IDs / candidate IDs','What was directly measured / recorded','Applicable confidence level','Evidence limitation'])
source_rows=[
 ['G59B-S001','NCBI RefSeq GCF_049532655.1, AFEATCC23270_v5.2 genomic.gff and cds_from_genomic.fna; local snapshot per G5X frozen baseline',2026,'','Reference genome/CDS','A. ferrooxidans ATCC 23270','Genome annotation','RU820_RS12885–RU820_RS13370','','Locus, coordinates, strand, protein_id, product in same assembly','Not applicable','Annotation identity only; product names do not establish biochemical activity.'],
 ['G59B-S002','B1 standardized old/new gene mapping.tsv; local snapshot read 2026-09-24',2026,'','Curated mapping dataset','A. ferrooxidans ATCC 23270','Exact RefSeq protein ID mapping','Requested 98 loci','','Old AFE locus and mapping status/method where exact protein_id is matched','Not applicable','Mapping source is a local project crosswalk; sequence/neighborhood evidence not independently realigned in this subtask.'],
 ['G59B-S003','D3-1-R whole-genome cross-reference and annotation evidence TSVs; local snapshot read 2026-09-24',2026,'','Cross-database annotation aggregation','A. ferrooxidans ATCC 23270','Database cross-reference inspection','Requested 98 loci','','Existing KEGG/BioCyc/Rhea/2016 pointers as recorded in D3 tables','Not applicable','Inherited database links are leads, not independent experimental validation; source accession/version rows must be revisited for formal edits.'],
 ['G59B-S004','iMC507 supplementary model workbook mmc1.xlsx (converted read-only local copy); sheet cells scanned for AFE identifiers',2016,'','Model supplementary table','A. ferrooxidans ATCC 23270','Model/GPR comparison','Requested AFE loci','','Model row text containing an exact legacy locus when present','As published in workbook; exact reaction score requires row-level review','Only an identifier-text scan was performed here; not a full audit of all 2016 reaction fields or culture bounds.']]
for row in source_rows: sources.append(row)
# GFF sequence is region-based; build row records
for i,(seq,start,end,strand,at) in enumerate(subset,2443):
 locus=at.get('locus_tag',''); pid,product=cds.get(locus,('',''))
 pos=f'{seq}:{start}-{end}({strand})'
 mr=map_by.get(locus,{})
 xr=xby.get(locus,{})
 afe=mr.get('旧AFE位点','') or xr.get('旧AFE locus','未确认')
 if not afe: afe='未确认'
 method=mr.get('映射方法','') or xr.get('B1映射方法','')
 status=mr.get('映射状态','') or xr.get('B1映射状态','') or '未确认'
 row2016='; '.join(model_byAFE.get(afe,[])) if afe in model_byAFE else '未在2016补表文本扫描中找到明确AFE标识；反应层复核未完成'
 # current D3 records are cross-reference leads. No reaction assumed from names.
 kegg='; '.join(filter(None,[xr.get('KEGG KO',''),xr.get('KEGG EC',''),xr.get('KEGG reaction','')]))
 bio=xr.get('BioCyc reaction ID','')
 rhe=xr.get('Rhea精确交叉引用','')
 exact=bool(kegg or bio or rhe)
 action='证据不足暂缓'
 status_review='身份与本地交叉引用已核；独立原始条目、结构域、论文及2016反应字段未完成'
 unresolved='需回查UniProt/InterPro/Pfam、KEGG/Rhea/BioCyc/BRENDA/TCDB原始条目，定位同株论文，并按2016补表反应行核对化学式/分数/GPR/培养边界；底物产物/辅因子/方向/区室/复合体未确认'
 if exact:
  unresolved='存在数据库精确交叉引用线索，但链式来源与反应化学未独立复核；需核对底物产物/辅因子/方向/区室/复合体及2016对应反应行'
 rid=f'G59B-{i:04d}'
 idx.append([i,locus,pid,pos,product,afe,method or '未核验',status,kegg,bio,rhe,row2016,'是' if exact else '未找到精确链接','未确认：仅扫描AFE标识' if 'AFE_' in row2016 else '未确认：未检出明确AFE文本','数据库有精确链接' if exact else '未找到','证据不足暂缓',rid,status_review,unresolved,'GCF_049532655.1 genomic.gff/CDS; G5X 2026-09-24','B1 old/new mapping.tsv; D3 cross-reference/annotation TSV 2026-09-24'])
 # Reaction-level template: gene audit record, explicitly no formal chemistry claim
 alpha.append(['','','','','','','','','','','','','','','','','gene audit - incomplete',rid,'G59',row2016 if 'AFE_' in row2016 else '', '',afe,locus,'','local database cross-reference inspection','只确认当前基因身份及本地数据库登记线索；本基因精确生化反应、同株直接证据、反应方向/区室和化学式尚未核实。','证据不足暂缓','不评分','G59B-S001;G59B-S002;G59B-S003;G59B-S004','未定义精确反应；不评分'])
# Source tracking sheet
notes=wb.create_sheet('Open items')
notes.append(['总序号','RU820 locus','优先未决项'])
for row in idx.iter_rows(min_row=2,values_only=True): notes.append([row[0],row[1],row[18]])
# basic style
for ws in wb.worksheets:
 ws.freeze_panes='A2'; ws.auto_filter.ref=ws.dimensions
 for c in ws[1]: c.font=Font(bold=True,color='FFFFFF'); c.fill=PatternFill('solid',fgColor='1F4E78'); c.alignment=Alignment(wrap_text=True,vertical='top')
 for col in ws.columns:
  letter=col[0].column_letter
  ws.column_dimensions[letter].width=min(52,max(14,max((len(str(c.value or '')) for c in col[:30]),default=14)*0.9))
 for row in ws.iter_rows(min_row=2):
  for c in row: c.alignment=Alignment(wrap_text=True,vertical='top')
wb.save(out)
print(out)
print('count',len(subset),'unique loci',len({x[4].get('locus_tag') for x in subset}),'first',subset[0][4].get('locus_tag'),'last',subset[-1][4].get('locus_tag'))
