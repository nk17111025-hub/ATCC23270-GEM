import urllib.request,urllib.parse,csv
from openpyxl import load_workbook
p='G5X_2929基因逐基因证据链/G59_work/part_b.xlsx'
w=load_workbook(p); idx=w['Gene index']; genes=list(idx.iter_rows(min_row=2,values_only=True)); ids=[r[2] for r in genes]
query=' OR '.join('xref:'+x for x in ids)
url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode({'query':query,'format':'tsv','fields':'accession,reviewed,id,protein_name,gene_names,organism_name,xref_refseq,xref_interpro,xref_pfam'})
req=urllib.request.Request(url,headers={'User-Agent':'ATCC23270-evidence-audit/1.0'})
data=urllib.request.urlopen(req,timeout=60).read().decode('utf-8')
rows=list(csv.DictReader(data.splitlines(),delimiter='\t'))
print('records',len(rows))
by={}
for r in rows:
 for pid in r.get('RefSeq','').split(';'):
  # refseq crossref sometimes WP with version/no version
  by.setdefault(pid.strip(),[]).append(r)
s=w.create_sheet('UniProt InterPro Pfam')
s.append(['Current WP','UniProt accession','review status','Entry name','protein name','gene names','organism/strain','RefSeq xref','InterPro accessions','Pfam accessions','query URL','Interpretation boundary'])
for row in genes:
 pid=row[2]; matches=by.get(pid,[])
 if not matches: matches=by.get(pid.split('.')[0],[])
 if not matches:
  s.append([pid,'未返回','未查到UniProt交叉引用','','','','','','','','https://rest.uniprot.org/uniprotkb/search?query='+urllib.parse.quote('xref:'+pid),'本次REST查询未回传关联记录；UniProt缺失/索引规则所致可能性未区分。'])
 else:
  # ATCC/DSM entry prioritized if exact; otherwise all results retained
  matches=sorted(matches,key=lambda r:('ATCC 23270' not in r.get('Organism',''),r.get('Reviewed','')!='reviewed'))
  for r in matches:
   s.append([pid,r.get('Entry'),r.get('Reviewed'),r.get('Entry Name'),r.get('Protein names'),r.get('Gene Names'),r.get('Organism'),r.get('RefSeq'),r.get('InterPro'),r.get('Pfam'),url,'UniProt注释交叉引用；其中InterPro/Pfam为数据库签名注释，未视为本株直接酶学证明。'])
# update summary status conservatively, clean duplicated 2016 strings
for i,row in enumerate(genes,2):
 pid=row[2]
 # retained annotation is explicit in UniProt sheet; only action remains based on row 2016 field correction below
 v=idx.cell(i,14).value or ''
 idx.cell(i,14).value='; '.join(dict.fromkeys(v.split('; ')))
w.save(p)
print('targeted WP',len(ids),'uniprot rows',s.max_row-1,'file_bytes',__import__('os').path.getsize(p))

