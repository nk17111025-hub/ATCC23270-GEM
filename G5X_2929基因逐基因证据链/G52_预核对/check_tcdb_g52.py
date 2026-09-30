import json
import re
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

DIR=Path(__file__).resolve().parent
ROOT=Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1')
data=json.loads((DIR/'g52_precheck.json').read_text(encoding='utf-8'))
decisions=json.loads((DIR/'g52_gene_decisions_draft.json').read_text(encoding='utf-8'))
byloc={r['当前locus']:r for r in data['records']}
selected={r['当前locus'] for r in decisions if r['决策类别']=='转运底物未定'}
selected.update({'RU820_RS01835','RU820_RS02640','RU820_RS02790','RU820_RS01515','RU820_RS01680'})
wp={byloc[l]['WP']:l for l in selected}
query={}
for record in SeqIO.parse(ROOT/'protein.faa','fasta'):
    if record.id in wp:query[wp[record.id]]=str(record.seq).replace('*','')
assert len(query)==len(selected),(len(query),len(selected),selected-set(query))

tcpath=DIR/'tcdb_2026_09_25.faa'
if not tcpath.exists():
    req=urllib.request.Request('https://www.tcdb.org/public/tcdb',headers={'User-Agent':'G52-evidence-review/1.0'})
    with urllib.request.urlopen(req,timeout=90) as response:
        tcpath.write_bytes(response.read())

def kmers(seq,k=5):
    return {seq[i:i+k] for i in range(len(seq)-k+1) if 'X' not in seq[i:i+k]}

refs=[]
for record in SeqIO.parse(tcpath,'fasta'):
    seq=str(record.seq).replace('*','').upper()
    if len(seq)>=70:
        refs.append((record.id,record.description,seq,kmers(seq)))

aligner=PairwiseAligner(mode='local')
aligner.substitution_matrix=substitution_matrices.load('BLOSUM62')
aligner.open_gap_score=-10
aligner.extend_gap_score=-0.5

def alignment_metrics(q,t):
    a=aligner.align(q,t)[0]
    qa,ta=a.aligned
    aligned=0;match=0
    for (q0,q1),(t0,t1) in zip(qa,ta):
        n=min(q1-q0,t1-t0)
        aligned+=n
        match+=sum(x==y for x,y in zip(q[q0:q0+n],t[t0:t0+n]))
    return {'identity':float(round(match/aligned,4)) if aligned else 0.0,'query_coverage':float(round(aligned/len(q),4)),'target_coverage':float(round(aligned/len(t),4)),'aligned_residues':int(aligned),'score':float(round(a.score,1))}

results={}
for loc,seq in query.items():
    qk=kmers(seq)
    pre=[]
    for rid,description,target,tk in refs:
        overlap=len(qk & tk)
        if overlap<4:continue
        jaccard=overlap/(len(qk)+len(tk)-overlap)
        pre.append((jaccard,rid,description,target))
    pre.sort(reverse=True)
    hits=[]
    for pre_score,rid,description,target in pre[:12]:
        m=alignment_metrics(seq,target)
        tcid=(rid.rsplit('|',1)[-1] if '|' in rid else '')
        hits.append({'tcdb_id':rid,'tcid':tcid,'description':description,'prefilter_jaccard':round(pre_score,4),**m})
    hits.sort(key=lambda x:(x['identity']*x['query_coverage'],x['score']),reverse=True)
    results[loc]={'wp':byloc[loc]['WP'],'product':byloc[loc]['当前product'],'query_length':len(seq),'hits':hits[:5]}

out={'source':'https://www.tcdb.org/public/tcdb','accessed_utc':datetime.now(timezone.utc).isoformat(),'reference_count':len(refs),'by_locus':results}
(DIR/'g52_tcdb_homology.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'query_loci':len(query),'references':len(refs),'strong_hits':sum(bool(v['hits']) and v['hits'][0]['identity']>=.3 and v['hits'][0]['query_coverage']>=.6 for v in results.values())}))
for loc,v in sorted(results.items()):
    h=v['hits'][0] if v['hits'] else None
    print(loc,v['product'][:40],h['tcid'] if h else '',h['identity'] if h else '',h['query_coverage'] if h else '',h['description'][:95] if h else '')
