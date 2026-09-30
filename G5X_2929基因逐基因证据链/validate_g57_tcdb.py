"""Screen G57 transport-related WPs against TCDB family sequences.

Family assignment is homology only. Substrate, coupling, direction and GPR are never inferred.
"""
import csv
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

root = Path(__file__).resolve().parents[1]
out = root / 'G5X_2929基因逐基因证据链/G57_审查输出'
tcdb_path = root / 'G5X_2929基因逐基因证据链/G52_预核对/tcdb_2026_09_25.faa'
protein_path = root / '01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/protein.faa'
with (out / 'G57_逐基因索引与状态.tsv').open(encoding='utf-8-sig',newline='') as f:
    genes = list(csv.DictReader(f,delimiter='\t'))

def prefixes(product):
    p = product.lower()
    if 'tonb-dependent receptor' in p:
        return ['1.B.14.']
    if 'carbohydrate porin' in p:
        return ['1.B.19.','1.B.23.']
    if p == 'porin':
        return ['1.B.']
    if 'tonb' in p and 'receptor' not in p:
        return ['1.A.30.','2.C.1.']
    if 'exbd' in p or 'exbb' in p or 'tolr' in p or 'tolq' in p or 'mota' in p:
        return ['1.A.30.']
    if 'rnd' in p or 'mmpl' in p or 'hlyd' in p:
        return ['2.A.6.']
    if 'eama' in p:
        return ['2.A.7.']
    if 'mfs' in p or 'sugar porter' in p:
        return ['2.A.1.']
    if 'nhaa' in p:
        return ['2.A.33.']
    if 'cation:proton antiporter' in p:
        return ['2.A.36.','2.A.37.']
    if 'kdp' in p:
        return ['3.A.3.7.']
    if 'p-type atpase' in p:
        return ['3.A.3.']
    if 'lpt' in p:
        return ['3.A.1.152.']
    if 'abc' in p:
        return ['3.A.1.']
    return []

selected = {g['WP']:(g,prefixes(g['当前 product'])) for g in genes if prefixes(g['当前 product'])}
sequences = {r.id:str(r.seq).replace('*','').upper() for r in SeqIO.parse(protein_path,'fasta') if r.id in selected}
assert set(sequences) == set(selected), (len(sequences),len(selected))

refs = []
for rec in SeqIO.parse(tcdb_path,'fasta'):
    tcid = rec.id.rsplit('|',1)[-1]
    seq = str(rec.seq).replace('*','').upper()
    if len(seq) >= 70:
        refs.append((tcid,rec.description,seq))

aligner = PairwiseAligner(mode='local')
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

def kmers(s,k=5):
    return {s[i:i+k] for i in range(len(s)-k+1) if 'X' not in s[i:i+k]}

def metrics(q,t):
    aln = aligner.align(q,t)[0]
    aligned = matches = 0
    for (q0,q1),(t0,t1) in zip(*aln.aligned):
        qpart,tpart = q[q0:q1],t[t0:t1]
        aligned += min(len(qpart),len(tpart))
        matches += sum(a == b for a,b in zip(qpart,tpart))
    return {'identity':round(matches/aligned,4) if aligned else 0,'query_coverage':round(aligned/len(q),4),'aligned_residues':aligned,'score':round(aln.score,1)}

rows=[]
for i,(wp,(g,expected)) in enumerate(selected.items(),1):
    q=sequences[wp]
    qk=kmers(q)
    candidates=[]
    for tcid,desc,t in refs:
        if any(tcid.startswith(prefix) for prefix in expected):
            tk=kmers(t)
            shared=len(qk & tk)
            if shared:
                candidates.append((shared,tcid,desc,t))
    candidates.sort(reverse=True)
    hits=[]
    for shared,tcid,desc,t in candidates[:20]:
        hits.append({'tcid':tcid,'description':desc,'shared_5mers':shared,**metrics(q,t)})
    hits.sort(key=lambda z:(z['identity']*z['query_coverage'],z['score']),reverse=True)
    best=hits[0] if hits else None
    supported=bool(best and best['identity']>=0.30 and best['query_coverage']>=0.60)
    rows.append({'总序号':g['总序号'],'RU820':g['RU820'],'WP':wp,'当前product':g['当前 product'],'查询TCDB大类':';'.join(expected),'候选序列数':len(candidates),'最佳TCDB ID':best['tcid'] if best else '', '最佳TCDB描述':best['description'] if best else '', '同一性':best['identity'] if best else '', '当前蛋白覆盖率':best['query_coverage'] if best else '', '比对残基数':best['aligned_residues'] if best else '', '家族同源支持':'是' if supported else '未达阈值','证据边界':'仅同源家族线索；不证明具体底物、运输方向、耦联、区室或旧 GPR'})
    if i % 10 == 0:
        print(f'{i}/{len(selected)}',flush=True)

with (out / 'G57_TCDB家族同源核验.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
(out / 'G57_TCDB家族同源核验说明.json').write_text(json.dumps({'source':'https://www.tcdb.org/public/tcdb','cached_reference':str(tcdb_path.relative_to(root)),'created_utc':datetime.now(timezone.utc).isoformat(),'selection':'当前 RefSeq product 预选 TCDB 大类；5-mer 选前20；BLOSUM62 局部比对','threshold':'identity >= 0.30 and query coverage >= 0.60','transport_properties_inferred':False,'queries':len(rows)},ensure_ascii=False,indent=2),encoding='utf-8')
print({'queries':len(rows),'supported':sum(r['家族同源支持']=='是' for r in rows)})
