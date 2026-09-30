from pathlib import Path
from difflib import SequenceMatcher

ROOT=Path(__file__).resolve().parents[2]
paths=list((ROOT/'01_原始数据/01_基因组与官方注释').glob('GCF_*_20260922/解压/ncbi_dataset/data/GCF_*/protein.faa'))
ids={'WP_012536195.1','WP_225487408.1','WP_012606555.1','WP_012536130.1'}
seq={}
for path in paths:
    key=None
    for line in path.open(encoding='utf-8'):
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in ids:
                seq[(path.parent.name,key)]=''
        elif key in ids:
            seq[(path.parent.name,key)]+=line.strip()
for key,value in seq.items():
    print(key,len(value))
a=seq[('GCF_000021485.1','WP_012536195.1')]
b=seq[('GCF_049532655.1','WP_225487408.1')]
match=SequenceMatcher(None,a,b,autojunk=False)
blocks=[(m.a,m.b,m.size) for m in match.get_matching_blocks() if m.size]
print('AFE_0545 old-new',{'old_aa':len(a),'new_aa':len(b),'sequence_match_ratio':match.ratio(),'matching_blocks':blocks[:10],'matched_residues':sum(x[2] for x in blocks),'shorter_coverage':sum(x[2] for x in blocks)/min(len(a),len(b))})
