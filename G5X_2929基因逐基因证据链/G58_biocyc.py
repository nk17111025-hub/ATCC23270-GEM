import csv
import json
from collections import defaultdict
from pathlib import Path

root=Path(__file__).resolve().parents[1]
out=Path(__file__).resolve().parent/"G58_审查输出"
data=json.loads((out/"G58_中间数据.json").read_text(encoding="utf-8"))
oldgff=root/"B1_最新基因组与旧AFE编号映射/标准化数据/GCF_000021485.1/genomic.gff"
relations=root/"C2_外部数据库与注释来源核验/C2-1_BioCyc与数据库版本核验/C2-1_基因_蛋白_反应_通路.tsv"
metacyc=root/"D3_全基因组候选功能发现/D3-1_初版候选筛选/MetaCyc/MetaCyc_BioCyc_ID精确映射.tsv"

afe2rs=defaultdict(set)
with oldgff.open(encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        x=line.split("\t")
        if len(x)<9 or x[2]!="gene": continue
        a={k:v for z in x[8].rstrip().split(";") if "=" in z for k,v in [z.split("=",1)]}
        for old in a.get("old_locus_tag","").split(","):
            if old.startswith("AFE_") and a.get("locus_tag","").startswith("AFE_RS"):
                afe2rs[old].add(a["locus_tag"])

gene2protein=defaultdict(set)
reaction2enzyme=defaultdict(set)
enzyme2enzrxn=defaultdict(set)
enzrxn2reaction=defaultdict(set)
reaction2pathway=defaultdict(set)
metabyid={}
with metacyc.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        metabyid[row["MetaCyc/Pathway Tools ID"]]=row
with relations.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        t,s,o=row["关系类型"],row["主体ID"],row["客体ID"]
        if t=="gene_to_protein": gene2protein[s].add(o)
        elif t=="reaction_to_enzyme": reaction2enzyme[o].add(s)
        elif t=="enzymatic_reaction_to_enzyme": enzyme2enzrxn[o].add(s)
        elif t=="enzymatic_reaction_to_reaction": enzrxn2reaction[s].add(o)
        elif t=="reaction_to_pathway": reaction2pathway[s].add(o)

rows=[]
for g in data["index"]:
    serial,locus,wp,afe=g[0],g[1],g[2],g[5]
    if not afe: continue
    for rs in sorted(afe2rs.get(afe,[])):
        for protein in sorted(gene2protein.get(rs,[])):
            reaction_ids=set(reaction2enzyme.get(protein,[]))
            erxs=enzyme2enzrxn.get(protein,set())
            for erx in erxs: reaction_ids.update(enzrxn2reaction.get(erx,set()))
            for rid in sorted(reaction_ids):
                matches=sorted(erx for erx in erxs if rid in enzrxn2reaction.get(erx,set()))
                meta=metabyid.get(rid,{})
                rows.append([serial,locus,wp,afe,rs,protein,protein,";".join(matches),rid,
                             ";".join(sorted(reaction2pathway.get(rid,[]))),
                             meta.get("BioCyc equation participants",""),meta.get("direction",""),
                             meta.get("EC",""),meta.get("BioCyc evidence code",""),
                             meta.get("Rhea ID exact cross-reference",""),meta.get("MetaCyc direct retrieval",""),
                             "BioCyc GCF_000021485 v30.0 Tier 3; legacy genome; computational gene association",
                             f"https://biocyc.org/GCF_000021485/NEW-IMAGE?type=REACTION&object={rid}"])
(out/"BioCyc_G58_证据链.json").write_text(json.dumps(rows,ensure_ascii=False),encoding="utf-8")
print("genes with chain",len({r[1] for r in rows}),"gene-reaction rows",len(rows),"unique reaction",len({r[8] for r in rows}))
