import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

out=Path(__file__).resolve().parent/"G58_审查输出"
payload=json.loads((out/"G58_中间数据.json").read_text(encoding="utf-8"))
rawfile=out/"TCDB_pfam_20260927.tsv"
if not rawfile.exists():
    req=urllib.request.Request("https://www.tcdb.org/public/pfam.tsv",headers={"User-Agent":"G58 evidence audit/1.0"})
    with urllib.request.urlopen(req,timeout=30) as r: rawfile.write_bytes(r.read())
by_pfam=defaultdict(set)
for line in rawfile.read_text(encoding="utf-8",errors="replace").splitlines():
    x=line.split("\t")
    if len(x)>=3:
        by_pfam[x[0]].add((x[1],x[2].strip()))
rows=[]
for g in payload["index"]:
    if not re.search(r"transport|permease|porin|efflux|antiporter|exporter",g[4],re.I):
        continue
    pfams=[x for x in g[27].split(";") if x]
    hits=sorted({(tc.split(".")[0]+"."+tc.split(".")[1]+"."+tc.split(".")[2],name) for pf in pfams for tc,name in by_pfam.get(pf,[]) if len(tc.split("."))>=3})
    rows.append([g[1],g[2],g[4],";".join(pfams),len(hits),";".join(f"{tc} {name}" for tc,name in hits[:15]),
                 "Pfam→TCDB family mapping only; no WP-specific transporter assignment",
                 "https://www.tcdb.org/public/pfam.tsv"])
(out/"TCDB_Pfam_线索.json").write_text(json.dumps(rows,ensure_ascii=False),encoding="utf-8")
print("transport rows",len(rows),"with family hints",sum(bool(r[4]) for r in rows))
