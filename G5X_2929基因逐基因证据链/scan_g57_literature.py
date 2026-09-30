"""Find candidate papers mentioning each verified historical AFE or current WP.

Hits are search leads only; they never count as functional validation.
"""
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path(__file__).resolve().parent / "G57_审查输出"
rows = list(csv.DictReader((BASE / "G57_逐基因索引与状态.tsv").open(encoding="utf-8-sig"), delimiter="\t"))
cache_path = BASE / "EuropePMC_gene_scan_20260927.json"
cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}

def scan(r):
    key = r["旧 AFE"] if r["旧 AFE"] != "未确认" else r["WP"]
    if key in cache:
        return key, cache[key]
    params = urllib.parse.urlencode({"query":key,"format":"json","pageSize":20,"resultType":"core"})
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"Codex-G57-review/1.0"})
            with urllib.request.urlopen(req, timeout=30) as f: data = json.load(f)
            result = {"query":key,"url":url,"hitCount":data.get("hitCount",0),"papers":[{"title":x.get("title", ""),"doi":x.get("doi", ""),"pmid":x.get("pmid", ""),"pmcid":x.get("pmcid", ""),"year":x.get("pubYear", ""),"isOpenAccess":x.get("isOpenAccess", "")} for x in data.get("resultList",{}).get("result",[])]}
            return key, result
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if attempt == 2:return key,{"query":key,"url":url,"error":repr(e)}
            time.sleep(1.5*(attempt+1))

todo = [r for r in rows if (r["旧 AFE"] if r["旧 AFE"] != "未确认" else r["WP"]) not in cache]
with ThreadPoolExecutor(max_workers=5) as pool:
    futures = [pool.submit(scan,r) for r in todo]
    for i, future in enumerate(as_completed(futures),1):
        key, result = future.result()
        cache[key] = result
        if i % 20 == 0 or i == len(futures):
            cache_path.write_text(json.dumps(cache,ensure_ascii=False,indent=2),encoding="utf-8")
            print(f"{i}/{len(futures)}", flush=True)

out=[]
for r in rows:
    key = r["旧 AFE"] if r["旧 AFE"] != "未确认" else r["WP"]
    c = cache[key]
    out.append({"总序号":r["总序号"],"RU820":r["RU820"],"WP":r["WP"],"查询ID":key,"命中数":c.get("hitCount", ""),"DOI候选":";".join(p["doi"] for p in c.get("papers",[]) if p["doi"]),"PMID候选":";".join(p["pmid"] for p in c.get("papers",[]) if p["pmid"]),"查询URL":c["url"],"状态":"检索失败" if "error" in c else "候选待逐篇核实；零命中不证明无论文"})
with (BASE / "G57_逐基因文献检索线索.tsv").open("w",encoding="utf-8-sig",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
    w.writeheader();w.writerows(out)
print({"genes":len(out),"errors":sum(x["状态"]=="检索失败" for x in out),"with_hits":sum(isinstance(x["命中数"],int) and x["命中数"]>0 for x in out)})
