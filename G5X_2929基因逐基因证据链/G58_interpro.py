import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

out = Path(__file__).resolve().parent / "G58_审查输出"
raw = json.loads((out / "UniProt_WP_查询原始.json").read_text(encoding="utf-8"))
cache_file = out / "InterPro_Pfam_直接查询.json"
cache = json.loads(cache_file.read_text(encoding="utf-8")) if cache_file.exists() else {}
todo = [(wp,items[0]["primaryAccession"]) for wp,items in raw.items() if items]

def fetch(item):
    wp,acc = item
    data={}
    for db in ("interpro","pfam"):
        url=f"https://www.ebi.ac.uk/interpro/api/entry/{db}/protein/uniprot/{acc}/?page_size=100"
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"G58 evidence audit/1.0"})
            with urllib.request.urlopen(req,timeout=25) as response:
                obj={} if response.status==204 else json.load(response)
            data[db]={"url":url,"count":obj.get("count"),"entries":[
                {"id":r["metadata"].get("accession"),"name":r["metadata"].get("name")}
                for r in obj.get("results",[])]}
        except Exception as e:
            data[db]={"url":url,"error":str(e)}
    return wp,acc,data

with ThreadPoolExecutor(max_workers=3) as pool:
    futures=[pool.submit(fetch,item) for item in todo if item[0] not in cache or any("error" in cache[item[0]].get(db,{}) for db in ("interpro","pfam"))]
    for n,future in enumerate(as_completed(futures),1):
        wp,acc,data=future.result()
        cache[wp]={"accession":acc,**data}
        if n%20==0 or n==len(futures):
            cache_file.write_text(json.dumps(cache,ensure_ascii=False),encoding="utf-8")
            print("InterPro/Pfam",n,"/",len(futures),flush=True)
print("done",len(cache))
