"""Search both current identifiers; hits remain leads, never functional proof."""
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

base = Path(__file__).resolve().parent / 'G57_审查输出'
with (base / 'G57_逐基因索引与状态.tsv').open(encoding='utf-8-sig', newline='') as f:
    genes = list(csv.DictReader(f, delimiter='\t'))
cache_path = base / 'EuropePMC_current_ids_20260927.json'
cache = json.loads(cache_path.read_text(encoding='utf-8')) if cache_path.exists() else {}

def scan(key):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query':key,'format':'json','pageSize':20,'resultType':'core'})
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent':'Codex-G57-review/1.0'})
            with urllib.request.urlopen(req, timeout=30) as f:
                data = json.load(f)
            papers = [{'title':r.get('title',''),'doi':r.get('doi',''),'pmid':r.get('pmid',''),'year':r.get('pubYear','')} for r in data.get('resultList',{}).get('result',[])]
            return key, {'query':key,'url':url,'hitCount':data.get('hitCount',0),'papers':papers}
        except (urllib.error.URLError,TimeoutError,json.JSONDecodeError) as e:
            if attempt == 2:
                return key, {'query':key,'url':url,'error':repr(e)}
            time.sleep(1.5*(attempt+1))

keys = list(dict.fromkeys(k for g in genes for k in (g['RU820'],g['WP'])))
todo = [k for k in keys if k not in cache]
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = [pool.submit(scan,k) for k in todo]
    for i, future in enumerate(as_completed(futures),1):
        k, value = future.result()
        cache[k] = value
        if i % 40 == 0 or i == len(todo):
            cache_path.write_text(json.dumps(cache,ensure_ascii=False,indent=2),encoding='utf-8')
            print(f'{i}/{len(todo)}', flush=True)

rows = []
for g in genes:
    ru, wp = cache[g['RU820']],cache[g['WP']]
    rows.append({'总序号':g['总序号'],'RU820':g['RU820'],'WP':g['WP'],'RU820命中':ru.get('hitCount',''),'WP命中':wp.get('hitCount',''),'RU820 DOI':';'.join(p['doi'] for p in ru.get('papers',[]) if p['doi']),'WP DOI':';'.join(p['doi'] for p in wp.get('papers',[]) if p['doi']),'RU820查询URL':ru['url'],'WP查询URL':wp['url'],'检索状态':'失败' if 'error' in ru or 'error' in wp else '完成；零命中不证明无论文'})
with (base / 'G57_当前编号双路文献检索.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    writer = csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
print({'genes':len(genes),'queries':len(keys),'errors':sum('error' in cache[k] for k in keys),'with_hits':sum(int(cache[k].get('hitCount',0))>0 for k in keys)})
