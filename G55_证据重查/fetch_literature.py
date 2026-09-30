import concurrent.futures
import json
import subprocess
import urllib.parse
from pathlib import Path

HERE=Path(__file__).parent
genes=json.loads((HERE/'g55_data.json').read_text(encoding='utf-8'))
terms=sorted({a for z in genes for a in z['afe'].split(';') if a.startswith('AFE_')})

def query(term):
    q=urllib.parse.urlencode({'query':term,'format':'json','pageSize':'25','resultType':'lite'})
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+q
    p=subprocess.run(['curl.exe','-f','-L','--silent','--show-error','--max-time','30',url],capture_output=True)
    if p.returncode:
        return term,{'error':p.stderr.decode(errors='replace')[:200],'url':url}
    try:
        x=json.loads(p.stdout)
        return term,{'hit_count':x.get('hitCount',0),'results':[{k:r.get(k) for k in ['id','source','pmid','pmcid','doi','title','authorString','journalTitle','pubYear','inPMC']} for r in x.get('resultList',{}).get('result',[])],'url':url}
    except Exception as e:
        return term,{'error':str(e),'url':url}

out={}
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    for i,(term,result) in enumerate(pool.map(query,terms),1):
        out[term]=result
        if i%50==0: print('done',i,'/',len(terms),flush=True)
(HERE/'literature_hits.json').write_text(json.dumps({'date':'2026-09-24','query_mode':'Europe PMC basic gene-locus full-text search','results':out},ensure_ascii=False,indent=2),encoding='utf-8')
print('terms',len(terms),'hit_terms',sum(v.get('hit_count',0)>0 for v in out.values()),'errors',sum('error' in v for v in out.values()))
