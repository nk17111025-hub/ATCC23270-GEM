import json,time,concurrent.futures,datetime
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from pathlib import Path
base=Path('outputs/G53')
genes=json.loads((base/'g53_data.json').read_text(encoding='utf-8'))['index']
terms=sorted({r[2] for r in genes if r[7]=='未确认'})
def scan(term):
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode({'query':term,'format':'json','pageSize':20,'resultType':'lite'})
    for attempt in range(2):
        try:
            with urlopen(Request(url,headers={'User-Agent':'G53-evidence-review/1.0'}),timeout=30) as resp:d=json.load(resp)
            return term,{'hitCount':d.get('hitCount',0),'results':[{'id':x.get('id'),'title':x.get('title'),'doi':x.get('doi'),'pubYear':x.get('pubYear'),'source':x.get('source')} for x in d.get('resultList',{}).get('result',[])],'url':url}
        except Exception as e:
            if attempt:return term,{'error':repr(e),'url':url}
            time.sleep(0.5)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    res=dict(pool.map(scan,terms))
(base/'epmc_wp_scan.json').write_text(json.dumps({'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queries':res},ensure_ascii=False),encoding='utf-8')
print('done',len(terms),'hits',sum(v.get('hitCount',0)>0 for v in res.values()),'errors',sum('error' in v for v in res.values()))
