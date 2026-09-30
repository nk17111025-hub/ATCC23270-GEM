import json,time,datetime
from urllib.request import Request,urlopen
from pathlib import Path
p=Path('outputs/G53/epmc_afe_scan.json')
d=json.loads(p.read_text(encoding='utf-8'))
bad=[k for k,v in d['queries'].items() if 'error' in v][:5]
for i,term in enumerate(bad,1):
    url=d['queries'][term]['url']
    try:
        with urlopen(Request(url,headers={'User-Agent':'G53-evidence-review/1.0'}),timeout=40) as resp:ans=json.load(resp)
        d['queries'][term]={'hitCount':ans.get('hitCount',0),'results':[{'id':x.get('id'),'title':x.get('title'),'doi':x.get('doi'),'pubYear':x.get('pubYear'),'source':x.get('source')} for x in ans.get('resultList',{}).get('result',[])],'url':url}
        print(i,term,'OK',ans.get('hitCount',0),flush=True)
    except Exception as e:print(i,term,'ERROR',repr(e),flush=True)
    p.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8')
    time.sleep(0.5)
d['retried_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
p.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8')
print('errors_remaining',sum('error' in v for v in d['queries'].values()))
