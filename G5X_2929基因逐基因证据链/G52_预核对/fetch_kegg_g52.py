import json
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

DIR=Path(__file__).resolve().parent
records=json.load((DIR/'g52_precheck.json').open(encoding='utf-8'))['records']
afe={r['旧AFE'] for r in records if r['旧AFE']}
out={'accessed_utc':datetime.now(timezone.utc).isoformat(),'source_info':'https://rest.kegg.jp/info/afr','by_afe':{a:{} for a in afe},'request_status':{}}
for target in ('ko','enzyme','pathway'):
    url=f'https://rest.kegg.jp/link/{target}/afr'
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'G52-evidence-review/1.0'})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req,timeout=45) as response:
                    body=response.read().decode('utf-8')
                break
            except (urllib.error.URLError,OSError):
                if attempt==2:
                    raise
                time.sleep(2*(attempt+1))
        out['request_status'][target]={'url':url,'lines':len(body.splitlines()),'status':'ok'}
        for line in body.splitlines():
            parts=line.split('\t')
            if len(parts)!=2:
                continue
            a=parts[0].removeprefix('afr:')
            if a in afe:
                out['by_afe'][a].setdefault(target,[]).append(parts[1])
    except (urllib.error.URLError,OSError) as exc:
        out['request_status'][target]={'url':url,'status':str(exc)}

kos=sorted({ko for values in out['by_afe'].values() for ko in values.get('ko',[])})
ko_reactions={ko:[] for ko in kos}
for index in range(0,len(kos),10):
    chunk=kos[index:index+10]
    url='https://rest.kegg.jp/link/reaction/' + '+'.join(chunk)
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'G52-evidence-review/1.0'})
            with urllib.request.urlopen(req,timeout=45) as response:
                body=response.read().decode('utf-8')
            for line in body.splitlines():
                parts=line.split('\t')
                if len(parts)==2 and parts[0] in ko_reactions:
                    ko_reactions[parts[0]].append(parts[1])
            break
        except (urllib.error.URLError,OSError):
            if attempt==2:
                out['request_status']['reaction_chunk_'+str(index//10)]={'url':url,'status':'failed after retries'}
            else:
                time.sleep(attempt+1)
out['ko_reactions']=ko_reactions

with (DIR/'g52_kegg.json').open('w',encoding='utf-8') as f:
    json.dump(out,f,ensure_ascii=False,indent=2)
print(json.dumps({'targets':len(afe),'status':out['request_status'],'KO':sum(bool(v.get('ko')) for v in out['by_afe'].values()),'enzyme':sum(bool(v.get('enzyme')) for v in out['by_afe'].values())},ensure_ascii=False))
