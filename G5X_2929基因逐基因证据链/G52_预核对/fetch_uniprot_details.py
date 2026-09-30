import concurrent.futures
import json
import time
import urllib.error
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

DIR = Path(__file__).resolve().parent
data = json.load((DIR/'g52_uniprot.json').open(encoding='utf-8'))
entries = {e['accession']:wp for wp,es in data['by_wp'].items() for e in es}
rawdir = DIR/'uniprot_raw'
rawdir.mkdir(exist_ok=True)


def fetch(item):
    accession,wp = item
    path=rawdir/f'{accession}.json'
    if path.exists():
        return accession,'cached'
    url=f'https://rest.uniprot.org/uniprotkb/{accession}.json'
    for attempt in range(3):
        try:
            request=urllib.request.Request(url,headers={'User-Agent':'G52-evidence-review/1.0','Accept':'application/json'})
            with urllib.request.urlopen(request,timeout=30) as response:
                obj=json.load(response)
            with path.open('w',encoding='utf-8') as f:
                json.dump(obj,f,ensure_ascii=False)
            return accession,'ok'
        except (urllib.error.URLError,TimeoutError,OSError) as exc:
            if attempt==2:
                return accession,'error:'+str(exc)
            time.sleep(2*(attempt+1))


outcomes=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    for item in pool.map(fetch,entries.items()):
        outcomes.append(item)
        if len(outcomes)%25==0 or len(outcomes)==len(entries):
            print('done',len(outcomes),'of',len(entries),Counter(x[1] for x in outcomes),flush=True)

summary={}
for accession,wp in entries.items():
    path=rawdir/f'{accession}.json'
    if not path.exists():
        continue
    j=json.load(path.open(encoding='utf-8'))
    citations=[]
    for r in j.get('references',[]):
        c=r.get('citation',{})
        citations.append({'title':c.get('title',''),'year':c.get('publicationDate',''),'pubmed':c.get('pubmedId',''),'doi':c.get('doi',''),'citation_type':c.get('citationType','')})
    reactions=[]
    for c in j.get('comments',[]):
        if c.get('commentType')=='CATALYTIC ACTIVITY':
            rx=c.get('reaction',{})
            reactions.append({'name':rx.get('name',''),'ec':rx.get('ecNumber',''),'cross_refs':rx.get('reactionCrossReferences',[]),'evidence':rx.get('evidences',[])})
    summary[wp]={'accession':accession,'url':f'https://www.uniprot.org/uniprotkb/{accession}/entry','entry_type':j.get('entryType',''),'organism':j.get('organism',{}).get('scientificName',''),'citations':citations,'reactions':reactions}

with (DIR/'g52_uniprot_details.json').open('w',encoding='utf-8') as f:
    json.dump({'accessed_utc':datetime.now(timezone.utc).isoformat(),'by_wp':summary,'fetch_outcomes':outcomes},f,ensure_ascii=False,indent=2)
print(json.dumps({'entries':len(entries),'details':len(summary),'errors':[x for x in outcomes if x[1].startswith('error')]},ensure_ascii=False))
