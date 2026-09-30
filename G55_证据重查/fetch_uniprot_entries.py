import concurrent.futures
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).parent
m=json.loads((HERE/'uniprot_idmap.json').read_text(encoding='utf-8'))
accs=sorted({r['to'] for r in m['results'] if isinstance(r['to'],str)})
raw=HERE/'uniprot_raw'
raw.mkdir(exist_ok=True)

def fetch(acc):
    target=raw/(acc+'.json')
    if target.exists():
        try:
            return acc,json.loads(target.read_text(encoding='utf-8')),None
        except Exception:
            target.unlink()
    url=f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    for attempt in range(3):
        p=subprocess.run(['curl.exe','-f','-L','--silent','--show-error','--max-time','30',url],capture_output=True)
        if p.returncode==0:
            try:
                data=json.loads(p.stdout)
                target.write_bytes(p.stdout)
                return acc,data,None
            except Exception as e:
                err=str(e)
        else:
            err=p.stderr.decode(errors='replace')[:200]
        time.sleep(2*(attempt+1))
    return acc,None,err

results={}
errors={}
with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool:
    for i,(acc,data,err) in enumerate(pool.map(fetch,accs),1):
        if data:
            results[acc]=data
        else:
            errors[acc]=err
        if i%50==0:
            print('done',i,'/',len(accs),'errors',len(errors),flush=True)
by_wp=defaultdict(list)
for r in m['results']:
    a=r['to']
    x=results.get(a)
    if not x or x.get('organism',{}).get('taxonId')!=243159:
        continue
    xrefs=x.get('uniProtKBCrossReferences',[])
    comments=x.get('comments',[])
    catalytic=[]
    for c in comments:
        if c.get('commentType')=='CATALYTIC ACTIVITY':
            catalytic.append(c.get('reaction',{}))
    by_wp[r['from']].append({
        'accession':a,'reviewed':x.get('entryType',''),'entry_version':x.get('entryAudit',{}).get('entryVersion'),
        'annotation_date':x.get('entryAudit',{}).get('lastAnnotationUpdateDate'),
        'protein_existence':x.get('proteinExistence',''),
        'name':x.get('proteinDescription',{}).get('recommendedName',{}).get('fullName',{}).get('value',''),
        'ec':[v['value'] for v in x.get('proteinDescription',{}).get('recommendedName',{}).get('ecNumbers',[])],
        'interpro':[v['id'] for v in xrefs if v.get('database')=='InterPro'],
        'pfam':[v['id'] for v in xrefs if v.get('database')=='Pfam'],
        'catalytic':catalytic,
        'url':'https://www.uniprot.org/uniprotkb/'+a+'/entry',
    })
out={'date':'2026-09-24','input_count':len(accs),'errors':errors,'by_wp':by_wp}
(HERE/'uniprot_g55.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print('entries',len(results),'target_wp',len(by_wp),'target_records',sum(len(v) for v in by_wp.values()),'errors',len(errors))
