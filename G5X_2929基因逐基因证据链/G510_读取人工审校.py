import csv
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'G510_UniProt_reviewed_JSON'
OUT.mkdir(exist_ok=True)
with (HERE/'G510_交付'/'G510_逐基因索引.tsv').open(encoding='utf-8-sig',newline='') as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
ids=sorted({a for r in rows if 'reviewed' in r['UniProt审校'].split(';') for a in r['UniProtKB'].split(';') if a})

def fetch(acc):
    target=OUT/f'{acc}.json'
    if not target.exists():
        req=urllib.request.Request(f'https://rest.uniprot.org/uniprotkb/{acc}.json',headers={'User-Agent':'G510 evidence audit/1.0'})
        target.write_bytes(urllib.request.urlopen(req,timeout=30).read())
    j=json.loads(target.read_text(encoding='utf-8'))
    rhea=[x['id'] for x in j.get('uniProtKBCrossReferences',[]) if x.get('database')=='Rhea']
    refs=[]
    for x in j.get('references',[]):
        c=x.get('citation',{})
        refs.append((c.get('title',''),c.get('id',''),c.get('citationType','')))
    return acc,rhea,refs

with ThreadPoolExecutor(max_workers=4) as pool:
    futures=[pool.submit(fetch,acc) for acc in ids]
    for future in as_completed(futures):
        try:
            acc,rhea,refs=future.result()
            print(acc,'Rhea',','.join(rhea),'refs',len(refs))
        except Exception as e:
            print('ERROR',str(e))
