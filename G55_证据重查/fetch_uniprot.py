import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
data = json.loads((HERE/'g55_data.json').read_text(encoding='utf-8'))
ids = [x['wp'] for x in data]
payload = urllib.parse.urlencode({'from':'RefSeq_Protein','to':'UniProtKB','ids':','.join(ids)}).encode()
req = urllib.request.Request('https://rest.uniprot.org/idmapping/run', data=payload)
with urllib.request.urlopen(req,timeout=60) as r:
    job = json.load(r)['jobId']
print('job',job,flush=True)
for i in range(30):
    time.sleep(3)
    with urllib.request.urlopen(f'https://rest.uniprot.org/idmapping/status/{job}',timeout=60) as r:
        status=json.load(r)
    if status.get('jobStatus') in ('NEW','RUNNING'):
        continue
    break
else:
    raise RuntimeError('UniProt ID mapping timed out')
results=[]
url=f'https://rest.uniprot.org/idmapping/results/{job}?format=json&size=500'
while url:
    with urllib.request.urlopen(url,timeout=90) as r:
        page=json.load(r)
        link=r.headers.get('Link','')
    results.extend(page.get('results',[]))
    url=link.split('<',1)[1].split('>',1)[0] if 'rel="next"' in link and '<' in link else None
out={'job':job,'date':'2026-09-24','endpoint':'https://rest.uniprot.org/idmapping','results':results}
(HERE/'uniprot_idmap.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print('mapped_rows',len(results),'mapped_input',len(set(x['from'] for x in results)),'unmapped',len(ids)-len(set(x['from'] for x in results)))
