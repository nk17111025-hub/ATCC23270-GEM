import csv,json,re,time,datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen

out=Path('outputs/G53')
genes=json.loads((out/'g53_data.json').read_text(encoding='utf-8'))['index']
found=set()
with (out/'uniprot_243159.tsv').open(encoding='utf-8') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        found.update(re.findall(r'WP_[0-9]+[.][0-9]+',row.get('RefSeq','')))
missing=[r[2] for r in genes if r[2] not in found]
fields='accession,id,reviewed,protein_name,xref_refseq,xref_interpro,xref_pfam,organism_name'
allrows=[]
for i in range(0,len(missing),10):
    group=missing[i:i+10]
    query=' OR '.join(f'xref:RefSeq-{wp}' for wp in group)
    url='https://rest.uniprot.org/uniprotkb/stream?'+urlencode({'query':query,'format':'tsv','fields':fields})
    try:
        with urlopen(Request(url,headers={'User-Agent':'G53-evidence-review/1.0'}),timeout=90) as resp:
            lines=resp.read().decode('utf-8').splitlines()
        allrows+=lines[1:]
        print(i+1,len(group),'rows',len(lines)-1,flush=True)
    except Exception as exc:
        print(i+1,'ERROR',repr(exc),flush=True)
    time.sleep(0.5)
header='Entry\tEntry Name\tReviewed\tProtein names\tRefSeq\tInterPro\tPfam\tOrganism'
(out/'uniprot_missing.tsv').write_text('\n'.join([header]+allrows)+'\n',encoding='utf-8')
print('missing',len(missing),'fetched_rows',len(allrows),'retrieved_utc',datetime.datetime.now(datetime.timezone.utc).isoformat())
