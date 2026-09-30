from urllib.parse import urlencode
from urllib.request import Request,urlopen
from pathlib import Path
import datetime

out=Path('outputs/G53/uniprot_243159.tsv')
query='organism_id:243159'
fields='accession,id,reviewed,protein_name,xref_refseq,xref_interpro,xref_pfam'
url='https://rest.uniprot.org/uniprotkb/stream?'+urlencode({'query':query,'format':'tsv','fields':fields})
req=Request(url,headers={'User-Agent':'G53-evidence-review/1.0'})
with urlopen(req,timeout=180) as resp:
    data=resp.read()
    modified=resp.headers.get('Last-Modified','')
out.write_bytes(data)
print('url',url)
print('bytes',len(data),'lines',data.count(b'\n'),'http_last_modified',modified,'retrieved_utc',datetime.datetime.now(datetime.timezone.utc).isoformat())
