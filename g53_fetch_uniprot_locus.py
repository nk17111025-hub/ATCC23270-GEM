from urllib.parse import urlencode
from urllib.request import Request,urlopen
from pathlib import Path
import datetime
out=Path('outputs/G53/uniprot_243159_locus.tsv')
url='https://rest.uniprot.org/uniprotkb/stream?'+urlencode({'query':'organism_id:243159','format':'tsv','fields':'accession,gene_oln,gene_primary,xref_refseq'})
with urlopen(Request(url,headers={'User-Agent':'G53-evidence-review/1.0'}),timeout=180) as resp:data=resp.read()
out.write_bytes(data)
print('bytes',len(data),'lines',data.count(b'\n'),'retrieved_utc',datetime.datetime.now(datetime.timezone.utc).isoformat())
