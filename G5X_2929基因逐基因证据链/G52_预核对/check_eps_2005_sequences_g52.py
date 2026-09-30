"""Map early ATCC 23270 EPS-precursor enzyme submissions to current WPs."""
import json
import urllib.request
from pathlib import Path
from Bio import SeqIO, pairwise2

HERE = Path(__file__).resolve().parent
ROOT = Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1')
proteins = list(SeqIO.parse(ROOT/'protein.faa', 'fasta'))
g52 = {r['WP']:r for r in json.loads((HERE/'g52_precheck.json').read_text(encoding='utf-8'))['records']}
out = {}
for accession in ['AY789511','AY789512']:
    fasta=HERE/f'{accession}_EPS_2005.fasta'
    if not fasta.exists():
        url=f'https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=false'
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
        fasta.write_bytes(urllib.request.urlopen(req,timeout=35).read())
    record=next(SeqIO.parse(fasta,'fasta'))
    aa=str(record.seq.translate(to_stop=True))
    hits=[]
    for p in proteins:
        target=str(p.seq).rstrip('*')
        score=pairwise2.align.localms(aa,target,2,-1,-5,-.5,score_only=True)
        hits.append((score,p.id,len(target)))
    hits.sort(reverse=True)
    out[accession]={
        'published_title':record.description,
        'nt_length':len(record.seq),
        'aa_length':len(aa),
        'top_hits':[{'score':score,'wp':wp,'current_aa_length':length,'g52_locus':g52.get(wp,{}).get('当前locus','')} for score,wp,length in hits[:5]]
    }
    if accession=='AY789511':
        current=next(x for x in SeqIO.parse(ROOT/'cds_from_genomic.fna','fasta') if 'RU820_RS02155' in x.description)
        assert str(record.seq)==str(current.seq)
        out[accession]['exact_current_cds_match']=True
        out[accession]['current_locus']='RU820_RS02155'
        out[accession]['current_wp']='WP_009567323.1'
    if accession=='AY789512':
        out[accession]['best_current_wp_outside_g52']=out[accession]['top_hits'][0]['wp']
(HERE/'g52_eps_2005_sequence_mapping.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False))
