"""Map published ATCC 23270 P15 and P16.2 sequences to current proteins."""
import json
import urllib.request
from pathlib import Path
from Bio import SeqIO, pairwise2

HERE=Path(__file__).resolve().parent
ROOT=Path('01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1')
proteins=list(SeqIO.parse(ROOT/'protein.faa','fasta'))
g52={r['WP']:r for r in json.loads((HERE/'g52_precheck.json').read_text(encoding='utf-8'))['records']}
out={}
for accession in ['AY863107','AY863108']:
    fasta=HERE/f'{accession}_rhodanese_2005.fasta'
    if not fasta.exists():
        req=urllib.request.Request(f'https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=false',headers={'User-Agent':'Mozilla/5.0'})
        fasta.write_bytes(urllib.request.urlopen(req,timeout=35).read())
    old=next(SeqIO.parse(fasta,'fasta'))
    translations=[]
    for strand,seq in [('+',old.seq),('-',old.seq.reverse_complement())]:
        for frame in range(3):
            aa=str(seq[frame:].translate(to_stop=False)).rstrip('*')
            if len(aa)>30 and aa.count('*')==0:
                translations.append((strand,frame,aa))
    best=[]
    for strand,frame,aa in translations:
        for p in proteins:
            target=str(p.seq).rstrip('*')
            aln=pairwise2.align.localms(aa,target,2,-1,-5,-.5,one_alignment_only=True,score_only=True)
            best.append((aln,p.id,strand,frame,len(aa),len(target),p.id in g52))
    best.sort(reverse=True)
    out[accession]={'nt_length':len(old.seq),'top_hits':[{'score':x[0],'wp':x[1],'strand':x[2],'frame':x[3],'query_aa':x[4],'target_aa':x[5],'in_g52':x[6],'g52_locus':g52.get(x[1],{}).get('当前locus','')} for x in best[:10]]}
    if accession=='AY863108':
        current=next(x for x in SeqIO.parse(ROOT/'cds_from_genomic.fna','fasta') if 'RU820_RS02245' in x.description)
        a,b=str(old.seq),str(current.seq)
        deletion_positions=[i+1 for i in range(len(a)) if a[:i]+a[i+1:]==b]
        assert len(a)==451 and len(b)==450 and deletion_positions==[361,362,363]
        out[accession]['current_locus']='RU820_RS02245'
        out[accession]['current_wp']='WP_009564831.1'
        out[accession]['current_cds_length_nt']=len(b)
        out[accession]['single_G_insertion_in_historical_allele_position_1_based']=deletion_positions
        out[accession]['all_other_nt_identical']=True
        out[accession]['boundary']='AY863108.1 is published P16.2 partial CDS; one G insertion relative to current CDS changes the C-terminal reading frame. The experimental P16.2 reaction is not direct assay of exact current WP sequence.'
(HERE/'g52_rhodanese_2005_sequence_mapping.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(out,ensure_ascii=False))
