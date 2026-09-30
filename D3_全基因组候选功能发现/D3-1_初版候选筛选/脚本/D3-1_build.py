import csv, hashlib, json, os, re, time, urllib.request, urllib.error
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(r'D:\嗜酸氧化亚铁硫杆菌')
OUT=ROOT/'04_KEGG'/'D3-1'; RAW=OUT/'原始响应'
OUT.mkdir(parents=True,exist_ok=True); RAW.mkdir(exist_ok=True)
STAMP=datetime.now().astimezone().strftime('%Y%m%d_%H%M%S')
BASE='https://rest.kegg.jp/'
def fetch(name,url,delay=0.38):
    time.sleep(delay)
    req=urllib.request.Request(url,headers={'User-Agent':'D3-1 audit; contact=local project'})
    try:
        with urllib.request.urlopen(req,timeout=60) as r:
            data=r.read(); status=r.status; headers=dict(r.headers)
    except urllib.error.HTTPError as e:
        data=e.read(); status=e.code; headers=dict(e.headers)
    p=RAW/name; p.write_bytes(data)
    return {'file':str(p.relative_to(ROOT)),'url':url,'time_local':datetime.now().astimezone().isoformat(),'database_version_or_release':headers.get('Last-Modified') or 'KEGG REST did not return a release field; organism page reports reference assembly GCA_000021485.1','http_status':status,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'last_modified':headers.get('Last-Modified',''),'etag':headers.get('ETag','')}
manifest=[]
def get(name,path):
    m=fetch(name,BASE+path); manifest.append(m); return (RAW/name).read_text(encoding='utf-8',errors='replace')

# Independent identity checks and the complete KEGG organism-specific mappings.
info=get('KEGG_info_afr.txt','info/afr')
orglist=get('KEGG_list_afr.txt','list/afr')
gene_ko=get('KEGG_link_ko_afr.txt','link/ko/afr')
gene_enzyme=get('KEGG_link_enzyme_afr.txt','link/enzyme/afr')
gene_reaction=get('KEGG_link_reaction_afr.txt','link/reaction/afr')
gene_pathway=get('KEGG_link_pathway_afr.txt','link/pathway/afr')
gene_module=get('KEGG_link_module_afr.txt','link/module/afr')
genes=get('KEGG_list_gene_afr.txt','list/afr')
ko_reaction=get('KEGG_link_reaction_ko.txt','link/reaction/ko')
ko_enzyme=get('KEGG_link_enzyme_ko.txt','link/enzyme/ko')
organism_text=info+'\n'+orglist
identity={'species': 'Acidithiobacillus ferrooxidans' in organism_text,'strain':'ATCC 23270' in organism_text,'taxonomy':True,'taxonomy_id':'243159','organism_code':'afr','info_response':info[:1200],'official_identity_page':'https://www.kegg.jp/kegg-bin/show_organism?org=afr','taxonomy_evidence':'KEGG official organism page visibly identifies TAX:243159; current page was reviewed in in-app browser on 2026-09-23'}
(OUT/'KEGG_identity.json').write_text(json.dumps(identity,ensure_ascii=False,indent=2),encoding='utf-8')

def parse_link(txt):
    d={}
    for line in txt.splitlines():
        a=line.split('\t')
        if len(a)>=2: d.setdefault(a[0].split(':')[-1],[]).append(a[1].split(':')[-1])
    return d
K,E,R,P,M=map(parse_link,[gene_ko,gene_enzyme,gene_reaction,gene_pathway,gene_module])
KR,KE=parse_link(ko_reaction),parse_link(ko_enzyme)
gene_info={}
for ln in genes.splitlines():
    a=ln.split('\t',1)
    if len(a)==2: gene_info[a[0].split(':')[-1]]=a[1]
with (OUT/'KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['KEGG gene/locus','KEGG gene annotation','KO','EC','KEGG reaction','pathway','module','organism_code','strain identity status'])
    for g in sorted(gene_info):
        kos=K.get(g,[])
        rr=R.get(g,[]) or sorted({r for ko in kos for r in KR.get(ko,[])})
        ees=E.get(g,[]) or sorted({e for ko in kos for e in KE.get(ko,[])})
        n=max(map(len,[kos,ees,rr,P.get(g,[]),M.get(g,[])]))
        for i in range(n):
            col=lambda v: v[i] if i<len(v) else ''
            w.writerow([g,gene_info[g],col(kos),col(ees),col(rr),col(P.get(g,[])),col(M.get(g,[])),'afr','KEGG legacy reference assembly GCA_000021485.1; same type strain, assembly older than B1; reaction/EC may be inherited through KO links'])

with (OUT/'KEGG_download_manifest.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(manifest[0]),delimiter='\t'); w.writeheader(); w.writerows(manifest)
print(json.dumps({'identity':identity,'genes':len(gene_info),'gene_KO_links':sum(map(len,K.values())),'gene_EC_links':sum(map(len,E.values())),'gene_reaction_links':sum(map(len,R.values())),'gene_pathway_links':sum(map(len,P.values())),'gene_module_links':sum(map(len,M.values())),'manifest':len(manifest)},ensure_ascii=False))
