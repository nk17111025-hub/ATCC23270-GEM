import csv, hashlib, json, re, time, urllib.request, urllib.parse, urllib.error
from datetime import datetime
from pathlib import Path

R=Path(r'D:\嗜酸氧化亚铁硫杆菌'); ST=datetime.now().astimezone().isoformat()
OUT=R/'09_差异比较'/'D3-1'; OUT.mkdir(parents=True,exist_ok=True)
dirs={x:R/x/'D3-1' for x in ['05_BRENDA','06_Rhea','07_MetaCyc']}
for p in dirs.values(): p.mkdir(parents=True,exist_ok=True)
def readtsv(path):
    with open(path,encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(path,cols,rows):
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def split(s): return [x.strip() for x in re.split(r'[;,|]',s or '') if x.strip()]
def fetchraw(folder,name,url):
    cached=folder/'原始响应'/name
    if cached.exists():
        data=cached.read_bytes(); return data.decode('utf-8',errors='replace'),{'文件':str(cached.relative_to(R)),'URL/API':url,'获取时间':datetime.fromtimestamp(cached.stat().st_mtime).astimezone().isoformat(),'数据库版本':'见响应内容/HTTP记录缺失于缓存复用','字节数':len(data),'SHA256':hashlib.sha256(data).hexdigest(),'HTTP状态':'200 (缓存复用；原HTTP状态见前次运行)'}
    req=urllib.request.Request(url,headers={'User-Agent':'D3-1 ATCC23270 audit; contact not supplied'})
    try:
        with urllib.request.urlopen(req,timeout=90) as z: data=z.read(); status=z.status; hs=dict(z.headers)
    except urllib.error.HTTPError as e: data=e.read(); status=e.code; hs=dict(e.headers)
    p=folder/'原始响应'; p.mkdir(exist_ok=True); f=p/name; f.write_bytes(data)
    return data.decode('utf-8',errors='replace'),{'文件':str(f.relative_to(R)),'URL/API':url,'获取时间':datetime.now().astimezone().isoformat(),'数据库版本':hs.get('Last-Modified','未提供'),'字节数':len(data),'SHA256':hashlib.sha256(data).hexdigest(),'HTTP状态':status}

# source tables
diff=readtsv(R/'08_模型基线'/'A0'/'2016_到2024_真实差异.tsv')
gpr=readtsv(R/'08_模型基线'/'A0'/'恢复GPR及原始证据定位.tsv')
mp=readtsv(R/'02_NCBI'/'标准化数据'/'旧新基因映射.tsv')
bio_genes=readtsv(R/'03_BioCyc'/'C2-1'/'标准化数据'/'genes.tsv')
bio_proteins=readtsv(R/'03_BioCyc'/'C2-1'/'标准化数据'/'proteins.tsv')
bio_rx=readtsv(R/'03_BioCyc'/'C2-1'/'标准化数据'/'reactions.tsv')
bio_enzyme_rx=readtsv(R/'03_BioCyc'/'C2-1'/'标准化数据'/'enzymatic-reactions.tsv')
bio_gene={x.get('对象ID',''):x for x in bio_genes}
mpold={x.get('旧AFE位点','').strip('"'):x for x in mp}
mpcur={x.get('当前GCF049位点','').strip('"'):x for x in mp}
biomap={}
for x in bio_genes:
    for alias in [x.get('对象ID',''),*split(x.get('同义词',''))]: biomap[alias]=x
bioprot={x.get('基因',''):x for x in bio_proteins}
bio_prot_by_id={x.get('对象ID',''):x for x in bio_proteins}
bioreact={x.get('对象ID',''):x for x in bio_rx}
bio_gene_rx={}
for er in bio_enzyme_rx:
    pr= bio_prot_by_id.get(er.get('酶',''),{}); gene=pr.get('基因','')
    for old,mm in mpold.items():
        if gene and gene==mm.get('旧RefSeq位点',''):
            bio_gene_rx.setdefault(old,set()).add(er.get('反应',''))

# Read KEGG raw-derived mapping; link KEGG locus ids back to the old AFE tags used by both papers.
kpath=R/'04_KEGG'/'D3-1'/'KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv'
krows=readtsv(kpath) if kpath.exists() else []
kby={}
for x in krows: kby.setdefault(x.get('KEGG gene/locus',''),[]).append(x)
def kegg_for(g):
    return kby.get(g,[])

# BRENDA is queried through the official public DSMZ SPARQL endpoint; retain strain names verbatim.
braw=R/'05_BRENDA'/'D3-1'/'原始响应'
def json_rows(name):
    p=braw/name
    if not p.exists(): return []
    try: return json.loads(p.read_text(encoding='utf-8'))['results']['bindings']
    except Exception: return []
bentries=json_rows('BRENDA_org_EC.json')
bcomp=json_rows('BRENDA_FeS_reaction_compounds.json')
bcof=json_rows('BRENDA_FeS_cofactors_references.json')
def bv(x,k): return x.get(k,{}).get('value','')
def short(uri): return uri.rsplit('/',1)[-1]
bren_by_ec={}
for x in bentries:
    ec=short(bv(x,'ec')); org=bv(x,'orgName'); oid=short(bv(x,'org'))
    if ec: bren_by_ec.setdefault(ec,[]).append((org,oid,bv(x,'ecName'),short(bv(x,'enzyme'))))
bren_detail={}
for x in [*bcomp,*bcof]:
    ec=short(bv(x,'ec')); org=bv(x,'orgName')
    if not ec or 'ATCC 23270' not in org: continue
    q=bren_detail.setdefault(ec,{'substrates':set(),'products':set(),'cofactors':set(),'references':set(),'reaction_ids':set()})
    role=bv(x,'roleClass').rsplit('/',1)[-1]
    if role=='Substrate': q['substrates'].add(bv(x,'compoundName'))
    if role=='Product': q['products'].add(bv(x,'compoundName'))
    for key,dst in [('cofName','cofactors'),('ref','references'),('rxn','reaction_ids')]:
        if bv(x,key): q[dst].add(short(bv(x,key)) if key!='cofName' else bv(x,key))
brout=[]
for ec,ents in sorted(bren_by_ec.items()):
    for org,oid,ecname,enzyme in ents:
        detail=bren_detail.get(ec,{}) if oid in ('16285','17182','22717') else {}
        same='同株标签明确' if oid in ('16285','17182') else ('BRENDA记录明确写明 ATCC 23270，但兼列 CCM 4253' if oid=='22717' else 'DSM 14882；作为同型株相关记录保留，未据此标为同株')
        brout.append({'EC号':ec,'酶名称':ecname,'BRENDA organism标签':org,'organism record ID':oid,'同株判断':same,'BRENDA酶ID':enzyme,'反应条目ID':';'.join(sorted(detail.get('reaction_ids',set()))),'底物/反应参与物':';'.join(sorted(detail.get('substrates',set()))),'产物/反应参与物':';'.join(sorted(detail.get('products',set()))),'辅因子':';'.join(sorted(detail.get('cofactors',set()))),'BRENDA文献ID':';'.join(sorted(detail.get('references',set()))),'动力学':('本轮SPARQL未导出动力学值；可按文献ID追溯' if detail else '本条未取得'),'证据说明':'BRENDA enzyme/EC 条目带菌株标签；反应参与物可能来自该酶类通用条目，需回溯文献确认是否为该菌株实验'})
write(dirs['05_BRENDA']/'BRENDA_酶学结果.tsv',list(brout[0]) if brout else ['EC'],brout)

# Rhea subset: direct MetaCyc IDs first, then ECs from BioCyc reactions; preserve every API response.
rdir=dirs['06_Rhea']; ec_set=sorted({ec.replace('EC-','').replace('EC:','') for x in bio_rx for ec in split(x.get('EC号','')) if re.match(r'^(?:EC[-:]?)?\d+\.\d+\.\d+\.[\d-]+$',ec.replace('EC-','').replace('EC:',''))})
queries=[]; rmanifest=[]
cols='rhea-id,equation,ec,chebi-id,reaction-xref(KEGG),reaction-xref(MetaCyc)'
for i in range(0,len(ec_set),8):
    batch=ec_set[i:i+8]; q=' OR '.join('ec:'+e for e in batch)
    url='https://www.rhea-db.org/rhea/?'+urllib.parse.urlencode({'query':q,'columns':cols,'format':'tsv','limit':1000})
    txt,meta=fetchraw(rdir,f'Rhea_EC_batch_{i//8+1:03d}.tsv',url); rmanifest.append(meta)
    lines=txt.splitlines()
    if lines:
        if not queries: queries.append(lines[0])
        queries.extend(lines[1:])
rrhea={}
if queries:
    lines=list(dict.fromkeys(queries))
    with (rdir/'Rhea_modelBioCyc_EC_crossrefs.tsv').open('w',encoding='utf-8-sig') as f: f.write('\n'.join(lines)+'\n')
    hdr=lines[0].split('\t')
    for ln in lines[1:]:
        a=ln.split('\t'); d=dict(zip(hdr,a))
        for ref in split(d.get('Cross-reference (MetaCyc)','')): rrhea.setdefault(ref.replace('MetaCyc:',''),[]).append(d)
        for ref in split(d.get('Cross-reference (KEGG)','')): rrhea.setdefault(ref.replace('KEGG:',''),[]).append(d)
write(rdir/'Rhea_API_download_manifest.tsv',list(rmanifest[0]) if rmanifest else ['文件','URL/API','获取时间','数据库版本','字节数','SHA256','HTTP状态'],rmanifest)

# MetaCyc IDs are already embedded as BioCyc/Pathway Tools reaction object IDs. Retain as ID-exact crosswalk only.
metarows=[]
for rx in bio_rx:
    rid=rx.get('对象ID',''); refs=rrhea.get(rid,[])
    metarows.append({'MetaCyc/Pathway Tools ID':rid,'BioCyc PGDB':rx.get('PGDB','GCF_000021485'),'BioCyc equation participants':f"{rx.get('左侧','')} => {rx.get('右侧','')}",'direction':rx.get('方向',''),'EC':rx.get('EC号',''),'pathway':rx.get('通路',''),'BioCyc evidence code':rx.get('证据代码',''),'Rhea ID exact cross-reference':';'.join(sorted({q.get('Reaction identifier','') for q in refs})),'MetaCyc direct retrieval':'not independently retrieved; ID inherited from BioCyc PGDB','interpretation':'MetaCyc pathway definition candidate; not ATCC23270 current-genome evidence'})
write(dirs['07_MetaCyc']/'MetaCyc_BioCyc_ID精确映射.tsv',list(metarows[0]) if metarows else ['MetaCyc/Pathway Tools ID'],metarows)

# Map model genes and reactions, never infer gene loss from locus-tag changes.
gprby={x['反应ID']:x for x in gpr}
master=[]; candidates=[]; conflicts=[]
for d in diff:
    rid=d.get('反应ID',''); gp=gprby.get(rid,{})
    genes=sorted(set(re.findall(r'AFE_\d{4}',(d.get('2016GPR','')+' '+d.get('2024GPR','')+' '+gp.get('恢复GPR','')))))
    biocand=[]; kec=[]; rhea=[]
    for old in genes:
        mm=mpold.get(old,{})
        bio_id=mm.get('旧RefSeq位点','')
        bg=biomap.get(bio_id,{}) or biomap.get(old,{})
        ke=kegg_for(old); kec.extend([x.get('EC','') for x in ke if x.get('EC')])
        for k in ke:
            for ref in [*split(k.get('KEGG reaction','')),*split(k.get('EC',''))]: rhea.extend(rrhea.get(ref,[]))
        bprot=bioprot.get(bio_id,{})
        brx=';'.join(sorted(bio_gene_rx.get(old,set())))
        biocand.append(f"{old} => current {mm.get('当前GCF049位点','未映射')}; NCBI={mm.get('当前产物','未映射')}; BioCyc gene product={bprot.get('名称','') or bg.get('名称','') or bg.get('对象ID','')}; BioCyc gene-associated reaction IDs (not reaction-equivalence mapping)={brx}; KEGG afr KO={';'.join(sorted({x.get('KO','') for x in ke if x.get('KO')}))}")
    rx=b in [] if False else None
    model_rxn=d.get('2024规范化反应式','') or d.get('2016规范化反应式','')
    matches=[x for x in bio_rx if x.get('对象ID','')==rid or (rid and rid in (x.get('外部ID','')+' '+x.get('对象ID','')))]
    matches=list({x.get('对象ID',''):x for x in matches}.values())
    for x in matches:
        kec.extend(ec.replace('EC-','').replace('EC:','') for ec in split(x.get('EC号','')) if ec)
    br='; '.join(f"EXACT ID {x.get('对象ID','')}[{x.get('方向','')};EC={x.get('EC号','')};evidence={x.get('证据代码','')};{x.get('左侧','')}=>{x.get('右侧','')}]" for x in matches[:5])
    meta_ids=';'.join(sorted({x.get('对象ID','') for x in matches}))
    refs=[]
    for x in matches: refs+=rrhea.get(x.get('对象ID',''),[])
    krefs=[]
    for old in genes:
        for k in kegg_for(old):
            for kr in split(k.get('KEGG reaction','')): krefs.extend(rrhea.get(kr,[]))
    refs=list({x.get('Reaction identifier',''):x for x in [*refs,*krefs]}.values())
    rhea_str='; '.join(f"{x.get('Reaction identifier','')} {x.get('Equation','')} [EC={x.get('EC number','')};KEGG={x.get('Cross-reference (KEGG)','')};MetaCyc={x.get('Cross-reference (MetaCyc)','')}]" for x in refs[:5])
    state='数据库一致性需人工判读'
    suggestion='暂不处理'; reason='模型基线反应保留；数据库命中不自动构成模型修改证据。'
    status=d.get('状态','')
    if genes and any(not mpold.get(x,{}).get('当前GCF049位点') for x in genes): suggestion='证据冲突'; reason='B1旧新基因映射有未决条目；不能据此判断缺失。'
    elif '差异' in status and status not in ('共同且无语义差异','共同且仅ID变化'):
        suggestion='反应定义修改候选'; reason='2016/2024模型字段存在差异，须结合 A0 证据和化学计量复核；不得自动采用。'
    if not genes and not gp.get('恢复GPR',''): suggestion='暂不处理'; reason='模型反应没有可恢复GPR或该条目可能为交换/汇总反应。'
    br_text=[]
    for ec in sorted(set(kec)):
        ents=bren_by_ec.get(ec,[]); direct=[z for z in ents if z[1] in ('16285','17182','22717')]
        if direct:
            de=bren_detail.get(ec,{})
            br_text.append(f"EC {ec}: ATCC23270-specific BRENDA organism record(s)={len(direct)}; substrates={';'.join(sorted(de.get('substrates',set())))}; products={';'.join(sorted(de.get('products',set())))}; cofactors={';'.join(sorted(de.get('cofactors',set())))}; BRENDA refs={';'.join(sorted(de.get('references',set())))}; reaction IDs={';'.join(sorted(de.get('reaction_ids',set())))}")
        elif ents: br_text.append(f"EC {ec}: only non-exact organism record(s): "+';'.join(sorted({z[0] for z in ents})))
    rhea_cofactors=sorted({m for x in refs for m in re.findall(r'\b(?:NADP?H?|FAD|FMN|ATP|ADP|GTP|GDP|ferredoxin|ubiquinone|ubiquinol|cytochrome|CoA|acetyl-CoA|thiamine diphosphate|pyridoxal 5.?phosphate)\b',x.get('Equation',''),re.I)})
    row={'旧模型ID':rid,'2016 GPR':gp.get('2016 Excel GPR',d.get('2016GPR','')),'2024 GPR':gp.get('2024 S1 GPR',d.get('2024GPR','')),'旧 AFE locus':';'.join(genes),'旧 RefSeq AFE_RS locus':';'.join(mpold.get(x,{}).get('旧RefSeq位点','') for x in genes),'当前 NCBI locus tag (GCF_049532655.1)':';'.join(mpold.get(x,{}).get('当前GCF049位点','') for x in genes),'NCBI 当前功能':'; '.join(mpold.get(x,{}).get('当前产物','') for x in genes),'BioCyc GCF_000021485':br,'KEGG afr':'; '.join(biocand),'BRENDA':' | '.join(br_text) if br_text else '此模型行未关联到已查询的ATCC23270-specific BRENDA enzyme record','Rhea':rhea_str,'MetaCyc':meta_ids,'EC':';'.join(sorted(set(kec))),'reaction equation (2024/2016)':model_rxn,'direction':d.get('2024上下界',''),'Rhea chemical comparison':'Rhea candidate equations/ChEBI IDs are shown alongside the model equation; exact compound mapping and stoichiometry comparison remain manual.','cofactor':'Rhea candidate participants: '+';'.join(rhea_cofactors) if rhea_cofactors else 'not resolved from an explicit reaction match','transport information':'model transport/boundary reaction; substrate specificity must be established independently','数据库一致性':state,'同株实验支持情况':'BRENDA organism-specific records only when exact ATCC 23270 culture label is present; BioCyc EV-COMP remains computational','建议处理':suggestion,'理由':reason,'A0差异状态':status}
    has_current=any(mpold.get(x,{}).get('当前GCF049位点') and mpold.get(x,{}).get('当前产物') for x in genes)
    has_kegg=any(kegg_for(x) for x in genes)
    has_direct_br=any(any(z[1] in ('16285','17182','22717') for z in bren_by_ec.get(ec,[])) for ec in set(kec))
    if has_direct_br and has_current and has_kegg and suggestion!='证据冲突': grade='B (organism-specific enzyme record + current genome annotation; reaction stoichiometry still requires review)'
    elif has_current and has_kegg: grade='C (current genome + KEGG functional annotation; no same-strain gene-level biochemical linkage established)'
    elif any('DSM 14882' in z[0] for ec in set(kec) for z in bren_by_ec.get(ec,[])): grade='D (related strain only)'
    else: grade='E (automated prediction or incomplete mapping only)'
    row['证据等级']=grade
    master.append(row)
    if suggestion!='暂不处理': candidates.append(row)
    if matches and (d.get('2016规范化反应式','') and model_rxn and not any(x.get('对象ID','') in rrhea for x in matches)): pass
    if suggestion=='证据冲突': conflicts.append(row)

cols=list(master[0]) if master else ['旧模型ID']
write(OUT/'D3-1_跨数据库主表.tsv',cols,master); write(OUT/'D3-1_候选修改.tsv',cols,candidates); write(OUT/'D3-1_数据库冲突.tsv',cols,conflicts)

# Focused thematic audit from IDs, GPR, annotations and model equations.
def audit_rows(kind,pattern):
    return [x for x in master if re.search(pattern,' '.join([x.get('旧模型ID',''),x.get('2016 GPR',''),x.get('2024 GPR',''),x.get('NCBI 当前功能',''),x.get('KEGG afr',''),x.get('BioCyc GCF_000021485','')]),re.I)]
write(OUT/'D3-1_中央碳代谢重点审计.tsv',cols,audit_rows('central',r'EMP|glycol|glucose|glc|pfk|PPP|pentose|pyruvat|acetyl|ACCOA|TCA|citrate|succ|malate|oxoglutar|oxoglutarate|CBB|ribulose|RuBisCO|AFE_1807|AFE_2841'))
transport_rows=audit_rows('transport',r'transport|permease|MFS|ABC transporter|sugar|glucose|organic acid|amino acid|symporter|antiporter')
for tr in transport_rows:
    textblob=' '.join([tr.get('NCBI 当前功能',''),tr.get('BioCyc GCF_000021485',''),tr.get('KEGG afr','')]).lower()
    if ('transporter' in textblob or 'transport' in textblob) and any(x in textblob for x in ('glucose','glc','sugar','succinate','citrate','amino acid')): tr['转运底物证据类型']='注释点名底物/底物类别；仍需同株功能验证'
    elif any(x in textblob for x in ('mfs','permease','abc transporter','transporter','transport')): tr['转运底物证据类型']='家族级或计算预测；未明确具体底物'
    else: tr['转运底物证据类型']='仅计算预测/间接推断；不能指定底物'
write(OUT/'D3-1_有机碳转运重点审计.tsv',list(transport_rows[0]) if transport_rows else [*cols,'转运底物证据类型'],transport_rows)
write(OUT/'D3-1_FeS_ETC重点审计.tsv',cols,audit_rows('resp',r'cyc2|rus|petI|petII|bc1|nuo|oxidase|ATP synthase|SQR|sulfide:quinone|tce|TetH|TQO|Sox|HDR|iron oxidation|sulfur oxidation|respiratory'))

# Coverage, evidence hierarchy, audit and per-source raw checksum catalog.
bio_ids={x.get('对象ID','') for x in bio_rx}; model_ids={x.get('旧模型ID','') for x in master}
coverage=[{'数据库':'KEGG','覆盖分母':'631 model reaction rows；按GPR中的旧AFE gene locus连接afr gene annotation','覆盖命中':sum(1 for x in master if x.get('KEGG afr')),'覆盖率':round(sum(1 for x in master if x.get('KEGG afr'))/max(1,len(master)),4),'说明':'afr是ATCC23270 KEGG条目；KEGG参考assembly为旧GCA_000021485.1。获取gene→KO/EC/pathway/module，并经KO→reaction补全gene→KO→reaction链；直接gene→reaction REST link为空。'},
{'数据库':'BioCyc','覆盖分母':len(master),'覆盖命中':sum(x.get('BioCyc GCF_000021485','').startswith('EXACT ID ') for x in master),'覆盖率':round(sum(x.get('BioCyc GCF_000021485','').startswith('EXACT ID ') for x in master)/max(1,len(master)),4),'说明':'仅 exact reaction ID / declared external ID 算模型反应覆盖；gene→reaction 映射作为酶基因候选，未视为反应等价。PGDB Tier 3、2021旧assembly，EV-COMP为计算预测。'},
{'数据库':'Rhea','覆盖分母':len(master),'覆盖命中':sum(bool(x.get('Rhea')) for x in master),'覆盖率':round(sum(bool(x.get('Rhea')) for x in master)/max(1,len(master)),4),'说明':'REST按BioCyc EC查询；KEGG reaction xref映射保留为化学反应候选，需比较计量、质子、电子和方向后才算确定映射。'},
{'数据库':'MetaCyc','覆盖分母':len(master),'覆盖命中':sum(bool(x.get('MetaCyc')) for x in master),'覆盖率':round(sum(bool(x.get('MetaCyc')) for x in master)/max(1,len(master)),4),'说明':'只计BioCyc/MetaCyc反应ID exact match；未独立读取MetaCyc条目。'},
{'数据库':'BRENDA','覆盖分母':'目标同株BRENDA SPARQL enzyme/EC records','覆盖命中':len({(short(bv(x,'enzyme')),short(bv(x,'ec'))) for x in bentries if short(bv(x,'org')) in ('16285','17182','22717')}),'覆盖率':'按同株 BRENDA EC record count；不是模型reaction覆盖率','说明':'已通过官方DSMZ BRENDA SPARQL prototype取得有效enzyme/EC记录和Fe/S反应底物/产物/辅因子/文献ID。BRENDA SOAP须注册，但不影响公开SPARQL抓取。DSM 14882单独保留，不标同株。'}]
write(OUT/'D3-1_数据库覆盖率.tsv',list(coverage[0]),coverage)

meta_file=R/'04_KEGG'/'D3-1'/'KEGG_download_manifest.tsv'
ke_meta=readtsv(meta_file) if meta_file.exists() else []
bmanifest=R/'05_BRENDA'/'D3-1'/'BRENDA_download_manifest.json'
if bmanifest.exists():
    for bm in json.loads(bmanifest.read_text(encoding='utf-8')):
        ke_meta.append({'文件':bm.get('file',''),'URL/API':bm.get('URL/API',''),'获取时间':bm.get('access_time',''),'数据库版本':bm.get('version',''),'字节数':bm.get('bytes',''),'SHA256':bm.get('SHA256',''),'HTTP状态':bm.get('HTTP_status','')})
for src,folder in [('KEGG',R/'04_KEGG'/'D3-1'/'原始响应'),('BRENDA',R/'05_BRENDA'/'D3-1'/'原始响应'),('Rhea',rdir/'原始响应')]:
    for f in folder.glob('*'):
        if f.is_file():
            rel=str(f.relative_to(R))
            if any(x.get('文件')==rel for x in ke_meta): continue
            endpoint={'KEGG':'https://rest.kegg.jp/','BRENDA':'https://sparql.dsmz.de/api/brenda','Rhea':'https://www.rhea-db.org/rhea/'}[src]
            version={'KEGG':'版本字段未由REST响应提供；同株旧参考assembly GCA_000021485.1','BRENDA':'BRENDA 2026.1 / DSMZ public SPARQL prototype','Rhea':'Rhea REST查询；release header按原始响应为准'}[src]
            ke_meta.append({'文件':rel,'URL/API':endpoint,'获取时间':datetime.fromtimestamp(f.stat().st_mtime).astimezone().isoformat(),'数据库版本':version,'字节数':f.stat().st_size,'SHA256':hashlib.sha256(f.read_bytes()).hexdigest(),'HTTP状态':'200；各最终查询状态见来源manifest；限量探索响应按原文件标注'})
manifest_cols=['文件','来源URL/API','获取时间','数据库版本/注释','字节数','SHA256','HTTP状态','ETag','补充备注']
manifest_norm=[]
for x in ke_meta:
    manifest_norm.append({'文件':x.get('文件',x.get('file','')),'来源URL/API':x.get('URL/API',x.get('url','')),'获取时间':x.get('获取时间',x.get('time_local',x.get('access_time',''))),'数据库版本/注释':x.get('数据库版本',x.get('database_version_or_release',x.get('version','未提供'))),'字节数':x.get('字节数',x.get('bytes','')),'SHA256':x.get('SHA256',x.get('sha256','')),'HTTP状态':x.get('HTTP状态',x.get('HTTP_status',x.get('http_status',''))),'ETag':x.get('etag',''),'补充备注':x.get('last_modified','')})
write(OUT/'D3-1_原始数据校验清单.tsv',manifest_cols,manifest_norm)

# counts comparing BioCyc old-locus gene symbols/function to current NCBI product annotations via old reference mapping
biodiff=0; bion=0
for old,mm in mpold.items():
    ref=mm.get('旧RefSeq位点',''); bg=biomap.get(ref,{}) or biomap.get(old,{})
    bprot=bioprot.get(ref,{})
    if bprot and mm.get('当前产物'):
        bion+=1
        if bprot.get('名称','').strip().lower()!=mm.get('当前产物','').strip().lower(): biodiff+=1
summary={'model_reactions':len(master),'candidates':len(candidates),'conflicts':len(conflicts),'bio_gene_product_comparable':bion,'bio_vs_ncbi_product_text_differences':biodiff,'bio_reactions':len(bio_rx),'rhea_hits':sum(bool(x.get('Rhea')) for x in master),'rhea_queries':len(rmanifest),'brenda_same_strain_EC_enzyme_records':sum(1 for x in bentries if short(bv(x,'org')) in ('16285','17182','22717')),'brenda_distinct_EC_total':len(bren_by_ec),'evidence_status':'BRENDA data acquired via public SPARQL; independent MetaCyc access limited'}
(OUT/'D3-1_metrics.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False))
