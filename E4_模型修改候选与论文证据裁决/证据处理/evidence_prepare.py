import json, os, re, urllib.request
from pathlib import Path
import openpyxl

root=Path(__file__).parent
text=(root/'ATCC23270-2026 最终模型修改清单.md').read_text(encoding='utf-8')
base=Path(os.environ['TEMP'])/'atcc2016read'/'mmc1.xlsx'
ws=openpyxl.load_workbook(base,read_only=True,data_only=True)['Table 1']
original=[[c.value for c in row] for row in ws.iter_rows(min_row=1,max_row=617,min_col=1,max_col=16)]
headers=original[1]
records=original[2:]
byrid={str(r[0]):r for r in records if r[0]}

groups={}
for part,body in re.findall(r'### (0[1-5]) [^\n]+\n(.*?)(?=\n### |\n\*\*专项结论)',text,re.S):
    for line in body.splitlines():
        if not line.startswith('|') or line.startswith('|---') or line.startswith('| 候选'):
            continue
        cells=[c.strip() for c in line.strip('|').split('|')]
        if len(cells)<5: continue
        cid=cells[0]
        if not (cid.startswith('B-') or cid in {'BDGK-GPR','BDGK-方向','GAP2','GLCP','SBPASE','CA2tpp','NA1tpp','ACCOAC','生物量反应','CTPS2','GLBRAN1','GLCS1','GLDBRAN2','MTRP'}): continue
        groups[cid]={'group':part,'object':cells[1], 'source_summary':cells[2], 'score_text':cells[3], 'boundary':' | '.join(cells[4:])}

evidence={'B-002','B-004','GLCP','SBPASE','B-014','B-015','B-016','B-018','B-027','B-028','B-029'}
keep={'BDGK-GPR','GAP2','CA2tpp','NA1tpp','B-012','ACCOAC','CTPS2','GLBRAN1','GLCS1','GLDBRAN2','MTRP'}
assert len(groups)==45,(len(groups),list(groups))
assert len(evidence)==11 and len(keep)==11

reaction={
 'B-001':'','B-001-MDH-GPR':'MDH','B-002':'TKT1','B-003':'PPC','B-004':'','BDGK-GPR':'BDGK','BDGK-方向':'BDGK','GAP2':'GAP2','GLCP':'GLCP','SBPASE':'SBPASE',
 'B-005':'','B-006':'','B-007':'','B-008':'','B-009':'','B-010':'','B-011':'','CA2tpp':'CA2tpp','NA1tpp':'NA1tpp',
 'B-012':'TSQOC','B-013':'','B-014':'','B-015':'','B-016':'','B-017':'','B-018':'',
 'B-019':'','B-020':'','B-021':'ATPS5rpp','B-022':'','B-023':'','B-024':'ATPS5rpp','B-025':'NADHI','B-026':'ATPS5rpp','B-027':'','B-028':'','B-029':'','B-030':'',
 'ACCOAC':'ACCOAC','生物量反应':'Afe_biomass_mc507_WT_139p0M','CTPS2':'CTPS2','GLBRAN1':'GLBRAN1','GLCS1':'GLCS1','GLDBRAN2':'GLDBRAN2','MTRP':'MTRP'
}
doi={
 'B-001':'10.3389/fmicb.2019.00592','B-001-MDH-GPR':'10.3389/fmicb.2019.00592','B-002':'10.1128/AEM.03057-12','B-003':'10.1186/1471-2164-9-597','B-004':'10.3389/fmicb.2016.01365','BDGK-GPR':'','BDGK-方向':'','GAP2':'','GLCP':'10.3389/fmicb.2016.01365','SBPASE':'10.3389/fmicb.2019.00592',
 'B-005':'10.1186/1471-2164-9-597','B-006':'10.1186/1471-2164-9-597','B-007':'10.1186/1471-2164-9-597','B-008':'10.1186/1471-2164-9-597','B-009':'10.1186/1471-2164-9-597','B-010':'10.3389/fmicb.2019.00592','B-011':'10.3389/fmicb.2019.00592','CA2tpp':'10.1186/1471-2164-9-597','NA1tpp':'',
 'B-012':'10.1128/JB.01472-13','B-013':'10.1186/1471-2164-9-597','B-014':'10.1007/s00253-014-5830-4','B-015':'10.1007/s00284-018-1453-9','B-016':'10.1007/s00284-018-1453-9','B-017':'10.1186/1471-2164-9-597','B-018':'10.3389/fmicb.2019.00592',
 'B-019':'10.1186/1471-2164-9-597','B-020':'10.1186/1471-2164-9-597','B-021':'','B-022':'10.1186/1471-2164-9-597','B-023':'10.1186/1471-2164-9-597','B-024':'','B-025':'','B-026':'','B-027':'10.1016/j.hydromet.2006.03.030','B-028':'10.3389/fmicb.2019.00592','B-029':'10.1128/AEM.03057-12','B-030':'10.1186/1471-2164-9-597',
 'ACCOAC':'10.1186/1471-2164-9-597','生物量反应':'','CTPS2':'10.1186/1471-2164-9-597','GLBRAN1':'10.1016/j.hydromet.2006.03.029','GLCS1':'10.1016/j.hydromet.2006.03.029','GLDBRAN2':'10.1016/j.hydromet.2006.03.029','MTRP':'10.1016/j.hydromet.2006.03.029'
}

def locus_for(cid):
    # Authoritative final-list tables are in sections III-V.
    pre=text.split('## 六、')[0]
    for line in pre.splitlines():
        if line.startswith('| '+cid+' /') or line.startswith('| '+cid+' |'):
            found=re.findall(r'RU820_RS\d+',line)
            return '; '.join(dict.fromkeys(found))
    return {'ACCOAC':'RU820_RS06765'}.get(cid,'')

def legacy_for(cid,object_text):
    pre=text.split('## 六、')[0]
    for line in pre.splitlines():
        if line.startswith('| '+cid+' /') or line.startswith('| '+cid+' |'):
            found=re.findall(r'AFE_\d+',line)
            if found:return '; '.join(dict.fromkeys(found))
    return {'BDGK-方向':'AFE_2841'}.get(cid,'; '.join(dict.fromkeys(re.findall(r'AFE_\d+',object_text))))

rows=[]
for cid,item in groups.items():
    score_match=re.search(r'(?<!\d)([1-4])(?!\d)',item['score_text'])
    score=int(score_match.group(1)) if score_match and '不评分' not in item['score_text'] else None
    core=byrid.get(reaction[cid],[None]*16)
    if cid=='B-014': core=[None]*16 # SDO assay is chemically distinct from legacy SULDO.
    action='Evidence update only' if cid in evidence else 'Keep' if cid in keep else 'Hold'
    if score==4: etype='biochemical characterization'
    elif score==3: etype='physiological evidence'
    elif score==2: etype='indirect evidence or sequence homology'
    elif score==1: etype='gap-filling / network functionality'
    else: etype='unscored — no defined reaction'
    old_loci=legacy_for(cid,item['object'])
    rows.append({'id':cid,'group':item['group'],'core':core,'linked_reaction':reaction[cid] or ('SULDO' if cid=='B-014' else ''),'object':item['object'].replace('**',''),'legacy_locus':old_loci,'locus':locus_for(cid),'doi':doi[cid],'evidence_type':etype,'boundary':item['boundary'].replace('**',''),'action':action,'score':score,'source_summary':item['source_summary'],'score_text':item['score_text']})

out={'headers':headers,'group_titles':[original[0][10],original[0][12],original[0][14]],'original':records,'candidates':rows}

source_meta={
'10.3389/fmicb.2019.00592':('Experimental proteomics','DSM 14882T = ATCC 23270','Comparative Fe(II) and pyrite biofilm proteomics','AFE_0660, AFE_1667, AFE_1799, AFE_3253, AFE_3112, AFE_2554, AFE_2836, AFE_2678, AFE_1462','Protein detection and relative abundance','2','Protein abundance does not establish specific catalysis, GPR membership, direction or stoichiometry.'),
'10.1128/AEM.03057-12':('Experimental transcriptomics / proteomics','ATCC 23270','S/O2 versus S/Fe3+ protein spots, microarrays and qPCR','AFE_1667, AFE_2553, AFE_2586, main ATP synthase cluster','Protein spots and RNA abundance under two conditions','2; 3 only for legacy core physiology','Main ATP synthase expression does not prove extra subunit replacement or H+/ATP ratio.'),
'10.1186/1471-2164-9-597':('Genome sequencing and computational annotation','ATCC 23270','Genome sequencing, annotation and pathway reconstruction','AFE_1883 and transport, redox and biosynthesis candidates','Genome sequence and inferred homology/pathways','2','No candidate-specific substrate, transport, enzyme kinetics or exact reaction was measured.'),
'10.3389/fmicb.2016.01365':('Experimental transcriptomics','ATCC 23270','Quorum-sensing analog exposure; RNA comparisons','AFE_2024, AFE_1799','RNA abundance changes','2','Transcript abundance does not prove 6PG oxidative decarboxylation or glycogen reaction stoichiometry.'),
'10.1007/s00253-014-5830-4':('Biochemical characterization / genetic perturbation','ATCC 23270; A. caldus MTH-04','Recombinant SDO activity, sdo knockout and overexpression','AFE_0269 / Sdo','GSH-dependent SDO activity and S0-growth phenotype','4 for measured GSH-dependent SDO','Does not validate legacy SULDO elemental-sulfur substrate, periplasmic location or reversible direction.'),
'10.1007/s00284-018-1453-9':('Experimental transcriptomics','ATCC 23270','Fe2+ versus S0 expression and sequence analysis','AFE_1428, AFE_1429','Condition-dependent RNA expression','2','No purified electron transfer, acceptor or physical complex measurement.'),
'10.1128/JB.01472-13':('Genetic perturbation and transcriptomics','ATCC 23270','tetH knockout / overexpression, growth and qPCR','tetH, AFE_0048/doxD1','tetH phenotype and doxD1 transcription','2 for AFE_0048/TSQOC association','tetH perturbation does not demonstrate AFE_0048 catalysis or q8 as TSQOC acceptor.'),
'10.1074/jbc.M115.657551':('Whole-cell physiological spectroscopy','ATCC 23270','Iron respiratory-chain redox kinetics','Main respiratory chain; NADHI context','Whole-cell cytochrome and rusticyanin redox changes','3 for core respiration','Does not assign extra NuoB AFE_2411 or quantify NADHI q8/H+ stoichiometry.'),
'10.1016/j.hydromet.2006.03.030':('Experimental transcriptomics','ATCC 23270','Fe2+ versus elemental sulfur microarrays','AFE_2553/HdrA and respiratory pathway genes','Condition-dependent transcript changes','2','No Hdr substrate, acceptor or coupling chemistry measured.'),
'10.1016/j.hydromet.2006.03.029':('Experimental transcriptomics','ATCC 23270','Fe(II) versus sulfur microarray and RT-PCR','glgA/AFE_2678, glgB1/AFE_2836, glgP2/AFE_0527','Gene transcription ratios','2','Does not determine glycogen chain length, proton balance or reaction direction.'),
'10.1016/j.meteno.2016.03.003':('Model reconstruction','ATCC 23270 model','Model reconstruction, fitting and supplementary reaction table','iMC507 model and biomass objective','Model predictions and fits to prior physiological data','1-4 original model scores; 3 legacy BOF','Model coefficients and growth fit are not direct strain-specific glycogen or proton-translocation measurements.'),
'10.1128/AEM.07230-11':('Genetic knockout / physiology','ATCC 23270','AFE_1807/pfkB knockout and glucose-culture phenotype','AFE_1807/pfkB','Mutant growth and glucose use','3 for pfkB-linked phenotype only','Does not support AFE_2841 BDGK GPR, BDGK reversibility or candidate transport proteins.'),
'10.1128/AEM.02251-12':('Biochemical characterization','ATCC 23270','AFE_0042 thiosulfate dehydrogenase assay','AFE_0042','Thiosulfate to tetrathionate with ferricyanide; tested quinone negative','4 for measured AFE_0042 assay','Cannot transfer to AFE_0048 or legacy q8-dependent TSQOC equation.'),
'10.3389/fmicb.2019.00603':('Experimental physiology / qPCR','ATCC 23270','CO2 uptake and selected qPCR; pathway annotation','AFE_3253/cbbF predicted','CO2 assimilation and selected transcript measurements','2 for cbbF/SBPASE inference','cbbF SBPase substrate and direction were not assayed.'),
'10.1186/1471-2164-10-394':('Transcriptomics and pathway modeling','ATCC 23270','Sulfur/iron qPCR and electron-transfer pathway reconstruction','AFE_2836 and Fe/S pathway genes','RNA changes and modeled pathway links','2','Expression and pathway model do not establish specific reaction chemistry.'),
'10.1074/mcp.M700042-MCP200':('Experimental proteomics','A. ferrooxidans; exact culture collection not established here','Periplasmic proteome identification','Sdo-like periplasmic context','Detected periplasmic proteins','2 at most for legacy SULDO','Protein detection does not verify elemental sulfur substrate, direction or AFE_0269-specific old reaction.'),
'10.1099/mic.0.26212-0':('Biochemical characterization','A. ferrooxidans R1 and other strains; not ATCC 23270','Persulfide/sulfane sulfur enzyme substrate tests','Sulfur dioxygenase activity in other strains','Persulfide sulfur as assay substrate','4 for assayed enzyme in tested strains only','R1 result cannot certify the ATCC 23270 legacy SULDO formula or locus.'),
'10.1016/j.resmic.2023.104115':('Model study','ATCC 23270-derived model','In silico compatible-solute flux analysis and model update','2024 model differences, especially glycogen and reaction bounds','Model revision and simulated flux','Not a new experimental score','No reported wet-lab validation of seven 05-group coefficient/direction changes.')
}
sources=[]
for idx,(d,meta) in enumerate(source_meta.items(),1):
    with urllib.request.urlopen('https://api.crossref.org/works/'+d,timeout=20) as response:
        cr=json.load(response)['message']
    authors=cr.get('author',[])
    author_text=', '.join((a.get('given','')+' '+a.get('family','')).strip() for a in authors)
    if d=='10.1128/AEM.03057-12': author_text='Héctor Osorio, Stefanie Mangold, Yann Denis, Ivan Ñancucheo, Mario Esparza, D. Barrie Johnson, Violaine Bonnefoy, Mark Dopson, David S. Holmes'
    year=cr.get('published',cr.get('issued',{})).get('date-parts',[[None]])[0][0]
    citation=f"{author_text} ({year}). {cr.get('title',[''])[0]}. {cr.get('container-title',[''])[0]} {cr.get('volume','')}{'('+cr.get('issue','')+')' if cr.get('issue') else ''}, {cr.get('page',cr.get('article-number',''))}. https://doi.org/{d}"
    linked=[r['id'] for r in rows if r['doi'].lower()==d.lower()]
    extra={'10.1128/AEM.03057-12':['B-021','B-026'], '10.1074/jbc.M115.657551':['B-025'], '10.1016/j.meteno.2016.03.003':['B-021','B-024','B-025','B-026','生物量反应','BDGK-GPR','GAP2','NA1tpp'], '10.1128/AEM.07230-11':['BDGK-GPR','BDGK-方向','B-006','B-007'], '10.1128/AEM.02251-12':['B-012'], '10.3389/fmicb.2019.00603':['SBPASE'], '10.1186/1471-2164-10-394':['B-025','GLBRAN1'], '10.1074/mcp.M700042-MCP200':['B-014'], '10.1099/mic.0.26212-0':['B-014'], '10.1016/j.resmic.2023.104115':['ACCOAC','CTPS2','GLBRAN1','GLCS1','GLDBRAN2','MTRP','生物量反应']}.get(d,[])
    linked=list(dict.fromkeys(linked+extra))
    sources.append({'id':f'S{idx:02d}','doi':d,'citation':citation,'year':year,'source_type':meta[0],'strain':meta[1],'experimental_type':meta[2],'genes':meta[3],'measured':meta[4],'applicable_confidence':meta[5],'limitation':meta[6],'related':'; '.join(linked)})
out['sources']=sources
(root/'evidence_data.json').write_text(json.dumps(out,ensure_ascii=False),encoding='utf-8')
print('original',len(records),'candidates',len(rows),'scores',{n:sum(r['score']==n for r in rows) for n in [4,3,2,1]},'unscored',sum(r['score'] is None for r in rows),'actions',{a:sum(r['action']==a for r in rows) for a in ['Evidence update only','Keep','Hold']})
print('missing reaction IDs',[r['id'] for r in rows if r['linked_reaction'] and r['linked_reaction'] not in byrid])
print('unique chosen DOIs',len(set(r['doi'].lower() for r in rows if r['doi'])))
print('source papers',len(sources),'original experimental',sum(s['source_type'] not in {'Genome sequencing and computational annotation','Model reconstruction','Model study'} for s in sources))
