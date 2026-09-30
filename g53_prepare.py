import csv,json,re,collections
from urllib.parse import unquote
from pathlib import Path
from openpyxl import load_workbook

R=Path('.')
def tsv(name):
    with open(name,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
genes=[]
for line in (R/'G53.txt').read_text(encoding='utf-8-sig').splitlines():
    m=re.match(r'^(\d{4})｜(RU820_RS\d+)｜(WP_[\d.]+)｜NZ_CP136162.1:(\d+)-(\d+)\(([+-])\)｜(.+)$',line)
    if m:genes.append(dict(zip(['seq','locus','wp','start','end','strand','product'],m.groups())))
assert len(genes)==293 and [int(g['seq']) for g in genes]==list(range(587,880))
ids={g['locus'] for g in genes}
gff={}
for line in (R/'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff').read_text(encoding='utf-8').splitlines():
    p=line.split('\t')
    if len(p)!=9 or p[0]!='NZ_CP136162.1' or p[2] not in ('gene','CDS'):continue
    a=dict(x.split('=',1) for x in p[8].split(';') if '=' in x)
    locus=a.get('locus_tag')
    if locus in ids:gff.setdefault(locus,{})[p[2]]=(p,a)
for g in genes:
    pair=gff.get(g['locus'],{})
    gene=pair.get('gene');cds=pair.get('CDS')
    g['gff_check']='一致' if gene and cds and gene[0][3:5]==[g['start'],g['end']] and gene[0][6]==g['strand'] and cds[1].get('protein_id')==g['wp'] and unquote(cds[1].get('product',''))==g['product'] else '差异待核'
    g['pseudogene']=bool(gene and ('pseudo' in gene[1] or 'pseudogene' in gene[1]))
cross={r['当前locus']:r for r in tsv(R/'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv') if r['当前locus'] in ids}
maps={}
for r in tsv(R/'B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv'):
    if r['当前GCF049位点'] in ids:maps.setdefault(r['当前GCF049位点'],[]).append(r)
alpha=load_workbook(R/'00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx',read_only=True,data_only=True)
original={}
for row in list(alpha.worksheets[0].values)[2:]:
    if row[16]!='2016 original':continue
    for afe in set(re.findall(r'AFE_\d+',str(row[7] or ''))):original.setdefault(afe,[]).append(list(row[:16]))
uniprot_by_wp={}
uniprot_by_acc={}
for fn in ['uniprot_243159.tsv','uniprot_missing.tsv']:
    p=R/'outputs/G53'/fn
    if not p.exists():continue
    with p.open(encoding='utf-8',newline='') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            uniprot_by_acc[row['Entry']]=row
            for wp in re.findall(r'WP_[0-9]+[.][0-9]+',row.get('RefSeq','')):
                uniprot_by_wp.setdefault(wp,row)
uniprot_extra={}
p=R/'outputs/G53/uniprot_243159_extra.tsv'
if p.exists():
    with p.open(encoding='utf-8',newline='') as f:
        uniprot_extra={r['Entry']:r for r in csv.DictReader(f,delimiter='\t')}
uniprot_by_afe={}
p=R/'outputs/G53/uniprot_243159_locus.tsv'
if p.exists():
    with p.open(encoding='utf-8',newline='') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            for afe in re.findall(r'AFE_[0-9]+',row.get('Gene Names (ordered locus)','')):
                uniprot_by_afe.setdefault(afe,row)
rhea_rows={}
p=R/'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/Rhea_全基因EC查询汇总.tsv'
with p.open(encoding='utf-8-sig',newline='') as f:
    rhea_rows={r['Reaction identifier']:r for r in csv.DictReader(f,delimiter='\t')}
bellenberg_pairs=set()
p=R/'outputs/G53/bellenberg2019_supp.txt'
if p.exists():
    bellenberg_pairs={tuple(x.split()) for x in re.findall(r'AFE_[0-9]+\s+WP_[0-9]+[.][0-9]+',p.read_text(encoding='utf-8'))}
epmc_afe={}
epmc_wp={}
for fn,target in [('epmc_afe_scan.json',epmc_afe),('epmc_wp_scan.json',epmc_wp)]:
    p=R/'outputs/G53'/fn
    if p.exists():
        target.update(json.loads(p.read_text(encoding='utf-8')).get('queries',{}))
index=[]
review=[]
sources=[
['G53-S01','NCBI RefSeq genomic.gff, GCF_049532655.1, AFEATCC23270_v5.2',2026,'','Official genome annotation','ATCC 23270','GFF gene/CDS identity and product','G53 all','G53 all','Locus, protein ID, coordinates, current product','Not a reaction score','Functional product labels are annotations; https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/'],
['G53-S02','B1 verified old/new gene mapping, local 旧新基因映射.tsv',2026,'','Sequence and protein ID mapping','ATCC 23270','Exact RefSeq protein ID / B1 mapping','G53 mapped subset','G53 all','Historical AFE↔RU820 links','Not a reaction score','One-to-many/absent mappings remain unresolved; local B1 file'],
['G53-S03','Campodonico et al. 2016, iMC507 original supplementary mmc1.xls and Alpha frozen copy',2016,'10.1016/j.meteno.2016.03.003','Model reconstruction','ATCC 23270 model','Reaction table and GPR','G53 mapped AFE genes','2016 reaction IDs','Original formula, confidence, GPR and media bounds','Original score only','Model scoring does not prove gene-specific catalysis; https://doi.org/10.1016/j.meteno.2016.03.003'],
['G53-S04','D3-1-R full-genome cross-database table, 2026-09-23',2026,'','Computational cross-reference integration','ATCC 23270 / legacy KEGG afr','KEGG, BioCyc, BRENDA, Rhea cross-references','G53 all','G53 all','Database identifiers and coverage flags','Usually indirect; not assigned automatically','Cross-references may inherit each other; local D3-1-R_全基因跨数据库交叉引用.tsv'],
['G53-S05','Valdés et al. 2008, Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications',2008,'10.1186/1471-2164-9-597','Genome annotation','ATCC 23270','Genome sequencing and computational pathway prediction','Selected G53 genes','Selected pathways','Annotation and pathway proposals','2 only for matched reactions','Genome prediction does not prove protein-specific reaction; https://doi.org/10.1186/1471-2164-9-597'],
['G53-S06','Bellenberg et al. 2019, Proteomics Reveal Enhanced Oxidative Stress Responses...',2019,'10.3389/fmicb.2019.00592','Proteomics','DSM 14882T = ATCC 23270','Fe(II)/pyrite biofilm proteomics','AFE_0660','MDH, B-001','AFE_0660 protein detection','Expression only','No AFE_0660 substrate, cofactor or product measurement; https://doi.org/10.3389/fmicb.2019.00592'],
['G53-S07','Osorio et al. 2013, Anaerobic Sulfur Metabolism Coupled to Dissimilatory Iron Reduction...',2013,'10.1128/AEM.03057-12','Transcriptomics/proteomics','ATCC 23270','Expression by growth condition','ATP synthase cluster','ATPS5rpp, B-021','Expression of main ATP synthase cluster','Expression only','No AFE_0758 epsilon-subunit replacement or H+/ATP stoichiometry; https://doi.org/10.1128/AEM.03057-12'],
['G53-S08','Quatrini et al. 2009, Extending the models for iron and sulfur oxidation...',2009,'10.1186/1471-2164-10-394','Pathway model / expression','ATCC 23270 and cited other strains','Literature synthesis and expression','AFE_0631-0634','CYTBO3','bo3 terminal oxidase is discussed','Indirect','Reported bo3 experimental basis includes ATCC 19859; strain attribution must stay explicit; https://doi.org/10.1186/1471-2164-10-394'],
['G53-S09','UniProtKB REST organism 243159 and exact RefSeq WP queries, retrieved 2026-09-24/25',2026,'','Protein database','ATCC 23270 and shared RefSeq proteins','Reviewed status, protein name, InterPro/Pfam, EC and cross-references','G53 exact WP subset','G53 all','UniProt accession, review state and family IDs','Annotation only','Old AFE-only links are marked separately and cannot prove current sequence equivalence; https://rest.uniprot.org/uniprotkb/'],
['G53-S10','Bellenberg et al. 2019, original supplementary Data Sheet 2',2019,'10.3389/fmicb.2019.00592','Proteomics supplementary data','DSM 14882T = ATCC 23270','Fe(II)/pyrite biofilm proteomics','31 G53 AFE/WP exact pairs','G53 selected genes','Selected proteins detected and pyrite/Fe abundance comparisons','Expression only','No isolated enzyme activity or reaction stoichiometry; https://www.frontiersin.org/api/v4/articles/434788/file/Data_Sheet_2.pdf/434788_supplementary-materials_datasheets_2_pdf/2'],
['G53-S11','Drobner et al. 1990, Hydrogen oxidation by Acidithiobacillus ferrooxidans ATCC 23270',1990,'10.1128/aem.56.9.2922-2923.1990','Whole-strain physiology','ATCC 23270','Hydrogen-grown culture and hydrogenase induction','No individual G53 locus assigned','HYD candidates','Strain-level H2 utilization','No exact reaction score','No Hup/Hox assignment, electron acceptor or proton coupling; https://doi.org/10.1128/aem.56.9.2922-2923.1990'],
['G53-S12','Almarcegui et al. 2014, ATCC 23270 copper-stress proteomics',2014,'10.1016/j.resmic.2014.07.005','Proteomics','ATCC 23270','Copper-stress protein abundance','AFE_0105 context, not AFE_0920 flux','CU2tpp/MN2tpp/MNTH/ZN2tpp','Stress-associated protein change','Expression only','No direct metal substrate/transport or H+ coupling measurement for AFE_0920; https://doi.org/10.1016/j.resmic.2014.07.005'],
['G53-S13','IUBMB EC 1.17.1.1 official nomenclature',2026,'','Enzyme nomenclature','General enzyme class','Curated reaction definition','No G53 protein experimentally assigned','AFE_0860 conflict','Defines CDP-deoxyglucose reductase chemistry','Not a gene-specific score','Shows D3 EC/KO association cannot identify 2Fe-2S protein function; https://iubmb.qmul.ac.uk/enzyme/EC1/17/1/1.html'],
['G53-S14','Europe PMC accession text search, retrieved 2026-09-25',2026,'','Literature search index','All strains; individually unverified','AFE or WP identifier query','G53 all','G53 all','Publication hits and DOI candidates','Screening only; no reaction score','Search coverage is limited to indexed text; zero hits does not prove no paper. Each hit needs strain/object/experiment review; https://www.ebi.ac.uk/europepmc/webservices/rest/search'],
['G53-S15','Rhea full-gene EC query table with ChEBI, KEGG and MetaCyc cross-references, 2026-09-23',2026,'','Reaction database','General reaction definitions','Generic chemical equation and cross-reference','82 G53 genes with D3 exact cross-reference','119 gene–Rhea pairs','Reaction equation, EC, ChEBI and external IDs','Not gene-specific; no automatic score','The D3 cross-reference is not a protein-specific assay or compartment assignment; https://www.rhea-db.org/'],
['G53-S16','Zeng et al. 2007, Expression, purification and characterization of IscU from Acidithiobacillus ferrooxidans',2007,'10.1007/s10529-007-9488-1','Original biochemical study','A. ferrooxidans; strain not stated in abstract','Recombinant IscU Fe-S assembly and mutagenesis','IscU, IscS and IscA','CYSDSS context','[2Fe-2S] assembly on IscU in vitro','Does not alone establish model CYSDSS exact step or four-gene AND','2016 CYSDSS PMID 17660945 is this IscU paper; https://pubmed.ncbi.nlm.nih.gov/17660945/'],
['G53-S17','Zeng et al. 2007, Expression, purification and characterization of IscS from Acidithiobacillus ferrooxidans',2007,'10.1007/s10529-007-9491-6','Original biochemical study','ATCC 23270 source gene; recombinant expression in E. coli','Purified IscS cysteine desulfurase assay','IscS','CYSDSS context','L-cysteine yielded L-alanine and elemental sulfur or H2S, depending on reducing agent','Different measured sulfur product from model-bound persulfide; PMID 17660944 is not the 2016 cited PMID; https://pubmed.ncbi.nlm.nih.gov/17660944/'],
['G53-S18','Suzuki and Knaff 2005, Glutamate synthase: structural, mechanistic and regulatory properties',2005,'10.1007/s11120-004-3478-0','Review','Multiple organisms, predominantly plants','Literature synthesis','No ATCC 23270 GltB/GltD assay','GLUSy','Review of glutamate synthase classes','Cannot support 2016 GLUSy score 4 as same-strain direct biochemistry; https://pubmed.ncbi.nlm.nih.gov/16143852/'],
['G53-S19','Selkov et al. 2000, Functional analysis of gapped microbial genomes: amino acid metabolism of Thiobacillus ferrooxidans',2000,'10.1073/pnas.97.7.3509','Computational genome study','ATCC 23270','Gapped genome annotation and pathway reconstruction','Predicted AroB context','DHQS','Sequence and pathway-based gene assignment','No purified DHQS enzyme measurement; https://pubmed.ncbi.nlm.nih.gov/10737802/'],
['G53-S20','Islam et al. 2020, A widely distributed hydrogenase oxidises atmospheric H2 during bacterial growth',2020,'10.1038/s41396-020-0713-4','Original physiology and qRT-PCR','DSM 14882 = ATCC 23270','AFE_0702 qRT-PCR, whole-cell gas chromatography and H2-supported growth','AFE_0702; group 2a [NiFe] hydrogenase','Hup candidate outside 2016 iMC507','AFE_0702 transcript and strain-level atmospheric H2 uptake','Supports H2 oxidation role, but no gene knockout, purified enzyme, electron acceptor, proton stoichiometry or exclusive attribution; other hydrogenases expressed; https://pmc.ncbi.nlm.nih.gov/articles/PMC7784904/'],
]
for g in genes:
    locus=g['locus'];c=cross.get(locus,{});mm=maps.get(locus,[])
    afes=sorted({r['旧AFE位点'] for r in mm if r['旧AFE位点']})
    map_status='已核验' if len(afes)==1 and all(r['映射状态']=='完全匹配' for r in mm) else ('歧义' if len(afes)>1 else '未确认')
    up=uniprot_by_wp.get(g['wp'])
    up_basis='exact WP'
    if not up and len(afes)==1 and afes[0] in uniprot_by_afe:
        old=uniprot_by_afe[afes[0]]
        up=uniprot_by_acc.get(old['Entry'])
        up_basis='legacy AFE only; WP mismatch'
    if not up:up_basis='未找到'
    ux=uniprot_extra.get(up['Entry'],{}) if up else {}
    proteomics=any((afe,g['wp']) in bellenberg_pairs for afe in afes)
    epmc=epmc_afe.get(afes[0],{}) if len(afes)==1 else epmc_wp.get(g['wp'],{})
    epmc_dois=sorted({x.get('doi','') for x in epmc.get('results',[]) if x.get('doi')})
    modelrows=[]
    for afe in afes:
        modelrows+=original.get(afe,[])
    seen=set();modelrows=[x for x in modelrows if not (x[0] in seen or seen.add(x[0]))]
    rhea_ids=sorted(set(re.findall(r'RHEA:[0-9]+',c.get('Rhea精确交叉引用',''))))
    if locus=='RU820_RS04125':rhea_ids=[]  # EC/KO assignment conflicts with the current 2Fe–2S annotation.
    keg=c.get('KEGG reaction','');bio=c.get('BioCyc reaction ID','')
    exact='候选，未核实反应等价' if (keg or bio) else '未找到'
    if modelrows: exact='2016 已有公式；基因专属等价未全面复核'
    if locus=='RU820_RS03065':exact='PHFT 为旧模型候选关联；两个旧 AFE 共用 WP，位点对应未确认'
    if locus in {'RU820_RS03380','RU820_RS03385'}:exact='KEGG R08034 与氢化酶化学不符；不得作候选反应'
    if locus=='RU820_RS04125':exact='KEGG K00523/EC 1.17.1.1 疑似错配；不作位点反应'
    nonmetabolic=('transcription','ribosomal','restriction','recombinase','integrase','transposase','cell division','pilus','toxin-antitoxin','addiction module','chaperone','elongation factor','chromosome segregation','integration host factor','regulatory protein','dna-binding','excinuclease','endonuclease','dna gyrase','dna repair','hsp20','stress protein')
    if modelrows:action='仅更新证据或编号（待确认）'
    elif not (keg or bio) and any(word in g['product'].lower() for word in nonmetabolic):action='无需修改（本轮未见代谢反应）'
    else:action='证据不足暂缓'
    issue=[]
    if map_status!='已核验':issue.append('旧 AFE 映射'+map_status)
    if map_status!='已核验' and c.get('旧AFE locus'):issue.append('D3 历史候选 '+c['旧AFE locus']+'；不得直接采用')
    if g['wp']=='WP_012606582.1':issue.append('WP_012606582.1 在本组出现两份转座酶拷贝；位点分别核查')
    if not modelrows and (keg or bio):issue.append('数据库反应与 2016 公式/区室未作精确等价核验')
    if not modelrows and not (keg or bio):issue.append('未找到精确反应')
    if c.get('映射冲突')=='是':issue.append('历史映射冲突')
    if g['gff_check']!='一致':issue.append('同版 GFF 身份差异')
    if locus=='RU820_RS03190':
        action='证据不足暂缓';issue.append('B-001 苹果酸酶及 MDH OR 均缺基因专属酶学')
    if locus=='RU820_RS03655':
        action='证据不足暂缓';issue.append('B-021 不能增入主 ATP 合酶 GPR')
    if locus=='RU820_RS04125':
        action='证据不足暂缓';issue.append('B-022 电子供体/受体不明')
    if locus=='RU820_RS04055':issue.append('CTPS2 方向缺本株双向实验证据')
    if locus=='RU820_RS03065':issue.append('PHFT 的 AFE_0635/AFE_3143 两旧位点共用 WP；不可单独指派 RU820_RS03065')
    if locus in {'RU820_RS03380','RU820_RS03385'}:issue.append('R08034 实为非氢化酶反应；Hup 电子受体和质子耦联未明')
    if locus=='RU820_RS03385':issue.append('Islam 2020 在同株测得 AFE_0702 转录和全细胞 H2 摄取；有其他氢化酶共同表达，未测基因专属催化或敲除')
    if locus=='RU820_RS03345':issue.append('FDH 的 NAD 受体与产物未由本株位点专属实验确定；保留旧模型假设')
    if locus=='RU820_RS04125':issue.append('K00523/EC 1.17.1.1 对应脱氧糖还原，与广义 2Fe–2S 注释冲突')
    if locus=='RU820_RS04425':issue.append('NRAMP 四种模型底物、OR GPR 和 MNTH 的 H+ 耦联未在本株测定')
    if locus in {'RU820_RS04505','RU820_RS04510','RU820_RS04515','RU820_RS04520'}:issue.append('HYD3pp 的 Q8 受体、AND GPR 及 2 H+ 跨膜计量未在本株精确验证')
    if locus in {'RU820_RS03045','RU820_RS03050','RU820_RS03055','RU820_RS03060'}:issue.append('CYTBO3 亚基身份及 1.8 H+ 系数待核；异株实验证据不能移用')
    if locus in {'RU820_RS03255','RU820_RS03260'}:issue.append('CYSDSS 原引 PMID 17660945 是 IscU 装配研究；IscS 活性研究 PMID 17660944 测得产物与模型式不同；四基因 AND 待 QC')
    if locus in {'RU820_RS03520','RU820_RS03525'}:issue.append('GLUSy 原引 PMID 16143852 是综述，不能支持原分数 4 的同株直接酶学')
    if locus=='RU820_RS03540':issue.append('DHQS 原引 PMID 10737802 为缺口基因组推断，无直接酶学')
    if up_basis=='legacy AFE only; WP mismatch':issue.append('UniProt 仅旧 AFE 对应，当前 WP 未精确对应')
    if up_basis=='未找到':issue.append('UniProt 未命中当前 WP')
    if epmc.get('hitCount',0):issue.append('Europe PMC DOI 命中仅为检索线索；逐篇菌株/实验待核')
    refs=['G53-S01','G53-S04']
    if mm:refs.append('G53-S02')
    if modelrows:refs.append('G53-S03')
    if locus=='RU820_RS03190':refs+=['G53-S05','G53-S06']
    if locus=='RU820_RS03655':refs.append('G53-S07')
    if locus in {'RU820_RS03045','RU820_RS03050','RU820_RS03055','RU820_RS03060'}:refs.append('G53-S08')
    if locus in {'RU820_RS03380','RU820_RS03385','RU820_RS04505','RU820_RS04510','RU820_RS04515','RU820_RS04520'}:refs.append('G53-S11')
    if locus=='RU820_RS04425':refs.append('G53-S12')
    if locus=='RU820_RS04125':refs.append('G53-S13')
    if up:refs.append('G53-S09')
    if locus in {'RU820_RS03255','RU820_RS03260'}:refs+=['G53-S16','G53-S17']
    if locus in {'RU820_RS03380','RU820_RS03385'}:refs.append('G53-S20')
    if locus in {'RU820_RS03520','RU820_RS03525'}:refs.append('G53-S18')
    if locus=='RU820_RS03540':refs.append('G53-S19')
    if proteomics:refs.append('G53-S10')
    refs.append('G53-S14')
    rid='G53-G'+g['seq']
    if locus=='RU820_RS04500':
        action='无需修改（本轮未见代谢反应）'
        issue.append('氢化酶成熟加工辅助蛋白；当前无可定义代谢反应')
    linked=';'.join(x[0] for x in modelrows)
    if locus=='RU820_RS03065':linked='PHFT（历史歧义；非已核验 GPR）'
    related=[rid]+[rid+'-'+str(x[0]) for x in modelrows]+[rid+'-'+x.replace(':','') for x in rhea_ids]
    if locus=='RU820_RS03065':related.append(rid+'-PHFT-ambiguous')
    if rhea_ids:refs.append('G53-S15')
    transporter=any(word in g['product'].lower() for word in ('transporter','permease','efflux','symporter','antiporter','channel','porin'))
    tcdb=ux.get('TCDB','') or ('未取得 WP 对应 TCDB ID' if transporter else '')
    index.append([int(g['seq']),locus,g['wp'],g['start'],g['end'],g['strand'],g['product'],';'.join(afes) or '未确认',map_status,g['gff_check'],c.get('KEGG KO',''),c.get('KEGG EC',''),keg,bio,c.get('Rhea精确交叉引用',''),c.get('Rhea仅EC候选',''),c.get('2016 GPR反应',''),exact,action,';'.join(related),linked,'; '.join(issue),';'.join(refs),up.get('Entry','') if up else '',up_basis,up.get('Reviewed','') if up else '',up.get('Protein names','') if up else '',up.get('InterPro','') if up else '',up.get('Pfam','') if up else '',tcdb,ux.get('EC number',''), '同株蛋白检出；不证明反应' if proteomics else '',epmc.get('hitCount',0),' ; '.join(epmc_dois),epmc.get('url',''),c.get('BRENDA EC记录','')])
    base=[None]*30
    gene_doi='10.1038/s41396-020-0713-4' if locus=='RU820_RS03385' else ('10.3389/fmicb.2019.00592' if proteomics else '')
    gene_boundary='直接证实：同版 gene/CDS 身份'+('与同株蛋白检出' if proteomics else '')+'；尚未证实：注释所指精确酶活与新 GPR。'
    if locus=='RU820_RS03385':gene_boundary='直接证实：同株 AFE_0702 转录及全细胞低浓度 H2 摄取。尚未证实：该蛋白独立催化、电子受体、膜侧、质子计量及新模型 GPR。'
    base[16:30]=['gene-review',rid,'G53',linked,'',';'.join(afes) or '未确认',locus,gene_doi,'NCBI/B1/D3/2016 cross-reference'+(' + same-strain proteomics' if proteomics else '')+(' + qRT-PCR/whole-cell physiology' if locus=='RU820_RS03385' else ''),gene_boundary,action,'不评分',';'.join(refs),'基因审查记录；不转用 2016 反应分数']
    review.append(base)
    for x in modelrows:
        row=x+[None]*14
        pmid=int(x[5]) if x[5] else None
        doi={10737802:'10.1073/pnas.97.7.3509',16143852:'10.1007/s11120-004-3478-0',17660945:'10.1007/s10529-007-9488-1',19077236:'10.1186/1471-2164-9-597',19703284:'10.1186/1471-2164-10-394'}.get(pmid,'')
        score=1 if x[0]=='XAND' else 2
        boundary='直接证实：2016 原表存在该反应/GPR；当前蛋白家族和数据库交叉引用为间接证据。尚未证实：本株该位点催化模型精确反应及其方向、区室、辅因子或耦联。'
        note=f'2016 原分数 {x[3]}；本轮反应分数 {score}，仅限旧模型化学式的注释/序列支持，不传递给 GPR 或计量参数。'
        if x[0]=='CYSDSS':
            boundary='直接证实：ATCC 23270 来源 IscS 重组蛋白可将 L-半胱氨酸转为 L-丙氨酸及单质硫/H2S；IscU 体外装配 Fe-S。尚未证实：模型所写 IscS 结合过硫化物产物及四个基因必须同时存在。'
            note='2016 原分数 4；原引 PMID 17660945 测 IscU，PMID 17660944 测 IscS 但产物不同。按模型精确 CYSDSS 式暂评 2；四基因 AND 另列 QC。'
        elif x[0]=='GLUSy':
            boundary='直接证实：2016 模型有 GLUSy；当前 GltB/GltD 注释支持复合体候选。尚未证实：ATCC 23270 纯化酶的 NADPH 依赖反应。'
            note='2016 原分数 4；PMID 16143852 为综述，未测本株 GltB/GltD；本轮对模型式暂评 2。'
        elif x[0]=='DHQS':note='2016 原分数 2；PMID 10737802 为基因组途径推断；本轮按注释/序列暂评 2。'
        elif x[0]=='CYTBO3':note='2016 原分数 3；异株/通路研究不能支持 ATCC 23270 的 1.8 H+ 精确计量；本轮按亚基注释暂评 2。'
        elif x[0]=='XAND':note='2016 原分数 1；模型式仅保留网络假设，当前无位点专属精确反应证据；本轮维持 1。'
        row[16:30]=['reaction-review',rid+'-'+str(x[0]),'G53',x[0],x[2],';'.join(afes),locus,doi,'2016 model + annotation'+(' + PMID recheck' if pmid else ''),boundary,action,score,';'.join(refs),note]
        review.append(row)
    if locus=='RU820_RS03065':
        phft=next((x for x in original.get('AFE_0635',[]) if x[0]=='PHFT'),None)
        if phft:
            row=list(phft)+[None]*14
            row[16:30]=['2016-ambiguous-link',rid+'-PHFT-ambiguous','G53','PHFT',phft[2],'AFE_0635 or AFE_3143',locus,'','2016 model + shared-WP ambiguity','直接证实：2016 PHFT 存在；尚未确认：此当前位点与两个旧 AFE 中哪一个对应。','证据不足暂缓','不重新评分','G53-S01;G53-S03;G53-S04','2016 原分数 2 仅指 PHFT；不传递给当前位点映射']
            review.append(row)
    for rhea_id in rhea_ids:
        definition=rhea_rows.get(rhea_id,{})
        row=[None]*30
        row[0]=rhea_id
        row[2]=definition.get('Equation','')
        row[4]=definition.get('EC number','')
        row[16:30]=['database-reaction-candidate',rid+'-'+rhea_id.replace(':',''),'G53',';'.join(x[0] for x in modelrows),'Rhea generic equation; no model compartment',';'.join(afes) or '未确认',locus,'','D3 Rhea cross-reference; generic chemistry','直接证实：Rhea 定义通用反应；尚未证实：该 WP 在本株催化此反应、区室、方向与 GPR。','证据不足暂缓','不评分','G53-S04;G53-S15','数据库候选；不得继承 2016 分数']
        review.append(row)
paper_judgment={
 '10.1186/1471-2164-10-394':('ATCC 23270 通路/表达推断；所引 bo3 酶学含异株','间接；不证明 1.8 H+ 计量'),
 '10.3389/fmicb.2017.01277':('RegB/RegA 对 cyoA 调控研究','转录调控；不证明 CYTBO3 精确计量'),
 '10.3389/fmicb.2019.00592':('同株 DSM 14882 蛋白组；原补表已核对 AFE/WP','蛋白检出；不证明催化'),
 '10.1007/s12223-013-0244-8':('研究菌株 DC 的金属相关基因测序/生信','异株；不转成本株直接证据'),
 '10.3389/fmicb.2016.01365':('同株群体感应扰动后的转录组','表达；不证明精确反应'),
 '10.1111/1462-2920.15163':('群体感应调控研究','调控/表达；不证明精确反应'),
 '10.1038/s41396-020-0713-4':('同株 DSM 14882；AFE_0702 qRT-PCR 与全细胞 H2 摄取','支持 H2 氧化功能方向；无基因敲除、电子受体和计量'),
 '10.1074/mcp.m700042-mcp200':('周质蛋白组研究','定位/蛋白检出；不证明酶活'),
 '10.1128/aem.03057-12':('同株硫代谢条件转录/蛋白组','表达；不证明蛋白精确反应'),
 '10.1128/jb.00791-12':('AtzD/巴比妥酸酶家族的序列与产物研究','家族/异源对象；需核对本株蛋白实测'),
 '10.1371/journal.pone.0112226':('A. ferrooxidans 可移动基因组毒素-抗毒素研究','非代谢功能线索'),
 '10.1128/aem.01545-08':('Metallosphaera sedula 铁/硫电子传递转录组','异种；不作 ATCC 23270 直接证据'),
}
paper_by_doi=collections.defaultdict(list)
for r in index:
    for doi in filter(None,(r[33] or '').split(' ; ')):
        paper_by_doi[doi].append(str(r[0])+' '+r[1])
assert set(paper_by_doi)==set(paper_judgment)
literature=[[doi,len(paper_by_doi[doi]),'; '.join(paper_by_doi[doi]),*paper_judgment[doi],'https://doi.org/'+doi] for doi in sorted(paper_by_doi)]
data={'index':index,'review':review,'sources':sources,'literature':literature,'count':len(genes)}
out=R/'outputs/G53'
out.mkdir(parents=True,exist_ok=True)
(out/'g53_data.json').write_text(json.dumps(data,ensure_ascii=False),encoding='utf-8')
print(json.dumps({'genes':len(index),'review_rows':len(review),'model_reaction_rows':sum(r[16]=='reaction-review' for r in review),'rhea_candidate_rows':sum(r[16]=='database-reaction-candidate' for r in review),'gff_consistent':sum(g['gff_check']=='一致' for g in genes),'actions':dict(collections.Counter(r[18] for r in index))},ensure_ascii=False))
