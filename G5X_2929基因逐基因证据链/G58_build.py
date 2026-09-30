import csv
import json
import re
import urllib.parse
import urllib.request
from collections import defaultdict, Counter
from datetime import date
from pathlib import Path
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "G58_审查输出"
OUT.mkdir(exist_ok=True)
DOC = Path(__file__).resolve().parent / "G58_原始清单.txt"
GFF = ROOT / "01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/genomic.gff"
MAP = ROOT / "B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv"
G5MAP = ROOT / "G5_全新基因识别与编号/G5_04_排除记录.tsv"
G5STATUS = ROOT / "G5_全新基因识别与编号/G5_02_全新基因身份核验.tsv"
ALPHA = ROOT / "00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx"
E4RESP = ROOT / "E4_模型修改候选与论文证据裁决/04_呼吸链与能量/04_呼吸与能量_逐候选裁决表.md"
KEGG = ROOT / "D3_全基因组候选功能发现/D3-1_初版候选筛选/KEGG/KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv"
D3 = ROOT / "D3_全基因组候选功能发现/D3-1_初版候选筛选/差异比较/D3-1_跨数据库主表.tsv"
BRENDA = ROOT / "D3_全基因组候选功能发现/D3-1_初版候选筛选/BRENDA/BRENDA_酶学结果.tsv"

def attr(s):
    return {k: urllib.parse.unquote(v) for x in s.split(";") if "=" in x for k,v in [x.split("=",1)]}

genes = []
for line in DOC.read_text(encoding="utf-8-sig").splitlines():
    m = re.match(r"^(\d{4})｜(RU820_RS\d+)｜(WP_\d+\.\d+)｜([^｜]+)｜(.+)$",line)
    if m:
        genes.append(dict(zip(["serial","locus","wp","coordinate","product"],m.groups())))
assert len(genes)==293 and [int(x["serial"]) for x in genes]==list(range(2052,2345))
assert len({x["locus"] for x in genes})==293

official = defaultdict(dict)
with GFF.open(encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.rstrip("\n").split("\t")
        if len(parts)<9 or parts[2] not in ("gene","CDS"): continue
        a=attr(parts[8]); locus=a.get("locus_tag")
        if locus:
            official[locus][parts[2]]=dict(seq=parts[0],start=parts[3],end=parts[4],strand=parts[6],**a)
for g in genes:
    x=official[g["locus"]]
    assert x["gene"]["gene_biotype"]=="protein_coding"
    assert x["CDS"]["protein_id"]==g["wp"]
    assert x["CDS"]["product"]==g["product"]
    assert g["coordinate"]==f'{x["gene"]["seq"]}:{x["gene"]["start"]}-{x["gene"]["end"]}({x["gene"]["strand"]})'

mapped=defaultdict(list)
with MAP.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        if row["当前GCF049位点"]:
            mapped[row["当前GCF049位点"]].append(row)
g5mapped={}
with G5MAP.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        g5mapped[row["当前locus"]]=row
g5status={}
with G5STATUS.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        g5status[row["当前locus"]]=row

kegg=defaultdict(list)
with KEGG.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        kegg[row["KEGG gene/locus"]].append(row)

wb=load_workbook(ALPHA,read_only=True,data_only=True)
rows=list(wb.worksheets[0].values)
header=[str(x).strip() if x is not None else "" for x in rows[1]]
original=[list(r[:30]) for r in rows[2:] if len(r)>16 and r[16]=="2016 original"]
extended=[list(r[:30]) for r in rows[2:] if len(r)>16 and r[16]!="2016 original" and any(x is not None for x in r[:30])]
by_afe=defaultdict(list)
for r in original:
    for afe in set(re.findall(r"AFE_\d{4}",str(r[7])+" "+str(r[8]))):
        by_afe[afe].append(r)
by_current=defaultdict(list)
for r in extended:
    if r[22]: by_current[str(r[22])].append(r)

def uniprot_batch():
    cache=OUT/"UniProt_WP_查询原始.json"
    if cache.exists(): return json.loads(cache.read_text(encoding="utf-8"))
    result={}
    for start in range(0,len(genes),15):
        chunk=genes[start:start+15]
        query="organism_id:243159 AND ("+" OR ".join("xref:refseq-"+g["wp"] for g in chunk)+")"
        url="https://rest.uniprot.org/uniprotkb/search?"+urllib.parse.urlencode({"query":query,"format":"json","size":"500"})
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"G58 evidence audit/1.0"})
            with urllib.request.urlopen(req,timeout=35) as f: data=json.load(f)
            for item in data.get("results",[]):
                for ref in item.get("uniProtKBCrossReferences",[]):
                    if ref.get("database")=="RefSeq":
                        if ref.get("id","").startswith("WP_"):
                            result.setdefault(ref["id"],[]).append(item)
            print("UniProt",start+len(chunk),len(data.get("results",[])),flush=True)
        except Exception as e:
            print("UniProt ERROR",start,str(e),flush=True)
    cache.write_text(json.dumps(result,ensure_ascii=False),encoding="utf-8")
    return result

up=uniprot_batch()
direct_file=OUT/"InterPro_Pfam_直接查询.json"
direct=json.loads(direct_file.read_text(encoding="utf-8")) if direct_file.exists() else {}
tcfile=OUT/"TCDB_Pfam_线索.json"
tcrows=json.loads(tcfile.read_text(encoding="utf-8")) if tcfile.exists() else []
tcbylocus={r[0]:r for r in tcrows}
biofile=OUT/"BioCyc_G58_证据链.json"
biorows=json.loads(biofile.read_text(encoding="utf-8")) if biofile.exists() else []
biobylocus=defaultdict(set)
for r in biorows: biobylocus[r[1]].add(r[8])
brenda_by_ec=defaultdict(list)
with BRENDA.open(encoding="utf-8-sig",newline="") as f:
    for row in csv.DictReader(f,delimiter="\t"):
        brenda_by_ec[row["EC号"]].append(row)
d3rows=[]
with D3.open(encoding="utf-8-sig",newline="") as f:
    reader=csv.DictReader(f,delimiter="\t")
    d3header=reader.fieldnames
    d3rows=list(reader)
g58_loci={g["locus"] for g in genes}
d3g58=[r for r in d3rows if r["当前 NCBI locus tag (GCF_049532655.1)"] in g58_loci]
index=[]
review=[]
rhea_templates=[]
brenda_links=[]
sources=[
 ["G58-S01","NCBI RefSeq GCF_049532655.1 genomic.gff","2026","","official genome annotation","ATCC 23270","genomic feature/CDS","all 293","all","CDS coordinates, protein_id and product; no biochemical measurement","","Automated product is not gene-specific reaction proof","https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/"],
 ["G58-S02","B1 verified legacy-to-current mapping","2026","","sequence/identity mapping","ATCC 23270","exact WP or documented mapping","mapped loci","2016 reaction linkage","Old AFE/current RU820 mapping; no activity measurement","","Ambiguous mappings retained unresolved",str(MAP)],
 ["G58-S03","2016 iMC507 original supplement, mmc1.xls; Alpha read-only extended copy","2016","","original model supplement","ATCC 23270","model reconstruction","model-linked AFE","2016 reaction IDs","Original reaction IDs, formulas, confidence, EC, PMID, GPR and condition bounds","original score only","Original score cannot validate a new GPR or chemical reaction",str(ALPHA)],
 ["G58-S04","UniProtKB REST accession records, taxon 243159, queried 2026-09-24","2026","","protein annotation database","ATCC 23270 reference proteome","mostly automated annotation","matched WP","protein domains/EC","Entry type, cross-references, InterPro/Pfam IDs when present","","Cross-referenced InterPro/Pfam IDs share this database evidence chain; not direct experiment","https://rest.uniprot.org/uniprotkb/search"],
 ["G58-S05","KEGG afr legacy gene table, local D3 snapshot","2026","","pathway database","ATCC 23270 legacy assembly","KO/EC/reaction annotation","verified AFE only","candidate reactions","Legacy gene to KO/EC/reaction reference","","KO/EC transfers cannot prove exact current protein reaction",str(KEGG)],
 ["G58-S06","Campodonico et al. Acidithiobacillus ferrooxidans comprehensive model-driven analysis. Metab Eng Commun 3:84-96","2016","10.1016/j.meteno.2016.03.003","primary model paper","ATCC 23270","model reconstruction","iMC507 genes","model baseline","Published model scope and simulation; not gene-specific enzyme assay","","Model prediction is not biochemical proof","https://doi.org/10.1016/j.meteno.2016.03.003"],
 ["G58-S07","E4 prior adjudication B-017; original cited paper evidence not found","2026","","prior project adjudication","ATCC 23270","literature/database review","AFE_2450 / RU820_RS11310","Sox candidate","Only computational Sox-linked annotation; no gene-specific catalysis","","Historical adjudication, not new primary evidence",str(ROOT/"E4_模型修改候选与论文证据裁决/交接文档_20260923/03_铁硫_逐候选裁决表.md")],
 ["G58-S08","Haskamp et al. The radical SAM protein HemW is a heme chaperone. J Biol Chem 293:2558-2572","2018","10.1074/jbc.RA117.000229","primary biochemical paper","E. coli and other heterologous proteins","purified protein, heme transfer, complementation","HemW/HemN, not ATCC 23270 WP_012537138.1","CPPPGO2 QC","E. coli HemW lacked HemN coproporphyrinogen III dehydrogenase activity and transferred heme","cannot score this strain","Different organism; sequence family and fourth cysteine need gene-specific validation","https://doi.org/10.1074/jbc.RA117.000229"],
 ["G58-S09","Vibrio cholerae FeoB hydrolyzes ATP and GTP in vitro. Metallomics 12:2065-2074","2020","10.1039/d0mt00195c","primary biochemical paper","Vibrio cholerae","purified FeoB NTPase assay","FeoB, not ATCC 23270 WP_012537240.1","FE2tpp QC","Vibrio FeoB hydrolyzed ATP and GTP in vitro","cannot score this strain","NTP hydrolysis does not establish transport coupling or FeoABC subunit necessity in ATCC 23270","https://doi.org/10.1039/d0mt00195c"],
 ["G58-S10","Valdes et al. Extending the models for iron and sulfur oxidation in the extreme acidophile Acidithiobacillus ferrooxidans. BMC Genomics 10:394","2009","10.1186/1471-2164-10-394","primary transcript study and pathway proposal","ATCC 23270","expression comparison, pathway inference","AFE_2550 and hdr cluster","SUCD QC","AFE_2550 annotated in an hdr-associated locus; transcript differences measured","expression only","Does not establish succinate/fumarate catalysis or exact Hdr chemistry","https://doi.org/10.1186/1471-2164-10-394"],
 ["G58-S11","Osorio et al. Anaerobic sulfur metabolism coupled to dissimilatory iron reduction. Appl Environ Microbiol","2013","10.1128/AEM.03057-12","primary physiology and omics paper","ATCC 23270","anaerobic culture, transcript/protein comparison","hdr/dsrE/tusA region","sulfur-pathway QC","Growth and expression/protein changes under sulfur/iron conditions","pathway-level indirect","No purified G58 enzyme reaction, membrane side or complex stoichiometry","https://doi.org/10.1128/AEM.03057-12"],
 ["G58-S12","G5 excluded-list verified original AFE mappings","2026","","sequence/coordinate/neighbor mapping","ATCC 23270","old/current CDS comparison","17 G58 loci","2016 identity","Current loci mapped to original AFE CDS by sequence and coordinate evidence","","Mapping evidence is not enzyme activity",str(G5MAP)],
 ["G58-S13","G5 novel/ambiguous identity status table","2026","","identity adjudication","ATCC 23270","old/current CDS and neighborhood comparison","29 G58 loci without accepted AFE","legacy identity","Distinguishes AFE-free records from seven ambiguous old AFE candidates","","Ambiguous AFE candidates cannot link 2016 GPR",str(G5STATUS)],
 ["G58-S14","D3 cross-database model comparison snapshot","2026","","derived cross-reference table","ATCC 23270 older and current assemblies","database ID comparison","12 G58 loci","13 old model rows","BioCyc/MetaCyc/Rhea/BRENDA summary as originally collected by D3","","Derivative records; each database chain must be verified before reaction change",str(D3)],
 ["G58-S15","EMBL-EBI InterPro API direct query, 2026-09-27","2026","","protein family database","UniProtKB mapped ATCC 23270 entries","InterPro and Pfam entry lookup","WP with matched UniProt accession","protein family IDs","Direct API entry IDs and names; no enzyme assay","","Family signatures do not establish precise substrate or model compartment","https://www.ebi.ac.uk/interpro/api/"],
 ["G58-S16","E4 respiratory and Hdr adjudication B-025 to B-029","2026","","prior project adjudication","ATCC 23270","paper/database evidence review","AFE_2411, AFE_2423, AFE_2553, AFE_2554, AFE_2586","NADHI and Hdr candidates","B-027 to B-029 have conditional expression/protein abundance support; no exact Hdr reaction","expression-level indirect","Prior decision reused as historical evidence, not new enzyme measurement",str(E4RESP)],
 ["G58-S17","Fe/S redox-coupled mercury transformation mediated by A. ferrooxidans ATCC 23270. Microorganisms 11:1028","2023","10.3390/microorganisms11041028","primary physiology/transcriptomics paper","ATCC 23270","Hg stress culture, solution/solid speciation, transcriptomics","mer-associated genes, exact current WP assignment not reported","G58-MERA","Hg speciation and expression changes under 1 mg/L Hg(II) and Fe/S regimes","indirect/3 at pathway scope only","Does not measure purified AFE_2481 MerA reaction, compartment or flux","https://doi.org/10.3390/microorganisms11041028"],
 ["G58-S18","Rhea RHEA:23856 mercury(II) reductase standard reaction","2026","","curated reaction database","general biochemical definition","balanced reaction ontology","EC 1.16.1.1","G58-MERA","Hg + NADP+ + H+ = Hg2+ + NADPH","2 for candidate chemistry with annotation","Does not prove ATCC 23270 WP activity or compartment","https://www.rhea-db.org/rhea/23856"],
 ["G58-S19","Quatrini et al. Insights into iron and sulfur energetic metabolism by microarray profiling. Hydrometallurgy 83:263-272","2006","10.1016/j.hydromet.2006.03.030","primary transcript study","ATCC 23270","Fe2+/S0 culture comparison; microarray with follow-up methods","Nuo/Hdr pathway loci","B-025;B-027;B-028;B-029","Condition-dependent transcript abundance","indirect only","No purified enzyme reaction, complex composition or proton count","https://doi.org/10.1016/j.hydromet.2006.03.030"],
 ["G58-S20","Bellenberg et al. Proteomics of A. ferrooxidans biofilm on pyrite. Front Microbiol 10:592","2019","10.3389/fmicb.2019.00592","primary proteomics study","DSM 14882 type strain = ATCC 23270","Fe liquid versus pyrite biofilm proteomics","AFE_2553, AFE_2554, AFE_2586 and main Nuo cluster","B-027;B-028;B-029;NADHI","Protein abundance differences; main Nuo cluster detected","indirect only","Proteomics does not establish exact catalytic reaction or Hdr subunit AND/OR","https://doi.org/10.3389/fmicb.2019.00592"],
 ["G58-S21","Li et al. Multicenter aerobic iron respiratory chain functions as an ensemble. J Biol Chem 290:18293-18303","2015","10.1074/jbc.M115.657551","primary physiology/spectroscopy paper","ATCC 23270","whole-cell optical redox kinetics","iron respiratory chain, not individual Nuo WP","NADHI","Whole-chain redox kinetics after Fe2+ addition","pathway-level indirect","Macroscopic rate is not Nuo kcat or direct H+/2e- stoichiometry","https://doi.org/10.1074/jbc.M115.657551"],
 ["G58-S22","TCDB public Pfam-to-family mapping, fetched 2026-09-27","2026","","transporter family mapping","all organisms","Pfam family cross-reference","24 G58 transporter-like products","TCDB family candidates","Family-level TC classes inferred from matching Pfam IDs","","No WP-specific TCDB assignment or transported substrate","https://www.tcdb.org/public/pfam.tsv"],
 ["G58-S23","BioCyc GCF_000021485 v30.0 gene-protein-enzyme-reaction snapshot from C2-1","2026","","Tier 3 uncurated PGDB","ATCC 23270 legacy GCF_000021485.1","computed gene/reaction relations","78 G58 current loci mapped via verified AFE","98 distinct reaction IDs","Gene→protein→enzymatic-reaction→reaction chain in official PGDB snapshot","","Old assembly and computational evidence; no independent current-genome enzyme measurement","https://biocyc.org/GCF_000021485/organism-summary"],
 ["G58-S24","D3 exact Pathway Tools/MetaCyc ID mapping from BioCyc PGDB","2026","","reaction ID mapping snapshot","same legacy PGDB","reaction participants, direction, EC, pathway evidence code","98 BioCyc reaction IDs","BioCyc chain rows","Equation participants and direction from PGDB IDs","","MetaCyc direct retrieval absent; one inherited source chain",str(ROOT/"D3_全基因组候选功能发现/D3-1_初版候选筛选/MetaCyc/MetaCyc_BioCyc_ID精确映射.tsv")],
 ["G58-S25","D3 BRENDA organism-tagged enzyme records snapshot","2026","","enzyme database snapshot","ATCC 23270 organism label","EC-based cross-reference","2 G58 loci matched by EC","MerA and AFE_2601 candidate","EC record exists with strain label; no gene-specific assay retrieved","","EC match is not a WP-specific experiment or exact reaction proof",str(BRENDA)],
]

for g in genes:
    locus=g["locus"]; wp=g["wp"]
    matches=mapped[locus]
    exact=[m for m in matches if m["映射状态"]=="完全匹配" and m["当前蛋白ID"]==wp]
    afe=exact[0]["旧AFE位点"] if len(exact)==1 else ""
    mapping_note=(exact[0]["映射方法"] if afe else ("; ".join(m["旧AFE位点"]+":"+m["映射状态"] for m in matches) or "B1 无行"))
    if not afe and locus in g5mapped and g5mapped[locus]["当前protein ID"]==wp:
        afe=g5mapped[locus]["旧AFE位点"]
        mapping_note="G5原始AFE序列/坐标核验："+g5mapped[locus]["重新识别依据"]
    if not afe and locus in g5status:
        mapping_note="G5结论："+g5status[locus]["本轮核验结论"]+"；候选="+g5status[locus]["旧AFE候选"]
    model=by_afe.get(afe,[]) if afe else []
    unik=up.get(wp,[])
    # Require exact RefSeq WP cross-reference; multiple UniProt records remain separate.
    uniacc=";".join(dict.fromkeys(x["primaryAccession"] for x in unik))
    reviewed=";".join(dict.fromkeys(x.get("entryType","") for x in unik))
    ipr=[]; pfam=[]; uec=[]
    for x in unik:
        for z in x.get("uniProtKBCrossReferences",[]):
            if z.get("database")=="InterPro": ipr.append(z.get("id",""))
            if z.get("database")=="Pfam": pfam.append(z.get("id",""))
        for z in x.get("proteinDescription",{}).get("recommendedName",{}).get("ecNumbers",[]):
            uec.append(z.get("value",""))
        for q in x.get("comments",[]):
            if q.get("commentType")=="CATALYTIC ACTIVITY":
                rx=q.get("reaction",{})
                rids=[z.get("id","") for z in rx.get("reactionCrossReferences",[]) if z.get("database")=="Rhea" and z.get("id","").startswith("RHEA:")]
                ecos=";".join(sorted({z.get("evidenceCode","") for z in rx.get("evidences",[])}))
                rhea_templates.append([locus,wp,x["primaryAccession"],";".join(rids),rx.get("name",""),rx.get("ecNumber",""),ecos,"UniProt automated catalytic activity; not direct ATCC 23270 assay","https://rest.uniprot.org/uniprotkb/"+x["primaryAccession"]])
    kres=kegg.get(afe,[]) if afe else []
    kid=";".join(dict.fromkeys(z["KO"] for z in kres if z["KO"]))
    kec=";".join(dict.fromkeys(z["EC"] for z in kres if z["EC"]))
    krxn=";".join(dict.fromkeys(z["KEGG reaction"] for z in kres if z["KEGG reaction"]))
    becs=sorted((set(uec)|set(kec.split(";")))&set(brenda_by_ec))
    for ec in becs:
        for b in brenda_by_ec[ec]:
            brenda_links.append([locus,wp,afe,ec,b["酶名称"],b["BRENDA organism标签"],b["同株判断"],b["底物/反应参与物"],b["产物/反应参与物"],b["辅因子"],b["动力学"],b["证据说明"],"EC/organism match only; no gene-specific linkage",str(BRENDA)])
    modelids=";".join(dict.fromkeys(str(r[0]) for r in model))
    prior=by_current.get(locus,[])
    priorids=";".join(str(r[17] or r[0] or "") for r in prior)
    special=locus=="RU820_RS11310"
    action="证据不足暂缓"
    boundary="当前 GFF 蛋白身份已核对；外部注释和2016模型关联均不能单独证明本株该蛋白对精确反应的催化。"
    nonmet=re.search(r"transcriptional regulator|transcription factor|response regulator|ribosomal protein|translation initiation factor|DNA helicase|DNA polymerase|recombinase|integrase|transposase|restriction|antitoxin|family toxin|molecular chaperone|heat-inducible transcriptional repressor",g["product"],re.I)
    if nonmet and not model and not krxn and not priorids and not biobylocus.get(locus) and not any(z[0]==locus for z in rhea_templates):
        action="无需修改（当前无代谢反应关联）"
        boundary="当前注释指向调控、遗传信息或蛋白稳态功能；未找到可归给该基因的精确代谢反应。"
    if special:
        boundary="E4 B-017已核查：Sox相关仅计算注释；未找到本株该蛋白的纯化催化或扰动证据，不能指定SoxY结合硫反应。"
    if locus=="RU820_RS10820":
        boundary="当前蛋白的UniProt为自动HemW注释且有InterPro IPR004559；异株HemW实验未见HemN氧化活性，本株CPPPGO2 OR分支待QC核验。"
    if locus in ("RU820_RS11605","RU820_RS11610","RU820_RS11615","RU820_RS11620"):
        boundary="2016 FE2tpp写ATP水解与四基因AND；FeoB的核苷酸选择、FeoC必需性和本株实际转运耦联未测定。"
    if locus=="RU820_RS11745":
        boundary="当前蛋白是HdrB家族而2016 SUCD把它单独指为琥珀酸反应GPR；同株文献只给hdr簇表达/通路推断。"
    if locus=="RU820_RS11440":
        boundary="UniProt/KEGG给MerA EC 1.16.1.1，Rhea定义Hg(II)+NADPH反应；同株汞暴露论文测得群体汞形态/表达，未逐基因测定本蛋白催化与区室。"
    if priorids in ("B-027","B-028","B-029"):
        action="仅更新证据或编号（沿用E4历史裁决）"
        boundary="E4记录同株条件性转录或蛋白丰度；具体Hdr反应、复合体成员和电子受体未直接测定。"
    unresolved=("具体反应化学式、辅因子、方向、区室与GPR待核验" if model or krxn else "未找到可核实的精确代谢反应")
    if special: unresolved="SoxY底物、催化亚基、电子受体与本株直接实验证据"
    if locus=="RU820_RS10820": unresolved="CPPPGO2中AFE_2346/HemW的催化资格；需纯化或遗传验证"
    if locus in ("RU820_RS11605","RU820_RS11610","RU820_RS11615","RU820_RS11620"): unresolved="FE2tpp中ATP计量、FeoC必需性及各成员运输贡献"
    if locus=="RU820_RS11745": unresolved="SUCD与HdrB的错误关联、正确琥珀酸酶及Hdr精确底物"
    if locus=="RU820_RS11440": unresolved="MerA同株纯化/遗传证据、胞内区室和汞交换边界；暂不新增模型反应"
    extra_source=(";G58-S12" if locus in g5mapped and afe==g5mapped[locus]["旧AFE位点"] else "")
    if not afe and locus in g5status: extra_source+=";G58-S13"
    if priorids in ("B-025","B-026","B-027","B-028","B-029"): extra_source+=";G58-S16"
    if priorids in ("B-025","B-027","B-028","B-029"): extra_source+=";G58-S19"
    if priorids in ("B-027","B-028","B-029"): extra_source+=";G58-S20"
    if modelids and "NADHI" in modelids: extra_source+=";G58-S20;G58-S21"
    if locus=="RU820_RS11440": extra_source+=";G58-S17;G58-S18"
    if wp in direct: extra_source+=";G58-S15"
    if locus in tcbylocus: extra_source+=";G58-S22"
    if locus in biobylocus: extra_source+=";G58-S23;G58-S24"
    if becs: extra_source+=";G58-S25"
    if any(d["当前 NCBI locus tag (GCF_049532655.1)"]==locus for d in d3g58): extra_source+=";G58-S14"
    extra_source+=(";G58-S08" if locus=="RU820_RS10820" else
                  ";G58-S09" if locus in ("RU820_RS11605","RU820_RS11610","RU820_RS11615","RU820_RS11620") else
                  ";G58-S10;G58-S11" if locus=="RU820_RS11745" else "")
    dr=direct.get(wp,{})
    dr_ipr=";".join(sorted({z["id"] for z in dr.get("interpro",{}).get("entries",[]) if z.get("id")}))
    dr_pfam=";".join(sorted({z["id"] for z in dr.get("pfam",{}).get("entries",[]) if z.get("id")}))
    direct_status=("两库成功" if dr and "error" not in dr.get("interpro",{}) and "error" not in dr.get("pfam",{}) else "查询失败或未匹配UniProt")
    index.append([
        g["serial"],locus,wp,g["coordinate"],g["product"],afe,mapping_note,
        uniacc,reviewed,";".join(sorted(set(ipr))),";".join(sorted(set(pfam))),";".join(sorted(set(uec))),
        kid,kec,krxn,modelids,priorids,
        "否",action,("G58-G"+g["serial"]),unresolved,
        "G58-S01;G58-S02;G58-S03;G58-S04;G58-S05"+(";G58-S07" if special else "")+extra_source,
        boundary,
        "https://www.ncbi.nlm.nih.gov/protein/"+wp,
        ("https://rest.uniprot.org/uniprotkb/"+unik[0]["primaryAccession"] if unik else ""),
        ("https://www.kegg.jp/entry/afr:"+afe if afe else ""),dr_ipr,dr_pfam,direct_status,
        ";".join(sorted({r for z in rhea_templates if z[0]==locus for r in z[3].split(";") if r})),
        ";".join(sorted({x.split(" ")[0] for x in tcbylocus.get(locus,[None,None,None,None,None,""])[5].split(";") if x})),
        ";".join(sorted(biobylocus.get(locus,set()))),";".join(becs)
    ])
    rowid="G58-G"+g["serial"]
    # One review record for every gene. Original model rows are separately copied below.
    primarydoi=("10.1128/AEM.03057-12" if priorids in ("B-027","B-029") else
                "10.3389/fmicb.2019.00592" if priorids=="B-028" else None)
    rec=[None]*30
    rec[16:30]=["G58 gene audit",rowid,"G58",modelids,None,afe,locus,primarydoi,
                 "official CDS; exact WP mapping; protein/database annotation; 2016 original rows",
                 boundary,action,None,
                 "G58-S01;G58-S02;G58-S03;G58-S04;G58-S05"+(";G58-S07" if special else "")+extra_source,
                 "不评分；未定义新的精确反应"]
    review.append(rec)
reaction_ids=defaultdict(set)
reaction_afe=defaultdict(set)
reaction_source={}
for g in index:
    locus,afe=g[1],g[5]
    for r in by_afe.get(afe,[]) if afe else []:
        rid=str(r[0]); reaction_ids[rid].add(locus); reaction_afe[rid].add(afe); reaction_source[rid]=r
for rid in sorted(reaction_source):
    r=reaction_source[rid]
    line=list(r[:16])+["G58 2016 reaction comparison","G58-R-"+rid,"G58",rid,
                        str(r[2] or ""),";".join(sorted(reaction_afe[rid])),
                        ";".join(sorted(reaction_ids[rid])),None,
                        "2016 original supplement and current exact WP mapping",
                        "2016原行仅是模型基线；本轮未完成独立的反应化学与实验重评。",
                        "保留2016原行；模型动作暂缓",None,
                        "G58-S01;G58-S02;G58-S03",
                        "2016原反应分数"+str(r[3])+"；2026未评分"]
    review.append(line)
mer=[g for g in index if g[1]=="RU820_RS11440"][0]
mer_review=[None]*30
mer_review[16:30]=["G58 candidate reaction","G58-MERA","G58",None,
                    "Hg(2+) + NADPH -> Hg(0) + NADP(+) + H(+); RHEA:23856 reverse orientation",
                    "AFE_2481","RU820_RS11440","10.3390/microorganisms11041028",
                    "UniProt/KEGG enzyme annotation; Rhea chemistry; same-strain Hg physiology/transcriptomics",
                    "Rhea给出平衡的通式，同株论文未对AFE_2481纯化、敲除或测定反应区室与通量。",
                    "证据不足暂缓；不新增SBML反应",2,"G58-S04;G58-S05;G58-S17;G58-S18",
                    "仅对候选通式的间接注释评分2；不转给区室、方向、GPR或Hg交换"]
review.append(mer_review)
mer[19]+=";G58-MERA"

payload={
 "index_header":["总序号","当前 locus","当前 WP","坐标","当前 NCBI product","旧 AFE","旧号映射依据/状态","UniProt accession","UniProt 审校状态","InterPro IDs via UniProt","Pfam IDs via UniProt","UniProt EC","KEGG KO","KEGG EC","KEGG reaction","2016 关联 Reaction ID","历史候选 ID","精确反应已核实","2026 动作","关联行 ID","未决点","Source ID","证据边界","NCBI 蛋白链接","UniProt 链接","KEGG 基因链接","InterPro direct IDs","Pfam direct IDs","直接查询状态","Rhea IDs via UniProt","TCDB family clue","BioCyc legacy reaction IDs","BRENDA EC clue"],
 "index":index,"review_header":header,"review":review,
 "sources_header":["Source ID","Full citation","Year","DOI","Source type","Strain","Experimental type","Genes / proteins studied","Related Reaction IDs / candidate IDs","What was directly measured","Applicable confidence level","Evidence limitation","Source URL"],"sources":sources,
 "d3_header":d3header,"d3":[[r[h] for h in d3header] for r in d3g58],
 "rhea_header":["Current locus","WP protein","UniProt accession","Rhea ID","Standard equation","EC","Evidence code","Evidence boundary","Original protein URL"],"rhea_templates":rhea_templates,
 "tcdb_header":["Current locus","WP protein","NCBI product","Pfam direct IDs","TCDB candidate family count","TCDB family hints (up to 15)","Evidence boundary","TCDB source URL"],"tcdb":tcrows,
 "biocyc_header":["Serial","Current locus","WP protein","Verified AFE","Legacy AFE_RS","PGDB protein","PGDB enzyme","Enzymatic-reaction frame","BioCyc reaction ID","Pathway ID","Equation participants","Direction","EC","PGDB evidence code","Rhea xref","MetaCyc direct retrieval","Evidence boundary","Original BioCyc URL"],"biocyc":biorows,
 "brenda_header":["Current locus","WP protein","Verified AFE","EC","Enzyme name","BRENDA organism","Strain status","Substrate","Product","Cofactor","Kinetics","Snapshot limitation","Evidence boundary","Snapshot path"],"brenda":brenda_links,
 "counts":dict(genes=len(index),unique_loci=len({x[1] for x in index}),mapped=sum(bool(x[5]) for x in index),uniprot_matched=sum(bool(x[7]) for x in index),direct_interpro=sum(bool(x[26]) for x in index),direct_pfam=sum(bool(x[27]) for x in index),model_linked_genes=sum(bool(x[15]) for x in index),model_reaction_rows=len(review)-len(index),historical_candidates=sum(bool(x[16]) for x in index),d3_rows=len(d3g58),uniprot_rhea_templates=len(rhea_templates),tcdb_transporter_rows=len(tcrows),tcdb_family_hints=sum(bool(x[4]) for x in tcrows),biocyc_genes=len(biobylocus),biocyc_chain_rows=len(biorows),biocyc_reactions=len({r[8] for r in biorows}),brenda_ec_links=len(brenda_links))
}
(OUT/"G58_中间数据.json").write_text(json.dumps(payload,ensure_ascii=False),encoding="utf-8")
print(payload["counts"])
