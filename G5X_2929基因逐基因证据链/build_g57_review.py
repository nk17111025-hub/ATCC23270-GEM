"""Build a traceable G57 working review from frozen local source exports.

This script deliberately leaves unsupported chemistry and scores blank.
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

import openpyxl


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "G5X_2929基因逐基因证据链" / "G57_审查输出"
OUT.mkdir(exist_ok=True)
GFF = ROOT / "B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff"
D3 = ROOT / "D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv"
D3E = ROOT / "D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_注释证据表.tsv"
RHEA = ROOT / "D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/Rhea_全基因EC查询汇总.tsv"
ALPHA = ROOT / "00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx"
UNIPROT = OUT / "UniProt_243159_20260924.tsv"
LIT = OUT / "EuropePMC_gene_scan_20260927.json"
CURRENT_LIT = OUT / "G57_当前编号双路文献检索.tsv"
TCDB_REVIEW = OUT / "G57_TCDB家族同源核验.tsv"
LIT_TRIAGE = OUT / "G57_文献命中逐篇判读.tsv"
PRIMER_AUDIT = OUT / "G57_论文引物身份核验.tsv"
EC_MODEL = OUT / "G57_同EC旧模型候选.tsv"


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(name, rows, fields):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def attrs(s):
    return {k: unquote(v) for k, v in (p.split("=", 1) for p in s.split(";") if "=" in p)}


genes, cds = [], defaultdict(list)
with GFF.open(encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9:
            continue
        a = attrs(p[8])
        locus = a.get("locus_tag", "")
        if p[2] == "gene" and a.get("gene_biotype") == "protein_coding":
            genes.append((locus, p[0], p[3], p[4], p[6], a))
        elif p[2] == "CDS":
            cds[locus].append((p, a))
assert len(genes) == 2929, len(genes)
part = genes[1758:2051]
assert len(part) == 293 and part[0][0] == "RU820_RS09280" and part[-1][0] == "RU820_RS10800"

cross = {r["当前locus"]: r for r in read_tsv(D3)}
uniprot = defaultdict(list)
if UNIPROT.exists():
    for r in read_tsv(UNIPROT):
        for wp in r["RefSeq"].split(";"):
            if wp.strip():uniprot[wp.strip()].append(r)
literature = __import__("json").loads(LIT.read_text(encoding="utf-8")) if LIT.exists() else {}
current_literature = {r['RU820']:r for r in read_tsv(CURRENT_LIT)} if CURRENT_LIT.exists() else {}
tcdb_review = {r['RU820']:r for r in read_tsv(TCDB_REVIEW)} if TCDB_REVIEW.exists() else {}
rhea_by_id = {r["Reaction identifier"]:r for r in read_tsv(RHEA)}
evidence = defaultdict(list)
for r in read_tsv(D3E):
    evidence[r["当前locus"]].append(r)

wb = openpyxl.load_workbook(ALPHA, read_only=True, data_only=True)
model_rows = []
historical = defaultdict(list)
for sheetrow in wb.worksheets[0].iter_rows(min_row=3, values_only=True):
    if sheetrow[16] == "2016 original":
        model_rows.append(sheetrow)
    elif sheetrow[16] == "2026 candidate" and sheetrow[22]:
        historical[str(sheetrow[22])].append(sheetrow)
model_by_id = {str(r[0]):r for r in model_rows}
model_by_afe = defaultdict(list)
for row in model_rows:
    for afe in set(re.findall(r"AFE_\d{4}", " ".join(str(row[i] or "") for i in (7, 8, 9)))):
        model_by_afe[afe].append(row)

source_specs = [
    ("G57-S01", "NCBI RefSeq GCF_049532655.1 genomic.gff, frozen 2026-09-22", "2026", "官方同株基因组", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/"),
    ("G57-S02", "B1 old/new locus mapping, GCF_000021485.1 to GCF_049532655.1", "2026", "旧新身份映射", ""),
    ("G57-S03", "D3-1-R cross-database accession table, 2026-09-23", "2026", "数据库交叉引用", ""),
    ("G57-S04", "Campodonico et al. 2016 iMC507 supplementary table mmc1.xls; Alpha preserved copy", "2016", "2016 原始模型", ""),
    ("G57-S05", "KEGG afr organism annotations in D3-1-R; old AFE gene namespace", "2026", "自动注释", "https://www.genome.jp/kegg-bin/show_organism?org=afr"),
    ("G57-S06", "BioCyc AFER243159 PGDB cross-references in D3-1-R", "2026", "路径工具预测", "https://biocyc.org/AFER243159/"),
    ("G57-S07", "Rhea EC cross-references in D3-1-R", "2026", "EC 到反应交叉引用", "https://www.rhea-db.org/"),
    ("G57-S08", "BRENDA EC records in D3-1-R", "2026", "EC 条目", "https://www.brenda-enzymes.org/"),
    ("G57-S09", "Vera M, Pagliai F, Guiliani N, Jerez CA. The chemolithoautotroph Acidithiobacillus ferrooxidans can survive under phosphate-limiting conditions by expressing a C-P lyase operon that allows it to grow on phosphonates. Applied and Environmental Microbiology 74:1829-1835.", "2008", "同株原始研究", "https://doi.org/10.1128/AEM.02101-07"),
    ("G57-S10", "Tsai YCC, Lapina MC, Bhushan S, Mueller-Cajar O. Identification and characterization of multiple rubisco activases in chemoautotrophic bacteria. Nature Communications 6:8883.", "2015", "同株基因重组蛋白实验", "https://doi.org/10.1038/ncomms9883"),
    ("G57-S11", "UniProtKB taxon 243159 stream, RefSeq cross-reference, InterPro and Pfam fields; retrieved 2026-09-24", "2026", "蛋白数据库交叉引用", "https://rest.uniprot.org/uniprotkb/stream?query=organism_id%3A243159&format=tsv&fields=accession%2Creviewed%2Cprotein_name%2Cgene_primary%2Cxref_refseq%2Cxref_interpro%2Cxref_pfam%2Cec"),
    ("G57-S12", "Alpha ATCC23270-2026 final evidence table, 45 historical candidate decisions, frozen local copy", "2026", "历史候选裁决", "https://drive.google.com/file/d/11KtlD3IIvaI4iwBX4bFVJQzscnTjviHv/view"),
]
sources = []
for sid, title, year, typ, url in source_specs:
    is_rub = sid == "G57-S10"
    sources.append({"Source ID":sid,"Full citation":title,"Year":year,"DOI":"10.1038/ncomms9883" if is_rub else "10.1128/AEM.02101-07" if sid == "G57-S09" else "","Source type":typ,"Strain":"ATCC 23270 原始 DNA；E. coli 表达重组蛋白" if is_rub else "ATCC 23270" if sid in {"G57-S01","G57-S02","G57-S04","G57-S09","G57-S11","G57-S12"} else "需逐条核实","Experimental type":"重组蛋白纯化；Rubisco 活性与 CbbQO ATPase/活化实验" if is_rub else "生长实验、polyP、共转录、宏阵列、RT-PCR" if sid == "G57-S09" else "无；数据库/模型记录","Genes / proteins studied":"AFE_2155 CbbM；AFE_2156 CbbQ2；AFE_2157 CbbO2" if is_rub else "phnG—phnM；其他成员见原文" if sid == "G57-S09" else "G57 相关记录见逐基因索引","Related Reaction IDs / candidate IDs":"G57-1891—G57-1893；RUBISCO" if is_rub else "G57-1999—G57-2005；PPTtpp" if sid == "G57-S09" else "见反应审查表","What was directly measured":"同株基因来源的重组 CbbM 羧化酶活性；CbbQ2/CbbO2 对受抑 Rubisco 的 ATP 依赖活化" if is_rub else "ATCC 23270 在甲基/乙基膦酸盐作为唯一磷源时生长；phnG—phnM 共转录；甲基膦酸盐诱导 phnG、phnJ 转录；polyP 变化" if sid == "G57-S09" else "无基因专属实验测量","Applicable confidence level":"CbbM 正向羧化反应可评 4；活化过程另评，不向逆向或 GPR 外推" if is_rub else "通路生理支持；不能为具体酶步骤直接赋 4" if sid == "G57-S09" else "不能单独达到 3/4","Evidence limitation":"重组蛋白实验不证明胞内区室、逆向通量或 CbbQO 是羧化催化亚基" if is_rub else "未测各 Phn 蛋白的催化化学式、运输、耦联或单独敲除；不能把全链生长当成单酶测定" if sid == "G57-S09" else "数据库注释及旧模型关联不能证明具体基因催化反应；与其他数据库共源时不重复计数","Original URL":url,"Local source":"G5X_2929基因逐基因证据链/G57_审查输出/UniProt_243159_20260924.tsv" if sid == "G57-S11" else str({"G57-S01":GFF,"G57-S02":ROOT / "B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv","G57-S03":D3,"G57-S04":ALPHA,"G57-S05":D3,"G57-S06":D3,"G57-S07":D3,"G57-S08":D3,"G57-S12":ALPHA}[sid].relative_to(ROOT)) if sid not in {"G57-S09","G57-S10"} else "期刊原文"})
sources.append({"Source ID":"G57-S13","Full citation":"Europe PMC REST search of 293 verified AFE/WP identifiers, 2026-09-27","Year":2026,"DOI":"","Source type":"逐基因文献检索线索","Strain":"逐篇待核","Experimental type":"检索，不是实验","Genes / proteins studied":"G57 全部 293 基因","Related Reaction IDs / candidate IDs":"逐基因文献检索线索.tsv","What was directly measured":"无","Applicable confidence level":"不能据检索命中赋分","Evidence limitation":"命中可能只在背景提及，零命中不能证明无论文；每篇需查菌株、实验对象和数据","Original URL":"https://europepmc.org/","Local source":"G5X_2929基因逐基因证据链/G57_审查输出/EuropePMC_gene_scan_20260927.json"})
sources.append({"Source ID":"G57-S14","Full citation":"Osorio H et al. Anaerobic sulfur metabolism coupled to dissimilatory iron reduction in the extremophile Acidithiobacillus ferrooxidans. Applied and Environmental Microbiology 79:2172–2181.","Year":2013,"DOI":"10.1128/AEM.03057-12","Source type":"同株原始研究","Strain":"ATCC 23270T","Experimental type":"S0/Fe3+ 厌氧与 S0/O2 有氧培养；转录组、蛋白组、Cu2+ 沉淀观察","Genes / proteins studied":"AFE_2177—2179；AFE_2181 SreABCD","Related Reaction IDs / candidate IDs":"G57-1910—G57-1913；SQRED1","What was directly measured":"Sre 亚基表达差异及厌氧条件下与 H2S 形成相符的 Cu2+ 沉淀","Applicable confidence level":"表达属间接证据；不能为 Sre 精确反应赋 4","Evidence limitation":"SreABCD 单酶催化、醌种类、膜侧、精确质子数未测；H2S 证据属全细胞生理推断","Original URL":"https://doi.org/10.1128/AEM.03057-12","Local source":"期刊原文"})
sources.append({"Source ID":"G57-S15","Full citation":"Rhea EC query summary, captured 2026-09-23; Rhea reaction identifiers, equations and ChEBI cross-references","Year":2026,"DOI":"","Source type":"反应数据库 EC 交叉引用","Strain":"不适用；反应通用","Experimental type":"数据库记录","Genes / proteins studied":"仅由 EC/数据库记录连接候选基因","Related Reaction IDs / candidate IDs":"逐基因索引 Rhea 字段","What was directly measured":"无本株基因专属测量","Applicable confidence level":"不能单独升至 3/4","Evidence limitation":"Rhea 精确方程不等于当前 WP 催化该方程；方向、区室和体内电子供体需另证","Original URL":"https://www.rhea-db.org/","Local source":str(RHEA.relative_to(ROOT))})
sources.append({"Source ID":"G57-S16","Full citation":"Valdés J et al. Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications. BMC Genomics 9:597.","Year":2008,"DOI":"10.1186/1471-2164-9-597","Source type":"同株基因组功能预测","Strain":"ATCC 23270","Experimental type":"基因组注释与同源性推断","Genes / proteins studied":"AFE_2149—2154 氢化酶簇；AFE_2267—2276 与 AFE_2300—2305 附近 TonB/Exb 系统","Related Reaction IDs / candidate IDs":"HYD4；FE32tpp","What was directly measured":"未测本组单基因催化或转运","Applicable confidence level":"功能与反应为预测，不能据此升至 3/4","Evidence limitation":"formate→H2 与 Fe3+ 外膜转运的确切底物、方向、复合体成员和化学式未经同株实验证实","Original URL":"https://doi.org/10.1186/1471-2164-9-597","Local source":"期刊原文"})
sources.append({"Source ID":"G57-S17","Full citation":"Drobner E, Huber H, Stetter KO. Thiobacillus ferrooxidans, a facultative hydrogen oxidizer. Applied and Environmental Microbiology 56:2922–2923.","Year":1990,"DOI":"10.1128/AEM.56.9.2922-2923.1990","Source type":"同株原始生理研究","Strain":"ATCC 23270","Experimental type":"H2 培养与氢化酶诱导","Genes / proteins studied":"未指认具体 RU820/AFE 基因","Related Reaction IDs / candidate IDs":"HYD4 生理背景","What was directly measured":"菌株使用 H2 生长及氢化酶诱导","Applicable confidence level":"菌株层级，不可赋给 HYD4 精确反应","Evidence limitation":"未证明 AFE_2149—2154 催化甲酸氧化产氢","Original URL":"https://doi.org/10.1128/AEM.56.9.2922-2923.1990","Local source":"期刊原文"})
sources.append({"Source ID":"G57-S18","Full citation":"IUBMB Enzyme Nomenclature EC 1.3.5.2, EC 1.1.1.49 and EC 6.4.1.2, accessed 2026-09-27","Year":2026,"DOI":"","Source type":"酶反应标准定义","Strain":"通用，不限菌株","Experimental type":"标准酶命名及反应式","Genes / proteins studied":"RU820_RS10450、RU820_RS09360 与 RU820_RS09540 的功能类别对照","Related Reaction IDs / candidate IDs":"DHORD；G6PDH2；ACCOAC；BIOC1","What was directly measured":"无本株实验","Applicable confidence level":"指出旧方程与酶类别、反应方向或复合体构成的冲突；不能单独确定本株体内方向","Evidence limitation":"DHORD 醌种类/膜侧、G6PDH2 通量方向、ACCOAC 复合体确切成员及 BIOC1 载体化学式均待本株核验","Original URL":"https://iubmb.qmul.ac.uk/enzyme/EC1/3/5/2.html ; https://iubmb.qmul.ac.uk/enzyme/EC1/1/1/49.html ; https://iubmb.qmul.ac.uk/enzyme/EC6/4/1/2.html","Local source":"IUBMB 官网"})
sources.append({"Source ID":"G57-S19","Full citation":"Ribeiro DA et al. The small heat shock proteins from Acidithiobacillus ferrooxidans: gene expression, phylogenetic analysis, and structural modeling. BMC Microbiology 11:259.","Year":2011,"DOI":"10.1186/1471-2180-11-259","Source type":"文献编号排除记录","Strain":"实验为 A. ferrooxidans LR；ATCC 23270 仅作为基因组参考","Experimental type":"qRT-PCR 引物定位核对","Genes / proteins studied":"文中 Afe_2172 Hsp20；不能归给当前 RU820_RS10045","Related Reaction IDs / candidate IDs":"G57-1905 文献排除","What was directly measured":"论文 Afe_2172 两引物精确落在当前基因组 796313/796390 的 Hsp20 区域","Applicable confidence level":"不支持当前 RU820_RS10045 的热休克蛋白功能或代谢反应","Evidence limitation":"论文编号与 G57 的旧 AFE_2172 字面相同，但实验对象对应另一 WP；不可转移","Original URL":"https://doi.org/10.1186/1471-2180-11-259","Local source":"期刊原文及本地 check_phn_primers.py 引物定位"})
sources.append({"Source ID":"G57-S20","Full citation":"Pang JJ, Shin JS, Li SY. The catalytic role of RuBisCO for in situ CO2 recycling in Escherichia coli. Frontiers in Bioengineering and Biotechnology 8:543807.","Year":2020,"DOI":"10.3389/fbioe.2020.543807","Source type":"同株来源基因的异源表达实验","Strain":"基因来源 ATCC 23270；实验宿主 E. coli","Experimental type":"异源表达 CbbM/CbbQ2/CbbO2；RuBisCO 活性及 CO2 再利用","Genes / proteins studied":"AFE_2155 CbbM；AFE_2156 CbbQ2；AFE_2157 CbbO2","Related Reaction IDs / candidate IDs":"G57-1891—G57-1893；RUBISCO","What was directly measured":"工程 E. coli 中同株来源 CbbM 与活化因子共同表达后的酶活及碳回收变化","Applicable confidence level":"支持 CbbM 羧化及 CbbQ2/CbbO2 活化身份；与 Tsai 2015 共享生物学对象","Evidence limitation":"宿主非 A. ferrooxidans；不能证实体内方向、区室或直接把宿主代谢通量移植入本株","Original URL":"https://doi.org/10.3389/fbioe.2020.543807","Local source":"期刊原文"})
sources.append({"Source ID":"G57-S21","Full citation":"Li XJ et al. A one-step route for the conversion of Cd waste into CdS quantum dots by Acidithiobacillus sp. via unique biosynthesis pathways. RSC Chemical Biology 6:281–294.","Year":2025,"DOI":"10.1039/d4cb00195h","Source type":"同株来源 WP 的重组蛋白与晶体结构实验","Strain":"ATCC 23270 来源 WP_012536942.1；异源表达宿主 E. coli","Experimental type":"蛋白纯化；半胱氨酸+Cd2+ 体系的 UV/荧光实验；X 射线晶体结构 PDB 9IQH","Genes / proteins studied":"WP_012536942.1＝RU820_RS09520，O-succinylhomoserine sulfhydrylase","Related Reaction IDs / candidate IDs":"G57-1804；HACT；半胱氨酸脱硫候选","What was directly measured":"纯化蛋白促进半胱氨酸/Cd2+ 体系 CdS 纳米粒子生成；结构含 PLP","Applicable confidence level":"蛋白身份与体系活性有直接证据；旧 HACT 精确反应不能升分","Evidence limitation":"未完整测定半胱氨酸→丙酮酸+NH3+H2S 的独立化学计量；未测 O-乙酰高丝氨酸+H2S→高半胱氨酸，且实验为 pH 8.0 体外体系","Original URL":"https://doi.org/10.1039/d4cb00195h","Local source":"期刊原文；PDB 9IQH"})
sources.append({"Source ID":"G57-S22","Full citation":"TCDB public protein sequence export cached 2026-09-25; G57 product-anchored family alignment, 2026-09-27","Year":2026,"DOI":"","Source type":"转运蛋白家族同源性","Strain":"跨物种参考序列；当前查询蛋白来自 ATCC 23270","Experimental type":"当前 product 预选大类；5-mer 预筛；BLOSUM62 局部比对","Genes / proteins studied":"G57 44 条具有可映射转运家族的当前 WP","Related Reaction IDs / candidate IDs":"G57_TCDB家族同源核验.tsv","What was directly measured":"无运输实测；仅氨基酸序列相似性","Applicable confidence level":"家族线索，不足以为精确转运反应评分","Evidence limitation":"同源性不证明底物、耦联、方向或复合体；未达阈值也不排除家族功能；参考库选择依当前 product","Original URL":"https://www.tcdb.org/public/tcdb","Local source":"G5X_2929基因逐基因证据链/G57_审查输出/G57_TCDB家族同源核验.tsv"})

index = []
review = []
alpha_headers = ["Reaction ID","Reaction Name","Reaction Formula","Confidence Level","EC Number","PMID","Subsystem","Gene-Reaction Association","Gene-Protein-Reaction Association","Protein-Reaction-Association","Fe2 lb","Fe2 ub","tetrathionate lb","tetrathionate ub","sulfur lb","sulfur ub","Record class","Candidate ID","Candidate group","2016 linked Reaction ID","2026 reaction object","Legacy AFE locus","Current locus","Primary evidence DOI","Evidence type","Evidence boundary","2026 action","2026 reaction score","Source ID","Score scope"]
reaction_rows = []

for n, (locus, seqid, start, end, strand, ga) in enumerate(part, 1759):
    cs = cds[locus]
    assert len(cs) == 1, (locus, len(cs))
    cp, ca = cs[0]
    wp, product = ca.get("protein_id", ""), ca.get("product", "")
    ur = uniprot.get(wp, [])
    assert len(ur) <= 1, (wp, len(ur))
    u = ur[0] if ur else {}
    assert wp.startswith("WP_") and cp[3:5] == [start, end]
    d = cross.get(locus, {})
    assert d.get("当前protein ID") == wp, (locus, wp, d.get("当前protein ID"))
    afe = d.get("旧AFE locus", "")
    map_status = d.get("B1映射状态", "")
    map_method = d.get("B1映射方法", "")
    if d.get("映射冲突") == "是" or map_status not in {"完全匹配", "高可信匹配"}:
        afe_claim = "未确认"
    else:
        afe_claim = afe or "未确认"
    linked = model_by_afe.get(afe_claim, []) if afe_claim != "未确认" else []
    historic = historical.get(locus, [])
    ids = list(dict.fromkeys(str(x[0]) for x in linked))
    erows = evidence[locus]
    dbids = defaultdict(list)
    for e in erows:
        if e["数据库"] not in {"NCBI", "B1"}:
            dbids[e["数据库"]].append(e["记录ID"])
    mentions = {"KEGG KO":d.get("KEGG KO", ""),"KEGG EC":d.get("KEGG EC", ""),"KEGG reaction":d.get("KEGG reaction", ""),"BioCyc reaction":d.get("BioCyc reaction ID", ""),"Rhea crossref":d.get("Rhea精确交叉引用", ""),"BRENDA":d.get("BRENDA EC记录", "")}
    source_ids = ["G57-S01", "G57-S03", "G57-S04"]
    if afe_claim != "未确认":source_ids.append("G57-S02")
    if mentions["KEGG KO"] or mentions["KEGG reaction"]:source_ids.append("G57-S05")
    if mentions["BioCyc reaction"]:source_ids.append("G57-S06")
    if mentions["Rhea crossref"]:source_ids.append("G57-S07")
    if mentions["BRENDA"]:source_ids.append("G57-S08")
    if 1999 <= n <= 2005:source_ids.append("G57-S09")
    if 1891 <= n <= 1893:source_ids.append("G57-S10")
    if u:source_ids.append("G57-S11")
    if historic:source_ids.append("G57-S12")
    if literature:source_ids.append("G57-S13")
    if 1910 <= n <= 1913:source_ids.append("G57-S14")
    if 1885 <= n <= 1890:source_ids.extend(["G57-S16","G57-S17"])
    if n in {1988,1989,1990,1991,1992,1994,1995,1996,1997,2019,2020,2021,2022,2023,2024}:source_ids.append("G57-S16")
    if n in {1774,1808,1985}:source_ids.append("G57-S18")
    if n == 1905:source_ids.append("G57-S19")
    if 1891 <= n <= 1893:source_ids.append("G57-S20")
    if n == 1804:source_ids.append("G57-S21")
    if locus in tcdb_review:source_ids.append("G57-S22")
    if mentions["Rhea crossref"]:source_ids.append("G57-S15")
    product_metabolic_hint = bool(re.search(r"decarboxylase|transferase|hydrolase|oxidoreductase|dehydrogenase|synthase|synthetase|kinase|epimerase|isomerase|reductase|oxygenase|permease|transporter|porin|receptor|antiporter|symporter|ATPase|deiminase",product,re.I) and not re.search(r"DNA|RNA|tRNA|ribosom|histidine kinase|Sec[DF]|protein translocase|secretion|endonuclease|protease|peptidase",product,re.I))
    metabolism_hint = bool(mentions["KEGG reaction"] or mentions["BioCyc reaction"] or mentions["Rhea crossref"] or d.get("旧模型未覆盖功能候选") == "是" or ids or product_metabolic_hint)
    action = "证据不足暂缓" if metabolism_hint else "无需修改（暂定）"
    if 1891 <= n <= 1893:action = "仅更新证据或编号"
    if n == 1773:action = "仅更新证据或编号"
    gaps = ["UniProt/InterPro/Pfam 基因专属记录及保守位点待核", "同株原始实验与精确反应待核"]
    if afe_claim == "未确认":gaps.append("旧 AFE 未确认")
    if metabolism_hint:gaps += ["反应化学计量、方向、区室与复合体待核", "2016 反应等价性待逐项复核"]
    if product_metabolic_hint and not ids and not any((mentions["KEGG reaction"],mentions["BioCyc reaction"],mentions["Rhea crossref"])):gaps.append("当前 product 提示可能有代谢/运输功能，但跨库未获得基因专属精确反应；不能暂判无需修改")
    if "transport" in product.lower() or "permease" in product.lower() or "antiporter" in product.lower() or "receptor" in product.lower():gaps.append("TCDB、底物和耦联机制待核")
    if 1999 <= n <= 2005:gaps.append("Vera 2008 早期 TIGR 称 AFE-811—821；论文 8 条 phnG—phnM 引物全部精确命中当前 RS10520—RS10550。原始论文支持同株共转录和膦酸盐条件生长；PPTtpp 被写为 ABC 转运且这些降解蛋白被并入 AND，需独立 QC；该反应三个条件的上下界均为 0")
    if n == 1891:gaps.append("旧 RUBISCO 反应的 GPR 将 AFE_2155 称 CbbL2；原始重组实验确认其为 form II CbbM；旧可逆方向中的逆向未获支持")
    if n in {1892,1893}:gaps.append("同株来源重组实验确认为 CbbQ2/CbbO2 活化复合体；不得并入 RUBISCO 羧化催化 GPR；活化反应化学式尚未模型化")
    if 1910 <= n <= 1913:gaps.append("Osorio 2013 同株厌氧硫/铁实验支持 SreABCD 表达及 H2S 生成推断；旧 SQRED1 将 SreABCD 并入硫化物:醌氧化 GPR，反应三个条件均关闭；底物方向、醌、膜侧未证，需独立 QC")
    if 1885 <= n <= 1890:gaps.append("旧 HYD4 写作甲酸氧化产 H2；Valdés 2008 仅以同源性推断该氢化酶簇和可能反应，Drobner 1990 的同株 H2 生长未鉴定本簇；精确反应与方向待核")
    if n in {1988,1989,1990,1991,1992,1994,1995,1996,1997,2019,2020,2021,2022,2023,2024}:gaps.append("旧 FE32tpp 为 Fe3+ 外膜转运，GPR 的 AND 分支含磷酸酶/调节/其他非典型转运成分；TonB/Exb 家族注释不能确定 Fe3+ 底物或全部必需亚基，需独立 QC")
    if n == 1774:gaps.append("旧 G6PDH2 被写成单向 6PGL+NADPH→G6P+NADP；EC 1.1.1.49 标准反应为 G6P+NADP→6PGL+NADPH，方向与底物异构体需独立 QC")
    if n == 1985:gaps.append("旧 DHORD 以 O2 为电子受体并产 H2O2；当前同版 product、UniProt/KEGG EC 1.3.5.2 及 IUBMB 标准定义均指向醌受体，须复核反应式与醌种类")
    if n == 2012:gaps.append("2016 PSD120 将 AFE_2292 关联到磷脂酰丝氨酸脱羧；2008 与当前 RefSeq 均将相同 WP_012537105.1 标为 TonB 外膜受体，旧 GPR 基因功能冲突，优先 QC")
    if n == 1808:gaps.append("旧 ACCOAC 将四个不同羧化酶组分写为 OR；IUBMB EC 6.4.1.2 说明分体酶需载体、biotin carboxylase 和 carboxytransferase 协同；旧 BIOC1 一分子羧化载体变成两分子未羧化载体，化学计量需 QC")
    if n == 1804:gaps.append("Li 2025 精确使用当前 WP_012536942.1 纯化并解析含 PLP 结构，半胱氨酸/Cd2+ 体系有 CdS 生成；未测旧 HACT 的 O-乙酰高丝氨酸反应或完整半胱氨酸脱硫化学式")
    if n == 2038:gaps.append("当前 product 是 alpha-D-glucose phosphate-specific phosphoglucomutase，UniProt/KEGG EC 5.4.2.2 和 Rhea:23536 均为 alpha-D-G1P↔alpha-D-G6P；旧 PGMT 使用 g1p-B/g6p-B，β 异构体关联待核；PGMT2 的 g1p 未标异构体")
    if n == 1905:gaps.append("Ribeiro 2011 文中 Afe_2172 Hsp20 的两条 qPCR 引物精确落在当前基因组 796313/796390 的 Hsp20 区域，不在当前 RS10045（旧 AFE_2172）；此文不得作为本基因功能证据")
    for h in historic:gaps.append(f"历史候选 {h[17]}：{h[26]}；边界：{h[25]}")
    rid = f"G57-{n}"
    boundary = "直接核实同版身份及模型记录；尚未证实本基因的精确酶反应"
    if 1999 <= n <= 2005:boundary = "同株生长与 phnG—phnM 共转录有原始实验；单个 Phn 酶反应和膦酸盐转运未获直接测定"
    if n == 1891:boundary = "同株基因来源重组 CbbM 羧化活性有直接实验；旧反应逆向与体内通量未证实"
    if n in {1892,1893}:boundary = "同株基因来源重组 CbbQ2/CbbO2 活化作用有直接实验；具体模型活化化学式与体内必需性未证实"
    if 1910 <= n <= 1913:boundary = "同株厌氧硫/铁培养中 SreABCD 表达有组学证据；该复合体的精确硫还原反应未直接测定"
    if 1885 <= n <= 1890:boundary = "同株 H2 生长有生理实验，本氢化酶簇的甲酸→CO2+H2 仅有同源性推断，未获催化实测"
    if n in {1988,1989,1990,1991,1992,1994,1995,1996,1997,2019,2020,2021,2022,2023,2024}:boundary = "TonB/Exb 等成分为基因组预测；FE32tpp 的 Fe3+ 底物、膜侧及旧 AND 分支成员未经单簇验证"
    if n == 1774:boundary = "该 WP 注释为 G6PDH；旧 G6PDH2 的单向反应与 IUBMB EC 1.1.1.49 标准方向相反，体内方向未测"
    if n == 1985:boundary = "当前 WP 的醌依赖型功能与 EC 1.3.5.2 一致；旧 DHORD 的 O2/H2O2 电子受体写法缺本株实验支持"
    if n == 2012:boundary = "2008 和 2026 RefSeq 同一 WP 均是 TonB 外膜受体注释；旧 PSD120 的磷脂脱羧 GPR 与蛋白类别冲突；缺单蛋白功能实验"
    if n == 1808:boundary = "当前 WP 为乙酰辅酶A羧化酶羧基转移 β 亚基；旧 ACCOAC 的四组分 OR GPR 与分体酶协同作用冲突，BIOC1 的 bcCP 载体计数需核"
    if n == 1804:boundary = "当前 WP 重组蛋白和 PLP 晶体结构有直接实验；半胱氨酸/Cd2+ 体系产 CdS；旧 HACT 的精确底物/产物未测"
    if n == 2038:boundary = "当前 WP 注释和 EC/Rhea 候选均指向 α-D 葡萄糖磷酸互变；旧 PGMT 的 β-D 底物与产物异构体未获基因专属证明"
    if n == 1905:boundary = "旧 AFE 编号文献命中经引物定位排除；当前 WP 为丝氨酸水解酶结构域蛋白，具体底物未证实"
    lit_key = afe_claim if afe_claim != "未确认" else wp
    lit = literature.get(lit_key, {})
    current_lit = current_literature.get(locus, {})
    tc = tcdb_review.get(locus, {})
    rhea_ids = list(dict.fromkeys(re.findall(r"RHEA:\d+", mentions["Rhea crossref"])))
    rhea_eq = "; ".join(f"{ri}: {rhea_by_id[ri]['Equation']}" for ri in rhea_ids if ri in rhea_by_id)
    base = {"总序号":n,"RU820":locus,"WP":wp,"染色体":seqid,"起点":int(start),"终点":int(end),"链向":strand,"当前 product":product,"当前 gene symbol":ga.get("gene", ""),"CDS 坐标核对":"通过","旧 AFE":afe_claim,"映射状态":map_status,"映射方法":map_method,"2016 GPR 反应":";".join(ids),"历史候选ID":";".join(str(h[17]) for h in historic),"历史候选动作":";".join(str(h[26]) for h in historic),"历史关联反应":";".join(str(h[19]) for h in historic if h[19]),"是否找到精确 reaction":"是；RUBISCO 正向羧化" if n == 1891 else "否；仅获得候选交叉引用" if metabolism_hint else "否","动作":action,"审查状态":"身份/映射/旧模型索引已核；完整证据链未完成","关联行 ID":rid,"未决点":"；".join(gaps),"UniProt accession":u.get("Entry", ""),"UniProt reviewed":u.get("Reviewed", ""),"UniProt protein name":u.get("Protein names", ""),"InterPro accession":u.get("InterPro", ""),"Pfam accession":u.get("Pfam", ""),"UniProt EC":u.get("EC number", ""),"KEGG KO":mentions["KEGG KO"],"KEGG EC":mentions["KEGG EC"],"KEGG reaction":mentions["KEGG reaction"],"BioCyc reaction":mentions["BioCyc reaction"],"Rhea crossref":mentions["Rhea crossref"],"Rhea 方程候选（EC级）":rhea_eq,"BRENDA":mentions["BRENDA"],"TCDB 家族候选":tc.get("最佳TCDB ID", ""),"TCDB 家族支持":tc.get("家族同源支持", ""),"TCDB 同一性":tc.get("同一性", ""),"TCDB 当前蛋白覆盖率":tc.get("当前蛋白覆盖率", ""),"EuropePMC 命中数":lit.get("hitCount", ""),"EuropePMC DOI线索":";".join(p["doi"] for p in lit.get("papers", []) if p.get("doi")),"EuropePMC 查询URL":lit.get("url", ""),"当前编号 RU820 命中":current_lit.get("RU820命中", ""),"当前编号 WP 命中":current_lit.get("WP命中", ""),"当前编号 DOI线索":current_lit.get("WP DOI", "") or current_lit.get("RU820 DOI", ""),"Source ID":";".join(source_ids),"NCBI protein URL":f"https://www.ncbi.nlm.nih.gov/protein/{wp}","UniProt URL":f"https://www.uniprot.org/uniprotkb/{u['Entry']}/entry" if u else "","证据边界":boundary}
    index.append(base)
    row = dict(zip(alpha_headers, [""]*30))
    row.update({"Record class":"2026 gene review","Candidate ID":rid,"Candidate group":"G57","2016 linked Reaction ID":";".join(ids),"2026 reaction object":"RUBISCO 既有羧化反应" if n == 1891 else "CbbQ2/CbbO2 Rubisco 活化，具体反应待定义" if n in {1892,1893} else f"当前注释：{product}；精确反应待核" if metabolism_hint else f"当前注释：{product}；未发现可确认的代谢反应","Legacy AFE locus":afe_claim,"Current locus":locus,"Primary evidence DOI":"10.1038/ncomms9883" if 1891 <= n <= 1893 else "10.1128/AEM.02101-07" if 1999 <= n <= 2005 else "","Evidence type":"重组蛋白活性；同版 RefSeq 身份；历史映射" if 1891 <= n <= 1893 else "同株生长、共转录；同版 RefSeq 身份；历史映射；2016 模型索引" if 1999 <= n <= 2005 else "同版 RefSeq 身份；历史映射；数据库预测；2016 模型索引","Evidence boundary":base["证据边界"],"2026 action":action,"2026 reaction score":"4（仅 RUBISCO 正向羧化）" if n == 1891 else "不评分","Source ID":base["Source ID"],"Score scope":"同株来源重组 CbbM 的精确正向羧化；不覆盖逆向、区室或其他同工酶" if n == 1891 else "活化反应未精确定义，不能赋分" if n in {1892,1893} else "未定义精确反应，不能赋分"})
    if n == 1773:row.update({"Primary evidence DOI":"10.3389/fmicb.2016.01365","Evidence type":"同株表达变化；历史候选 B-004","Evidence boundary":"转录变化支持表达；6PG 氧化脱羧的 Ru5P、CO2 和 NAD(P) 产物未直接测定","2026 reaction object":"6PG 氧化脱羧为候选；旧 PGDH 2 是另一脱水反应","Score scope":"原 B-004 候选分数 2 仅作历史记录；本轮未确定精确反应，不转分"})
    if 1910 <= n <= 1913:row.update({"Primary evidence DOI":"10.1128/AEM.03057-12","Evidence type":"同株硫/铁培养、表达组学及 H2S 生理线索","Evidence boundary":boundary,"2026 reaction object":"SreABCD 硫还原候选；精确反应待定义","Score scope":"旧 SQRED1 分数 3 不能转给 Sre 硫还原候选"})
    if 1885 <= n <= 1890:row.update({"Primary evidence DOI":"10.1186/1471-2164-9-597;10.1128/AEM.56.9.2922-2923.1990","Evidence type":"同株基因组同源性推断；同株 H2 生长生理（未指认本基因）","Evidence boundary":boundary,"2026 reaction object":"旧 HYD4 甲酸氧化产 H2 假说；精确反应待核","Score scope":"2016 HYD4 原分 2 保留历史；同株 H2 生长不提升该具体反应评分"})
    if n in {1988,1989,1990,1991,1992,1994,1995,1996,1997,2019,2020,2021,2022,2023,2024}:row.update({"Primary evidence DOI":"10.1186/1471-2164-9-597","Evidence type":"同株基因组转运系统推断；2016 模型 GPR","Evidence boundary":boundary,"2026 reaction object":"FE32tpp 旧 Fe3+ 外膜转运，复合体成员待核","Score scope":"旧 FE32tpp 分数 2 不确认每个 AND 成员或确切 Fe3+ 底物"})
    if n == 1774:row.update({"Evidence type":"当前 WP 产品与 IUBMB 标准酶反应对照","Evidence boundary":boundary,"2026 reaction object":"G6PDH2 旧单向方程待校正；候选标准式 G6P+NADP→6PGL+NADPH","Score scope":"旧分 2 不证明逆向单向约束；具体体内方向待核"})
    if n == 1985:row.update({"Evidence type":"同版 WP 产品；UniProt/KEGG EC 1.3.5.2；IUBMB 标准反应","Evidence boundary":boundary,"2026 reaction object":"DHORD 旧 O2/H2O2 方程待校正；候选醌→醌醇式","Score scope":"旧分 2 不支持氧作为直接电子受体；本株醌种类未定"})
    if n == 2012:row.update({"Evidence type":"同一 WP 的 2008/2026 RefSeq 一致标注 TonB 受体；2016 PSD120 GPR","Evidence boundary":boundary,"2026 reaction object":"PSD120 与 RU820_RS10585 的 GPR 关联待撤销核验","Score scope":"旧分 2 不可转为该 WP 的磷脂脱羧证据"})
    if n == 1808:row.update({"Evidence type":"当前 WP 羧基转移 β 亚基；IUBMB EC 6.4.1.2；旧 ACCOAC/BIOC1 对照","Evidence boundary":boundary,"2026 reaction object":"乙酰辅酶A羧化酶旧 GPR 与 BIOC1 载体计数待修订核验","Score scope":"旧 ACCOAC 2 分与 BIOC1 3 分不能证明四组分 OR 或载体复制"})
    if n == 1804:row.update({"Primary evidence DOI":"10.1039/d4cb00195h","Evidence type":"同株来源 WP 的重组蛋白纯化、PLP 晶体结构与体外 CdS 生成体系","Evidence boundary":boundary,"2026 reaction object":"O-succinylhomoserine sulfhydrylase；半胱氨酸脱硫体系活性；旧 HACT 化学式待核","Score scope":"可确认蛋白身份/体系活性，未确认独立、完整的代谢反应式；旧 HACT 分 2 不升分"})
    if n == 2038:row.update({"Evidence type":"当前 WP product；UniProt/KEGG EC 5.4.2.2；Rhea:23536 反应式","Evidence boundary":boundary,"2026 reaction object":"α-D-G1P↔α-D-G6P 候选；旧 PGMT β 异构体待核","Score scope":"数据库标准方程与当前注释一致，仍无同株精确酶测定；旧分数不自动改变"})
    if n == 1905:row.update({"Evidence type":"当前同版 RefSeq；同号文献引物定位排除","Evidence boundary":boundary,"2026 reaction object":"未确定代谢反应；同号 Hsp20 论文不能归属此 WP","Score scope":"无可接受精确反应；不得转移论文实验"})
    review.append(row)
    for j, original in enumerate(linked, 1):
        rr = {h:original[k] if original[k] is not None else "" for k,h in enumerate(alpha_headers[:16])}
        rr.update({"Record class":"2026 reaction review - existing","Candidate ID":f"{rid}-R{j:02d}","Candidate group":"G57","2016 linked Reaction ID":original[0],"2026 reaction object":"2016 既有反应；基因对应关系待审","Legacy AFE locus":afe_claim,"Current locus":locus,"Primary evidence DOI":"10.1128/AEM.02101-07" if 1999 <= n <= 2005 else "","Evidence type":"2016 原始补表及 B1 精确映射；功能仍待独立核实","Evidence boundary":boundary if 1999 <= n <= 2005 else "直接核实 2016 反应字段和历史基因对应；尚未证实当前蛋白催化此精确反应","2026 action":"独立 QC 问题；证据不足暂缓" if 1999 <= n <= 2005 and original[0] == "PPTtpp" else "证据不足暂缓","2026 reaction score":"不评分","Source ID":"G57-S01;G57-S02;G57-S04;G57-S09" if 1999 <= n <= 2005 else "G57-S01;G57-S02;G57-S04","Score scope":"2016 原分保留在 D 列；本轮尚无足够基因专属反应证据"})
        if n == 1891 and original[0] == "RUBISCO":rr.update({"2026 reaction object":"h2o[c] + co2[c] + rb15bp[c] -> 2 h[c] + 2 3pg[c]（旧式正向）","Primary evidence DOI":"10.1038/ncomms9883","Evidence type":"同株基因来源的重组 CbbM 羧化活性实验；2016 旧反应","Evidence boundary":boundary,"2026 action":"仅更新证据或编号：AFE_2155 蛋白标签改 CbbM；不修改反应式","2026 reaction score":"4（仅正向羧化）","Source ID":"G57-S01;G57-S02;G57-S04;G57-S10","Score scope":"反应正向羧化；原 2016 分数为 3；逆向及其他同工酶不继承 4"})
        if 1910 <= n <= 1913 and original[0] == "SQRED1":rr.update({"Primary evidence DOI":"10.1128/AEM.03057-12","Evidence type":"同株厌氧硫/铁表达与生理线索；2016 原反应","Evidence boundary":boundary,"2026 action":"独立 QC 问题；证据不足暂缓","2026 reaction score":"不评分（Sre 候选）","Source ID":"G57-S01;G57-S02;G57-S04;G57-S14","Score scope":"D 列 3 是旧 SQRED1 原分；不证明 SreABCD 催化该氧化反应"})
        if 1885 <= n <= 1890 and original[0] == "HYD4":rr.update({"Primary evidence DOI":"10.1186/1471-2164-9-597;10.1128/AEM.56.9.2922-2923.1990","Evidence type":"氢化酶簇同源性预测；同株 H2 生长未指认本基因","Evidence boundary":boundary,"2026 action":"独立 QC 问题；证据不足暂缓","Source ID":"G57-S01;G57-S02;G57-S04;G57-S16;G57-S17","Score scope":"旧分 2 对应历史甲酸产氢假说；不据 H2 生长升分"})
        if n in {1988,1989,1990,1991,1992,1994,1995,1996,1997,2019,2020,2021,2022,2023,2024} and original[0] == "FE32tpp":rr.update({"Primary evidence DOI":"10.1186/1471-2164-9-597","Evidence type":"同株基因组 TonB/Exb 系统预测；旧 GPR 对照","Evidence boundary":boundary,"2026 action":"独立 QC 问题；证据不足暂缓","Source ID":"G57-S01;G57-S02;G57-S04;G57-S16","Score scope":"旧分 2 未确认 Fe3+ 底物、转运方向或每一 AND 成员"})
        if n in {1774,1808,1985,2012}:rr.update({"Evidence boundary":boundary,"2026 action":"独立 QC 问题；反应/基因关联待修订","Source ID":";".join(source_ids),"Score scope":"2016 D 列原分保留；其方程或基因关联不能由现有证据直接确认"})
        if n == 1808:rr.update({"Evidence type":"当前羧基转移 β 亚基；IUBMB 标准机制；旧模型复合体/载体计数对照","2026 reaction object":"ACCOAC 四组分 OR 与 BIOC1 载体计数冲突，待校正"})
        if n == 1804 and original[0] == "HACT":rr.update({"Primary evidence DOI":"10.1039/d4cb00195h","Evidence type":"精确 WP 的重组蛋白/PLP 结构及半胱氨酸-Cd2+ 体系活性","Evidence boundary":boundary,"2026 action":"独立 QC 问题；旧 HACT 底物/方向待核","Source ID":";".join(source_ids),"Score scope":"旧分 2 不因 CdS 体系实验提升；未测旧式 O-乙酰高丝氨酸反应"})
        if n == 2038 and original[0] in {"PGMT","PGMT2"}:rr.update({"Evidence type":"当前 alpha-D 特异 product；UniProt/KEGG EC 5.4.2.2；Rhea:23536","Evidence boundary":boundary,"2026 action":"独立 QC 问题；异构体与基因关联待核","Source ID":";".join(source_ids),"Score scope":"旧反应分数不证明 β-D 异构体被此 WP 催化" if original[0] == "PGMT" else "旧 g1p 未明确异构体；与 α-D 标准式需逐项核对"})
        if n == 1774:rr.update({"Evidence type":"当前 G6PDH 产品与 IUBMB EC 1.1.1.49 方向对照","2026 reaction object":"G6PDH2 方向冲突：旧式仅允许 6PGL→G6P；标准氧化方向为 G6P→6PGL"})
        if n == 1985:rr.update({"Evidence type":"当前醌依赖型 WP；UniProt/KEGG EC 1.3.5.2；IUBMB 标准式","2026 reaction object":"DHORD 受体冲突：旧 O2/H2O2；候选醌/醌醇"})
        if n == 2012:rr.update({"Evidence type":"2008 与 2026 同一 WP 均为 TonB 外膜受体；2016 旧 GPR","2026 reaction object":"PSD120 旧 GPR 与蛋白类别冲突；待核对真正脱羧酶"})
        reaction_rows.append(rr)
    if n == 1791:
        old = model_by_id["ATPS5rpp"]
        rr = {h:old[k] if old[k] is not None else "" for k,h in enumerate(alpha_headers[:16])}
        rr.update({"Record class":"2026 reaction review - historical candidate","Candidate ID":f"{rid}-B024","Candidate group":"G57","2016 linked Reaction ID":"ATPS5rpp","2026 reaction object":"历史 B-024 拟让 AFE_2047 替代 γ；该基因不在 2016 原 GPR","Legacy AFE locus":afe_claim,"Current locus":locus,"Evidence type":"历史候选 B-024；2016 原反应字段","Evidence boundary":"旧 ATP 合酶整体功能有生理证据；AFE_2047 独立替代 γ 与 5 H+/ATP 未被证明","2026 action":"证据不足暂缓；保留旧反应与原 GPR","2026 reaction score":"不评分（候选 GPR）","Source ID":"G57-S01;G57-S02;G57-S04;G57-S12","Score scope":"D 列 3 是旧主反应原分，不转移至 AFE_2047 OR"})
        reaction_rows.append(rr)

index_fields = list(index[0])
write_tsv("G57_逐基因索引与状态.tsv", index, index_fields)
write_tsv("G57_Alpha逐基因审查.tsv", review, alpha_headers)
write_tsv("G57_Alpha既有反应逐项对照.tsv", reaction_rows, alpha_headers)
write_tsv("G57_来源表.tsv", sources, list(sources[0]))
chain = []
for r in index:
    n = r["总序号"]
    transport = bool(re.search(r"transport|permease|antiporter|porter|receptor|porin|channel|efflux",r["当前 product"],re.I))
    ids = [f"KEGG KO {r['KEGG KO']}" if r["KEGG KO"] else "", f"KEGG R {r['KEGG reaction']}" if r["KEGG reaction"] else "", f"BioCyc {r['BioCyc reaction']}" if r["BioCyc reaction"] else "", f"Rhea {r['Rhea crossref']}" if r["Rhea crossref"] else ""]
    direct = "；同株来源实验已核，边界见逐基因审查" if n in {1804,1891,1892,1893,1910,1911,1912,1913,1999,2000,2001,2002,2003,2004,2005} else "；旧编号同名论文已用引物排除" if n == 1905 else "；命中 DOI 初判见文献表，本基因专属实验仍待核"
    chain.append({"总序号":n,"RU820":r["RU820"],"WP":r["WP"],"1 身份":"通过：同版 gene/CDS/WP/坐标/链向/product 一致","2 历史映射":f"{r['旧 AFE']}；{r['映射状态']}；{r['映射方法']}","3 蛋白功能":f"UniProt {r['UniProt accession']} ({r['UniProt reviewed']})；InterPro {r['InterPro accession']}；Pfam {r['Pfam accession']}；仍须保守位点/功能边界核验" if r["UniProt accession"] else "同株 UniProt RefSeq 交叉引用未命中；仅有 NCBI product，功能未独立核实","4 反应通路":"；".join(x for x in ids if x) + "；数据库链接不算本基因催化验证" if any(ids) else "D3-1-R 未提供基因专属精确反应；仍需功能检索","4 TCDB":f"家族候选 {r['TCDB 家族候选']}；阈值支持 {r['TCDB 家族支持']}；同一性 {r['TCDB 同一性']}、覆盖率 {r['TCDB 当前蛋白覆盖率']}；底物/耦联/方向未证" if r['TCDB 家族候选'] else "转运类：TCDB 家族未获可靠同源支持；底物/耦联未证" if transport else "非明显转运类；不适用","5 原始论文":f"旧编号命中 {r['EuropePMC 命中数']}；当前 RU/WP 命中 {r['当前编号 RU820 命中']}/{r['当前编号 WP 命中']}{direct}","6 化学与定位":"RUBISCO 正向羧化有重组蛋白证据；逆向和体内区室未证" if n == 1891 else "Rhea EC 方程仅作候选；基因关联、方向、区室、原子/电荷与辅因子未逐项确证" if r["Rhea crossref"] else "未定义可接受的基因专属精确反应；化学字段暂缓","7 2016 对照":f"原 GPR：{r['2016 GPR 反应']}；完整 A—P 见反应表" if r["2016 GPR 反应"] else f"历史候选：{r['历史关联反应']}；原 A—P 见反应表" if r["历史关联反应"] else "未命中原 GPR；不能据此判旧模型无同化学反应","8 动作裁决":r["动作"],"9 复核":"身份与表间引用已复核；功能/化学证据链未全部验收","Source ID":r["Source ID"]})
write_tsv("G57_九项证据链状态.tsv",chain,list(chain[0]))
missing = [r for r in index if r["动作"].startswith("证据不足")]
write_tsv("G57_冲突缺证与待验证.tsv", [{"总序号":r["总序号"],"RU820":r["RU820"],"WP":r["WP"],"2016反应":r["2016 GPR 反应"],"待验证":r["未决点"],"初步动作":r["动作"]} for r in index], ["总序号","RU820","WP","2016反应","待验证","初步动作"])
payload = {"index":index,"chain":chain,"gene_review":review,"reaction_review":reaction_rows,"sources":sources,"literature_triage":read_tsv(LIT_TRIAGE) if LIT_TRIAGE.exists() else [],"primer_audit":read_tsv(PRIMER_AUDIT) if PRIMER_AUDIT.exists() else [],"ec_model_candidates":read_tsv(EC_MODEL) if EC_MODEL.exists() else [],"current_id_search":read_tsv(CURRENT_LIT) if CURRENT_LIT.exists() else [],"tcdb_review":read_tsv(TCDB_REVIEW) if TCDB_REVIEW.exists() else [],"counts":dict(Counter(r["动作"] for r in index)),"date":dt.date.today().isoformat()}
(OUT / "g57_data.json").write_text(json.dumps(payload,ensure_ascii=False,default=str),encoding="utf-8")
counts = Counter(r["动作"] for r in index)
report = f"""# G57 阶段审查与验收记录

日期：{dt.date.today().isoformat()}。范围：1759—2051，293 个当前蛋白编码基因。

## 已完成的核对

- 同版 GCF_049532655.1 GFF 中逐条核对 gene/CDS、RU820、WP、坐标、链向和 product；293 条连续，RU820 和 WP 均唯一。
- 与 B1 旧新映射及 D3-1-R 跨库记录逐基因连接；有歧义的旧 AFE 标为“未确认”。
- 逐个查询 Alpha 中保留的 2016 原始反应行，并将完整 A—P 字段复制到对照表。历史反应分数原样保留；本轮未给候选反应打分。
- 每个基因有一条 Alpha Q—AD 逐基因审查行；既有反应另有 {len(reaction_rows)} 条逐项对照行。来源表 {len(sources)} 条来源链，均标明自动注释与模型的证据限制。
- UniProt 同株 taxon 243159 全量导出后以同版 WP 精确连接，本组有 {sum(bool(uniprot.get(cds[x[0]][0][1].get('protein_id',''))) for x in part)} 条获得 UniProt accession，并记录 reviewed、InterPro、Pfam、EC 及原始记录链接；没有命中的不推断为无蛋白。
- Europe PMC 对 293 个已确认 AFE 或当前 WP 逐条检索，{sum(bool(literature.get((cross.get(x[0],{}).get('旧AFE locus') or cds[x[0]][0][1].get('protein_id','')),{}).get('hitCount',0)) for x in part)} 条有检索命中。命中仅作论文线索；零命中不证明无文献，未读原文不得升分。
- 对 293 个基因的 RU820 与 WP 当前编号另做 586 次 Europe PMC 检索，0 次错误、1 个 WP 命中；并入旧编号线索后共 17 个 DOI，逐篇初判异株、异种、异源表达、组学和编号误命中的使用边界。编号检索仍不等于完整文献穷尽。
- 当前无旧 GPR 的基因与 2016 原模型按完整 EC 号交叉，仅得到 5 个同 EC 候选；同 EC 不能证明同反应或遗漏 GPR，已单独列为待核线索。
- 44 个当前 product 可锚定转运家族的 WP 与 TCDB 参考序列做家族内同源比对，41 条达到预设的同一性和覆盖率阈值，3 条未达阈值。该结果只支持家族候选，不支持具体底物、转运方向或能量耦联；另有 2 个含运输/ATPase 字样但无适当家族锚点，未强行比对。

## 已核实的重点冲突

1999—2005（RU820_RS10520—RS10550）对应旧 AFE_2278—2284，即 phnG—phnM。Vera 等 2008 年原文使用早期 TIGR 编号 AFE-811—821；其表 1 的 8 条 phnG—phnM 引物全部精确匹配当前基因组 RU820_RS10520—RS10550 的对应区域，完成了论文对象与当前身份的序列核验。该研究在 ATCC 23270 测得甲基/乙基膦酸盐作为唯一磷源时的生长、phnG—phnM 共转录和甲基膦酸盐条件下的转录诱导（DOI: 10.1128/AEM.02101-07）。2016 PPTtpp 把这些降解系统蛋白与后续两个蛋白全列为 ABC 转运的 AND 条件，反应三个培养条件上下界又均为 0。这是需要单独复核的模型 QC 与 GPR 问题；论文未直接测定单酶反应或转运耦联，现不提出正式替代反应。

1891—1893（RU820_RS09975—RS09985）经 B1 映射到 AFE_2155—2157。Tsai 等 2015 年用 ATCC 23270 来源基因表达重组 CbbM、CbbQ2、CbbO2，测得 Rubisco 羧化活性及 CbbQO 活化作用（DOI: 10.1038/ncomms9883）。旧 RUBISCO 将 AFE_2155 蛋白称为 CbbL2，应更正为 CbbM；其正向羧化由 2016 分数 3 提议升至 4，逆向不继承该评分。CbbQ2/CbbO2 是活化因子，不写入羧化反应 GPR；当前 RU820_RS09985 的 NorD product 与同株来源实验的 CbbO2 身份冲突，建议更新证据注释。

1910—1913（RU820_RS10070—RS10085）经 B1 映射到 AFE_2177—2179/2181。Osorio 等 2013 年在 ATCC 23270T 的 S0/Fe3+ 厌氧与 S0/O2 有氧培养中测得 SreABCD 表达差异，并观察与 H2S 生成相符的生理现象（DOI: 10.1128/AEM.03057-12）。旧 SQRED1 将该组并入“硫化物:醌氧化”GPR，且三个条件上下界均为 0。原论文没有逐酶确认氧化方向、醌种类与膜侧；该关联列为独立 QC 问题，暂不更改原反应或新增硫还原反应。

1885—1890（RU820_RS09945—RS09970）的旧 HYD4 把甲酸氧化与产 H2 连在同一反应。Valdés 等 2008 年以基因组同源性提出该氢化酶簇和可能的甲酸产氢联系；Drobner 等 1990 年测得 ATCC 23270 使用 H2 生长，但未指认本簇。旧 2 分保留为历史记录，具体甲酸产氢反应与方向仍须实测（DOI: 10.1186/1471-2164-9-597；10.1128/AEM.56.9.2922-2923.1990）。

1988—1997 及 2019—2024 命中的旧 FE32tpp 是 Fe3+ 外膜转运模型反应；其 GPR 分支混入磷酸酶、调节蛋白等非典型转运成分。Valdés 等 2008 年提供 TonB/Exb 家族层面的基因组预测，未测 Fe3+ 专属底物或所有 AND 成员的必要性（DOI: 10.1186/1471-2164-9-597）。该 GPR 列为独立 QC 问题，暂不改写。

另有三处明确的旧模型一致性冲突。1774（RU820_RS09360，旧 AFE_2025）的 G6PDH2 被限定在 6PGL→G6P，而 IUBMB EC 1.1.1.49 的标准氧化方向为 G6P→6PGL；1985（RU820_RS10450，旧 AFE_2264）的 DHORD 写作 O2→H2O2，而当前 product、UniProt/KEGG 的 EC 1.3.5.2 与 IUBMB 标准式都指向醌受体；2012（RU820_RS10585，旧 AFE_2292）的 PSD120 将同一 WP_012537105.1 外膜 TonB 受体指认为磷脂酰丝氨酸脱羧酶，旧版和新版 RefSeq 注释均显示类别冲突。这三项进入独立 QC，保留原式作历史对照，具体修订需化学式和功能复核。

1808（RU820_RS09540，旧 AFE_2067）的当前产物是乙酰辅酶A羧化酶羧基转移 β 亚基。旧 ACCOAC 把四个不同组分写成 OR，旧 BIOC1 又把一分子羧化载体转成两分子未羧化载体。IUBMB EC 6.4.1.2 的分体酶说明包含载体、biotin carboxylase 和 carboxytransferase 的协作。这两条旧反应的复合体逻辑与载体计数列为独立 QC，原 2/3 分不转移至修订候选。

1804（RU820_RS09520，WP_012536942.1）的当前编号补查命中 Li 等 2025 年研究（DOI: 10.1039/d4cb00195h）。作者将 ATCC 23270 来源蛋白在 E. coli 表达、纯化，获得含 PLP 的晶体结构 PDB 9IQH，并在半胱氨酸/Cd2+ 体外体系中观察到 CdS 生成。该证据支持蛋白身份及体系活性。旧 HACT 写作 O-乙酰高丝氨酸与 H2S 形成高半胱氨酸，论文未测此反应；半胱氨酸脱硫的完整化学计量也未逐项测得。本轮补充实验来源，旧反应继续 QC 暂缓。

2038（RU820_RS10725，WP_012537123.1）的当前 product 写明 α-D 葡萄糖磷酸特异，UniProt/KEGG EC 5.4.2.2 与 Rhea:23536 候选方程都是 α-D-G1P↔α-D-G6P。旧 PGMT 使用 g1p-B/g6p-B，另一个 PGMT2 使用未标异构体的 g1p 与 g6p-A；两条旧式的异构体及基因对应关系需独立核验，现不更改原模型分数。

1905（RU820_RS10045，旧 AFE_2172）的 Europe PMC 检索命中 Ribeiro 等 2011 年 Hsp20 论文。论文实验使用 A. ferrooxidans LR，文中“Afe_2172”的两条 qPCR 引物在当前 ATCC 23270 基因组精确落于 796313 与 796390 的另一 Hsp20 区域，不落在本组 RS10045。因此该文献的热休克实验不能归给 G57-1905；旧 AFE 字面相同不等于当前 WP 相同（DOI: 10.1186/1471-2180-11-259）。本轮只作证据排除，不审查范围外基因。

Pang 等 2020 年又将 ATCC 23270 来源 CbbM/CbbQ2/CbbO2 异源表达到 E. coli，观察到 Rubisco 活性与活化因子共同表达后的碳回收变化（DOI: 10.3389/fbioe.2020.543807）。这进一步支持 1891—1893 的蛋白身份，但实验宿主是 E. coli，不能据此确认 ATCC 23270 体内反应方向或区室。

## 与既有 45 项裁决重叠

本组命中 B-004（RU820_RS09355）、B-006（RU820_RS10385）、B-007（RU820_RS10670）、B-024（RU820_RS09450）。B-004 延续“仅证据更新”，没有新增 6PG 氧化脱羧反应；其旧候选 2 分不直接转成本轮精确反应分数。B-006、B-007、B-024 继续暂缓。B-024 的 ATPS5rpp 原 A—P 字段已在反应对照表另列；旧分数 3 属旧主反应，不证明 AFE_2047 可替代 γ 亚基。没有发现新证据足以推翻这 4 项历史动作。

TCDB 家族比对为 B-006 的糖孔蛋白找到 OprB 家族候选（1.B.19.1.2，同一性 0.3333，当前蛋白覆盖率 0.8606），为 B-007 的糖转运蛋白找到 MFS 家族候选（2.A.1.1.63，同一性 0.4123，覆盖率 0.9712）。跨物种最相似蛋白的具体底物不可转给本株，两个历史候选的暂缓裁决保持不变。

## 当前动作统计

{chr(10).join('- ' + k + '：' + str(v) for k,v in counts.items())}

统计时将 31 个仅凭当前 product 已提示酶或转运功能、但跨库未找到基因专属精确反应的基因改列“证据不足暂缓”。它们不能因未命中旧 GPR 而暂判无需修改。全部动作仍为阶段性结论。

## 验收结论

**尚未通过总说明的全量证据链验收。** 全部 293 个基因的身份已核对，但 UniProt、InterPro、Pfam、TCDB、原始论文及精确化学反应没有逐条完成。2016 反应行虽然逐项复制了原始字段，基因到精确反应的等价关系仍未全部验证。因此本批没有正式新增或修改 reaction/GPR，也没有把自动注释当成同株实验证据。

后续逐基因优先核查现有模型关联与代谢候选；完成底物、产物、原子/电荷、辅因子、方向、区室、GPR 及同株实验后再把暂缓项升级为具体动作。无代谢提示的“无需修改（暂定）”同样未通过完整功能审查。
"""
(OUT / "G57_阶段验收报告.md").write_text(report, encoding="utf-8")
print({"genes":len(index),"unique_locus":len(set(x["RU820"] for x in index)),"unique_wp":len(set(x["WP"] for x in index)),"reaction_rows":len(reaction_rows),"counts":dict(counts),"output":str(OUT)})
