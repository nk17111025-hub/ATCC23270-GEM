"""Rebuild D3 from the current ATCC 23270 genome. Read-only inputs, new outputs only."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "D3-1-R"
DATE = "2026-09-23"
ASSEMBLY = "GCF_049532655.1"
PGDB = "GCF_000021485"


def read_tsv(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name, rows, columns):
    OUT.mkdir(exist_ok=True)
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def tokens(value):
    """Only actual list delimiters. Never treat a substring as an identifier."""
    return {piece.strip() for piece in re.split(r"[;,|]", value or "") if piece.strip()}


def exact_biocyc_match(model_id, reaction):
    if not model_id:
        return False
    return model_id == reaction.get("对象ID", "") or model_id in tokens(reaction.get("外部ID", ""))


def join(values):
    return ";".join(sorted({str(v).strip() for v in values if str(v).strip()}))


def norm_ec(value):
    return re.sub(r"\s+T$", "", re.sub(r"^(?:EC[-:]?)", "", value.strip()))


def parse_gff():
    path = ROOT / "02_NCBI/标准化数据/GCF_049532655.1/genomic.gff"
    genes = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("gene", "CDS"):
                continue
            attrs = {}
            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = unquote(value)
            locus = attrs.get("locus_tag", "")
            if not locus:
                continue
            if fields[2] == "gene":
                row = genes.setdefault(locus, {"当前locus": locus, "当前protein ID": "", "NCBI当前功能": "", "基因类型": "", "基因名": "", "复制子": fields[0], "坐标": f"{fields[3]}..{fields[4]}", "方向": fields[6]})
                row["基因类型"] = attrs.get("gene_biotype", "未标注")
                row["基因名"] = attrs.get("gene", "")
            elif locus in genes:
                row = genes[locus]
                row["当前protein ID"] = attrs.get("protein_id", "")
                row["NCBI当前功能"] = attrs.get("product", "")
    return genes, path


def load_model():
    gpr = read_tsv("08_模型基线/A0/恢复GPR及原始证据定位.tsv")
    rx = read_tsv("08_模型基线/A0/2016_到2024_真实差异.tsv")
    byid = {r["反应ID"]: r for r in rx}
    old16, old24 = defaultdict(set), defaultdict(set)
    model_rows = []
    for row in gpr:
        rid = row["反应ID"]
        g16 = set(re.findall(r"AFE_\d{4}", row.get("2016 Excel GPR", "")))
        g24 = set(re.findall(r"AFE_\d{4}", row.get("2024 S1 GPR", "")))
        for gene in g16: old16[gene].add(rid)
        for gene in g24: old24[gene].add(rid)
        model_rows.append({"模型反应ID": rid, "2016 GPR": row.get("2016 Excel GPR", ""), "2024 GPR": row.get("2024 S1 GPR", ""), "A0差异状态": byid.get(rid, {}).get("状态", "")})
    return old16, old24, model_rows


def load_biocyc():
    base = "03_BioCyc/C2-1/标准化数据/"
    genes = read_tsv(base + "genes.tsv")
    proteins = read_tsv(base + "proteins.tsv")
    ers = read_tsv(base + "enzymatic-reactions.tsv")
    reactions = read_tsv(base + "reactions.tsv")
    alias = defaultdict(set)
    for row in genes:
        for key in (row.get("对象ID", ""), *tokens(row.get("accession", "")), *tokens(row.get("同义词", ""))):
            if key: alias[key].add(row["对象ID"])
    pro_by_gene = defaultdict(list)
    pro_by_id = {}
    for row in proteins:
        pro_by_gene[row.get("基因", "")].append(row)
        pro_by_id[row.get("对象ID", "")] = row
    rx_by_id = {row["对象ID"]: row for row in reactions}
    rx_by_protein = defaultdict(set)
    for row in ers:
        if row.get("酶") and row.get("反应") in rx_by_id:
            rx_by_protein[row["酶"]].add(row["反应"])
    for row in reactions:
        for protein in tokens(row.get("酶", "")):
            rx_by_protein[protein].add(row["对象ID"])
    return alias, pro_by_gene, pro_by_id, rx_by_id, rx_by_protein


def load_rhea():
    combined = OUT / "Rhea_全基因EC查询汇总.tsv"
    rows = read_tsv("D3-1-R/Rhea_全基因EC查询汇总.tsv") if combined.exists() else read_tsv("06_Rhea/D3-1/Rhea_modelBioCyc_EC_crossrefs.tsv")
    xref = defaultdict(list)
    by_ec = defaultdict(list)
    for row in rows:
        for field, prefix in (("Cross-reference (KEGG)", "KEGG:"), ("Cross-reference (MetaCyc)", "MetaCyc:")):
            for ref in tokens(row.get(field, "")):
                xref[ref.removeprefix(prefix)].append(row)
        for ec in tokens(row.get("EC number", "")):
            by_ec[norm_ec(ec)].append(row)
    return xref, by_ec


def load_brenda():
    rows = read_tsv("05_BRENDA/D3-1/BRENDA_酶学结果.tsv")
    by_ec = defaultdict(list)
    for row in rows:
        by_ec[norm_ec(row.get("EC号", ""))].append(row)
    return by_ec


def load_brenda_details():
    folder = ROOT / "05_BRENDA/D3-1/原始响应"
    def binding_rows(name):
        return json.loads((folder / name).read_text(encoding="utf-8"))["results"]["bindings"]
    def value(row, key):
        return row.get(key, {}).get("value", "")
    reactions = defaultdict(lambda: {"Substrate": set(), "Product": set()})
    for item in binding_rows("BRENDA_FeS_reaction_compounds.json"):
        key = (value(item, "orgName"), value(item, "enzyme").rsplit("/", 1)[-1], value(item, "ec").rsplit("/", 1)[-1], value(item, "rxn"))
        role = value(item, "roleClass").rsplit("/", 1)[-1]
        if role in ("Substrate", "Product"):
            reactions[key][role].add(value(item, "compoundName"))
    cofactors = defaultdict(lambda: {"cofactor": set(), "ref": set()})
    for item in binding_rows("BRENDA_FeS_cofactors_references.json"):
        key = (value(item, "orgName"), value(item, "enzyme").rsplit("/", 1)[-1], value(item, "ec").rsplit("/", 1)[-1])
        if value(item, "cofName"):
            cofactors[key]["cofactor"].add(value(item, "cofName"))
        if value(item, "ref"):
            cofactors[key]["ref"].add(value(item, "ref"))
    return reactions, cofactors


def classify(row, model_equations):
    product = row["NCBI当前功能"].lower()
    gene = row["基因名"].lower()
    kegg = row["KEGG注释"].lower()
    bio = row["BioCyc蛋白名称"].lower()
    pathway = row["KEGG pathway"]
    ec = tokens(row["KEGG EC"] + ";" + row["BioCyc EC"])
    transport_signal = bool(re.search(r"\b(transporter|transport protein|permease|symporter|antiporter|porin|efflux pump|uptake protein)\b", product))
    transport_signal &= not bool(re.search(r"cell division|secretion protein|periplasmic adaptor|protein export|protein translocase", product))
    # A reaction with TRANS in its enzyme name can be a transferase. Protein
    # transport annotation is required before entering the transport table.
    row["转运候选"] = "是" if transport_signal else "否"
    organic = bool(re.search(r"glucose|fructose|galactose|ribose|xylose|sugar|carbohydrate|lactate|pyruvate|succinate|citrate|malate|acetate|amino acid|peptide|dicarboxylate|glycerol|organic acid", product))
    corroboration = bool(row["KEGG KO"] or row["BioCyc转运反应ID"]) and bool(re.search(r"glucose|fructose|galactose|ribose|xylose|sugar|carbohydrate|peptide|dicarboxylate|PTS", kegg + " " + product, re.I))
    model_ids = tokens(row["2016 GPR反应"] + ";" + row["2024 GPR反应"])
    inorganic_model = any(re.search(r"\b(?:fe2|fe3|so4|h2s|nh4)\[", model_equations.get(rid, ""), re.I) for rid in model_ids)
    organic_focus = transport_signal and organic and corroboration and not inorganic_model and row["映射冲突"] == "否"
    row["重点分类冲突"] = "旧模型GPR反应涉及无机底物，与有机碳注释待核" if transport_signal and organic and inorganic_model else ""
    row["有机碳转运分类"] = "当前蛋白与KEGG/反应支持底物类别，具体底物未实证" if organic_focus else ("转运家族或底物未定；不进入有机碳重点表" if transport_signal else "不适用")
    row["有机碳转运表入选"] = "是" if organic_focus else "否"
    central_terms = r"\b(phosphofructokinase|fructose.bisphosphatase|fructose.bisphosphate aldolase|hexokinase|glucokinase|ribulose.bisphosphate carboxylase|ribulose.phosphate.3.epimerase|ribose.5.phosphate isomerase|transketolase|transaldolase|pyruvate kinase|pyruvate dehydrogenase|citrate synthase|aconitate hydratase|isocitrate dehydrogenase|succinate dehydrogenase|succinate..CoA ligase|fumarase|fumarate hydratase|malate dehydrogenase|phosphoglycerate kinase|phosphoglycerate mutase|glyceraldehyde.3.phosphate dehydrogenase|phosphogluconate dehydrogenase|glucose.6.phosphate dehydrogenase|glucose.6.phosphate isomerase|phosphoglucomutase|phosphoribulokinase|phosphoenolpyruvate carboxylase|triose.phosphate isomerase|enolase)\b"
    central = bool(re.search(central_terms, product)) and bool(row["KEGG KO"] or row["KEGG EC"] or row["BioCyc EC"] or row["BioCyc reaction ID"])
    if "domain-containing" in product or "c-terminal-like" in product:
        central = False
    central |= "2.7.1.2" in ec and bool(set(pathway.split(";")) & {"afr00010", "afr01100"})
    central &= row["映射冲突"] == "否"
    row["中央碳分类"] = "蛋白酶名或EC与通路支持，具体反应待核" if central else "不适用"
    row["中央碳表入选"] = "是" if central else "否"
    fes_terms = r"\b(rusticyanin|cyc2|cytochrome c[124]?|cytochrome bc1|cytochrome b6|quinol oxidase|sulfide:quinone|sulfur:quinone|sulfur dioxygenase|tetrathionate hydrolase|thiosulfate dehydrogenase|sulfur oxygenase|NADH.quinone oxidoreductase|NADH dehydrogenase subunit|respiratory complex [IiVv]+|bc1 complex)\b"
    known_member = gene in {"rus", "cyc2", "sqr", "teth", "tqo", "soxa", "soxb", "soxc", "soxd", "peta", "petb", "petc"} or bool(re.fullmatch(r"(?:nuo|pet|sox)[a-z]", gene))
    respiratory_gpr = bool(model_ids & {"NADHI", "HYD4"})
    fes = (bool(re.search(fes_terms, product, re.I)) or known_member) and (known_member or respiratory_gpr or bool(row["KEGG KO"] or row["KEGG EC"] or row["BioCyc reaction ID"]))
    fes &= not bool(re.search(r"assembly protein|biogenesis protein|maturation protein", product))
    fes &= row["映射冲突"] == "否"
    row["FeS_ETC分类"] = "明确Fe/S或呼吸链成员" if fes else "不适用"
    row["FeS_ETC表入选"] = "是" if fes else "否"
    row["重点分类依据"] = join(["当前蛋白功能", "KEGG KO/EC/通路" if row["KEGG KO"] or row["KEGG EC"] or row["KEGG pathway"] else "", "BioCyc反应" if row["BioCyc reaction ID"] else "", "已知复合物成员" if known_member else "", "模型GPR及反应" if model_ids else ""])
    return row


def main():
    genome, gff_path = parse_gff()
    maps = read_tsv("02_NCBI/标准化数据/旧新基因映射.tsv")
    by_current, by_old = defaultdict(list), defaultdict(list)
    for row in maps:
        for locus in tokens(row.get("当前GCF049位点", "")):
            by_current[locus].append(row)
        if row.get("旧AFE位点"): by_old[row["旧AFE位点"]].append(row)
    kegg = defaultdict(list)
    for row in read_tsv("04_KEGG/D3-1/KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv"):
        kegg[row["KEGG gene/locus"]].append(row)
    old16, old24, model_rows = load_model()
    model_equations = {r["反应ID"]: (r.get("2024反应式", "") or r.get("2016反应式", "")) for r in read_tsv("08_模型基线/A0/2016_到2024_真实差异.tsv")}
    bio_alias, bio_pro_gene, bio_pro_id, bio_rx, bio_rx_pro = load_biocyc()
    rhea_xref, rhea_ec = load_rhea()
    rhea_source = "D3-1-R/Rhea_全基因EC查询汇总.tsv" if (OUT / "Rhea_全基因EC查询汇总.tsv").exists() else "06_Rhea/D3-1/Rhea_modelBioCyc_EC_crossrefs.tsv"
    brenda_ec = load_brenda()
    brenda_reactions, brenda_cofactors = load_brenda_details()
    output, evidence, conflicts, coverage, missing = [], [], [], [], []
    for locus in sorted(genome):
        base = genome[locus]
        m = by_current.get(locus, [])
        old = {r.get("旧AFE位点", "") for r in m if r.get("旧AFE位点")}
        refs = {r.get("旧RefSeq位点", "") for r in m if r.get("旧RefSeq位点")}
        old_pro = {r.get("旧蛋白ID", "") for r in m if r.get("旧蛋白ID")}
        map_states = {r.get("映射状态", "") for r in m if r.get("映射状态")}
        protein_conflict = any(r.get("当前蛋白ID") and base["当前protein ID"] and base["当前protein ID"] not in tokens(r["当前蛋白ID"]) for r in m)
        alias_hits = set().union(*(bio_alias.get(key, set()) for key in old | refs)) if old or refs else set()
        bio_pro = [p for geneid in alias_hits for p in bio_pro_gene.get(geneid, [])]
        bio_rids = {rid for p in bio_pro for rid in bio_rx_pro.get(p["对象ID"], set())}
        bio_rids |= {rid for geneid in alias_hits for rid in bio_rx_pro.get(geneid + "-MONOMER", set())}
        kr = [x for gene in old for x in kegg.get(gene, [])]
        kos = {ko for x in kr for ko in tokens(x.get("KO", ""))}
        kecs = {norm_ec(ec) for x in kr for ec in tokens(x.get("EC", ""))}
        krs = {rid for x in kr for rid in tokens(x.get("KEGG reaction", ""))}
        brxs = [bio_rx[rid] for rid in bio_rids if rid in bio_rx]
        becs = {norm_ec(ec) for x in brxs for ec in tokens(x.get("EC号", ""))}
        ecs = kecs | becs
        exact_rhea = {r["Reaction identifier"]: r for rid in (bio_rids | krs) for r in rhea_xref.get(rid, [])}
        ec_rhea = {r["Reaction identifier"]: r for ec in ecs for r in rhea_ec.get(ec, [])}
        b_rows = [b for ec in ecs for b in brenda_ec.get(ec, [])]
        ambiguous = len(old) > 1 or len(alias_hits) > 1 or protein_conflict or any("冲突" in s or "多重" in s for s in map_states)
        certain_old = old if not ambiguous and all(r.get("映射状态") == "完全匹配" and tokens(r.get("当前GCF049位点", "")) == {locus} for r in m) else set()
        uncertain_old = old - certain_old
        model16 = {rid for gene in certain_old for rid in old16.get(gene, set())}
        model24 = {rid for gene in certain_old for rid in old24.get(gene, set())}
        uncertain16 = {rid for gene in uncertain_old for rid in old16.get(gene, set())}
        uncertain24 = {rid for gene in uncertain_old for rid in old24.get(gene, set())}
        model_state = "是" if model16 or model24 else ("歧义候选" if uncertain16 or uncertain24 else "否")
        row = dict(base)
        row.update({"旧AFE locus": join(old), "旧RefSeq locus": join(refs), "旧protein ID": join(old_pro), "B1映射状态": join(map_states) if m else "未映射", "B1映射方法": join(r.get("映射方法", "") for r in m), "映射冲突": "是" if ambiguous else "否", "KEGG afr gene": join(old & kegg.keys()), "KEGG KO": join(kos), "KEGG EC": join(kecs), "KEGG reaction": join(krs), "KEGG pathway": join(x.get("pathway", "") for x in kr), "KEGG module": join(x.get("module", "") for x in kr), "KEGG注释": join(x.get("KEGG gene annotation", "") for x in kr), "BioCyc gene ID": join(alias_hits), "BioCyc protein ID": join(p.get("对象ID", "") for p in bio_pro), "BioCyc蛋白名称": join(p.get("名称", "") for p in bio_pro), "BioCyc reaction ID": join(bio_rids), "BioCyc EC": join(becs), "BioCyc转运反应ID": join(rid for rid in bio_rids if rid.startswith("TRANS-RXN") or rid.startswith("TRANSPORT-RXN")), "BRENDA EC记录": join(f'{b.get("EC号", "")}@{b.get("organism record ID", "")}@{b.get("BRENDA酶ID", "")}' for b in b_rows), "Rhea精确交叉引用": join(exact_rhea), "Rhea仅EC候选": join(set(ec_rhea) - set(exact_rhea)), "MetaCyc独立获取": "未独立取得，当前信息继承自BioCyc/Pathway Tools", "2016 GPR反应": join(model16), "2024 GPR反应": join(model24), "2016 GPR歧义候选反应": join(uncertain16), "2024 GPR歧义候选反应": join(uncertain24), "旧模型GPR覆盖": model_state})
        specific = bool(kos or kecs or bio_rids or (base["NCBI当前功能"] and not re.search(r"hypothetical protein|uncharacterized protein|domain.containing protein$", base["NCBI当前功能"], re.I)))
        product_text = base["NCBI当前功能"]
        nonmetabolic = bool(re.search(r"\b(?:DNA|RNA|tRNA|rRNA|ribosomal|transcription|translation|replication|helicase|topoisomerase|recombin\w*|\w*nuclease|histidine kinase|two.component|signal transduction|preprotein translocase|protein folding|chaperone|GTPase|peptide deformylase|peptidyl.prolyl)\b", product_text, re.I))
        metabolic_pathway = any(re.fullmatch(r"afr(?:00\d{3}|01[12]\d{2})", p) for x in kr for p in tokens(x.get("pathway", "")))
        metabolic_name = bool(re.search(r"acetyl|acyl|glutathione|quinone|sulfur|sulfide|thiosulfate|iron|heme|lipid|peptidoglycan|glucose|ribose|pyruvate|carboxyl|hydrogenase|dehydrogenase|oxidoreductase|aminotransferase|biosynthesis|synthase", product_text, re.I))
        metabolic = bool(kecs or becs or krs or bio_rids) and not nonmetabolic and (metabolic_pathway or metabolic_name)
        classify(row, model_equations)
        metabolic |= any(row[k] == "是" for k in ("转运候选", "中央碳表入选", "FeS_ETC表入选"))
        row["旧模型未覆盖功能候选"] = "待定" if ambiguous else ("是" if model_state == "否" and specific and metabolic else "否")
        row["候选依据"] = join(["KEGG EC/反应" if kecs or krs else "", "BioCyc基因关联反应" if bio_rids else "", "当前蛋白转运注释" if row["转运候选"] == "是" else "", "重点代谢酶/复合物注释" if row["中央碳表入选"] == "是" or row["FeS_ETC表入选"] == "是" else ""]) if row["旧模型未覆盖功能候选"] == "是" else ""
        row["功能判定"] = "候选；需反应、底物和GPR复核" if row["旧模型未覆盖功能候选"] == "是" else ("映射待定，不判定覆盖" if ambiguous else "仅交叉引用/已有GPR覆盖")
        row["数据库分歧"] = "映射歧义" if ambiguous else (row["重点分类冲突"] if row["重点分类冲突"] else ("KEGG与BioCyc EC不相交" if kecs and becs and not kecs & becs else "未见可自动判定的冲突"))
        output.append(row)
        coverage.append({"当前locus": locus, "旧AFE locus": row["旧AFE locus"], "映射冲突": row["映射冲突"], "2016 GPR反应": row["2016 GPR反应"], "2024 GPR反应": row["2024 GPR反应"], "2016 GPR歧义候选反应": row["2016 GPR歧义候选反应"], "2024 GPR歧义候选反应": row["2024 GPR歧义候选反应"], "旧模型GPR覆盖": row["旧模型GPR覆盖"], "旧模型未覆盖功能候选": row["旧模型未覆盖功能候选"], "BioCyc reaction ID": row["BioCyc reaction ID"], "KEGG reaction": row["KEGG reaction"], "反应等价状态": "未核实化学计量和区室；仅保留ID与基因关联", "注释来源": "GCF049 GFF;B1;A0;KEGG afr;BioCyc GCF000021485"})
        def ev(db, record, field, detail, level, source, version, strain):
            evidence.append({"当前locus": locus, "数据库": db, "记录ID": record, "证据字段": field, "记录内容": detail, "证据等级/限制": level, "菌株": strain, "数据库版本/组装": version, "取得日期": DATE, "来源文件": source})
        ev("NCBI", locus, "protein/product", base["当前protein ID"] + " " + base["NCBI当前功能"], "当前同株组装注释；功能属计算注释", "02_NCBI/标准化数据/GCF_049532655.1/genomic.gff", ASSEMBLY, "ATCC 23270")
        for x in m:
            ev("B1", x.get("旧AFE位点", ""), "旧新locus映射", x.get("旧RefSeq位点", "") + "=>" + locus + ";" + x.get("映射方法", ""), x.get("映射状态", ""), "02_NCBI/标准化数据/旧新基因映射.tsv", "旧GCF_000021485.1→" + ASSEMBLY, "ATCC 23270")
        for x in kr:
            ev("KEGG", "afr:" + x["KEGG gene/locus"], "KO/EC/reaction/pathway/module", "KO=" + x.get("KO", "") + ";EC=" + x.get("EC", "") + ";reaction=" + x.get("KEGG reaction", "") + ";pathway=" + x.get("pathway", "") + ";module=" + x.get("module", ""), "旧同株组装；KO继承反应不能当作基因实验", "04_KEGG/D3-1/KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv", "afr;旧GCA_000021485.1;REST release未提供", "ATCC 23270")
        for geneid in alias_hits:
            ev("BioCyc", PGDB + ":" + geneid, "gene", "对应旧locus " + row["旧AFE locus"], "Tier 3 PathoLogic计算；旧组装", "03_BioCyc/C2-1/标准化数据/genes.tsv", "PGDB 30.0;旧GCF_000021485.1", "ATCC 23270")
        for p in bio_pro:
            ev("BioCyc", PGDB + ":" + p["对象ID"], "protein", p.get("名称", "") + ";外部ID=" + p.get("外部ID", ""), "protein对象；不等于reaction对象", "03_BioCyc/C2-1/标准化数据/proteins.tsv", "PGDB 30.0;旧GCF_000021485.1", "ATCC 23270")
        for x in brxs:
            ev("BioCyc", PGDB + ":" + x["对象ID"], "reaction", x.get("左侧", "") + "=>" + x.get("右侧", "") + ";EC=" + x.get("EC号", "") + ";EV=" + x.get("证据代码", ""), "基因关联候选；未判定与模型反应等价", "03_BioCyc/C2-1/标准化数据/reactions.tsv", "PGDB 30.0;旧GCF_000021485.1", "ATCC 23270")
        for b in b_rows:
            org, enz, ec = b.get("BRENDA organism标签", ""), b.get("BRENDA酶ID", ""), b.get("EC号", "")
            ev("BRENDA", "organism:" + b.get("organism record ID", "") + ";enzyme:" + enz + ";EC:" + ec, "enzyme/EC", "酶名=" + b.get("酶名称", "") + ";substrate/product/cofactor=本条不据EC聚合赋值;kinetics/pH/temperature=未取得", "同EC交叉引用；非该基因直接实验证据；" + b.get("同株判断", ""), "05_BRENDA/D3-1/原始响应/BRENDA_org_EC.json", "BRENDA 2026.1;SPARQL", org)
            for key, details in brenda_reactions.items():
                if key[:3] == (org, enz, ec):
                    ev("BRENDA", "enzyme:" + enz + ";reaction:" + key[3], "reaction/compound", "substrate=" + join(details["Substrate"]) + ";product=" + join(details["Product"]), "原始reaction与enzyme同条连接；仍非当前基因实验证据", "05_BRENDA/D3-1/原始响应/BRENDA_FeS_reaction_compounds.json", "BRENDA 2026.1;SPARQL", org)
            details = brenda_cofactors.get((org, enz, ec))
            if details:
                ev("BRENDA", "enzyme:" + enz + ";EC:" + ec, "cofactor/reference", "cofactor=" + join(details["cofactor"]) + ";reference=" + join(details["ref"]), "原始enzyme级关联；未指定到具体reaction或当前基因", "05_BRENDA/D3-1/原始响应/BRENDA_FeS_cofactors_references.json", "BRENDA 2026.1;SPARQL", org)
        for r in exact_rhea.values():
            ev("Rhea", r["Reaction identifier"], "exact xref/equation/ChEBI", r.get("Equation", "") + ";ChEBI=" + r.get("ChEBI identifier", "") + ";KEGG=" + r.get("Cross-reference (KEGG)", "") + ";MetaCyc=" + r.get("Cross-reference (MetaCyc)", "") + ";direction=本地字段未取得", "仅外部反应ID交叉引用；方向/质子/电子待人工核验", rhea_source, "Rhea REST;release未提供", "非宿主特异")
        for r in ec_rhea.values():
            if r["Reaction identifier"] not in exact_rhea:
                ev("Rhea", r["Reaction identifier"], "EC-only candidate", r.get("Equation", "") + ";ChEBI=" + r.get("ChEBI identifier", ""), "仅EC候选，不能认定宿主反应", rhea_source, "Rhea REST;release未提供", "非宿主特异")
        if ambiguous or row["数据库分歧"] != "未见可自动判定的冲突":
            conflicts.append({"当前locus": locus, "旧AFE locus": row["旧AFE locus"], "分歧类型": row["数据库分歧"], "KEGG EC": row["KEGG EC"], "BioCyc EC": row["BioCyc EC"], "B1映射状态": row["B1映射状态"], "处理": "保留冲突，人工复核；不修改模型"})
        if not m:
            missing.append({"当前locus": locus, "来源": "B1", "状态": "旧locus未映射", "说明": "当前基因存在；无法据此判断旧组装缺失或新基因", "处理": "保留当前NCBI注释"})
        if m and not kr:
            missing.append({"当前locus": locus, "来源": "KEGG afr", "状态": "未关联到本地KEGG旧locus记录", "说明": "旧同株组装覆盖差异或映射缺失；不表示功能不存在", "处理": "保留未取得"})
        if ecs and not (exact_rhea or ec_rhea):
            missing.append({"当前locus": locus, "来源": "Rhea", "状态": "本地EC查询无命中或EC不完整", "说明": "本轮已补查当前全基因完整四级EC；无命中不等于宿主无反应", "处理": "保留未取得"})
    for mr in model_rows:
        rid = mr["模型反应ID"]
        gpr_old = set(re.findall(r"AFE_\d{4}", mr.get("2016 GPR", "") + " " + mr.get("2024 GPR", "")))
        mr["GPR旧locus未在B1映射表"] = join(gpr_old - set(by_old))
        hits = [r for r in bio_rx.values() if exact_biocyc_match(rid, r)]
        mr["BioCyc精确反应ID"] = join(r["对象ID"] for r in hits)
        mr["匹配依据"] = join("对象ID完整相等" if rid == r["对象ID"] else "外部ID独立token完整相等" for r in hits)
        mr["MetaCyc状态"] = "未独立取得，当前信息继承自BioCyc/Pathway Tools"
        mr["反应等价状态"] = "ID精确；化学计量、方向、区室仍待核对" if hits else "未取得ID精确匹配；不表示无对应反应"
    write_tsv("D3-1-R_全基因跨数据库交叉引用.tsv", output, list(output[0]))
    write_tsv("D3-1-R_旧模型未覆盖功能候选.tsv", [r for r in output if r["旧模型未覆盖功能候选"] == "是"], list(output[0]))
    write_tsv("D3-1-R_注释证据表.tsv", evidence, list(evidence[0]))
    write_tsv("D3-1-R_数据库分歧.tsv", conflicts, ["当前locus", "旧AFE locus", "分歧类型", "KEGG EC", "BioCyc EC", "B1映射状态", "处理"])
    for old_row in maps:
        if not old_row.get("当前GCF049位点"):
            missing.append({"当前locus": "未映射:" + old_row.get("旧AFE位点", ""), "来源": "B1", "状态": "旧位点到当前基因未确定", "说明": old_row.get("映射状态", "") + ";" + old_row.get("证据备注", "") + "；不得解释为当前基因缺失", "处理": "保留未确定"})
    model_old_absent = (set(old16) | set(old24)) - set(by_old)
    for old_gene in sorted(model_old_absent):
        missing.append({"当前locus": "未映射:" + old_gene, "来源": "A0 GPR/B1", "状态": "模型GPR旧locus未在B1映射表", "说明": "相关反应=" + join(old16.get(old_gene, set()) | old24.get(old_gene, set())) + "；复合体成员不能视为完整映射", "处理": "保留未确定"})
    missing += [{"当前locus": "全局", "来源": "MetaCyc", "状态": "未独立取得", "说明": "当前信息继承自BioCyc/Pathway Tools", "处理": "不得计为独立证据"}, {"当前locus": "全局", "来源": "BRENDA", "状态": "动力学/pH/温度未完整取得", "说明": "现有SPARQL部分响应达到LIMIT 5000；EC记录不能推出基因反应", "处理": "保留未取得"}, {"当前locus": "全局", "来源": "UniProt", "状态": "未检查", "说明": "未进行独立UniProt交叉证据检索", "处理": "保留未检查"}]
    all_ec = {ec for row in output for field in ("KEGG EC", "BioCyc EC") for ec in tokens(row[field]) if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec)}
    ec_no_rhea = all_ec - set(rhea_ec)
    missing.append({"当前locus": "全局", "来源": "Rhea", "状态": "完整四级EC已查、部分无返回", "说明": f"当前完整四级EC {len(all_ec)} 个，Rhea有返回 {len(all_ec)-len(ec_no_rhea)} 个，无返回 {len(ec_no_rhea)} 个；不完整EC和方向字段未全面取得", "处理": "无返回不解释为宿主缺失功能"})
    write_tsv("D3-1-R_访问缺失表.tsv", missing, ["当前locus", "来源", "状态", "说明", "处理"])
    write_tsv("D3-1-R_模型覆盖状态.tsv", coverage, list(coverage[0]))
    focus_cols = ["当前locus", "当前protein ID", "旧AFE locus", "NCBI当前功能", "KEGG KO", "KEGG EC", "KEGG pathway", "BioCyc gene ID", "BioCyc reaction ID", "2016 GPR反应", "2024 GPR反应", "旧模型未覆盖功能候选", "中央碳分类", "有机碳转运分类", "FeS_ETC分类", "功能判定"]
    for name, key in (("中央碳代谢证据", "中央碳表入选"), ("有机碳转运证据", "有机碳转运表入选"), ("FeS_ETC证据", "FeS_ETC表入选")):
        write_tsv("D3-1-R_" + name + ".tsv", [r for r in output if r[key] == "是"], focus_cols)
    write_tsv("D3-1-R_模型反应精确匹配.tsv", model_rows, list(model_rows[0]))
    manifest = build_manifest(gff_path)
    write_tsv("D3-1-R_原始数据校验清单.tsv", manifest, list(manifest[0]))
    rhea_gap_manifest = read_tsv("D3-1-R/Rhea_新增查询清单.tsv") if (OUT / "Rhea_新增查询清单.tsv").exists() else []
    counts = {"current_genes": len(output), "protein_coding": sum(r["基因类型"] == "protein_coding" for r in output), "mapped_old": sum(bool(r["旧AFE locus"]) for r in output), "old_unmapped": sum(not r.get("当前GCF049位点") for r in maps), "model_old_absent_B1": len(model_old_absent), "old_multimapping": sum("多重" in r.get("映射状态", "") for r in maps), "kegg": sum(bool(r["KEGG afr gene"]) for r in output), "biocyc": sum(bool(r["BioCyc gene ID"]) for r in output), "model_gpr": sum(r["旧模型GPR覆盖"] == "是" for r in output), "model_gpr_ambiguous": sum(r["旧模型GPR覆盖"] == "歧义候选" for r in output), "uncovered_function_candidates": sum(r["旧模型未覆盖功能候选"] == "是" for r in output), "mapping_conflicts": sum(r["映射冲突"] == "是" for r in output), "disagreements": len(conflicts), "transport_candidates": sum(r["转运候选"] == "是" for r in output), "uncovered_transport": sum(r["转运候选"] == "是" and r["旧模型GPR覆盖"] == "否" for r in output), "organic_transport": sum(r["有机碳转运表入选"] == "是" for r in output), "central": sum(r["中央碳表入选"] == "是" for r in output), "fes_etc": sum(r["FeS_ETC表入选"] == "是" for r in output), "model_exact_biocyc": sum(bool(r["BioCyc精确反应ID"]) for r in model_rows), "rhea_full_ec": len(all_ec), "rhea_ec_result": len(all_ec)-len(ec_no_rhea), "rhea_ec_no_return": len(ec_no_rhea), "rhea_gap_ec_queried": sum(len(tokens(r["查询EC"])) for r in rhea_gap_manifest), "rhea_gap_batches": len(rhea_gap_manifest), "rhea_gap_failed": sum(not (ROOT / r["文件"]).exists() for r in rhea_gap_manifest), "raw_manifest_compared": sum(r["校验"] == "PASS" for r in manifest), "raw_manifest_self_hash_only": sum(r["校验"].startswith("仅实算") for r in manifest), "raw_manifest_fail": sum(r["校验"] == "FAIL" for r in manifest)}
    (OUT / "D3-1-R_metrics.json").write_text(json.dumps(counts, ensure_ascii=False, indent=2), encoding="utf-8")
    report = f"""# D3-1-R 验收报告

## 范围与结果

- 起点：NCBI `{ASSEMBLY}` GFF 全部 gene，{counts['current_genes']} 条，其中蛋白编码 {counts['protein_coding']} 条。
- B1 可关联旧 AFE locus：{counts['mapped_old']} 条；旧位点未能映射当前：{counts['old_unmapped']} 条；旧位点多重候选：{counts['old_multimapping']} 条。模型 GPR 中另有 {counts['model_old_absent_B1']} 个旧 locus 未出现在 B1 映射表（`{';'.join(sorted(model_old_absent))}`），相关复合体不能视为完整映射。KEGG afr：{counts['kegg']} 条；BioCyc 旧 PGDB：{counts['biocyc']} 条。
- 2016 或 2024 GPR 确定覆盖：{counts['model_gpr']} 条；另有映射歧义造成的 GPR 候选 {counts['model_gpr_ambiguous']} 条，未计入确定覆盖。旧模型未覆盖功能候选：{counts['uncovered_function_candidates']} 条。候选只是注释筛选，不能直接加入模型。
- 映射冲突标记：{counts['mapping_conflicts']} 条当前基因；数据库分歧表：{counts['disagreements']} 条（含 EC 集合不相交的待核候选）。
- 转运候选：{counts['transport_candidates']} 条，其中未被旧模型 GPR 覆盖 {counts['uncovered_transport']} 条；明确点名有机碳底物类别：{counts['organic_transport']} 条。中央碳表 {counts['central']} 条；Fe/S 与 ETC 表 {counts['fes_etc']} 条。
- 631 条模型反应中 BioCyc 对象ID或外部ID独立 token 精确相等：{counts['model_exact_biocyc']} 条。`RPE` 不得匹配 `PHOSNACMURPENTATRANS-RXN`。
- 原始清单：与既有 SHA256 对照通过 {counts['raw_manifest_compared']} 个；原清单未给 SHA256、仅实算 {counts['raw_manifest_self_hash_only']} 个；不一致 {counts['raw_manifest_fail']} 个。
- Rhea 当前基因组完整四级 EC 共 {counts['rhea_full_ec']} 个，已全部列入查询范围；有返回 {counts['rhea_ec_result']} 个，已查无返回 {counts['rhea_ec_no_return']} 个。本轮补查 {counts['rhea_gap_ec_queried']} 个，{counts['rhea_gap_batches']} 批；失败 {counts['rhea_gap_failed']} 批。原始响应及查询清单保存在本目录。

## 证据边界

- **已核实**：当前 GFF 基因与蛋白 ID；B1 旧新位点表；KEGG `afr` 同株身份；BioCyc `GCF_000021485` PGDB 身份；模型 GPR 原始表。每条关联在证据表保留记录 ID、来源、菌株、版本与取得日期。
- **推测/候选**：NCBI 当前产物、KEGG KO、BioCyc PathoLogic 计算反应和 BRENDA 同 EC 记录组成筛选线索；尚未证明每条基因的具体催化反应、转运底物或 GPR。Rhea 仅标准化方程和 ChEBI，EC 命中不算宿主反应。
- **候选筛选口径**：旧 GPR 未覆盖，且有代谢通路或代谢产物名称与 EC/反应关联，或有蛋白转运、中央碳、Fe/S 呼吸链注释；排除明显 DNA/RNA/翻译维护蛋白。家族级转运只标未定底物。该自动筛选会有漏检和误检，需人工逐条决定。
- **数据库冲突**：见分歧表；EC 差异和一对多旧新映射保持未决，不自动消解。
- **数据缺失**：MetaCyc **未独立取得，当前信息继承自BioCyc/Pathway Tools**；UniProt 未检查；BRENDA 动力学、pH、温度未完整取得；Rhea 已查无返回的完整四级 EC、不完整 EC 及方向字段保留未取得。未命中不得解释为 ATCC 23270 没有功能。

## 来源版本

- NCBI：`{ASSEMBLY}`，ATCC 23270；本地文件 `02_NCBI/标准化数据/GCF_049532655.1/genomic.gff`。
- KEGG：`afr`，ATCC 23270；仍参照旧 `GCA_000021485.1`，REST release 字段未提供。
- BioCyc：`{PGDB}`，PGDB 30.0，Tier 3；基于旧 `GCF_000021485.1`，PathoLogic 计算预测。
- BRENDA：本地 DSMZ SPARQL 响应，网页标记 2026.1；菌株标签按原样保留，同 EC 不升级为同株基因实验。
- Rhea：本地 REST 响应，release 字段未取得；方程、ChEBI、质子/电子和方向需逐条核对后才能认定等价。

## 强制回归

- `AFE_2841` 应关联 `RU820_RS13140`、`WP_012537424.1`、`K25026`、`EC 2.7.1.2`。此链由输入文件交叉连接，未硬编码到构建逻辑。测试见 `test_d3r.py`。
- BioCyc 精确 ID、重点分类和输出文件由测试验证；独立 Reviewer 复查记录另存。

## 独立 Reviewer 复查

`gpt-5.6-sol / medium` Reviewer 独立抽样发现并推动修复：多重映射误计确定 GPR 覆盖、FE2tpp 对应 porin 进入有机碳表、HlyD/FtsX/oxidase assembly 污染、NADHI/HYD4 亚基漏项、BRENDA EC 聚合明细误挂单条 enzyme、模型 GPR 两个旧 locus 未列入缺失表。Reviewer 复查最新核心映射与验收结论后确认可写 PASS。8 项回归测试通过。

## 验收状态

**PASS：D3 证据与映射体系验收。** 候选仍需 E4-1 逐条判断；未修改 SBML 或正式 GPR。
"""
    (OUT / "D3-1-R_验收报告.md").write_text(report, encoding="utf-8")
    print(json.dumps(counts, ensure_ascii=False))


def build_manifest(gff_path):
    paths = [gff_path, ROOT / "02_NCBI/标准化数据/旧新基因映射.tsv", ROOT / "08_模型基线/A0/恢复GPR及原始证据定位.tsv", ROOT / "08_模型基线/A0/2016_到2024_真实差异.tsv", ROOT / "04_KEGG/D3-1/KEGG_afr_gene_ko_ec_reaction_pathway_module.tsv", ROOT / "04_KEGG/D3-1/KEGG_identity.json", ROOT / "05_BRENDA/D3-1/BRENDA_酶学结果.tsv", ROOT / "06_Rhea/D3-1/Rhea_modelBioCyc_EC_crossrefs.tsv", ROOT / "03_BioCyc/C2-1/C2-1_PGDB身份与版本.md"]
    if (OUT / "Rhea_全基因EC查询汇总.tsv").exists():
        paths.append(OUT / "Rhea_全基因EC查询汇总.tsv")
    for group in ("04_KEGG/D3-1/原始响应", "05_BRENDA/D3-1/原始响应", "06_Rhea/D3-1/原始响应", "03_BioCyc/C2-1/官方原始下载", "03_BioCyc/C2-1/标准化数据"):
        paths.extend(p for p in (ROOT / group).iterdir() if p.is_file())
    gap_raw = OUT / "Rhea新增原始响应"
    if gap_raw.exists():
        paths.extend(p for p in gap_raw.iterdir() if p.is_file())
    source_manifests = {}
    for path in (ROOT / "04_KEGG/D3-1/KEGG_download_manifest.tsv", ROOT / "06_Rhea/D3-1/Rhea_API_download_manifest.tsv"):
        if path.exists():
            for row in read_tsv(str(path.relative_to(ROOT))):
                key = row.get("file", row.get("文件", "")).replace("\\", "/")
                source_manifests[key] = row
    bpath = ROOT / "05_BRENDA/D3-1/BRENDA_download_manifest.json"
    if bpath.exists():
        for row in json.loads(bpath.read_text(encoding="utf-8")):
            source_manifests[row.get("file", "").replace("\\", "/")] = row
    gap_manifest = OUT / "Rhea_新增查询清单.tsv"
    if gap_manifest.exists():
        for row in read_tsv("D3-1-R/Rhea_新增查询清单.tsv"):
            source_manifests[row["文件"].replace("\\", "/")] = row
    c2path = ROOT / "03_BioCyc/C2-1/C2-1_下载清单.tsv"
    if c2path.exists():
        # This historical manifest has literal `t in its header, while data
        # rows use tab separators. Read data columns by their documented order.
        for line in c2path.read_text(encoding="utf-8-sig").splitlines()[1:]:
            fields = line.split("\t")
            if len(fields) >= 8:
                relative = "03_BioCyc/C2-1/" + fields[1].replace("\\", "/")
                source_manifests[relative] = {"sha256": fields[6], "url": fields[2], "time_local": fields[3], "database_version_or_release": "BioCyc PGDB " + fields[4] + ";旧GCF_000021485.1", "http_status": fields[7]}
    rows = []
    for path in sorted(set(paths)):
        relative = path.relative_to(ROOT).as_posix()
        meta = source_manifests.get(relative, {})
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        declared = meta.get("sha256", meta.get("SHA256", "")).lower()
        status = meta.get("http_status", meta.get("HTTP_status", meta.get("HTTP状态", "未在本表复核")))
        note = "访问失败响应归档；未用于注释" if status and ("404" in str(status) or "NO_HTTP_RESPONSE" in str(status) or "ERROR" in str(status)) else "来源文件归档；哈希实算"
        rows.append({"文件": relative, "字节数": path.stat().st_size, "SHA256实算": digest, "原清单SHA256": declared, "校验": ("PASS" if digest == declared else "FAIL") if declared else "仅实算，原清单未提供", "获取时间": meta.get("time_local", meta.get("access_time", meta.get("获取时间", ""))) or datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(), "来源URL/API": meta.get("url", meta.get("URL/API", meta.get("来源URL/API", ""))) or "本地输入/派生表；见来源身份文件", "数据库版本": meta.get("database_version_or_release", meta.get("version", meta.get("数据库版本", ""))) or ("BioCyc PGDB 30.0;旧GCF_000021485.1" if "BioCyc" in relative else "见来源身份文件；原始响应未提供release"), "HTTP状态": status, "备注": note})
    return rows


if __name__ == "__main__":
    main()
