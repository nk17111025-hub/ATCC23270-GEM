"""Assemble G51 working tables without treating copied annotations as experiments."""
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = json.loads((HERE / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
gene_hits = json.loads((HERE / "europe_pmc_gene_search.json").read_text(encoding="utf-8"))
current_hits = json.loads((HERE / "europe_pmc_current_id_search.json").read_text(encoding="utf-8"))
function_hits = json.loads((HERE / "europe_pmc_function_search.json").read_text(encoding="utf-8"))
curated = json.loads((HERE / "curated_findings.json").read_text(encoding="utf-8"))
doi_sources = json.loads((HERE / "G51_doi_source_ids.json").read_text(encoding="utf-8"))
balance = {r["reaction_id"]: r for r in json.loads((HERE / "2016_G51_balance.json").read_text(encoding="utf-8"))}
reconciled = {r["locus"]: r for r in json.loads((HERE / "G51_sequence_reconciled_2016.json").read_text(encoding="utf-8"))}

def semi(value):
    return "; ".join(str(x) for x in value if x)


def uniprot_name(record):
    value = record.get("name", [])
    if isinstance(value, str):
        return value
    for description in value:
        if not isinstance(description, dict):
            continue
        if description.get("fullName", {}).get("value"):
            return description["fullName"]["value"]
        for key in ("recommendedName", "submissionNames", "alternativeNames"):
            value = description.get(key)
            if isinstance(value, dict):
                return value.get("fullName", {}).get("value", "")
    return ""


def protein_database_detail(row):
    out = []
    for record in row["uniprot"]:
        name = uniprot_name(record)
        ecs = semi(record.get("ec", []))
        kind = "审校" if "reviewed (Swiss" in record["entry_type"] else "未审校"
        out.append(f"UniProt {record['accession']} {kind}: {name or '未给具体蛋白名'}; EC={ecs or '无'}; protein existence={record['protein_existence']}")
    if not out:
        out.append("未取得同 WP 的 UniProt 记录")
    domains = sorted({f"{x['db']}:{x['id']}" for record in row["uniprot"] for x in record["crossrefs"] if x["db"] in {"InterPro", "Pfam"}})
    return semi(out) + "; domains=" + (semi(domains) if domains else "未取得")


def metabolic_database_detail(row):
    cross = row["d3_crossref"]
    parts = [f"KEGG KO={cross.get('KEGG KO') or '无'}; EC={cross.get('KEGG EC') or '无'}; reaction={cross.get('KEGG reaction') or '无'}",
             f"BioCyc EC={cross.get('BioCyc EC') or '无'}; reaction={cross.get('BioCyc reaction ID') or '无'}",
             f"Rhea={semi(x['Reaction identifier'] + ': ' + x['Equation'] for x in row['rhea_definitions']) or '无'}",
             f"BRENDA={cross.get('BRENDA EC记录') or '无'}"]
    return "; ".join(parts)


def ec_disagreement(row):
    cross = row["d3_crossref"]
    kegg = set(filter(None, cross.get("KEGG EC", "").split(";")))
    biocyc = set(filter(None, cross.get("BioCyc EC", "").split(";")))
    uniprot = {ec for record in row["uniprot"] for ec in record.get("ec", [])}
    if ((kegg and biocyc and not kegg & biocyc) or
        (uniprot and kegg and not uniprot & kegg) or
        (uniprot and biocyc and not uniprot & biocyc)):
        return f"KEGG={semi(sorted(kegg)) or '无'}; BioCyc={semi(sorted(biocyc)) or '无'}; UniProt={semi(sorted(uniprot)) or '无'}；需核底物和EC沿革"
    return "未见跨来源 EC 完全无交集（不代表精确反应已验证）"


index = []
alpha = []
for row in rows:
    locus = row["locus"]
    old = row["old_afe"]
    uni = row["uniprot"]
    cross = row["d3_crossref"]
    found = curated.get(locus, {})
    recovered = reconciled.get(locus, {})
    model_rows = row["model_2016"] + recovered.get("2016_reactions", [])
    model_ids = list(dict.fromkeys(r["Reaction ID"] for r in model_rows))
    old_display = semi(old) if old else ((recovered.get("old_afe_candidate") or "未确认") + "（序列追踪；见状态）" if recovered.get("old_afe_candidate") else "未确认")
    rhea_ids = [r["Reaction identifier"] for r in row["rhea_definitions"]]
    interpro = sorted({r["id"] for u in uni for r in u["crossrefs"] if r["db"] == "InterPro"})
    pfam = sorted({r["id"] for u in uni for r in u["crossrefs"] if r["db"] == "Pfam"})
    uni_accessions = [u["accession"] for u in uni]
    reviewed = any("reviewed (Swiss" in u["entry_type"] for u in uni)
    gene_search = gene_hits[locus]
    function_search = function_hits[locus]
    candidate_papers = list(dict.fromkeys((r.get("doi") or "").lower() for h in gene_search for r in h["results"] if r.get("doi")))
    functional_papers = list(dict.fromkeys((r.get("doi") or "").lower() for r in function_search["results"] if r.get("doi")))
    functional_titles = list(dict.fromkeys(r.get("title", "") for r in function_search["results"] if r.get("title")))
    hits_total = sum(h.get("hit_count", 0) for h in gene_search)
    current_count = current_hits[locus].get("hit_count", 0) + current_hits[row["wp"]].get("hit_count", 0)
    literature_check = (f"旧 AFE/未映射时当前编号精确检索命中 {hits_total} 条；当前 RU820/WP 精确检索命中 {current_count} 条；功能词检索命中 {function_search.get('hit_count', 0)} 条。"
                        f"功能词前列标题={semi(functional_titles[:3]) or '无'}。命中须核论文中的序列、菌株及实验；零命中也不能证明不存在论文。")
    direct = found.get("direct_evidence", literature_check + "尚无可指向此 WP 精确反应的已核定原始实验。")
    if found:
        boundary = found["boundary"]
    elif model_ids:
        boundary = f"2016 已关联 {semi(model_ids)}；原表反应和旧 GPR 已定位，当前 WP 的精确底物/区室/方向仍缺直接实验。"
    elif rhea_ids or cross.get("KEGG reaction") or cross.get("BioCyc reaction ID"):
        boundary = "代谢数据库给出反应候选，但这些交叉引用未独立证明本 WP 的底物、区室、方向或体内通量。"
    elif "hypothetical" in row["product"].lower():
        boundary = "当前官方注释为 hypothetical；未取得足以确定底物和产物的本蛋白证据，不能定义精确代谢反应。"
    else:
        boundary = "当前蛋白家族/名称可追溯，但没有可核定本 WP 的精确代谢反应；不以基因存在推断通量。"
    if found:
        action = found["action"]
    elif model_ids:
        action = f"核对并更新旧反应 {semi(model_ids)} 的当前编号；以 Alpha 逐反应评分为准，实验缺口保留"
    elif rhea_ids or cross.get("KEGG reaction") or cross.get("BioCyc reaction ID"):
        action = "记录数据库反应候选；当前证据尚不足以确认新增精确式或 GPR，暂不写入模型"
    else:
        action = "证据不足暂缓；无已确认的模型新增或修改"
    reaction_state = "已确认旧模型式；基因催化待复核" if model_ids else ("数据库给出反应定义；基因链接未独立验证" if rhea_ids else "未找到可确认的精确反应")
    sources = ["G51-S001"]
    if row["map_rows"]:
        sources.append("G51-S002")
    if model_rows or row["model_2016_same_ec"]:
        sources.append("G51-S003")
    if recovered:
        sources.append("G51-S012")
    if locus in {"RU820_RS00280", "RU820_RS00315", "RU820_RS01200"}:
        sources.append("G51-S014")
    if locus == "RU820_RS01485":
        sources.append("G51-S015")
    if uni:
        sources.append("G51-S004")
    if interpro or pfam:
        sources.append("G51-S005")
    if cross.get("KEGG KO") or cross.get("KEGG EC") or cross.get("KEGG reaction"):
        sources.append("G51-S006")
    if cross.get("BioCyc reaction ID") or cross.get("BioCyc EC"):
        sources.append("G51-S007")
    if rhea_ids:
        sources.append("G51-S008")
    if cross.get("BRENDA EC记录"):
        sources.append("G51-S009")
    sources.append("G51-S010")
    sources.append("G51-S013")
    if re.search(r"transport|permease|efflux|porin|channel|sulp|aquaporin", row["product"], re.I):
        sources.append("G51-S011")
    if found.get("evidence_doi"):
        sources.extend(doi_sources[d.strip()] for d in found["evidence_doi"].split(";") if d.strip())
    status = "原文已核并裁决；未决点见本行" if found else "逐基因检索完成；同编号原始实验未见；精确反应证据不足"
    idx = {
        "总序号": row["ordinal"], "当前locus": locus, "当前WP": row["wp"],
        "坐标": f"{row['chromosome']}:{row['start']}-{row['end']}({row['strand']})",
        "当前product": row["product"], "CDS同版一致": row["cds_match_gene"],
        "旧AFE": old_display, "AFE映射依据": semi(sorted({m["映射方法"] + "/" + m["映射状态"] for m in row["map_rows"]})) + ("; " + recovered["status"] if recovered else ""),
        "UniProt": semi(uni_accessions), "UniProt审校": "有" if reviewed else "无/未取得",
        "InterPro": semi(interpro), "Pfam": semi(pfam),
        "KEGG KO": cross.get("KEGG KO", ""), "KEGG EC": cross.get("KEGG EC", ""), "KEGG reaction": cross.get("KEGG reaction", ""),
        "BioCyc reaction": cross.get("BioCyc reaction ID", ""), "BioCyc EC": cross.get("BioCyc EC", ""),
        "Rhea反应定义": semi(rhea_ids), "BRENDA记录": cross.get("BRENDA EC记录", ""),
        "2016直接关联反应": semi(model_ids), "2016同EC比较对象": semi([r["Reaction ID"] for r in row["model_2016_same_ec"] if r["Reaction ID"] not in model_ids]),
        "精确反应状态": reaction_state, "原始论文DOI": found.get("evidence_doi", ""),
        "蛋白功能数据库逐项": protein_database_detail(row), "代谢数据库逐项": metabolic_database_detail(row),
        "EC来源分歧": ec_disagreement(row),
        "原始论文检索判断": literature_check,
        "原始论文直接事实": direct, "证据边界": boundary, "本轮动作": action,
        "反应评分": "逐反应见 Alpha 行；未核定反应不评分", "审查状态": status,
        "文献ID查询命中": sum(h.get("hit_count", 0) for h in gene_search), "文献功能词查询命中": function_search.get("hit_count", 0),
        "当前RU/WP查询命中": current_count,
        "文献候选DOI": semi(list(dict.fromkeys(candidate_papers + functional_papers))[:20]),
        "文献检索链接": semi([h["url"] for h in gene_search]) + ("; " + function_search["url"] if function_search.get("url") else ""),
        "来源ID": semi(sources), "关联行ID": f"G51-{row['ordinal']:04d}",
        "未决点": boundary if not found else found["boundary"],
    }
    index.append(idx)
    alpha.append({
        "Reaction ID": "", "Reaction Name": "", "Reaction Formula": "", "Confidence Level": "", "EC Number": "", "PMID": "", "Subsystem": "",
        "Gene-Reaction Association": "", "Gene-Protein-Reaction Association": "", "Protein-Reaction-Association": "",
        "Fe2 lb": "", "Fe2 ub": "", "tetrathionate lb": "", "tetrathionate ub": "", "sulfur lb": "", "sulfur ub": "",
        "Record class": "2026 gene review (working)", "Candidate ID": f"G51-{row['ordinal']:04d}", "Candidate group": "G51", "2016 linked Reaction ID": semi(model_ids),
        "2026 reaction object": reaction_state + "; Rhea=" + semi(rhea_ids), "Legacy AFE locus": old_display, "Current locus": locus,
        "Primary evidence DOI": found.get("evidence_doi", ""), "Evidence type": "RefSeq; UniProt/InterPro/Pfam; KEGG/BioCyc/Rhea/BRENDA; literature index; 2016 supplement",
        "Evidence boundary": boundary, "2026 action": action, "2026 reaction score": "不评分（逐反应待核）", "Source ID": semi(sources), "Score scope": "逐基因索引不赋反应分数",
    })
    for j, model in enumerate(model_rows, 1):
        score = "待重评"
        scope = "2016 精确反应；2026 分数待原始证据与式子核对"
        if model["Reaction ID"] == "4THASE1":
            score = 3
            scope = "同株 TetH 纯化、遗传及产物类别支持；精确各产物 1:1:1:1 系数未直接测定"
        elif model["Reaction ID"] == "4THASE2":
            score = 1
            scope = "旧模型备选精确式；引文中的元素硫产物和 thiosulfate 比例与本式不符"
        elif model["Reaction ID"] == "TSQOC":
            score = 1
            scope = "旧 q8 依赖精确式：当前无基因专属 q8 催化证据；旧建模假设"
        elif model["Reaction ID"] in {"UAAD", "KDOCT", "KDOCT3", "LDAS", "LIPS"}:
            score = 2
            scope = "旧引文 PMID 15044494 实测 GnnA/GnnB，不是本行酶；当前仅注释/同源支持，本式特设底物待核"
        elif model["Reaction ID"] == "BIOC1":
            score = 1
            scope = "旧引文 PMID 238956 仅为脂多糖分析；本式载体代谢物缺化学定义且 GPR 需拆审"
        elif model["Reaction ID"] == "ACCOAC":
            score = 2
            scope = "当前亚基注释支持总体羧化；旧四亚基 OR GPR 不符合复合酶关系"
        elif model["Reaction ID"] == "GHMT2":
            score = 2
            scope = "UniProt EC 2.1.2.1 与 Rhea 15481 支持 SHMT 式；同株直接酶学未见"
        elif model["Reaction ID"] == "GHMT3":
            score = 1
            scope = "旧式为 NAD 依赖甘氨酸裂解系统整体反应；当前 SHMT 单基因 GPR 错配"
        elif model["Reaction ID"] == "AIRC":
            score = 1
            scope = "旧 PurE/PurK 合并式省略 ATP 与 N5-CAIR；作为精确反应只可视为模型假设"
        elif model["Reaction ID"] == "SPMS":
            score = 2
            scope = "亚精胺合成酶同源支持反应；旧 AND 含上游 SAM 脱羧酶，GPR 需修正"
        elif model["Reaction ID"] == "OPHHX":
            score = 1
            scope = "旧羟化精确式：UbiB 羟化 GPR 不受当前证据支持"
        elif model["Reaction ID"] == "HCO3E":
            score = 2
            scope = "旧水合精确式：功能注释支持；同株全细胞 CO2 数据不证明单酶方向"
        elif model["Reaction ID"] in {"CU2tpp", "ZN2tpp"}:
            score = 1
            scope = "AFE_0105 的 Cu/Zn 精确跨膜反应未有本蛋白实验证据；2008 原文仅推断 Mn/Fe"
        elif model["Reaction ID"] in {"MN2tpp", "MNTH"}:
            score = 2
            scope = "AFE_0105 的 Mn/Fe 转运仅基因组同源推断；离子耦联与转运方向未实测"
        elif model["Reaction ID"] == "HSERTA":
            score = 1
            scope = "旧乙酰 CoA 精确式与同 WP 的 UniProt 琥珀酰供体注释冲突；本蛋白供体未实测"
        elif model["Reaction ID"] == "PGLT":
            score = 1
            scope = "旧聚合反应由脂质 II 聚合酶承担；当前 MurG 证据只支持脂质 I 到 II 单体生成，GPR 错配"
        elif model["Reaction ID"] == "UCMAT":
            score = 2
            scope = "MurG EC 2.4.1.227 与 KEGG/BioCyc/Rhea 同向支持脂质 II 生成；同株直接酶学未见"
        elif model["Reaction ID"] == "UAMAGS":
            score = 2
            scope = "当前 MurD 与旧 AFE_0208 蛋白除旧 N 端2 aa 外其余446 aa完全一致；2016 式与酶名相符，无同株直接酶学"
        elif model["Reaction ID"] == "S7PI":
            score = 1
            scope = "当前仅覆盖旧 AFE_0141 C 端约46%，现行产物注释为 heptose-1-phosphate adenylyltransferase；旧异构酶 GPR 错配"
        elif model["Reaction ID"] == "BIOC2":
            score = 1
            scope = "旧式消耗 2 个 bcCP 却只生成 1 个 bcCPco2；自定义载体无分子式，精确式待重建"
        elif model["Reaction ID"] == "PGM1":
            score = 1
            scope = "当前仅组氨酸磷酸酶/磷酸甘油酸变位酶家族域；本 WP 的 3PG/2PG 底物尚无专属证据"
        elif model["Reaction ID"] == "SO42tpp":
            score = 1
            scope = "SulP 家族提示阴离子转运；当前 WP 的硫酸盐底物与外向内净通量未被实验证实"
        elif model["Reaction ID"] in {"CPPPGO2", "FCLT", "UAPGR"}:
            score = 2
            scope = "同 WP 当前酶名/KEGG/BioCyc 支持相应化学转化；旧 EC 编号与现行 EC 不同，需更新 EC；同株直接实验未见"
        elif model["EC Number"] and model["EC Number"] in (row["d3_crossref"].get("KEGG EC", "") + ";" + row["d3_crossref"].get("BioCyc EC", "")).split(";"):
            score = 2
            scope = "当前蛋白注释及 KEGG/BioCyc EC 与旧式相符；同株该 WP 的精确底物酶学未确认，仅注释级证据"
        else:
            score = 1
            scope = "旧模型精确式缺同 WP 专属反应实证，当前数据库亦未给出同 EC 精确支持；暂列模型假设"
        alpha.append({
            "Reaction ID": model["Reaction ID"], "Reaction Name": model["Reaction Name"], "Reaction Formula": model["Reaction Formula "],
            "Confidence Level": model["Confidence Level"], "EC Number": model["EC Number"], "PMID": model["PMID"], "Subsystem": model["Subsystem"],
            "Gene-Reaction Association": model["Gene-Reaction Association"], "Gene-Protein-Reaction Association": model["Gene-Protein-Reaction Association"], "Protein-Reaction-Association": model["Protein-Reaction-Association"],
            "Fe2 lb": model["lb"], "Fe2 ub": model["ub"], "tetrathionate lb": model["lb.1"], "tetrathionate ub": model["ub.1"], "sulfur lb": model["lb.2"], "sulfur ub": model["ub.2"],
            "Record class": "2016 linked reaction review (working)", "Candidate ID": f"G51-{row['ordinal']:04d}-R{j:02d}", "Candidate group": "G51", "2016 linked Reaction ID": model["Reaction ID"],
            "2026 reaction object": model["Reaction Formula "], "Legacy AFE locus": old_display, "Current locus": locus,
            "Primary evidence DOI": found.get("evidence_doi", ""), "Evidence type": "2016 source row; current identity/database; curated original article where listed",
            "Evidence boundary": boundary + "; 2016 balance=" + balance[model["Reaction ID"]]["status"], "2026 action": action,
            "2026 reaction score": score, "Source ID": semi(sources), "Score scope": scope,
        })

assert len(index) == 293 and len({r["当前locus"] for r in index}) == 293
assert len(alpha) == 293 + sum(len(r["model_2016"]) + len(reconciled.get(r["locus"], {}).get("2016_reactions", [])) for r in rows)
output = {"index": index, "alpha": alpha, "counts": {"genes": len(index), "alpha_rows": len(alpha), "linked_reaction_rows": len(alpha)-len(index), "curated_loci": len(curated), "review_status": dict(Counter(r["审查状态"] for r in index))}}
(HERE / "G51_working_tables.json").write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
print(output["counts"])
