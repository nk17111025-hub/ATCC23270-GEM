"""Assemble G51 working tables without treating copied annotations as experiments."""
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = json.loads((HERE / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
gene_hits = json.loads((HERE / "europe_pmc_gene_search.json").read_text(encoding="utf-8"))
function_hits = json.loads((HERE / "europe_pmc_function_search.json").read_text(encoding="utf-8"))
curated = json.loads((HERE / "curated_findings.json").read_text(encoding="utf-8"))
doi_sources = json.loads((HERE / "G51_doi_source_ids.json").read_text(encoding="utf-8"))
balance = {r["reaction_id"]: r for r in json.loads((HERE / "2016_G51_balance.json").read_text(encoding="utf-8"))}

def semi(value):
    return "; ".join(str(x) for x in value if x)


index = []
alpha = []
for row in rows:
    locus = row["locus"]
    old = row["old_afe"]
    uni = row["uniprot"]
    cross = row["d3_crossref"]
    found = curated.get(locus, {})
    model_ids = list(dict.fromkeys(r["Reaction ID"] for r in row["model_2016"]))
    rhea_ids = [r["Reaction identifier"] for r in row["rhea_definitions"]]
    interpro = sorted({r["id"] for u in uni for r in u["crossrefs"] if r["db"] == "InterPro"})
    pfam = sorted({r["id"] for u in uni for r in u["crossrefs"] if r["db"] == "Pfam"})
    uni_accessions = [u["accession"] for u in uni]
    reviewed = any("reviewed (Swiss" in u["entry_type"] for u in uni)
    gene_search = gene_hits[locus]
    function_search = function_hits[locus]
    candidate_papers = list(dict.fromkeys((r.get("doi") or "").lower() for h in gene_search for r in h["results"] if r.get("doi")))
    direct = found.get("direct_evidence", "本轮文献索引尚未确认该基因的同株精确反应实验。")
    boundary = found.get("boundary", "同版身份、家族/EC 和旧模型关联可核；底物、方向、区室及本基因直接催化证据仍需逐项裁决。")
    if found:
        action = found["action"]
    elif model_ids:
        action = "初核：旧反应已关联；仅更新证据或编号待逐反应裁决"
    else:
        action = "证据不足暂缓；无已确认的模型新增或修改"
    reaction_state = "已确认旧模型式；基因催化待复核" if model_ids else ("数据库给出反应定义；基因链接未独立验证" if rhea_ids else "未找到可确认的精确反应")
    sources = ["G51-S001", "G51-S002", "G51-S003", "G51-S004", "G51-S005", "G51-S006", "G51-S007", "G51-S008", "G51-S009", "G51-S010"]
    if found.get("evidence_doi"):
        sources.extend(doi_sources[d.strip()] for d in found["evidence_doi"].split(";") if d.strip())
    status = "初核；逐反应人工裁决未完成"
    idx = {
        "总序号": row["ordinal"], "当前locus": locus, "当前WP": row["wp"],
        "坐标": f"{row['chromosome']}:{row['start']}-{row['end']}({row['strand']})",
        "当前product": row["product"], "CDS同版一致": row["cds_match_gene"],
        "旧AFE": semi(old) if old else "未确认", "AFE映射依据": semi(sorted({m["映射方法"] + "/" + m["映射状态"] for m in row["map_rows"]})),
        "UniProt": semi(uni_accessions), "UniProt审校": "有" if reviewed else "无/未取得",
        "InterPro": semi(interpro), "Pfam": semi(pfam),
        "KEGG KO": cross.get("KEGG KO", ""), "KEGG EC": cross.get("KEGG EC", ""), "KEGG reaction": cross.get("KEGG reaction", ""),
        "BioCyc reaction": cross.get("BioCyc reaction ID", ""), "BioCyc EC": cross.get("BioCyc EC", ""),
        "Rhea反应定义": semi(rhea_ids), "BRENDA记录": cross.get("BRENDA EC记录", ""),
        "2016直接关联反应": semi(model_ids), "2016同EC比较对象": semi([r["Reaction ID"] for r in row["model_2016_same_ec"] if r["Reaction ID"] not in model_ids]),
        "精确反应状态": reaction_state, "原始论文DOI": found.get("evidence_doi", ""),
        "原始论文直接事实": direct, "证据边界": boundary, "本轮动作": action,
        "反应评分": "逐反应见 Alpha 行；未核定反应不评分", "审查状态": status,
        "文献ID查询命中": sum(h.get("hit_count", 0) for h in gene_search), "文献功能词查询命中": function_search.get("hit_count", 0),
        "文献候选DOI": semi(candidate_papers[:20]),
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
        "2026 reaction object": reaction_state + "; Rhea=" + semi(rhea_ids), "Legacy AFE locus": semi(old) if old else "未确认", "Current locus": locus,
        "Primary evidence DOI": found.get("evidence_doi", ""), "Evidence type": "RefSeq; UniProt/InterPro/Pfam; KEGG/BioCyc/Rhea/BRENDA; literature index; 2016 supplement",
        "Evidence boundary": boundary, "2026 action": action, "2026 reaction score": "不评分（逐反应待核）", "Source ID": semi(sources), "Score scope": "逐基因索引不赋反应分数",
    })
    for j, model in enumerate(row["model_2016"], 1):
        score = "待重评"
        scope = "2016 精确反应；2026 分数待原始证据与式子核对"
        if model["Reaction ID"] == "TSQOC":
            score = 1
            scope = "旧 q8 依赖精确式：当前无基因专属 q8 催化证据；旧建模假设"
        elif model["Reaction ID"] == "OPHHX":
            score = 1
            scope = "旧羟化精确式：UbiB 羟化 GPR 不受当前证据支持"
        elif model["Reaction ID"] == "HCO3E":
            score = 2
            scope = "旧水合精确式：功能注释支持；同株全细胞 CO2 数据不证明单酶方向"
        alpha.append({
            "Reaction ID": model["Reaction ID"], "Reaction Name": model["Reaction Name"], "Reaction Formula": model["Reaction Formula "],
            "Confidence Level": model["Confidence Level"], "EC Number": model["EC Number"], "PMID": model["PMID"], "Subsystem": model["Subsystem"],
            "Gene-Reaction Association": model["Gene-Reaction Association"], "Gene-Protein-Reaction Association": model["Gene-Protein-Reaction Association"], "Protein-Reaction-Association": model["Protein-Reaction-Association"],
            "Fe2 lb": model["lb"], "Fe2 ub": model["ub"], "tetrathionate lb": model["lb.1"], "tetrathionate ub": model["ub.1"], "sulfur lb": model["lb.2"], "sulfur ub": model["ub.2"],
            "Record class": "2016 linked reaction review (working)", "Candidate ID": f"G51-{row['ordinal']:04d}-R{j:02d}", "Candidate group": "G51", "2016 linked Reaction ID": model["Reaction ID"],
            "2026 reaction object": model["Reaction Formula "], "Legacy AFE locus": semi(old) if old else "未确认", "Current locus": locus,
            "Primary evidence DOI": found.get("evidence_doi", ""), "Evidence type": "2016 source row; current identity/database; curated original article where listed",
            "Evidence boundary": boundary + "; 2016 balance=" + balance[model["Reaction ID"]]["status"], "2026 action": action,
            "2026 reaction score": score, "Source ID": semi(sources), "Score scope": scope,
        })

assert len(index) == 293 and len({r["当前locus"] for r in index}) == 293
assert len(alpha) == 293 + sum(len(r["model_2016"]) for r in rows)
output = {"index": index, "alpha": alpha, "counts": {"genes": len(index), "alpha_rows": len(alpha), "linked_reaction_rows": len(alpha)-len(index), "curated_loci": len(curated), "review_status": dict(Counter(r["审查状态"] for r in index))}}
(HERE / "G51_working_tables.json").write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
print(output["counts"])
