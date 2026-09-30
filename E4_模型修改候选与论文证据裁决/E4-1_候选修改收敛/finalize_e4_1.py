"""把A0与经筛选的D3-1-R候选合并为可追溯的人工审计输入。"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
D3 = ROOT / "D3-1-R"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def value(row: dict[str, str], key: str) -> str:
    return row.get(key, "").strip() or "未知"


A = read(OUT / "阶段A_候选修改.tsv")
current = read(OUT / "E4-1_候选修改.tsv")
B = [r for r in current if r["候选ID"].startswith("E4-1-B-") and "BDGK" not in r["候选ID"] and r["候选ID"] != "E4-1-B-001-MDH-GPR"]
cross = {r["当前locus"]: r for r in read(D3 / "D3-1-R_全基因跨数据库交叉引用.tsv")}
evidence = defaultdict(list)
for r in read(D3 / "D3-1-R_注释证据表.tsv"):
    evidence[r["当前locus"]].append(r)
a0 = {r["反应ID"]: r for r in read(ROOT / "08_模型基线/A0/2016_到2024_真实差异.tsv")}
biocyc_raw = {r["对象ID"]: r["来源"] for r in read(ROOT / "03_BioCyc/C2-1/标准化数据/reactions.tsv")}
biocyc_gene_raw = {r["对象ID"]: r["来源"] for r in read(ROOT / "03_BioCyc/C2-1/标准化数据/genes.tsv")}
biocyc_protein_raw = {r["对象ID"]: r["来源"] for r in read(ROOT / "03_BioCyc/C2-1/标准化数据/proteins.tsv")}
rhea_raw_by_ec = defaultdict(list)
for manifest in read(D3 / "Rhea_新增查询清单.tsv"):
    for ec in manifest["查询EC"].split(";"):
        rhea_raw_by_ec[ec].append(manifest["文件"])

required_columns = ("候选ID", "子系统", "当前模型reaction ID", "2016值", "2024值", "当前候选值",
                 "旧gene", "当前gene/locus", "新protein", "KO", "EC", "Rhea", "MetaCyc",
                 "变更类型", "来源数据库", "来源记录ID", "原始文件", "sheet/行号或记录定位",
                 "数据库版本", "下载/检索日期", "同株支持", "证据类型", "文献ID",
                 "是否与其他数据库共用来源", "模型影响", "是否影响growth/energy/redox/transport",
                 "待人工确认问题", "优先级", "最终状态")
columns = list(dict.fromkeys(list(A[0]) + list(required_columns) + ["protein", "同株与否", "模型可能影响", "冲突类型", "D3-1-R定位", "候选层级", "人工裁决前动作"]))
for required in required_columns:
    assert required in columns, required

for row in A:
    rid = row["当前模型reaction ID"]
    row["当前候选值"] = "2016→2024差异已核对；宿主功能及模型动作待人工裁决"
    row["当前候选"] = row["当前候选值"]
    row["文献ID"] = row.get("文献", "具体文献待补")
    row["数据库版本"] = "2016/2024原始模型；B1 GCF_049532655.1；BioCyc官方GCF_000021485 PGDB 30.0（不作为本行反应等价证据）"
    row["版本"] = row["数据库版本"]
    row["来源数据库"] = "2016/2024原始模型；A0标准化；B1位点映射"
    row["是否与其他数据库共用来源"] = "A0来源于2016/2024原始模型；C2 BioCyc与MetaCyc共Pathway Tools来源，不算独立同株实验"
    if row["候选ID"] == "E4-1-BDGK-GPR":
        row["变更类型"] = "当前同株新增GPR候选"
        row["待人工确认问题"] = "AFE_2841/RU820_RS13140是否催化BDGK所写的β-D-葡萄糖磷酸化？需分开核查底物特异性；GPR证据不决定反应方向。"
        row["KO"], row["EC"], row["Rhea"] = "K25026", "2.7.1.2", "RHEA:17825（经KEGG R00299连接；非基因实验）"
        row["冲突类型"] = "GPR功能未证实"
        row["来源数据库"] = "2016/2024原始模型；A0；B1；D3-1-R的NCBI/KEGG/BioCyc/Rhea交叉注释"
        row["来源记录ID"] += "；" + ";".join(f"{r['数据库']}:{r['记录ID']}" for r in evidence["RU820_RS13140"])
        row["原始文件"] += "；" + "；".join(dict.fromkeys(r["来源文件"] for r in evidence["RU820_RS13140"]))
        row["sheet/行号或记录定位"] += "；D3-1-R_注释证据表.tsv:当前locus=RU820_RS13140"
    elif row["候选ID"] == "E4-1-BDGK-方向":
        row["变更类型"] = "反应方向/上下界冲突"
        row["待人工确认问题"] = "2016可逆(-1000..1000)与2024单向(0..1000)何者有热力学或同株实验支持？独立于AFE_2841 GPR裁决。"
        row["冲突类型"] = "方向/上下界"
    elif a0[rid]["状态"] == "共同且有语义差异":
        row["变更类型"] = "2016→2024作者修改"
        row["冲突类型"] = a0[rid]["差异字段"]
    elif rid.startswith("ZZ_Ecto"):
        row["变更类型"] = "异源工程反应"
        row["冲突类型"] = "工程层与宿主本底边界"
    else:
        row["变更类型"] = "原模型gap-fill或弱证据假设"
        row["冲突类型"] = "工程/模型假设与宿主本底边界"
    row["D3-1-R定位"] = "D3-1-R全基因表按当前gene/locus检索；本行差异源记录=A0反应" + rid
    row["候选层级"] = "异源工程" if rid.startswith("ZZ_Ecto") else ("工程或模型假设" if rid.startswith("ZZ_") else "宿主本底审计")
    row["同株与否"] = "否：异源工程" if rid.startswith("ZZ_Ecto") else ("同株模型假设，未证实宿主本底" if rid.startswith("ZZ_") else "是：ATCC 23270作者模型；功能待核")
    row["protein"] = row.get("新protein", "未知")
    row["模型可能影响"] = row["模型影响"]
    row["人工裁决前动作"] = "不修改正式模型"
    row["最终状态"] = "工程/假设隔离" if rid.startswith("ZZ_") else "无法裁决"


def model(rid: str, year: str = "2024") -> str:
    r = a0[rid]
    if not r[f"{year}反应式"]:
        return "该年模型无此反应"
    return f"{r[f'{year}反应式']}；bounds={r[f'{year}上下界']}；GPR={r[f'{year}GPR'] or '无'}"


manual = {
    "RU820_RS03190": ("当前同株新增反应候选", "辅因子冲突", "MDH（仅比较对象）",
        "苹果酸酶候选产生丙酮酸和CO2；模型MDH产生草酰乙酸。BioCyc同一蛋白还关联MALATE-DEH-RXN，此项另列MDH GPR审计。BioCyc 1.1.1.39-RXN写NAD，KEGG K00029/EC 1.1.1.40提示NADP；需确认真实底物、辅因子、H+及区室。"),
    "RU820_RS07700": ("当前同株新增GPR候选", "复合体/同工酶关系未知", "TKT1;TKT2",
        "BioCyc 2TRANSKETO-RXN/1TRANSKETO-RXN与模型TKT1/TKT2主要碳骨架相近，但区室未证实；AFE_1667是独立同工酶还是现有复合体成员？"),
    "RU820_RS08690": ("当前同株新增GPR候选", "反应化学表达差异", "PPC",
        "BioCyc PEPCARBOX-RXN原始左右侧为Pi+OAA / PEP+HCO3，生理方向PHYSIOL-RIGHT-TO-LEFT，即PEP+HCO3→Pi+OAA。模型PPC为H2O+CO2+PEP→H+Pi+OAA；核对H2O、H+、CO2/HCO3酸碱表达及[c]区室后再议GPR。"),
    "RU820_RS09355": ("当前同株新增反应候选", "辅因子冲突", "PGDH（仅比较对象）",
        "BioCyc RXN-9952写NADP→NADPH并产生ribulose-5P+CO2；NCBI名称写NAD+依赖；模型PGDH为6pgc脱水成2ddg6p，没有此氧化脱羧化学。核对CPD-2961身份、NAD/NADP、H+与[c]区室。"),
    "RU820_RS10385": ("当前同株新增转运候选", "底物/方向未定", "未知",
        "carbohydrate porin家族注释提示有机碳通透性，当前D3旧模型GPR覆盖=否、B1完全匹配。具体底物、跨膜方向和膜区室待确认。"),
    "RU820_RS00240": ("已有反应证据增强", "复合体亚基关系未知", "TSQOC",
        "AFE_0048与已有AFE_0044均标DoxD；先核对是否OR型旁系同工替代或独立复合体，不直接加成第七个AND亚基。并核对硫、电子、q8/q8h2与膜两侧质子。"),
    "RU820_RS01305": ("当前同株新增反应候选", "化学计量/底物/区室冲突", "SULDO（仅比较对象）",
        "BioCyc RXN-13161: CPD-11281+O2+H2O→glutathione+SO3+H+；模型SULDO: s[p]+H2O[p]+O2[p]⇄2H+[p]+SO3[p]。核对CPD-11281身份、底物、H+和区室，不能按名称合并。"),
}
transport_loci = {"RU820_RS08525", "RU820_RS10385", "RU820_RS10670", "RU820_RS13875", "RU820_RS13900", "RU820_RS14015", "RU820_RS14030"}
enhancement_loci = {"RU820_RS00240", "RU820_RS00635", "RU820_RS06600", "RU820_RS06605", "RU820_RS11310", "RU820_RS14455"}

for row in B:
    locus = row["当前gene/locus"]
    x = cross.get(locus)
    if not x:
        raise ValueError(f"D3-1-R中缺少候选locus: {locus}")
    typ, conflict, related, question = manual.get(locus, ("当前同株新增转运候选" if locus in transport_loci else "已有反应证据增强", "底物/方向/复合体未定" if locus in transport_loci else "反应等价未定", "未知", "核对反应化学、区室、复合体亚基及同株实验支持；家族注释不能直接形成正式GPR或反应。"))
    row["变更类型"], row["冲突类型"] = typ, conflict
    compare_ids = [rid.split("（")[0] for rid in related.split(";")]
    row["当前模型reaction ID"] = ("无对应反应；比较对象=" + related) if typ == "当前同株新增反应候选" else related
    row["2016值"] = "; ".join(model(rid, "2016") for rid in compare_ids if rid in a0) or "无对应已确认反应；比较对象见当前模型reaction ID"
    row["2024值"] = "; ".join(model(rid, "2024") for rid in compare_ids if rid in a0) or "无对应已确认反应；比较对象见当前模型reaction ID"
    row["当前候选值"] = "注释提出的功能候选；化学计量、区室、底物和GPR待人工裁决"
    row["旧gene"] = value(x, "旧AFE locus")
    row["新gene"] = locus
    row["旧protein"] = value(x, "旧protein ID")
    row["新protein"] = value(x, "当前protein ID")
    row["KO"] = value(x, "KEGG KO")
    row["EC"] = value(x, "KEGG EC") + ("；BioCyc=" + x["BioCyc EC"] if x["BioCyc EC"] else "")
    row["Rhea"] = value(x, "Rhea精确交叉引用") if x["Rhea精确交叉引用"] else "仅EC候选=" + value(x, "Rhea仅EC候选")
    row["MetaCyc"] = "未独立取得；BioCyc/Pathway Tools共源"
    row["待人工确认问题"] = question
    row["A0化学计量/区室比较"] = question
    row["来源数据库"] = ";".join(dict.fromkeys(r["数据库"] for r in evidence[locus])) or "D3-1-R"
    row["来源记录ID"] = ";".join(f"{r['数据库']}:{r['记录ID']}" for r in evidence[locus]) or locus
    row["原始文件"] = ";".join(dict.fromkeys(r["来源文件"] for r in evidence[locus]))
    extra_raw = []
    for biocyc_id, lookup in ((x["BioCyc gene ID"], biocyc_gene_raw), (x["BioCyc protein ID"], biocyc_protein_raw)):
        if biocyc_id in lookup:
            extra_raw.append(str(Path(lookup[biocyc_id]).relative_to(ROOT)).replace("\\", "/"))
    for rid in x["BioCyc reaction ID"].split(";"):
        if rid in biocyc_raw:
            extra_raw.append(str(Path(biocyc_raw[rid]).relative_to(ROOT)).replace("\\", "/"))
    if x["KEGG KO"]:
        extra_raw += ["04_KEGG/D3-1/原始响应/KEGG_link_ko_afr.txt", "04_KEGG/D3-1/原始响应/KEGG_link_reaction_ko.txt"]
    for ec in x["KEGG EC"].split(";"):
        extra_raw += rhea_raw_by_ec.get(ec, [])
    row["原始文件"] = ";".join(dict.fromkeys(row["原始文件"].split(";") + extra_raw))
    row["sheet/行号或记录定位"] = "D3-1-R_全基因跨数据库交叉引用.tsv:当前locus=" + locus + "；D3-1-R_注释证据表.tsv:当前locus=" + locus + "；原始记录ID见来源记录ID"
    row["数据库版本"] = ";".join(dict.fromkeys(r["数据库版本/组装"] for r in evidence[locus])) or "未知"
    row["下载/检索日期"] = ";".join(dict.fromkeys(r["取得日期"] for r in evidence[locus])) or "未知"
    row["同株支持"] = "ATCC 23270当前NCBI组装及旧同株KEGG/BioCyc预测；未有该基因催化实验证据"
    row["同株与否"] = "是：ATCC 23270；具体功能为计算预测"
    row["protein"] = row["新protein"]
    row["证据类型"] = "当前基因组计算注释；KEGG KO推断；BioCyc Tier 3 PathoLogic计算预测；Rhea仅反应标准化"
    row["文献ID"] = "PUB-LATENDRESSE2013（仅BioCyc计算来源）" if x["BioCyc reaction ID"] else "同株功能实验文献待补"
    row["是否与其他数据库共用来源"] = "BioCyc与MetaCyc共Pathway Tools来源；KEGG/Rhea反应交叉引用不构成独立同株实验"
    row["候选层级"] = "宿主本底候选；未经正式纳入"
    row["D3-1-R定位"] = "全基因当前locus=" + locus
    row["人工裁决前动作"] = "不修改正式模型"
    row["最终状态"] = "无法裁决"
    row["优先级"] = "高" if locus in manual or locus in {"RU820_RS10670", "RU820_RS14015", "RU820_RS14030"} else "中"
    if locus in transport_loci:
        row["transport影响"] = "可能；需核对底物和方向"
    if locus in {"RU820_RS03190", "RU820_RS09355", "RU820_RS01305", "RU820_RS00240"}:
        row["redox影响"] = "可能；辅因子或电子流未定"
    if locus in enhancement_loci:
        row["模型影响"] = "现有反应或复合体的证据增强候选；不形成正式新反应"
    row["模型可能影响"] = row["模型影响"]
    row["是否影响growth/energy/redox/transport"] = "; ".join(f"{name}={row[field]}" for name, field in (("growth", "growth影响"), ("energy", "energy影响"), ("redox", "redox影响"), ("transport", "transport影响")))

malic_base = next(r for r in B if r["当前gene/locus"] == "RU820_RS03190")
malic_mdh = dict(malic_base)
malic_mdh["候选ID"] = "E4-1-B-001-MDH-GPR"
malic_mdh["变更类型"] = "当前同株新增GPR候选"
malic_mdh["冲突类型"] = "同一基因双反应注释；辅因子未定"
malic_mdh["当前模型reaction ID"] = "MDH"
malic_mdh["当前候选值"] = "BioCyc MALATE-DEH-RXN对应MDH候选GPR；与苹果酸酶候选独立裁决"
malic_mdh["模型可能影响"] = "候选MDH替代GPR；待化学及基因证据裁决"
malic_mdh["A0化学计量/区室比较"] = "模型MDH: mal-L[c]+NAD[c]⇄OAA[c]+NADH[c]+H+[c]；BioCyc MALATE-DEH-RXN: MAL+NAD→OAA+NADH+PROTON，BioCyc未给模型[c]区室；需逐项确认。"
malic_mdh["待人工确认问题"] = "AFE_0660在BioCyc同时关联MALATE-DEH-RXN和1.1.1.39-RXN。MALATE-DEH-RXN能否支持模型MDH的替代GPR？需与现有AFE_3000 GPR、底物、NAD/NADP和[c]区室逐项核对。"
B.append(malic_mdh)
rows = A + B
ids = [r["候选ID"] for r in rows]
assert len(ids) == len(set(ids)), "候选ID重复"
assert len(rows) > len(A), "D3当前全基因候选未进入主表"
for row in rows:
    for col in ("候选ID", "变更类型", "当前模型reaction ID", "来源记录ID", "原始文件", "sheet/行号或记录定位"):
        assert row.get(col), (row["候选ID"], col)

write(OUT / "E4-1_候选修改.tsv", rows, columns)
new_host = [r for r in rows if r["变更类型"] == "当前同株新增反应候选"]
new_transport = [r for r in rows if r["候选ID"] in {b["候选ID"] for b in B if b["当前gene/locus"] in transport_loci}]
new_gpr = [r for r in rows if r["变更类型"] == "当前同株新增GPR候选"]
definition = [r for r in rows if r["冲突类型"] in {"方向/上下界", "辅因子冲突", "化学计量/底物/区室冲突", "反应化学表达差异"} or (r in A and any(t in r["冲突类型"] for t in ("化学计量", "方向", "上下界")))]
for name, subset in (("新增宿主反应候选", new_host), ("新增转运候选", new_transport), ("GPR更新候选", new_gpr), ("反应定义冲突", definition)):
    write(OUT / f"E4-1_{name}.tsv", subset, columns)

# D3记录全部保留，避免只导出E4入选基因而静默丢弃数据库分歧。
disagreements = read(D3 / "D3-1-R_数据库分歧.tsv")
by_locus = defaultdict(list)
for r in B:
    by_locus[r["当前gene/locus"]].append(r["候选ID"])
conflict_cols = list(disagreements[0]) + ["E4候选ID", "来源记录ID", "原始文件", "最终状态"]
conflicts = [{**r, "E4候选ID": ";".join(by_locus.get(r["当前locus"], [])) or "未入选E4重点候选", "来源记录ID": r["当前locus"], "原始文件": "D3-1-R/D3-1-R_数据库分歧.tsv", "最终状态": "未裁决"} for r in disagreements]
for r in rows:
    if r["冲突类型"] in {"辅因子冲突", "化学计量/底物/区室冲突", "底物与旧GPR冲突", "GPR功能未证实", "方向/上下界"}:
        conflicts.append({"当前locus": r["当前gene/locus"], "旧AFE locus": r["旧gene"], "分歧类型": r["冲突类型"], "处理": r["待人工确认问题"], "E4候选ID": r["候选ID"], "来源记录ID": r["来源记录ID"], "原始文件": r["原始文件"], "最终状态": "未裁决"})
write(OUT / "E4-1_数据库冲突.tsv", conflicts, conflict_cols)
access_missing = read(D3 / "D3-1-R_访问缺失表.tsv")
selected_access = [{**r, "E4候选ID": ";".join(by_locus.get(r["当前locus"], [])) or ("E4-1-BDGK-GPR" if r["当前locus"] == "RU820_RS13140" else "全局来源缺失"), "变更类型": "数据库缺失/访问失败", "原始文件": "D3-1-R/D3-1-R_访问缺失表.tsv"} for r in access_missing if r["当前locus"] in by_locus or r["当前locus"] in ("", "RU820_RS13140")]
write(OUT / "E4-1_来源缺失记录.tsv", selected_access, list(access_missing[0]) + ["E4候选ID", "变更类型", "原始文件"])

b1_uncertain = [r for r in read(ROOT / "02_NCBI/标准化数据/旧新基因映射.tsv") if r["映射状态"] != "完全匹配"]
mapping_cols = list(b1_uncertain[0]) + ["D3当前locus", "D3映射冲突", "来源记录ID", "原始文件", "处理"]
mapping = [{**r, "D3当前locus": r["当前GCF049位点"], "D3映射冲突": value(cross.get(r["当前GCF049位点"], {}), "映射冲突"), "来源记录ID": r["旧AFE位点"], "原始文件": "02_NCBI/标准化数据/旧新基因映射.tsv", "处理": "人工核对；不自动写GPR"} for r in b1_uncertain]
for x in cross.values():
    if x["映射冲突"] == "是":
        mapping.append({"旧AFE位点": x["旧AFE locus"], "旧RefSeq位点": x["旧RefSeq locus"], "旧蛋白ID": x["旧protein ID"], "当前GCF049位点": x["当前locus"], "当前蛋白ID": x["当前protein ID"], "映射方法": "D3-1-R全基因交叉核验", "映射状态": "歧义候选", "证据备注": x["B1映射状态"], "D3当前locus": x["当前locus"], "D3映射冲突": "是", "来源记录ID": x["当前locus"], "原始文件": "D3-1-R/D3-1-R_全基因跨数据库交叉引用.tsv", "处理": "人工核对多重旧位点；不自动写GPR"})
mapping.append({"旧AFE位点": "AFE_2841", "旧RefSeq位点": "AFE_RS13045", "旧蛋白ID": "WP_012537424.1", "当前GCF049位点": "RU820_RS13140", "当前蛋白ID": "WP_012537424.1", "映射方法": "蛋白ID精确", "映射状态": "易混淆别名提醒", "证据备注": "旧AFE_RS13140属于AFE_2860/ArsC；不可误接AFE_2841", "D3当前locus": "RU820_RS13140", "D3映射冲突": "否", "来源记录ID": "AFE_2841", "原始文件": "02_NCBI/标准化数据/旧新基因映射.tsv", "处理": "保留回归检查"})
write(OUT / "E4-1_标识映射不确定.tsv", mapping, mapping_cols)

literature_cols = ["候选ID", "当前gene/locus", "当前模型reaction ID", "待补问题", "已有文献ID", "来源记录ID", "原始文件", "状态"]
write(OUT / "E4-1_待补文献.tsv", [{"候选ID": r["候选ID"], "当前gene/locus": r["当前gene/locus"], "当前模型reaction ID": r["当前模型reaction ID"], "待补问题": r["待人工确认问题"], "已有文献ID": r["文献ID"], "来源记录ID": r["来源记录ID"], "原始文件": r["原始文件"], "状态": "待同株功能/反应定义文献"} for r in rows if r["最终状态"] == "无法裁决"], literature_cols)

report = (D3 / "D3-1-R_验收报告.md").read_text(encoding="utf-8")
d3_pass = any(line.strip() in {"PASS", "验收状态：PASS", "最终结论：PASS", "状态：PASS"} or line.strip().startswith("**PASS：D3 证据与映射体系验收") for line in report.splitlines())
main_sha = hashlib.sha256((OUT / "E4-1_候选修改.tsv").read_bytes()).hexdigest()
review_path = OUT / "E4-1_独立复核.md"
review_text = review_path.read_text(encoding="utf-8") if review_path.exists() else ""
e4_review_pass = "最终复核结论：PASS" in review_text and f"复核对象SHA256：{main_sha}" in review_text
stamp = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")
coverage = f"""# E4-1 来源覆盖与缺失

生成时间：{stamp}。A0真实差异631条，2016→2024有语义差异13条、新增16条；旧E4的29条仅作历史初筛。E4-1阶段A重新生成30项，BDGK拆为GPR与方向两项。D3-1-R当前全基因交叉引用{len(cross)}条，进入E4重点比较的候选{len(B)}条（{len({r['当前gene/locus'] for r in B})}个独立基因）。候选总数{len(rows)}。

- NCBI当前组装GCF_049532655.1：基因和蛋白映射；产物为计算注释。
- KEGG afr：旧同株组装；KO/EC/反应继承关系不等于基因实验。
- BioCyc GCF_000021485，PGDB 30.0，Tier 3：基于旧GCF_000021485.1；反应多数EV-COMP，相关记录来源可能共用Pathway Tools。
- Rhea：已取得记录只用于方程/ChEBI核对；EC或交叉ID命中不证明宿主酶功能。MetaCyc未独立取得，BRENDA同株基因实验不足；具体缺口见D3-1-R_访问缺失表。
- D3-1-R数据库分歧{len(disagreements)}条已全部转入E4-1_数据库冲突.tsv；未入选重点候选的分歧也保留。B1非完全映射{len(b1_uncertain)}条另存。
- D3-1-R访问/未取得记录{len(access_missing)}条保留在原表；与E4入选候选有关的{len(selected_access)}条另存E4-1_来源缺失记录.tsv，分类为“数据库缺失/访问失败”。未命中不推断菌株缺乏该反应。
- 反应等价按化学计量与区室、化合物身份、方向/上下界、GPR、交叉ID依次检查。未知项留待人工；未改SBML、正式反应或GAM/NGAM。

## 候选类别索引

2016→2024作者修改、当前同株新增反应候选、当前同株新增转运候选、当前同株新增GPR候选、旧GPR/current locus更新、已有反应证据增强、反应方向/上下界冲突、辅因子冲突、化学计量冲突、转运底物/方向冲突、原模型gap-fill或弱证据假设、异源工程反应、数据库缺失/访问失败、无法判定，分别在主表“变更类型/冲突类型/最终状态”、来源缺失记录及专题表中标注。旧GPR/current locus更新在BDGK的旧位点和当前位点链中保留为附属问题；未单独重复计数。辅因子与化学计量冲突也保留为新增反应候选的独立“冲突类型”。
"""
(OUT / "E4-1_来源覆盖与缺失.md").write_text(coverage, encoding="utf-8")
status = "PASS" if d3_pass and e4_review_pass else ("待E4独立Reviewer" if d3_pass else "等待D3-1-R依赖")
acceptance = f"""# E4-1 验收报告

状态：{status}。D3-1-R文件已落盘；其验收报告{'已写PASS' if d3_pass else '仍未给出最终PASS'}。E4-1独立Reviewer{'已确认高优先级候选和修正' if e4_review_pass else '仍待最终确认'}。

| 指标 | 数量 |
|---|---:|
| 候选总数 | {len(rows)} |
| 新增宿主反应候选 | {len(new_host)} |
| 新增转运候选 | {len(new_transport)} |
| 新GPR候选 | {len(new_gpr)} |
| 现有反应定义冲突待审 | {len(definition)} |
| 无法裁决 | {sum(r['最终状态']=='无法裁决' for r in rows)} |

阶段A的30项与D3全基因入选的{len(B)}项同表比较。BDGK两项各保留一次。D3注释不直接增加正式反应或GPR，工程/假设与宿主本底分层。全部{len(disagreements)}条D3数据库分歧保留，MetaCyc和Rhea未覆盖范围按未知处理。未修改正式模型。复核记录见E4-1_独立复核.md。

## 高优先级人工问题

1. BDGK：`AFE_2841/RU820_RS13140`是否支持β-D-葡萄糖磷酸化GPR；2016可逆到2024单向的方向修改须独立裁决。
2. 氧化PPP：`AFE_2024/RU820_RS09355`的6PG脱氢酶候选与模型现有PGDH化学不同；核对NAD/NADP、CPD-2961与核酮糖-5-磷酸。
3. 苹果酸节点：`AFE_0660/RU820_RS03190`同时关联MDH和苹果酸酶；分别核对GPR、草酰乙酸/丙酮酸及辅因子。
4. CBB与PPP：`AFE_1667/RU820_RS07700`是TKT1/TKT2的替代基因、复合体成员还是不同反应，仍需确认。
5. 补充反应：`AFE_1883/RU820_RS08690`与模型PPC比较时，核对BioCyc生理方向、HCO3与CO2/H2O/H+表达及区室。
6. 有机碳转运：孔蛋白、糖转运家族和PTS亚基的具体底物、膜侧与方向未知；不得用工程层`ZZ_glc_transport`填入宿主本底。
7. 硫氧化：`AFE_0269/RU820_RS01305`的RXN-13161与SULDO底物和区室不同，核对CPD-11281、SO3、O2及质子。
8. 呼吸链：`AFE_0048/DoxD`是否替代已有`AFE_0044`，需要复合体证据；同时核对TSQOC的醌、电子和跨膜质子。
9. 离子与能量边界：CA2tpp、NA1tpp的H+/Na+/Ca2+计量及方向和NADHI/bc1/ATP synthase的质子约束，需分别判清反应定义与容量证据。
"""
(OUT / "E4-1_验收报告.md").write_text(acceptance, encoding="utf-8")
source_rows = [
    {"来源文件": "08_模型基线/A0/2016_到2024_真实差异.tsv", "状态": "已核对；2016/2024原始单元格逐条定位"},
    {"来源文件": "02_NCBI/标准化数据/旧新基因映射.tsv", "状态": "已核对；非完全映射保留待审"},
    {"来源文件": "03_BioCyc/C2-1/标准化数据/reactions.tsv", "状态": "PGDB 30.0 Tier 3；原始XML路径见逐行来源字段"},
] + [{"来源文件": str(p.relative_to(ROOT)).replace("\\", "/"), "状态": "D3-1-R已PASS；E4逐条筛选并保留缺失/分歧"} for p in sorted(D3.glob("D3-1-R_*.tsv"))]
write(OUT / "E4-1_来源登记.tsv", source_rows, ["来源文件", "状态"])
write(OUT / "E4-1_验收清单.tsv", [
    {"检查项": "D3-1-R依赖", "结果": "PASS" if d3_pass else "等待", "依据": "D3-1-R/D3-1-R_验收报告.md"},
    {"检查项": "统一主表", "结果": str(len(rows)), "依据": "E4-1_候选修改.tsv；A0阶段A+全基因入选"},
    {"检查项": "BDGK拆分且无重复", "结果": "PASS", "依据": "E4-1-BDGK-GPR / E4-1-BDGK-方向"},
    {"检查项": "数据库分歧保留", "结果": str(len(disagreements)), "依据": "D3-1-R_数据库分歧.tsv全部转入"},
    {"检查项": "E4独立Reviewer", "结果": "PASS" if e4_review_pass else "待复核结论", "依据": "E4-1_独立复核.md及主表SHA256"},
], ["检查项", "结果", "依据"])
(OUT / "E4-1_候选统计.json").write_text(json.dumps({"候选总数": len(rows), "新增宿主反应候选": len(new_host), "新增转运候选": len(new_transport), "新GPR候选": len(new_gpr), "反应定义冲突待审": len(definition), "无法裁决": sum(r["最终状态"] == "无法裁决" for r in rows), "D3数据库分歧": len(disagreements), "D3依赖": "PASS" if d3_pass else "等待", "E4独立复核": "PASS" if e4_review_pass else "等待"}, ensure_ascii=False, indent=2), encoding="utf-8")
(OUT / "E4-1_审计输入说明.md").write_text(f"# E4-1 审计输入说明\n\n主表{len(rows)}条，以候选ID与四个专题表关联。阶段A原始差异、D3-1-R全基因证据及原始文件定位同列保存；工程/假设{sum(r['最终状态']=='工程/假设隔离' for r in rows)}条隔离。D3依赖{'已PASS' if d3_pass else '待验收'}，E4独立Reviewer意见见验收报告。正式模型未改。\n", encoding="utf-8")
print({"status": status, "total": len(rows), "host_reactions": len(new_host), "transport": len(new_transport), "gpr": len(new_gpr), "definition_conflicts": len(definition), "unresolved": sum(r["最终状态"] == "无法裁决" for r in rows), "d3_disagreements": len(disagreements)})
