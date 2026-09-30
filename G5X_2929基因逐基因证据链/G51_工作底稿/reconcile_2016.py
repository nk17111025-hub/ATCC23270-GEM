"""Link high-similarity legacy proteins to 2016 rows, retaining partial matches as conflicts."""
import json
import pandas as pd
from pathlib import Path

here = Path(__file__).resolve().parent
root = here.parents[1]
alignment = json.loads((here / "G51_unmapped_sequence_alignment.json").read_text(encoding="utf-8"))
table = pd.read_excel(root / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls", sheet_name="Table 1", dtype=str)
keys = ["Reaction ID", "Reaction Name", "Reaction Formula ", "Confidence Level", "EC Number", "PMID", "Subsystem",
        "Gene-Reaction Association", "Gene-Protein-Reaction Association", "Protein-Reaction-Association",
        "lb", "ub", "lb.1", "ub.1", "lb.2", "ub.2"]

out = []
for gene in alignment:
    top = gene.get("candidates", [{}])[0]
    status = "无可信 AFE 序列映射"
    if top.get("old_afe") and top["identity_fraction"] >= .99 and min(top["current_coverage"],top["old_coverage"]) >= .9:
        status = "高一致序列映射（仍以当前 WP 为主）"
    elif top.get("old_afe") and top["identity_fraction"] >= .99 and top["current_coverage"] >= .9:
        status = "仅当前蛋白覆盖的部分序列映射；不可继承旧全长蛋白功能"
    if gene["locus"] == "RU820_RS01505":
        status = "同一 WP 对应旧 AFE_0311/AFE_0325 两个旁系位点；仍不能唯一映射"
    old = top.get("old_afe", "") if status.startswith(("高一致", "仅当前蛋白")) else ""
    model = []
    if old:
        for _, row in table.iterrows():
            gpr = str(row.iloc[7])
            if old not in gpr:
                continue
            model.append({key: "" if pd.isna(value) else value for key, value in zip(keys,row.iloc[:16])})
    out.append({"locus": gene["locus"], "current_wp": gene["current_wp"], "status": status,
                "old_afe_candidate": old, "alignment": top, "2016_reactions": model})
(here / "G51_sequence_reconciled_2016.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
for gene in out:
    print(gene["locus"], gene["status"], gene["old_afe_candidate"], [r["Reaction ID"] for r in gene["2016_reactions"]])
