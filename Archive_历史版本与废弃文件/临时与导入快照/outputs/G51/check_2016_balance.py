"""Check atom and charge balance of unique 2016 reactions linked to G51."""
import json
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MODEL = ROOT / "A0_复现2016与2024模型/A0_主任务/脚本与环境/中间/mmc1.xlsx"
rows = json.loads((HERE / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
reactions = {r["Reaction ID"]: r for g in rows for r in g["model_2016"]}
met = pd.read_excel(MODEL, sheet_name="Table 2", header=0).fillna("")
metabolites = {str(r["Metabolite Abbreviation"]).strip(): (str(r["Charged Formula"]).strip(), r["Charge"]) for _, r in met.iterrows()}


def elements(formula):
    if not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
        return None
    out = defaultdict(int)
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        out[element] += int(number or 1)
    return out


def side(text):
    out = []
    for term in text.split(" + "):
        match = re.fullmatch(r"(?:(\d+(?:\.\d+)?) )?(.+)", term.strip())
        if not match:
            raise ValueError(term)
        out.append((float(match.group(1) or 1), match.group(2).strip()))
    return out


results = []
for rid, row in reactions.items():
    formula = row["Reaction Formula "].strip()
    parts = re.split(r"\s+(?:<=>|->)\s+", formula)
    missing = []
    unknown = []
    delta = defaultdict(float)
    charge = 0.0
    if len(parts) != 2:
        unknown.append("equation syntax")
    else:
        for sign, text in [(-1, parts[0]), (1, parts[1])]:
            for coeff, met_id in side(text):
                if met_id not in metabolites:
                    missing.append(met_id)
                    continue
                chem, q = metabolites[met_id]
                elems = elements(chem)
                if elems is None:
                    unknown.append(met_id + ":" + chem)
                    continue
                for element, count in elems.items():
                    delta[element] += sign * coeff * count
                try:
                    charge += sign * coeff * float(q)
                except Exception:
                    unknown.append(met_id + ":charge=" + str(q))
    discrepancy = {k: round(v, 8) for k, v in delta.items() if abs(v) > 1e-8}
    if missing or unknown:
        status = "无法完整核算"
    elif discrepancy or abs(charge) > 1e-8:
        status = "不平衡"
    else:
        status = "平衡"
    results.append({"reaction_id": rid, "formula": formula, "status": status, "element_delta": discrepancy, "charge_delta": round(charge, 8), "missing_metabolites": missing, "unknown_formula": unknown})

(HERE / "2016_G51_balance.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
from collections import Counter
print("reactions", len(results), Counter(r["status"] for r in results))
for r in results:
    if r["status"] != "平衡":
        print(r["reaction_id"], r["status"], r["element_delta"], r["charge_delta"], r["missing_metabolites"][:3], r["unknown_formula"][:3])
