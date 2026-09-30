# -*- coding: utf-8 -*-
"""12-combination physiological regression for updatedv3.1.

Conditions: FIM / TTM / TSM.
Perturbations: WT / no_donor / no_carbon / glucose_closed.

Runs the current authoritative model (models/current/updatedv3.1.xml) and compares
objective values against the frozen validated regression targets with tolerance 1e-9.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import scenarios  # noqa: E402


MODEL = REPO / "models" / "current" / "updatedv3.1.xml"
CONFIG = REPO / "config" / "scenarios.json"
OUT = REPO / "tests" / "regression_v3.1.tsv"

# Frozen validated regression targets (validation targets only; do not infer scenarios).
TARGETS = {
    "FIM_WT": 0.673529948005047,
    "FIM_no_donor": None,   # infeasible
    "FIM_no_carbon": 0.7302586898902611,
    "FIM_glucose_closed": 0.052076387102693945,
    "TTM_WT": 2.969254977069056,
    "TTM_no_donor": 0.6735299480050473,
    "TTM_no_carbon": 3.0134966113352473,
    "TTM_glucose_closed": 0.05207638710271767,
    "TSM_WT": 2.969254977069056,
    "TSM_no_donor": 0.6735299480050473,
    "TSM_no_carbon": 3.0134966113352473,
    "TSM_glucose_closed": 0.05207638710271767,
}
TOL = 1e-9


def main() -> int:
    model = scenarios.load_model(MODEL)
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    rows = []
    failures = []
    for cond in ("FIM", "TTM", "TSM"):
        for pert in ("WT", "no_donor", "no_carbon", "glucose_closed"):
            status, value, _ = scenarios.optimize(model, cond, pert, config)
            key = f"{cond}_{pert}"
            target = TARGETS[key]
            if target is None:
                ok = value is None
                value_str = "infeasible" if value is None else str(value)
            else:
                ok = value is not None and abs(value - target) <= TOL
                value_str = "" if value is None else repr(value)
            rows.append([cond, pert, status, value_str, "" if target is None else repr(target), "PASS" if ok else "FAIL"])
            if not ok:
                failures.append(key)

    header = ["condition", "test", "status", "objective_value", "target", "result"]
    with OUT.open("w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")

    print(f"REGRESSION {len(rows) - len(failures)}/{len(rows)} PASS")
    if failures:
        print("FAILURES:", ", ".join(failures))
        return 1
    print("WROTE", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
