# -*- coding: utf-8 -*-
"""Framework-common helpers: structured solver status, element/charge balance, IO."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def write_tsv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for row in rows:
            w.writerow(row)


def write_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def classify_lp_status(res) -> str:
    """Map a scipy linprog OptimizeResult to a structured status string."""
    if res is None:
        return "OTHER_ERROR"
    if res.success:
        return "OPTIMAL"
    st = res.status
    msg = (res.message or "").lower()
    if st == 2 or "infeasible" in msg:
        return "INFEASIBLE"
    if st == 3 or "unbounded" in msg:
        return "UNBOUNDED"
    if "time" in msg or st == 1:
        return "TIME_LIMIT"
    if "numerical" in msg or st == 4:
        return "NUMERICAL_ERROR"
    return "OTHER_ERROR"


def solve_lp(c, A_eq, b_eq, bounds, method="highs", A_ub=None, b_ub=None):
    """Return (status, res). Never converts INFEASIBLE to a numerical zero."""
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method=method,
                  A_ub=A_ub, b_ub=b_ub)
    return classify_lp_status(res), res


def formula_to_elements(formula: str):
    if not formula:
        return {}
    f = re.sub(r"[\(\)\[\]]", "", formula.strip())
    out = {}
    for m in re.finditer(r"([A-Z][a-z]?)([0-9.]*)", f):
        out[m.group(1)] = out.get(m.group(1), 0.0) + float(m.group(2) or 1)
    return out


def balance_metabolites(stoich, fchg):
    """Full element + charge balance. fchg: {met: (formula, charge_float_or_None)}.

    Returns (classification, {element: delta}, charge_delta).
    classification in BALANCED / UNBALANCED / BALANCE_UNKNOWN.
    """
    total = {}
    charge = 0.0
    unknown = False
    for m, c in stoich.items():
        fchg_m = fchg.get(m)
        if fchg_m is None or not fchg_m[0]:
            unknown = True
            continue
        formula, ch = fchg_m
        for e, n in formula_to_elements(formula).items():
            total[e] = total.get(e, 0.0) + c * n
        if ch is None:
            unknown = True
        else:
            charge += c * float(ch)
    if unknown:
        return "BALANCE_UNKNOWN", {}, None
    residual = {e: round(v, 6) for e, v in total.items() if abs(v) > 1e-6}
    charge_delta = round(charge, 6)
    if not residual and abs(charge_delta) < 1e-6:
        return "BALANCED", {}, charge_delta
    return "UNBALANCED", residual, charge_delta
