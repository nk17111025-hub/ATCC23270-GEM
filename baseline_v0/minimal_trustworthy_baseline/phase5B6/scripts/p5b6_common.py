# -*- coding: utf-8 -*-
"""Phase 5B-6 shared library.

Wraps the historical in-memory v2 model (phase4a_lib + build_glucose_model)
with the exact Phase 5B-5 R0 discovery constraints plus MILP / QC helpers.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
HIST = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B5"
if str(HIST) not in sys.path:
    sys.path.insert(0, str(HIST))

import phase4a_lib as lib  # noqa: E402
from build_glucose_model import build_memory  # noqa: E402

OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B6"

V2_T0_SHA256 = "97e31dd7b87852e8f639b72408e6e9ea05713dc0a48e99bc77edf087f0919e96"
V2_TH_SHA256 = "2a9c86d309d3f1e53e461f5a5df9cad65de974718db95ddc1bbe0945d917be39"

COND = "FIM"
BIOMASS_THRESHOLD = 1e-6
FEAS_TOL = 1e-7

# v2 correction set (reaction -> (lb, ub)) applied on top of build_memory
V2_CORRECTIONS = {
    "NADHI": (-1000.0, 1000.0),
}

BLANK = {
    "table1_lb_fe2": "", "table1_ub_fe2": "",
    "table1_lb_ttton": "", "table1_ub_ttton": "",
    "table1_lb_tsul": "", "table1_ub_tsul": "",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_v2(transport: str = "T0") -> Tuple[dict, dict]:
    metas, rxns = build_memory(transport)
    for rid, (lb, ub) in V2_CORRECTIONS.items():
        if rid in rxns:
            rxns[rid]["lb"] = lb
            rxns[rid]["ub"] = ub
    return metas, rxns


def setup(metas, rxns, cond, override=None):
    """LEGACY setup (v2 era). Applies full Table-1 internal+exchange bounds.

    Deprecated for v3+ work. Use v3_runtime.setup_v3, which applies only environmental
    (exchange) condition bounds and keeps internal reaction base bounds authoritative.
    """
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN:
            continue
        rx = rxns[r]
        lv = rx["table1_" + lbc]
        uv = rx["table1_" + ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv)
            ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl
            ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def fba(metas, rxns, cond, override=None):
    S, ml, rl, ri, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rl))
    c[ri[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return {r: res.x[ri[r]] for r in rl}


def max_biomass(metas, rxns, cond, override=None) -> float:
    v = fba(metas, rxns, cond, override)
    if v is None:
        return 0.0
    return v[lib.BIO_EX]


def reference_state(metas, rxns, cond):
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    v = fba(metas, rxns, cond, ov)
    if v is None:
        raise RuntimeError("reference state infeasible")
    donor = lib.DONOR_EX[cond]
    return {
        "biomass": v[lib.BIO_EX],
        "donor": v[donor],
        "o2": v["Ex_o2[e]"],
        "co2": v["Ex_h2co3[e]"],
        "rubisco": v["RUBISCO"] + v["RUBISCOX"],
    }


def r0_override(cond, ref, glucose=(-5.0, 0.0)):
    donor = lib.DONOR_EX[cond]
    return {
        "Ex_glc-B[e]": glucose,
        "Ex_h2co3[e]": (-1000.0, 0.0),
        donor: (ref["donor"], 0.0),
        "Ex_o2[e]": (ref["o2"], 0.0),
        "RUBISCO": (0.0, 0.0),
        "RUBISCOX": (0.0, 0.0),
    }


def add_reaction(rxns, rid, name, reversible, stoich, ec="", subsystem="candidate"):
    rx2 = dict(rxns)
    rx2[rid] = {
        "sbml_id": "R_" + rid,
        "name": name,
        "reversible": reversible,
        "lb": -1000.0 if reversible else 0.0,
        "ub": 1000.0,
        "stoich": dict(stoich),
        "confidence": "",
        "ec": ec,
        "pmid": "",
        "subsystem": subsystem,
        "gpr": "",
        "gpr2": "",
        "protein": "",
        **BLANK,
    }
    return rx2


def equation_str(stoich):
    lhs, rhs = [], []
    for met, c in stoich.items():
        if c < 0:
            lhs.append((met, -c))
        elif c > 0:
            rhs.append((met, c))

    def fmt(t):
        met, c = t
        cs = "" if abs(c - 1.0) < 1e-12 else (str(int(c)) if abs(c - round(c)) < 1e-12 else str(round(c, 6)))
        return (cs + " " + met).strip()

    return " + ".join(fmt(t) for t in sorted(lhs)) + " -> " + " + ".join(fmt(t) for t in sorted(rhs))


def stoich_key(stoich) -> tuple:
    return tuple(sorted(stoich.items()))


def reaction_in_model(rxns, stoich) -> bool:
    target = stoich_key(stoich)
    for r, rx in rxns.items():
        if stoich_key(rx["stoich"]) == target:
            return True
    return False


def reaction_equivalent_in_model(rxns, stoich) -> Optional[str]:
    """Return existing reaction id with identical stoichiometry, else None."""
    target = stoich_key(stoich)
    for r, rx in rxns.items():
        if stoich_key(rx["stoich"]) == target:
            return r
    return None


def write_tsv(path: Path, header: List[str], rows: List[List]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for row in rows:
            w.writerow(row)


def read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], rows[1:]
