# -*- coding: utf-8 -*-
"""MILP GapFill from v3 (v3_runtime.setup_v3). Native-R0 + large-universe search."""
from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

import p5b6_common as C
import phase4a_lib as lib
import v3_runtime as V3
from universe_curated import curated_dict

OUT = C.OUT
M = 1000.0
EPS = C.BIOMASS_THRESHOLD


def parse_equation(eq):
    stoich = {}
    try:
        lhs, rhs = eq.split("->")
    except ValueError:
        return None
    for side, sign in ((lhs, -1.0), (rhs, 1.0)):
        for term in side.split("+"):
            term = term.strip()
            if not term:
                continue
            p = term.split()
            if len(p) == 1:
                coef, met = 1.0, p[0]
            elif len(p) == 2:
                try:
                    coef = float(p[0])
                except ValueError:
                    coef, met = 1.0, term
                else:
                    met = p[1]
            else:
                return None
            stoich[met] = stoich.get(met, 0.0) + sign * coef
    return {m: c for m, c in stoich.items() if abs(c) > 1e-12}


def load_candidates():
    cands = []
    cur = curated_dict()
    for rid, d in cur.items():
        cands.append({"id": rid, "name": d["name"], "reversible": d["reversible"],
                      "stoich": d["stoich"], "ec": d["ec"], "source": d["source_db"], "tier": "U1"})
    u1 = OUT / "04_universe_clean" / "U1_clean.tsv"
    if u1.exists():
        _, rows = C.read_tsv(u1)
        for r in rows:
            s = parse_equation(r[3])
            if not s:
                continue
            cands.append({"id": "BC_" + r[0], "name": r[0], "reversible": r[4].strip().lower() == "true",
                          "stoich": s, "ec": r[8], "source": r[1], "tier": "U1"})
    return cands


def gapfill_v3(metas, rxns, cond, override, candidates, max_time=600):
    S_nat, met_list, rxn_list, rxn_idx, lb_nat, ub_nat = V3.setup_v3(metas, rxns, cond, override)
    met_set = set(met_list)
    met_idx = {m: i for i, m in enumerate(met_list)}
    valid = [d for d in candidates if all(m in met_set for m in d["stoich"])]
    n_nat = len(rxn_list)
    n_cand = len(valid)
    n_met = len(met_list)
    V = n_nat + n_cand
    n_var = V + n_cand

    S = np.zeros((n_met, n_nat + n_cand))
    S[:, :n_nat] = S_nat
    for j, d in enumerate(valid):
        for m, c in d["stoich"].items():
            S[met_idx[m], n_nat + j] = c

    c = np.zeros(n_var)
    c[V:] = 1.0
    integrality = np.zeros(n_var, dtype=int)
    integrality[V:] = 1
    lb = np.zeros(n_var)
    ub = np.zeros(n_var)
    lb[:n_nat] = lb_nat
    ub[:n_nat] = ub_nat
    for j, d in enumerate(valid):
        lb[n_nat + j] = -M if d["reversible"] else 0.0
        ub[n_nat + j] = M
        lb[V + j] = 0.0
        ub[V + j] = 1.0

    A_eq = np.zeros((n_met, n_var))
    A_eq[:, :n_nat + n_cand] = S
    n_link = n_cand + sum(1 for d in valid if d["reversible"])
    n_ineq = n_link + 1
    A_ineq = np.zeros((n_ineq, n_var))
    b_lb = np.full(n_ineq, -np.inf)
    b_ub = np.zeros(n_ineq)
    row = 0
    for j, d in enumerate(valid):
        A_ineq[row, n_nat + j] = 1.0
        A_ineq[row, V + j] = -M
        row += 1
        if d["reversible"]:
            A_ineq[row, n_nat + j] = -1.0
            A_ineq[row, V + j] = -M
            row += 1
    bio_j = rxn_idx[lib.BIO_EX]
    A_ineq[row, bio_j] = -1.0
    b_ub[row] = -EPS
    cons = [LinearConstraint(A_eq, np.zeros(n_met), np.zeros(n_met)),
            LinearConstraint(A_ineq, b_lb, b_ub)]
    bounds = Bounds(lb, ub)
    t0 = time.time()
    res = milp(c=c, integrality=integrality, bounds=bounds, constraints=cons,
               options={"time_limit": max_time, "mip_rel_gap": 1e-6, "presolve": True,
                        "primal_feasibility_tolerance": 1e-8, "dual_feasibility_tolerance": 1e-8,
                        "mip_feasibility_tolerance": 1e-8})
    elapsed = time.time() - t0
    if res.x is None:
        return {"status": res.message, "k": None, "added": [], "biomass": 0.0,
                "time": elapsed, "n_cand": n_cand}
    y = np.round(res.x[V:])
    added = [valid[j]["id"] for j in range(n_cand) if y[j] > 0.5]
    return {"status": res.message, "k": int(y.sum()), "added": added,
            "biomass": float(res.x[bio_j]), "time": elapsed, "n_cand": n_cand}


def main():
    metas, rxns = V3.load_v3("T0")
    # add the two native redox reactions as BASELINE (not selectable heterologous candidates)
    for rid, stoich, ec, gpr in [
        ("NDH2_NATIVE", {"nadh[c]": -1, "h[c]": -1, "q8[c]": -1, "nad[c]": 1, "q8h2[c]": 1},
         "1.6.5.9", "RU820_RS08555"),
        ("NOX_NATIVE", {"nadh[c]": -1, "h[c]": -1, "o2[c]": -0.5, "nad[c]": 1, "h2o[c]": 1},
         "1.6.3.4", "RU820_RS08315"),
    ]:
        rxns[rid] = {"sbml_id": "R_" + rid, "name": rid, "reversible": False, "lb": 0.0,
                     "ub": 1000.0, "stoich": dict(stoich), "confidence": "", "ec": ec,
                     "pmid": "", "subsystem": "Oxidative Phosphorylation", "gpr": gpr,
                     "gpr2": gpr, "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                     "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "",
                     "table1_ub_tsul": ""}
    cond = C.COND
    # native R0
    S, ml, rl, ri, lb, ub = V3.setup_v3(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0)})
    # reference for donor/O2 limits
    ref_ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    vref = V3.max_biomass_v3(metas, rxns, cond, ref_ov)
    donor = lib.DONOR_EX[cond]
    ref = {"donor": vref[donor], "o2": vref["Ex_o2[e]"]}
    ov = {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
          donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
          "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}
    b_native = V3.biomass_v3(metas, rxns, cond, ov)
    candidates = load_candidates()
    res = gapfill_v3(metas, rxns, cond, ov, candidates)
    res["native_R0_biomass"] = b_native
    print(json.dumps(res, ensure_ascii=False, indent=2))
    (OUT / "05_R0_search" / "phase5B6_v3_search_summary.tsv").write_text(
        "tier\tn_candidate\tsolver\tstatus\tmin_additions\tbiomass\tnote\n"
        f"U1\t{res['n_cand']}\tscipy.optimize.milp (HiGHS)\t{res['status']}\t{res['k']}\t{res['biomass']}\tv3-native\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
