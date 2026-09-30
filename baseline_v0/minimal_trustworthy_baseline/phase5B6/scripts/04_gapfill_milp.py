# -*- coding: utf-8 -*-
"""MILP GapFill: minimum heterologous additions for Rubisco-independent (R0) growth."""
from __future__ import annotations

import json
import sys
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

import p5b6_common as C
from universe_curated import curated_dict

OUT = C.OUT
SOLVER = "scipy.optimize.milp (HiGHS)"
M = 1000.0
EPS = C.BIOMASS_THRESHOLD

# MODEL_REPAIR (section 47): the v2 direction corrections (NADHI/GAPD1/GAPD2/PGK/G6PDH2
# opened) are stored in the base bounds but are silently overridden by the Table-1 FIM
# condition bounds (lb_fe2=0). Apply them explicitly as overrides so the promoted v2 state
# is honored under the R0/FIM discovery condition.
V2_REPAIR = {
    "NADHI": (-1000.0, 1000.0),
    "GAPD1": (-1000.0, 1000.0),
    "GAPD2": (-1000.0, 1000.0),
    "PGK": (-1000.0, 1000.0),
    "G6PDH2": (-1000.0, 1000.0),
}


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
            parts = term.split()
            if len(parts) == 1:
                coef, met = 1.0, parts[0]
            elif len(parts) == 2:
                try:
                    coef = float(parts[0])
                except ValueError:
                    coef, met = 1.0, term
                else:
                    met = parts[1]
            else:
                return None
            stoich[met] = stoich.get(met, 0.0) + sign * coef
    return {m: c for m, c in stoich.items() if abs(c) > 1e-12}


def load_candidates():
    cands = []
    cur = curated_dict()
    for rid, d in cur.items():
        cands.append({
            "id": rid, "name": d["name"], "reversible": d["reversible"],
            "stoich": d["stoich"], "ec": d["ec"], "source": d["source_db"], "tier": "U1",
        })
    u1_path = OUT / "04_universe_clean" / "U1_clean.tsv"
    if u1_path.exists():
        _, rows = C.read_tsv(u1_path)
        for r in rows:
            eq = r[3]
            rev = r[4].strip().lower() == "true"
            stoich = parse_equation(eq)
            if not stoich:
                continue
            cands.append({
                "id": "BC_" + r[0], "name": r[0], "reversible": rev,
                "stoich": stoich, "ec": r[8], "source": r[1], "tier": "U1",
            })
    return cands


def gapfill(metas, rxns, cond, override, candidates, max_time=600):
    import phase4a_lib as lib
    override = dict(override)
    override.update(V2_REPAIR)
    S_nat, met_list, rxn_list, rxn_idx, lb_nat, ub_nat = C.setup(metas, rxns, cond, override)
    met_set = set(met_list)
    met_idx = {m: i for i, m in enumerate(met_list)}

    valid = [d for d in candidates if all(m in met_set for m in d["stoich"])]
    n_nat = len(rxn_list)
    n_cand = len(valid)
    n_met = len(met_list)
    V = n_nat + n_cand            # flux variable count
    n_var = V + n_cand            # + binary y

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

    cons = [
        LinearConstraint(A_eq, np.zeros(n_met), np.zeros(n_met)),
        LinearConstraint(A_ineq, b_lb, b_ub),
    ]
    bounds = Bounds(lb, ub)

    t0 = time.time()
    opts = {
        "time_limit": max_time,
        "mip_rel_gap": 1e-6,
        "presolve": True,
        # Tighten feasibility tolerances so that the biomass threshold (1e-6)
        # is strictly enforced (default HiGHS MIP tolerance of 1e-6 is ambiguous).
        "primal_feasibility_tolerance": 1e-8,
        "dual_feasibility_tolerance": 1e-8,
        "mip_feasibility_tolerance": 1e-8,
    }
    res = milp(c=c, integrality=integrality, bounds=bounds, constraints=cons, options=opts)
    elapsed = time.time() - t0

    if res.x is None:
        return {"status": res.message, "k": None, "added": [], "biomass": 0.0,
                "time": elapsed, "n_cand": n_cand, "n_nat": n_nat, "optimal": False}
    y = np.round(res.x[V:])
    added = [valid[j]["id"] for j in range(n_cand) if y[j] > 0.5]
    v = res.x[:V]
    biomass = v[bio_j]
    return {"status": res.message, "k": int(y.sum()), "added": added,
            "biomass": float(biomass), "time": elapsed, "n_cand": n_cand,
            "n_nat": n_nat, "optimal": bool(res.success), "fun": float(res.fun),
            "mip_gap": float(res.mip_gap) if getattr(res, "mip_gap", None) is not None else None}


def main():
    metas, rxns = C.build_v2("T0")
    cond = C.COND
    ref = C.reference_state(metas, rxns, cond)
    override = C.r0_override(cond, ref)
    candidates = load_candidates()
    print(f"candidates loaded: {len(candidates)}", flush=True)
    res = gapfill(metas, rxns, cond, override, candidates, max_time=600)
    print(json.dumps(res, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
