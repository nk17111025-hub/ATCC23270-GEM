# -*- coding: utf-8 -*-
"""MILP GapFill supporting NEW intermediate metabolites + known-answer unit tests."""
from __future__ import annotations

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

M = 1000.0


def gapfill_milp(native_metas, native_rxns, native_bounds, candidates, biomass_rxn,
                 biomass_threshold=1e-6, max_time=300):
    """Minimize number of added candidate reactions for biomass >= threshold.

    native_rxns: dict rid -> {'stoich', 'lb', 'ub'}
    candidates: list of dict {'id','stoich','reversible'}
    native_bounds: dict rid -> (lb, ub) (already condition-applied)

    New candidate metabolites are added to the stoichiometric matrix automatically.
    Returns (status, k, added_ids, biomass).
    """
    # union of metabolites
    met_set = set(native_metas)
    for d in candidates:
        met_set.update(d["stoich"].keys())
    met_list = sorted(met_set)
    met_idx = {m: i for i, m in enumerate(met_list)}

    native_ids = list(native_rxns)
    cand_ids = [d["id"] for d in candidates]
    n_nat = len(native_ids)
    n_cand = len(cand_ids)
    n_met = len(met_list)
    V = n_nat + n_cand
    n_var = V + n_cand

    S = np.zeros((n_met, V))
    for j, r in enumerate(native_ids):
        for m, c in native_rxns[r]["stoich"].items():
            S[met_idx[m], j] = c
    for j, d in enumerate(candidates):
        for m, c in d["stoich"].items():
            S[met_idx[m], n_nat + j] = c

    c = np.zeros(n_var)
    c[V:] = 1.0
    integrality = np.zeros(n_var, dtype=int)
    integrality[V:] = 1
    lb = np.zeros(n_var); ub = np.zeros(n_var)
    for j, r in enumerate(native_ids):
        bl, bu = native_bounds[r]
        lb[j] = bl; ub[j] = bu
    for j, d in enumerate(candidates):
        lb[n_nat + j] = -M if d["reversible"] else 0.0
        ub[n_nat + j] = M
        lb[V + j] = 0.0; ub[V + j] = 1.0

    A_eq = np.zeros((n_met, n_var))
    A_eq[:, :V] = S
    n_link = n_cand + sum(1 for d in candidates if d["reversible"])
    n_ineq = n_link + 1
    A_ineq = np.zeros((n_ineq, n_var))
    b_lb = np.full(n_ineq, -np.inf)
    b_ub = np.zeros(n_ineq)
    row = 0
    for j, d in enumerate(candidates):
        A_ineq[row, n_nat + j] = 1.0; A_ineq[row, V + j] = -M; row += 1
        if d["reversible"]:
            A_ineq[row, n_nat + j] = -1.0; A_ineq[row, V + j] = -M; row += 1
    bio_j = native_ids.index(biomass_rxn) if biomass_rxn in native_ids else None
    if bio_j is None:
        return "OTHER_ERROR", None, [], 0.0
    A_ineq[row, bio_j] = -1.0
    b_ub[row] = -biomass_threshold
    cons = [LinearConstraint(A_eq, np.zeros(n_met), np.zeros(n_met)),
            LinearConstraint(A_ineq, b_lb, b_ub)]
    res = milp(c=c, integrality=integrality, bounds=Bounds(lb, ub), constraints=cons,
               options={"time_limit": max_time, "mip_rel_gap": 1e-6, "presolve": True,
                        "primal_feasibility_tolerance": 1e-8, "dual_feasibility_tolerance": 1e-8,
                        "mip_feasibility_tolerance": 1e-8})
    if res.x is None:
        return "INFEASIBLE", None, [], 0.0
    y = np.round(res.x[V:])
    added = [cand_ids[j] for j in range(n_cand) if y[j] > 0.5]
    return ("OPTIMAL" if res.success else "TIME_LIMIT", int(y.sum()), added,
            float(res.x[bio_j]))


# ---- known-answer unit tests ----


def test1_one_reaction():
    metas = {"A", "B"}
    rxns = {"bio": {"stoich": {"B": -1}, "lb": 0, "ub": 1000},
            "sinkA": {"stoich": {"A": 1}, "lb": 0, "ub": 1000}}
    bounds = {"bio": (0, 1000), "sinkA": (0, 1000)}
    cands = [{"id": "C1", "stoich": {"A": -1, "B": 1}, "reversible": False}]
    st, k, added, bio = gapfill_milp(metas, rxns, bounds, cands, "bio")
    return st == "OPTIMAL" and k == 1, st, k


def test2_new_intermediate():
    metas = {"A", "B"}
    rxns = {"bio": {"stoich": {"B": -1}, "lb": 0, "ub": 1000},
            "sinkA": {"stoich": {"A": 1}, "lb": 0, "ub": 1000}}
    bounds = {"bio": (0, 1000), "sinkA": (0, 1000)}
    cands = [{"id": "C1", "stoich": {"A": -1, "NEW_X": 1}, "reversible": False},
             {"id": "C2", "stoich": {"NEW_X": -1, "B": 1}, "reversible": False}]
    st, k, added, bio = gapfill_milp(metas, rxns, bounds, cands, "bio")
    return st == "OPTIMAL" and k == 2, st, k


def test3_impossible():
    metas = {"A", "B"}
    rxns = {"bio": {"stoich": {"B": -1}, "lb": 0, "ub": 1000},
            "sinkA": {"stoich": {"A": 1}, "lb": 0, "ub": 1000}}
    bounds = {"bio": (0, 1000), "sinkA": (0, 1000)}
    cands = [{"id": "C1", "stoich": {"A": -1, "C": 1}, "reversible": False}]
    st, k, added, bio = gapfill_milp(metas, rxns, bounds, cands, "bio")
    return st == "INFEASIBLE", st, k


def test4_duplicate_chemistry():
    metas = {"A", "B"}
    rxns = {"bio": {"stoich": {"B": -1}, "lb": 0, "ub": 1000},
            "sinkA": {"stoich": {"A": 1}, "lb": 0, "ub": 1000}}
    bounds = {"bio": (0, 1000), "sinkA": (0, 1000)}
    # identical chemistry, two aliases -> should still cost 1
    cands = [{"id": "C1", "stoich": {"A": -1, "B": 1}, "reversible": False},
             {"id": "C1_dup", "stoich": {"A": -1, "B": 1}, "reversible": False}]
    st, k, added, bio = gapfill_milp(metas, rxns, bounds, cands, "bio")
    return st == "OPTIMAL" and k == 1, st, k


def test5_reversible():
    metas = {"A", "B"}
    rxns = {"bio": {"stoich": {"B": -1}, "lb": 0, "ub": 1000},
            "sinkA": {"stoich": {"A": 1}, "lb": 0, "ub": 1000}}
    bounds = {"bio": (0, 1000), "sinkA": (0, 1000)}
    cands = [{"id": "C1", "stoich": {"A": -1, "B": 1}, "reversible": True}]
    st, k, added, bio = gapfill_milp(metas, rxns, bounds, cands, "bio")
    return st == "OPTIMAL" and k == 1, st, k


def test6_right_to_left():
    # source reaction marked RIGHT-TO-LEFT: B -> A means the usable direction is A <- B? no.
    # RIGHT-TO-LEFT means the physiological direction is right side -> left side.
    # Given listed (left: B, right: A), normalized stoich = A -> B (multiply by -1).
    from universe import normalize_direction, DIRECTION_MAP
    mode, rev = normalize_direction("PHYSIOL-RIGHT-TO-LEFT")
    # listed left=B, right=A  => raw stoich {B:-1, A:+1}; reversed => {B:+1, A:-1}
    raw = {"B": -1.0, "A": 1.0}
    if mode == "rev":
        norm = {m: -c for m, c in raw.items()}
    else:
        norm = dict(raw)
    ok = (norm == {"B": 1.0, "A": -1.0}) and rev is False
    return ok, mode, rev


def run_all_unit_tests():
    tests = [("test1_one_reaction", test1_one_reaction),
             ("test2_new_intermediate", test2_new_intermediate),
             ("test3_impossible", test3_impossible),
             ("test4_duplicate_chemistry", test4_duplicate_chemistry),
             ("test5_reversible", test5_reversible),
             ("test6_right_to_left", test6_right_to_left)]
    results = []
    for name, fn in tests:
        try:
            ok, *rest = fn()
            results.append((name, "PASS" if ok else "FAIL", rest))
        except Exception as e:
            results.append((name, "ERROR", str(e)))
    return results
