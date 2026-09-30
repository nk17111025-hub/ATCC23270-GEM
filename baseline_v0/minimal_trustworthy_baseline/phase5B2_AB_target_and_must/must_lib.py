# -*- coding: utf-8 -*-
"""Shared LP helpers for Phase 5B-2A/2B (FVA-based MUST analysis)."""
import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib

TOL = 1e-6


def setup(metas, rxns, cond, override):
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


def lp(metas, rxns, cond, c_obj, override, biomass_min=None):
    """Minimize c_obj (dict rxn->weight) subject to S v = 0, bounds, optional biomass >= min."""
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    n = len(rxn_list)
    c = np.zeros(n)
    for r, w in c_obj.items():
        c[rxn_idx[r]] = w
    kwargs = dict(A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if biomass_min is not None:
        A_ub = np.zeros((1, n))
        A_ub[0, rxn_idx[lib.BIO_EX]] = -1.0
        kwargs["A_ub"] = A_ub
        kwargs["b_ub"] = np.array([-biomass_min])
    res = linprog(c, **kwargs)
    if not res.success:
        return None
    return res.x


def max_biomass(metas, rxns, cond, override):
    x = lp(metas, rxns, cond, {lib.BIO_EX: -1.0}, override)
    if x is None:
        return None
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    return x[rxn_idx[lib.BIO_EX]]


def fva_range(metas, rxns, cond, override, biomass_min, reactions):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    n = len(rxn_list)
    A_ub = np.zeros((1, n))
    A_ub[0, rxn_idx[lib.BIO_EX]] = -1.0
    b_ub = np.array([-biomass_min])
    out = {}
    for r in reactions:
        lo = hi = None
        for sign, tag in [(1.0, "lo"), (-1.0, "hi")]:
            c = np.zeros(n)
            c[rxn_idx[r]] = sign
            res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                          bounds=list(zip(lb, ub)), method="highs")
            if res.success:
                if tag == "lo":
                    lo = res.x[rxn_idx[r]]
                else:
                    hi = res.x[rxn_idx[r]]
        out[r] = (lo, hi)
    return out
