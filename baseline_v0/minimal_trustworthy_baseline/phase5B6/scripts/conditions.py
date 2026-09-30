# -*- coding: utf-8 -*-
"""Single authoritative condition API for v3."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib  # historical metadata only (Table-1 exchange bounds)
import v3_model as M
import fw_common as FW

GLUCOSE_MIN_UPTAKE = 1e-4
GLUCOSE_LB = -5.0
GLUCOSE_UB = -GLUCOSE_MIN_UPTAKE

COND_EXCHANGE_COLS = {"FIM": ("lb_fe2", "ub_fe2"), "TTM": ("lb_ttton", "ub_ttton"),
                      "TSM": ("lb_tsul", "ub_tsul")}
DONOR_EX = {"FIM": "Ex_fe2[e]", "TTM": "Ex_ttton[e]", "TSM": "Ex_tsul[e]"}


def load_table1_exchanges():
    """Return {reaction: {cond: (lb, ub)}} for EXCHANGE reactions only."""
    t1 = lib._read_table1()
    out = {}
    for r, meta in t1.items():
        if not r.startswith("Ex_"):
            continue
        row = {}
        for cond, (lbc, ubc) in COND_EXCHANGE_COLS.items():
            lv = meta.get(lbc, "")
            uv = meta.get(ubc, "")
            if lv != "" and uv != "":
                row[cond] = (float(lv), float(uv))
        if row:
            out[r] = row
    return out


def build_matrix(metas, rxns):
    met_list = sorted(metas)
    rxn_list = list(rxns)
    met_idx = {m: i for i, m in enumerate(met_list)}
    S = np.zeros((len(met_list), len(rxn_list)))
    for j, r in enumerate(rxn_list):
        for m, c in rxns[r]["stoich"].items():
            if m in met_idx:
                S[met_idx[m], j] = c
    return S, met_list, rxn_list, met_idx


def apply_condition(metas, rxns, cond, overrides=None, exchange_overrides=None):
    """Base bounds from XML + Table-1 EXCHANGE-only condition bounds + overrides."""
    S, met_list, rxn_list, met_idx = build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]
    # environmental exchanges from Table-1 (only exchange reactions)
    for r, condmap in load_table1_exchanges().items():
        if r in rxn_idx and cond in condmap:
            el, eu = condmap[cond]
            lb[rxn_idx[r]] = el
            ub[rxn_idx[r]] = eu
    if exchange_overrides:
        for r, (ol, ou) in exchange_overrides.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    if overrides:
        for r, (ol, ou) in overrides.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def fba_status(metas, rxns, cond, overrides=None, exchange_overrides=None, objective=None):
    S, ml, rl, ri, lb, ub = apply_condition(metas, rxns, cond, overrides, exchange_overrides)
    obj = objective or M.BIOMASS_EX
    c = np.zeros(len(rl))
    if obj in ri:
        c[ri[obj]] = -1.0
    status, res = FW.solve_lp(c, S, np.zeros(len(ml)), list(zip(lb, ub)))
    if status != "OPTIMAL":
        return status, None, (S, ml, rl, ri, lb, ub)
    return "OPTIMAL", {r: res.x[ri[r]] for r in rl}, (S, ml, rl, ri, lb, ub)


def pfba_status(metas, rxns, cond, overrides=None, exchange_overrides=None, objective=None):
    S, ml, rl, ri, lb, ub = apply_condition(metas, rxns, cond, overrides, exchange_overrides)
    obj = objective or M.BIOMASS_EX
    n = len(rl)
    c = np.zeros(n)
    if obj in ri:
        c[ri[obj]] = -1.0
    status, res = FW.solve_lp(c, S, np.zeros(len(ml)), list(zip(lb, ub)))
    if status != "OPTIMAL":
        return status, None
    opt = -res.fun
    eps = 1e-6 * (1 + abs(opt))
    c2 = np.zeros(2 * n)
    c2[n:] = 1.0
    A_eq = np.zeros((len(ml), 2 * n))
    A_eq[:, :n] = S
    A_ub = np.zeros((2 * n + 1, 2 * n))
    b_ub = np.zeros(2 * n + 1)
    for i in range(n):
        A_ub[i, i] = 1.0
        A_ub[i, n + i] = -1.0
        A_ub[n + i, i] = -1.0
        A_ub[n + i, n + i] = -1.0
    A_ub[2 * n, :n] = c
    b_ub[2 * n] = -(opt - eps)
    bnd = list(zip(lb, ub)) + [(0.0, None)] * n
    st2, res2 = FW.solve_lp(c2, A_eq, np.zeros(len(ml)), bnd, A_ub=A_ub, b_ub=b_ub)
    if st2 != "OPTIMAL":
        return st2, None
    return "OPTIMAL", {r: res2.x[ri[r]] for r in rl}


def wt_reference(metas, rxns, cond):
    """Deterministic WT reference (glucose closed, CO2 fixed), via parsimonious flux.

    pFBA (maximize biomass, then minimize L1 flux) gives a deterministic donor/O2 cap
    instead of an arbitrary LP alternate optimum.
    """
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    status, v = pfba_status(metas, rxns, cond, overrides=ov)
    if status != "OPTIMAL":
        return None
    biomass = v[M.BIOMASS_EX]
    donor = DONOR_EX[cond]
    return {"biomass": biomass, "donor": v[donor], "o2": v["Ex_o2[e]"],
            "co2": v["Ex_h2co3[e]"], "flux": v}


def wt_donor_o2_fva(metas, rxns, cond):
    """FVA of donor/O2 at fixed WT biomass (deterministic reference check)."""
    S, ml, rl, ri, lb, ub = apply_condition(metas, rxns, cond,
                                            {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    c = np.zeros(len(rl)); c[ri[M.BIOMASS_EX]] = -1.0
    st, res = FW.solve_lp(c, S, np.zeros(len(ml)), list(zip(lb, ub)))
    if st != "OPTIMAL":
        return None
    opt = res.x[ri[M.BIOMASS_EX]]
    A_ub = np.zeros((1, len(rl))); A_ub[0, ri[M.BIOMASS_EX]] = -1.0
    b_ub = np.array([-0.999 * opt])
    out = {}
    for r in [DONOR_EX[cond], "Ex_o2[e]"]:
        j = ri[r]
        vals = []
        for sgn in (1.0, -1.0):
            cc = np.zeros(len(rl)); cc[j] = sgn
            rr = linprog(cc, A_eq=S, b_eq=np.zeros(len(ml)), A_ub=A_ub, b_ub=b_ub,
                         bounds=list(zip(lb, ub)), method="highs")
            if rr.success:
                vals.append(rr.x[j])
        out[r] = (min(vals), max(vals)) if vals else (None, None)
    return out


def r0_overrides(cond, ref, co2_bidirectional=False):
    """Return the authoritative R0 override dict (no hand-built dictionaries elsewhere)."""
    donor = DONOR_EX[cond]
    co2 = (-1000.0, 1000.0) if co2_bidirectional else (0.0, 1000.0)
    return {
        "Ex_glc-B[e]": (GLUCOSE_LB, GLUCOSE_UB),
        "Ex_h2co3[e]": co2,
        donor: (ref["donor"], 0.0),
        "Ex_o2[e]": (ref["o2"], 0.0),
        "RUBISCO": (0.0, 0.0),
        "RUBISCOX": (0.0, 0.0),
    }


def r0_strict(metas, rxns, cond="FIM"):
    ref = wt_reference(metas, rxns, cond)
    if ref is None:
        return None, None
    return ref, r0_overrides(cond, ref, co2_bidirectional=False)


def r0_co2_bidir(metas, rxns, cond="FIM"):
    ref = wt_reference(metas, rxns, cond)
    if ref is None:
        return None, None
    return ref, r0_overrides(cond, ref, co2_bidirectional=True)
