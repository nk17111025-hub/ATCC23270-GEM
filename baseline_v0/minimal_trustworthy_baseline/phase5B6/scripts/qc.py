# -*- coding: utf-8 -*-
"""Regression + artifact QC for v3. INFEASIBLE is NEVER counted as PASS for zero tests."""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
import v3_model as M
import conditions as C
import fw_common as FW

WT_REF = 0.0520764
TOL = 1e-7


def zero_test_result(status, value, tol=TOL):
    """Hard rule: only OPTIMAL with |value| <= tol is PASS. INFEASIBLE/None -> FAIL."""
    if status != "OPTIMAL":
        return "FAIL"
    if value is None:
        return "FAIL"
    return "PASS" if abs(value) <= tol else "FAIL"


def find_atp_maintenance(rxns):
    """Return the reaction id of the NGAM/ATP maintenance reaction, or None."""
    for r, rx in rxns.items():
        s = rx["stoich"]
        if s.get("atp[c]", 0) < 0 and s.get("h2o[c]", 0) < 0 and s.get("adp[c]", 0) > 0:
            if rx["lb"] > 0:
                return r
    # fallback: known id
    if "ATPM" in rxns:
        return "ATPM"
    return None


def closed_zero_background(metas, rxns, cond="FIM", disable_maintenance=True):
    """Build (metas, rxns, exchange_overrides, override) for a CLOSED zero-demand background.

    Closes all carbon/donor uptake, sets biomass export lb to 0, and (for loop tests)
    disables the ATP maintenance demand so the zero-flux state is feasible.
    """
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0),
          "Ex_fe2[e]": (0.0, 0.0), "Ex_ttton[e]": (0.0, 0.0), "Ex_tsul[e]": (0.0, 0.0),
          "Ex_s[e]": (0.0, 0.0), "Ex_h2s[e]": (0.0, 0.0), "Ex_bio[e]": (0.0, 1000.0)}
    overrides = {}
    if disable_maintenance:
        atm = find_atp_maintenance(rxns)
        if atm:
            overrides[atm] = (0.0, rxns[atm]["ub"])
    return ov, overrides


def closed_zero_background_feasible(metas, rxns, cond="FIM"):
    """Return (status, objective). Must be OPTIMAL (zero-flux feasible)."""
    ov, overrides = closed_zero_background(metas, rxns, cond)
    S, ml, rl, ri, lb, ub = C.apply_condition(metas, rxns, cond, exchange_overrides=ov, overrides=overrides)
    c = np.zeros(len(rl))
    status, res = FW.solve_lp(c, S, np.zeros(len(ml)), list(zip(lb, ub)))
    return status, (res.fun if res is not None else None)


def free_loop_test(metas, rxns, cond, label, demand_met=None, maximize_rxn=None, pmf=False):
    """Run one free-loop diagnostic on the closed zero-demand background.

    Returns (test, solver_status, maximum_flux, tolerance, result).
    """
    ov, overrides = closed_zero_background(metas, rxns, cond)
    rx2 = dict(rxns)
    sid = "DIAG_" + label
    if demand_met is not None:
        stoich = {demand_met: -1.0}
    elif maximize_rxn is not None:
        stoich = None  # maximize an existing reaction
    else:
        stoich = {}
    if stoich is not None:
        rx2[sid] = {"name": sid, "reversible": False, "lb": 0.0, "ub": 1000.0,
                    "stoich": stoich, "objective_coefficient": 0.0}
    S, ml, rl, ri, lb, ub = C.apply_condition(metas, rx2, cond, exchange_overrides=ov, overrides=overrides)
    c = np.zeros(len(rl))
    if maximize_rxn is not None and maximize_rxn in ri:
        c[ri[maximize_rxn]] = -1.0
    elif stoich is not None:
        c[ri[sid]] = -1.0
    status, res = FW.solve_lp(c, S, np.zeros(len(ml)), list(zip(lb, ub)))
    if status != "OPTIMAL":
        val = None
    elif maximize_rxn is not None and maximize_rxn in ri:
        val = res.x[ri[maximize_rxn]]
    elif stoich is not None:
        val = res.x[ri[sid]]
    else:
        val = 0.0
    return [label, status, val, TOL, zero_test_result(status, val)]


def free_loop_tests(metas, rxns, cond="FIM"):
    # PMF generation is produced by the proton-pumping ETC (e.g. terminal oxidase CYTBD),
    # which requires an electron donor and is therefore zero in the closed background.
    # A synthetic "proton export" would measure gross proton balance cycling, not PMF
    # production, so test 6 maximizes the native proton-pumping ETC instead.
    rows = [
        free_loop_test(metas, rxns, cond, "ATP", demand_met="atp[c]"),
        free_loop_test(metas, rxns, cond, "NADH", demand_met="nadh[c]"),
        free_loop_test(metas, rxns, cond, "NADPH", demand_met="nadph[c]"),
        free_loop_test(metas, rxns, cond, "q8h2_quinol", demand_met="q8h2[c]"),
        free_loop_test(metas, rxns, cond, "q8_quinone", demand_met="q8[c]"),
        free_loop_test(metas, rxns, cond, "PMF_generation", maximize_rxn="CYTBD"),
        free_loop_test(metas, rxns, cond, "ATP_synthase_loop", maximize_rxn="ATPS5rpp"),
    ]
    return rows


def phenotype_tests(metas, rxns, cond="FIM"):
    """Phenotype tests are separate from free-loop tests; report exact solver status."""
    rows = []
    # no-carbon growth
    st, v, _ = C.fba_status(metas, rxns, cond, exchange_overrides={
        "Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
    rows.append(["no_carbon_growth", st, (v[M.BIOMASS_EX] if v else None)])
    # no-donor growth
    st2, v2, _ = C.fba_status(metas, rxns, cond, exchange_overrides={
        "Ex_fe2[e]": (0.0, 0.0), "Ex_ttton[e]": (0.0, 0.0), "Ex_tsul[e]": (0.0, 0.0)})
    rows.append(["no_donor_growth", st2, (v2[M.BIOMASS_EX] if v2 else None)])
    return rows


def basic_smoke_regression(metas, rxns):
    """Four-row basic smoke panel (FIM/TTM/TSM WT + glucose-closed). NOT the 2016 panel."""
    rows = []
    for cnd in ["FIM", "TTM", "TSM"]:
        ref = C.wt_reference(metas, rxns, cnd)
        if ref is None:
            rows.append([f"{cnd}_WT_biomass", "INFEASIBLE", "N/A", str(WT_REF), "FAIL"])
            continue
        b = ref["biomass"]
        ok = abs(b - WT_REF) < 1e-4
        rows.append([f"{cnd}_WT_biomass", "OPTIMAL", f"{b:.6g}", str(WT_REF),
                     "PASS" if ok else "FAIL"])
    st, v, _ = C.fba_status(metas, rxns, "FIM", exchange_overrides={
        "Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    rows.append(["FIM_glucose_closed", st, f"{v[M.BIOMASS_EX]:.6g}" if v else "N/A",
                 str(WT_REF), "PASS" if v and abs(v[M.BIOMASS_EX] - WT_REF) < 1e-4 else "FAIL"])
    return rows
