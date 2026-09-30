# -*- coding: utf-8 -*-
"""2016 published-model reproduction using the current v3 XML + full Table-1 semantics."""
from __future__ import annotations

import sys

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib
import xlrd

MMC1 = C.ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
COND_COL = {"FIM": ("lb_fe2", "ub_fe2"), "TTM": ("lb_ttton", "ub_ttton"),
            "TSM": ("lb_tsul", "ub_tsul")}
DONOR_RXN = {"FIM": "Ex_fe2[e]", "TTM": "Ex_ttton[e]", "TSM": "Ex_tsul[e]"}


def read_table1():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 1")
    rows = []
    for r in range(2, s.nrows):
        rid = str(s.cell_value(r, 0)).strip()
        if not rid:
            continue
        def cv(c):
            v = s.cell_value(r, c)
            if v == "":
                return ""
            if isinstance(v, float) and v.is_integer():
                return int(v)
            return v
        rows.append({"reaction_id": rid, "lb_fe2": cv(10), "ub_fe2": cv(11),
                     "lb_ttton": cv(12), "ub_ttton": cv(13),
                     "lb_tsul": cv(14), "ub_tsul": cv(15)})
    book.release_resources()
    return rows


def read_table5():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s5 = book.sheet_by_name("Table 5")
    points = []
    def isnum(v):
        return isinstance(v, (int, float)) and not isinstance(v, bool)
    for r in range(3, s5.nrows):
        def cv5(c):
            v = s5.cell_value(r, c)
            return v if v != "" else None
        fe2_mu, fe2_up, fe2_o2, fe2_co2 = cv5(0), cv5(1), cv5(2), cv5(3)
        s4_mu, s4_up, s4_co2 = cv5(5), cv5(6), cv5(7)
        s2_mu, s2_up, s2_co2 = cv5(9), cv5(10), cv5(11)
        if isnum(fe2_mu):
            points.append(dict(idx=len(points) + 1, condition="FIM", donor="fe2",
                               mu=fe2_mu, donor_uptake=fe2_up, o2=fe2_o2, co2=fe2_co2))
        if isnum(s4_mu):
            points.append(dict(idx=len(points) + 1, condition="TTM", donor="s4o6",
                               mu=s4_mu, donor_uptake=s4_up, o2=None, co2=s4_co2))
        if isnum(s2_mu):
            points.append(dict(idx=len(points) + 1, condition="TSM", donor="s2o3",
                               mu=s2_mu, donor_uptake=s2_up, o2=None, co2=s2_co2))
    book.release_resources()
    return points


def apply_table1_full(metas, rxns, cond, rows, biomass_keep_sbml=True):
    """Apply FULL Table-1 condition bounds (internal + exchange) to a fresh bound dict.

    Reproduces the Phase 2 benchmark semantics: the biomass reaction keeps its SBML
    0/1000 bound, all other reactions get the published Table-1 condition bounds.
    """
    lbc, ubc = COND_COL[cond]
    t1 = {r["reaction_id"]: r for r in rows}
    bounds = {}
    for r, rx in rxns.items():
        if r == lib.BIOMASS_RXN and biomass_keep_sbml:
            bounds[r] = (rx["lb"], rx["ub"])
            continue
        meta = t1.get(r)
        if meta is not None and meta[lbc] != "" and meta[ubc] != "":
            bounds[r] = (float(meta[lbc]), float(meta[ubc]))
        else:
            bounds[r] = (rx["lb"], rx["ub"])
    return bounds


def run_scenario(metas, rxns, cond, bounds, co2):
    """Maximize biomass with CO2 fixed at -co2, using given bounds. Returns flux dict."""
    met_list = sorted(metas)
    rxn_list = list(rxns)
    met_idx = {m: i for i, m in enumerate(met_list)}
    S = np.zeros((len(met_list), len(rxn_list)))
    for j, r in enumerate(rxn_list):
        for m, c in rxns[r]["stoich"].items():
            if m in met_idx:
                S[met_idx[m], j] = c
    lb = np.zeros(len(rxn_list)); ub = np.zeros(len(rxn_list))
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    for j, r in enumerate(rxn_list):
        bl, bu = bounds[r]
        lb[j], ub[j] = bl, bu
    # close glucose (not in 2016 scenarios)
    if "Ex_glc-B[e]" in rxn_idx:
        lb[rxn_idx["Ex_glc-B[e]"]] = 0.0
        ub[rxn_idx["Ex_glc-B[e]"]] = 0.0
    # fix CO2
    if "Ex_h2co3[e]" in rxn_idx:
        lb[rxn_idx["Ex_h2co3[e]"]] = -float(co2)
        ub[rxn_idx["Ex_h2co3[e]"]] = -float(co2)
    c = np.zeros(len(rxn_list))
    c[rxn_idx["Ex_bio[e]"]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None, res.message
    return {r: res.x[ri] for r, ri in rxn_idx.items()}, "optimal"


def pearson(x, y):
    x = np.array(x, float); y = np.array(y, float)
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def r2_1_minus_sse_sst(x, y):
    x = np.array(x, float); y = np.array(y, float)
    sse = float(((x - y) ** 2).sum())
    sst = float(((x - x.mean()) ** 2).sum())
    if sst == 0:
        return 0.0
    return 1.0 - sse / sst


def linfit(x, y):
    x = np.array(x, float); y = np.array(y, float)
    if len(x) < 2:
        return 0.0, 0.0
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def run_2016(metas, rxns):
    rows = read_table1()
    points = read_table5()
    scenarios = []
    for p in points:
        cond = p["condition"]
        bounds = apply_table1_full(metas, rxns, cond, rows)
        v, status = run_scenario(metas, rxns, cond, bounds, p["co2"])
        if v is None:
            scenarios.append(dict(idx=p["idx"], condition=cond, donor_type=p["donor"],
                                  exp_growth=p["mu"], growth=None, exp_donor=p["donor_uptake"],
                                  donor=None, exp_o2=p["o2"], o2=None, status=status))
            continue
        growth = v["Ex_bio[e]"]
        donor_flux = abs(v[DONOR_RXN[cond]])
        o2_flux = abs(v["Ex_o2[e]"]) if p["o2"] is not None else None
        scenarios.append(dict(idx=p["idx"], condition=cond, donor_type=p["donor"],
                              exp_growth=p["mu"], growth=growth, exp_donor=p["donor_uptake"],
                              donor=donor_flux, exp_o2=p["o2"], o2=o2_flux, status=status))
    return scenarios


def validation_series(scenarios):
    series = []
    def add(name, kind, pred_key, exp_key, cond_filter):
        sub = [s for s in scenarios if s["condition"] == cond_filter]
        if not sub:
            return
        exp = [s[exp_key] for s in sub]
        pred = [s[pred_key] for s in sub]
        slope, intercept = linfit(exp, pred)
        r2c = r2_1_minus_sse_sst(exp, pred)
        series.append(dict(series=name, kind=kind, n=len(sub), slope=slope,
                           r2_1_minus_sse_sst=r2c))
    add("growth_FIM", "growth", "growth", "exp_growth", "FIM")
    add("growth_TTM", "growth", "growth", "exp_growth", "TTM")
    add("growth_TSM", "growth", "growth", "exp_growth", "TSM")
    add("donor_FIM_fe2", "donor_uptake", "donor", "exp_donor", "FIM")
    add("donor_TTM_tetrathionate", "donor_uptake", "donor", "exp_donor", "TTM")
    add("donor_TSM_thiosulfate", "donor_uptake", "donor", "exp_donor", "TSM")
    add("o2_FIM", "o2_uptake", "o2", "exp_o2", "FIM")
    return series
