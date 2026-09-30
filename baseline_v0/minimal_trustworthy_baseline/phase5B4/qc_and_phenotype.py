# -*- coding: utf-8 -*-
"""Phase 5B-4: QC + native phenotype + precursor reachability + Rubisco dependence."""
import csv
import json
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B4")
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}


def build_v2(t):
    metas, rxns = build_memory(t)
    rxns["NADHI"]["lb"] = -1000.0  # open respiratory NADH-oxidizing direction
    return metas, rxns


def setup(metas, rxns, cond, override):
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n); ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]; lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN: continue
        rx = rxns[r]
        lv = rx["table1_"+lbc]; uv = rx["table1_"+ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv); ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx: lb[rxn_idx[r]] = pl; ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx: lb[rxn_idx[r]] = ol; ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def max_biomass(metas, rxns, cond, override):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def free_energy(metas, rxns, cond, drain):
    rx2 = dict(rxns)
    ov = {r: (0.0, 0.0) for r in rx2 if r.startswith("Ex_")}
    ov.update({"MACPD": (0.0, 0.0), "ACOATA": (0.0, 1000.0), "Htpp": (0.0, 0.0)})
    if drain == "ATPM":
        ov["ATPM"] = (0.0, 1000.0)
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective="ATPM", override=ov, with_table1=False)
    else:
        d = "DM_drain"
        rx2[d] = {"sbml_id": d, "name": "drain", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {drain: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                  "gpr": "", "gpr2": "", "protein": "", **BLANK}
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective=d, override=ov, with_table1=False)
    return 0.0 if fl is None else float(obj)


def ref_state(metas, rxns, cond):
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, ov)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    v = res.x
    return {"biomass": v[rxn_idx[lib.BIO_EX]], "donor": v[rxn_idx[lib.DONOR_EX[cond]]],
            "o2": v[rxn_idx["Ex_o2[e]"]], "co2": v[rxn_idx["Ex_h2co3[e]"]],
            "rubisco": v[rxn_idx["RUBISCO"]] + v[rxn_idx["RUBISCOX"]]}


def demand(metas, rxns, cond, met, override):
    rx2 = dict(rxns)
    rx2["DM_demand"] = {"sbml_id": "DM_demand", "name": "demand", "reversible": False, "lb": 0.0,
                        "ub": 1000.0, "stoich": {met: -1.0}, "confidence": "", "ec": "", "pmid": "",
                        "subsystem": "hypo", "gpr": "", "gpr2": "", "protein": "", **BLANK}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rx2, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx["DM_demand"]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx["DM_demand"]]


CENTRAL = ["g6p-B[c]", "f6p-B[c]", "r5p[c]", "ru5p-D[c]", "e4p[c]", "g3p[c]", "dhap[c]", "13dpg[c]",
           "3pg[c]", "2pg[c]", "pep[c]", "pyr[c]", "accoa[c]", "cit[c]", "icit[c]", "akg[c]",
           "succoa[c]", "succ[c]", "fum[c]", "mal-L[c]", "oaa[c]", "ser-L[c]", "gly[c]", "ala-L[c]",
           "asp-L[c]", "glu-L[c]", "gln-L[c]"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns = build_v2("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns, cond)
    print("ref biomass", round(ref["biomass"], 7), "donor", round(ref["donor"], 4), "o2", round(ref["o2"], 4))

    # QC: free energy + baseline regression
    qc = []
    for c in ["FIM", "TTM", "TSM"]:
        qc.append({"condition": c, "test": "free_ATP", "value": free_energy(metas, rxns, c, "ATPM")})
        qc.append({"condition": c, "test": "free_NADH", "value": free_energy(metas, rxns, c, "nadh[c]")})
        qc.append({"condition": c, "test": "free_NADPH", "value": free_energy(metas, rxns, c, "nadph[c]")})
        b = max_biomass(metas, rxns, c, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
        qc.append({"condition": c, "test": "baseline_glucose_closed_biomass", "value": round(b, 7) if b is not None else "infeasible"})
        b2 = max_biomass(metas, rxns, c, {"Ex_glc-B[e]": (0.0, 0.0), lib.DONOR_EX[c]: (0.0, 0.0)})
        qc.append({"condition": c, "test": "no_donor_biomass", "value": round(b2, 7) if b2 is not None else "infeasible"})
        b3 = max_biomass(metas, rxns, c, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
        qc.append({"condition": c, "test": "no_carbon_biomass", "value": round(b3, 7) if b3 is not None else "infeasible"})
    with (OUT / "phase5B4_artifact_tests.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "test", "value"], delimiter="\t")
        w.writeheader()
        for r in qc:
            w.writerow(r)
    with (OUT / "phase5B4_baseline_regression.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "test", "value"], delimiter="\t")
        w.writeheader()
        for r in qc:
            w.writerow(r)
    print("=== QC ===")
    for r in qc:
        if "biomass" in r["test"] or "free" in r["test"]:
            print(" ", r["condition"], r["test"], r["value"])

    # native glucose phenotype (Phase 5B-2C constraints): donor<=WT, O2<=WT, CO2 free, glucose
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0), donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, ov)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    bio = res.x[rxn_idx[lib.BIO_EX]]
    # pFBA clean
    n = len(rxn_list)
    c2 = np.zeros(2*n); c2[n:] = 1.0
    A_eq = np.zeros((len(met_list), 2*n)); A_eq[:, :n] = S
    A_ub = np.zeros((2*n+1, 2*n)); b_ub = np.zeros(2*n+1)
    for i in range(n):
        A_ub[i, i]=1; A_ub[i, n+i]=-1; A_ub[n+i, i]=-1; A_ub[n+i, n+i]=-1
    A_ub[2*n, rxn_idx[lib.BIO_EX]] = -1; b_ub[2*n] = -(bio - 1e-6)
    res2 = linprog(c2, A_eq=A_eq, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                   bounds=list(zip(lb, ub))+[(0,None)]*n, method="highs")
    v = res2.x[:n]
    print("=== native glucose phenotype (v2, FIM) ===")
    print("max biomass", round(bio, 6), "vs WT", round(ref["biomass"], 6))
    flux_keys = ["Ex_glc-B[e]", "Ex_h2co3[e]", donor, "Ex_o2[e]", "BDGK", "PGI1", "PFK", "FBA",
                 "GAPD1", "GAPD2", "PGK", "PGM1", "ENO", "PYK", "PDH", "G6PDH2", "PGL", "PGDH", "DDGPA",
                 "TKT1", "TKT2", "TALA", "RPI", "RPE", "PRUK", "RUBISCO", "RUBISCOX", "PPC", "CS",
                 "ICDHyr", "MDH", "NADHI", "NADTRHD", "ATPS5rpp", "CYTAA31", "CYT2"]
    flux_rows = []
    for r in flux_keys:
        if r in rxn_idx:
            flux_rows.append({"reaction": r, "flux": round(v[rxn_idx[r]], 6)})
    with (OUT / "phase5B4_native_flux_map.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction", "flux"], delimiter="\t")
        w.writeheader()
        for r in flux_rows:
            w.writerow(r)
    for r in flux_rows:
        if abs(r["flux"]) > 1e-4:
            print(f"  {r['reaction']:12s} {r['flux']:+.5f}")

    # precursor reachability (Rubisco=0, glucose, relaxed)
    reach = []
    ov0 = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0), donor: (-1000.0, 0.0),
           "Ex_o2[e]": (-1000.0, 0.0), "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}
    for met in CENTRAL:
        dv = demand(metas, rxns, cond, met, ov0)
        reach.append({"metabolite": met, "name": metas[met]["name"],
                      "producible": "yes" if dv is not None and dv > 1e-6 else "no",
                      "max_flux": round(dv, 6) if dv is not None else "infeasible"})
    with (OUT / "phase5B4_precursor_reachability.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["metabolite", "name", "producible", "max_flux"], delimiter="\t")
        w.writeheader()
        for r in reach:
            w.writerow(r)
    blocked = [r["metabolite"] for r in reach if r["producible"] == "no"]
    prod = [r["metabolite"] for r in reach if r["producible"] == "yes"]
    print("=== precursor reachability (Rubisco=0, glucose) ===")
    print("producible:", prod)
    print("blocked:", blocked)

    # Rubisco dependence
    rub = []
    for eps in [None, 1e-6, 1e-7, 1e-8, 0.0]:
        ovr = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0), donor: (ref["donor"], 0.0),
               "Ex_o2[e]": (ref["o2"], 0.0)}
        if eps is None:
            ovr["RUBISCO"] = (0.0, 1000.0); ovr["RUBISCOX"] = (0.0, 1000.0)
        else:
            ovr["RUBISCO"] = (0.0, eps); ovr["RUBISCOX"] = (0.0, eps)
        b = max_biomass(metas, rxns, cond, ovr)
        rub.append({"rubisco_cap": "unrestricted" if eps is None else eps,
                    "biomass": round(b, 7) if b is not None else "infeasible"})
    with (OUT / "phase5B4_rubisco_dependency.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rubisco_cap", "biomass"], delimiter="\t")
        w.writeheader()
        for r in rub:
            w.writerow(r)
    print("=== Rubisco dependence (FIM, donor/O2 <= WT) ===")
    for r in rub:
        print(" ", r["rubisco_cap"], "->", r["biomass"])


if __name__ == "__main__":
    main()
