#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression / QC for minimal_trustworthy_baseline (post HCO3E patch)."""
import csv
import datetime as dt
import json
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import pfba, flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
WORK = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "working" / "minimal_trustworthy_baseline_working.xml"
MMC1 = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
TTTON = "Ex_ttton_LSQBKT_e_RSQBKT_"
TSUL = "Ex_tsul_LSQBKT_e_RSQBKT_"

def enc(s):
    return s.replace("[", "_LSQBKT_").replace("]", "_RSQBKT_").replace("-", "_DASH_")

def read_table1():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 1")
    rows = []
    for r in range(2, s.nrows):
        rid = s.cell_value(r, 0)
        if rid == "" or rid is None:
            continue
        def cv(c):
            v = s.cell_value(r, c)
            if v == "":
                return ""
            if isinstance(v, float) and v.is_integer():
                return int(v)
            return v
        rows.append(dict(reaction_id=str(rid), subsystem=str(cv(6)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
    book.release_resources()
    return rows

def load_working():
    raw = WORK.read_bytes()
    return _sbml_to_model(libsbml.readSBMLFromString(raw.decode("utf-8", "replace")))

def apply_operational_and_patch(m, cond, rows):
    ids = {r.id for r in m.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids or cid == "Afe_biomass_mc507_WT_139p0M":
            continue
        rx = m.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))
    # approved patch on top of operational bounds
    m.reactions.get_by_id("MACPD").bounds = (0.0, 0.0)
    m.reactions.get_by_id("ACOATA").bounds = (0.0, 1000.0)

def set_objective(m):
    m.objective = m.reactions.get_by_id(BIO)

def donor_rxn(cond):
    return {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}[cond]

def cond_col(cond):
    return {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}[cond]

def run(m, cond, rows, hco3=2.0, hco3e_off=False, donor_cap=None, o2_cap=None):
    mm = m.copy()
    apply_operational_and_patch(mm, cond_col(cond), rows)
    if hco3e_off:
        mm.reactions.get_by_id("HCO3E").bounds = (0.0, 0.0)
    if hco3 is None:
        mm.reactions.get_by_id(H2CO3).bounds = (0.0, 0.0)
    else:
        mm.reactions.get_by_id(H2CO3).bounds = (-float(hco3), -float(hco3))
    dr = donor_rxn(cond)
    if donor_cap is not None:
        mm.reactions.get_by_id(dr).bounds = (-float(donor_cap), 1000.0)
    if o2_cap is not None:
        mm.reactions.get_by_id(O2).bounds = (-float(o2_cap), 1000.0)
    set_objective(mm)
    try:
        sol = pfba(mm)
    except Exception:
        sol = None
    return mm, sol

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = load_working()
    rows = read_table1()
    conds = ["FIM", "TTM", "TSM"]

    results = {}
    growth_rows = []
    for cond in conds:
        # WT
        mm, sol = run(model, cond, rows)
        wt = float(sol.fluxes[BIO]) if sol is not None else None
        # HCO3E KO
        _, sol_ko = run(model, cond, rows, hco3e_off=True)
        ko = float(sol_ko.fluxes[BIO]) if sol_ko is not None else "infeasible"
        # bypass check at HCO3E KO: ACCOAC, MACPD, ACOATA
        ac = float(sol_ko.fluxes["ACCOAC"]) if sol_ko is not None else "n/a"
        macpd = float(sol_ko.fluxes["MACPD"]) if sol_ko is not None else "n/a"
        acoata = float(sol_ko.fluxes["ACOATA"]) if sol_ko is not None else "n/a"
        growth_rows.append(dict(condition=cond, WT_growth=wt, HCO3E_KO_growth=ko, HCO3E_KO_ACCOAC=ac, HCO3E_KO_MACPD=macpd, HCO3E_KO_ACOATA=acoata))
        results[cond] = dict(WT_growth=wt, HCO3E_KO_growth=ko)
        # FVA on WT for key reactions
        if sol is not None:
            fva = flux_variability_analysis(mm, reaction_list=["MACPD", "ACOATA", "RUBISCO", "Htpp"], fraction_of_optimum=1.0)
            for rid in ["MACPD", "ACOATA", "RUBISCO", "Htpp"]:
                results[cond]["FVA_" + rid] = [float(fva.loc[rid, "minimum"]), float(fva.loc[rid, "maximum"])]

    # EGC test (published bounds -> Htpp=0; and raw Htpp for comparison)
    egc = {}
    for cond in conds:
        mm, _ = run(model, cond, rows)
        # close all exchanges, maximize ATPM
        for r in mm.reactions:
            if r.id.startswith("Ex_"):
                r.bounds = (0.0, 0.0)
        atpm = mm.reactions.get_by_id("ATPM")
        atpm.bounds = (0.0, 1000.0)
        mm.objective = atpm
        sol = mm.optimize()
        egc[cond + "_published_htpp0"] = sol.objective_value
        # raw Htpp reversible test (temp)
        mm2 = mm.copy()
        mm2.reactions.get_by_id("Htpp").bounds = (-1000.0, 1000.0)
        sol2 = mm2.optimize()
        egc[cond + "_raw_htpp_reversible"] = sol2.objective_value

    # zero-carbon and zero-donor
    zerocarb = {}
    zerodonor = {}
    for cond in conds:
        _, s1 = run(model, cond, rows, hco3=None)
        zerocarb[cond] = float(s1.fluxes[BIO]) if s1 is not None else "infeasible"
        dr = donor_rxn(cond)
        mm, _ = run(model, cond, rows)
        mm.reactions.get_by_id(dr).bounds = (0.0, 0.0)
        set_objective(mm)
        try:
            s2 = pfba(mm)
            zerodonor[cond] = float(s2.fluxes[BIO])
        except Exception:
            zerodonor[cond] = "infeasible"

    summary = dict(
        generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        growth=growth_rows,
        egc=egc,
        zero_carbon=zerocarb,
        zero_donor=zerodonor,
        results=results,
    )
    (OUT / "working" / "regression_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    with (OUT / "working" / "regression_growth.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(growth_rows[0].keys()), delimiter="\t")
        w.writeheader()
        for d in growth_rows:
            w.writerow(d)

    print("REGRESSION done")
    for g in growth_rows:
        print(g)
    print("EGC", egc)
    print("zero_carbon", zerocarb)
    print("zero_donor", zerodonor)
    for cond in conds:
        print(cond, "FVA", results[cond].get("FVA_MACPD"), results[cond].get("FVA_ACOATA"), results[cond].get("FVA_Htpp"))

if __name__ == "__main__":
    main()
