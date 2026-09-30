#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 3: physiological sanity checks for the frozen 2016 iMC507 baseline."""
import csv
import datetime as dt
import json
import re
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import flux_variability_analysis, pfba

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "qc" / "physiology"
LOGS = ROOT / "baseline_v0" / "logs"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
FE3 = "Ex_fe3_LSQBKT_e_RSQBKT_"
TTTON = "Ex_ttton_LSQBKT_e_RSQBKT_"
TSUL = "Ex_tsul_LSQBKT_e_RSQBKT_"
SO4 = "Ex_so4_LSQBKT_e_RSQBKT_"
H_EX = "Ex_h_LSQBKT_e_RSQBKT_"

PHENO_RXNS = ["ATPM", "ATPS5rpp", "RUBISCO", "RUBISCOX", "HCO3tpp", "HCO3E",
              "CYT1", "CYT2", "CYTA2", "CYTAA31", "CYTAA32", "CYTBC1", "CYTRED",
              "CYTBO3", "CYTBD", "NADHI", "TSQOC", "SQRED1", "SULDO", "SCCR",
              "4THASE1", "4THASE2"]

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

def read_carbon_bases():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 2")
    out = set()
    for r in range(1, s.nrows):
        abbr = s.cell_value(r, 0)
        formula = s.cell_value(r, 2)
        if abbr == "" or abbr is None:
            continue
        f = str(formula)
        if "C" in f:
            base = str(abbr).split("[")[0]
            out.add(base)
    book.release_resources()
    return out

def load_model():
    raw = FROZEN_SBML.read_bytes()
    return _sbml_to_model(libsbml.readSBMLFromString(raw.decode("utf-8", "replace")))

def apply_operational(model, cond, rows):
    ids = {r.id for r in model.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids or cid == "Afe_biomass_mc507_WT_139p0M":
            continue
        rx = model.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))

def set_objective(model):
    model.objective = model.reactions.get_by_id(BIO)

def donor_rxn(cond):
    return {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}[cond]

def cond_col(cond):
    return {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}[cond]

def max_sv(model, sol):
    S = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    order = [r.id for r in model.reactions]
    v = [float(sol.fluxes.get(rid, 0.0)) for rid in order]
    return float(np.max(np.abs(S @ v)))

def phenotype(model, rows, cond, h2co3=2.0, knockouts=(), close_ex=()):
    m = model.copy()
    apply_operational(m, cond_col(cond), rows)
    for rid in knockouts:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    for rid in close_ex:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    if h2co3 is None:
        m.reactions.get_by_id(H2CO3).bounds = (0.0, 0.0)
    else:
        m.reactions.get_by_id(H2CO3).bounds = (-float(h2co3), -float(h2co3))
    set_objective(m)
    sol = m.optimize()
    return m, sol

def record(m, sol, cond):
    dr = donor_rxn(cond)
    d = dict(condition=cond, status=sol.status, growth=sol.objective_value)
    d["h2co3"] = sol.fluxes[H2CO3]
    d["o2"] = sol.fluxes[O2]
    d["donor"] = sol.fluxes[dr]
    d["fe2"] = sol.fluxes[FE2]
    d["fe3"] = sol.fluxes[FE3]
    d["tsul"] = sol.fluxes[TSUL]
    d["ttton"] = sol.fluxes[TTTON]
    d["so4"] = sol.fluxes[SO4]
    d["h_ex"] = sol.fluxes[H_EX]
    for rid in PHENO_RXNS:
        d[rid] = sol.fluxes[rid] if rid in m.reactions else ""
    return d

def carbon_exchange_uptake(m, sol, carbon_bases, tol=1e-7):
    out = []
    for r in m.reactions:
        if not r.id.startswith("Ex_"):
            continue
        f = float(sol.fluxes[r.id])
        if f < -tol:
            met = next(iter(r.metabolites))
            base = met.id[:-2].replace("_DASH_", "-")
            is_carbon = base in carbon_bases
            out.append(dict(reaction=r.id, metabolite=base, flux=f, carbon_containing=is_carbon))
    return out

def no_input_atpm_max(m):
    c = m.copy()
    for r in c.reactions:
        if r.id.startswith("Ex_"):
            r.bounds = (0.0, 0.0)
    atpm = c.reactions.get_by_id("ATPM")
    atpm.bounds = (0.0, 1000.0)
    c.objective = atpm
    sol = c.optimize()
    return sol.objective_value

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = load_model()
    rows = read_table1()
    carbon_bases = read_carbon_bases()
    conds = ["FIM", "TTM", "TSM"]

    # A. positive control phenotypes
    pos_rows = []
    fva_rows = []
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0)
        pos_rows.append(dict(method="FBA", **record(m, sol, cond)))
        p = pfba(m)
        p_rec = record(m, p, cond)
        p_rec["growth"] = float(p.fluxes[BIO])
        pos_rows.append(dict(method="pFBA", **p_rec))
        # FVA for key reactions
        key = [donor_rxn(cond), O2, "RUBISCO", "ATPS5rpp", "Htpp"]
        fva = flux_variability_analysis(m, reaction_list=key, fraction_of_optimum=1.0)
        for rid in key:
            fva_rows.append(dict(condition=cond, reaction=rid, minimum=float(fva.loc[rid, "minimum"]), maximum=float(fva.loc[rid, "maximum"])))
    with (OUT / "positive_control_phenotypes.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["method", "condition", "status", "growth", "h2co3", "o2", "donor", "fe2", "fe3", "tsul", "ttton", "so4", "h_ex"] + PHENO_RXNS
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in pos_rows:
            w.writerow(d)
    with (OUT / "positive_control_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "minimum", "maximum"], delimiter="\t")
        w.writeheader()
        for d in fva_rows:
            w.writerow(d)

    # B. carbon dependency
    carb_rows = []
    for cond in conds:
        # normal
        m, sol = phenotype(model, rows, cond, h2co3=2.0)
        carb_rows.append(dict(condition=cond, h2co3="normal", growth=sol.objective_value, carbon_exchange_uptake=""))
        # disabled
        m2, sol2 = phenotype(model, rows, cond, h2co3=None)
        carb_rows.append(dict(condition=cond, h2co3="disabled", growth=sol2.objective_value, carbon_exchange_uptake=""))
        up = carbon_exchange_uptake(m2, sol2, carbon_bases)
        for u in up:
            carb_rows.append(dict(condition=cond, h2co3="disabled", growth="", carbon_exchange_uptake=f"{u['reaction']} {u['metabolite']} flux={u['flux']} carbon={u['carbon_containing']}"))
    with (OUT / "carbon_dependency.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "h2co3", "growth", "carbon_exchange_uptake"], delimiter="\t")
        w.writeheader()
        for d in carb_rows:
            w.writerow(d)

    # C. Calvin dependency
    cal_rows = []
    for cond in conds:
        for name, ko in [("baseline", ()), ("RUBISCO_disabled", ("RUBISCO", "RUBISCOX")), ("HCO3E_disabled", ("HCO3E",)), ("HCO3tpp_disabled", ("HCO3tpp",))]:
            m, sol = phenotype(model, rows, cond, h2co3=2.0, knockouts=ko)
            cal_rows.append(dict(condition=cond, knockout=name, growth=sol.objective_value, status=sol.status))
    with (OUT / "calvin_dependency.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "knockout", "status", "growth"], delimiter="\t")
        w.writeheader()
        for d in cal_rows:
            w.writerow(d)

    # D. electron donor dependency
    donor_rows = []
    for cond in conds:
        dr = donor_rxn(cond)
        m, sol = phenotype(model, rows, cond, h2co3=2.0, close_ex=(dr,))
        donor_rows.append(dict(condition=cond, scenario="donor_disabled", growth=sol.objective_value, note=""))
        # alternative donor check: nonzero reduced-donor exchange uptake
        alts = []
        for rid in [FE2, TTTON, TSUL, "Ex_h2_LSQBKT_e_RSQBKT_", "Ex_h2s_LSQBKT_e_RSQBKT_", "Ex_s_LSQBKT_e_RSQBKT_"]:
            f = float(sol.fluxes[rid])
            if f < -1e-7:
                alts.append(f"{rid}={f}")
        if alts:
            donor_rows[-1]["note"] = "alt donor uptake: " + ";".join(alts)
    # stricter: all reduced donors disabled
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0, close_ex=(FE2, TTTON, TSUL, "Ex_h2_LSQBKT_e_RSQBKT_", "Ex_h2s_LSQBKT_e_RSQBKT_", "Ex_s_LSQBKT_e_RSQBKT_"))
        donor_rows.append(dict(condition=cond, scenario="all_reduced_donors_disabled", growth=sol.objective_value, note=""))
    with (OUT / "electron_donor_dependency.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "scenario", "growth", "note"], delimiter="\t")
        w.writeheader()
        for d in donor_rows:
            w.writerow(d)

    # E. electron acceptor dependency
    acc_rows = []
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0, close_ex=(O2,))
        # find alternative acceptor uptake (fe3)
        fe3f = float(sol.fluxes[FE3])
        acc_rows.append(dict(condition=cond, scenario="O2_removed_only", growth=sol.objective_value, fe3_uptake=fe3f, note=("Fe3 imported" if fe3f < -1e-7 else "")))
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0, close_ex=(O2, FE3))
        acc_rows.append(dict(condition=cond, scenario="no_terminal_acceptor", growth=sol.objective_value, fe3_uptake=0.0, note=""))
    with (OUT / "electron_acceptor_dependency.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "scenario", "growth", "fe3_uptake", "note"], delimiter="\t")
        w.writeheader()
        for d in acc_rows:
            w.writerow(d)

    # F. redox product directions (FBA at h2co3=2)
    redox_rows = []
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0)
        dr = donor_rxn(cond)
        redox_rows.append(dict(condition=cond, growth=sol.objective_value, fe2=sol.fluxes[FE2], fe3=sol.fluxes[FE3], o2=sol.fluxes[O2], h2co3=sol.fluxes[H2CO3], donor=sol.fluxes[dr], tsul=sol.fluxes[TSUL], ttton=sol.fluxes[TTTON], so4=sol.fluxes[SO4]))
    with (OUT / "redox_product_directions.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "growth", "fe2", "fe3", "o2", "h2co3", "donor", "tsul", "ttton", "so4"], delimiter="\t")
        w.writeheader()
        for d in redox_rows:
            w.writerow(d)

    # G. energy sanity
    energy_rows = []
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0)
        htpp = float(sol.fluxes["Htpp"])
        free = no_input_atpm_max(m)
        energy_rows.append(dict(condition=cond, Htpp_flux=htpp, free_ATP_with_all_ex_closed=free, note=("ISSUE" if free > 1e-7 else "OK")))
    with (OUT / "energy_sanity.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "Htpp_flux", "free_ATP_with_all_ex_closed", "note"], delimiter="\t")
        w.writeheader()
        for d in energy_rows:
            w.writerow(d)

    # H. exchange leakage inventory (for positive-growth scenarios)
    leak_rows = []
    mineral = {"k", "na1", "mg2", "ca2", "mn2", "zn2", "cu2", "pi", "so4", "so4aa", "n2", "nh4", "polypi"}
    for cond in conds:
        m, sol = phenotype(model, rows, cond, h2co3=2.0)
        for r in m.reactions:
            if not r.id.startswith("Ex_"):
                continue
            f = float(sol.fluxes[r.id])
            if abs(f) < 1e-7:
                continue
            met = next(iter(r.metabolites))
            base = met.id[:-2].replace("_DASH_", "-")
            is_carbon = base in carbon_bases
            if base == "h2co3":
                grp = "intended_carbon"
            elif base in ("fe2", "ttton", "tsul"):
                grp = "intended_donor"
            elif base == "o2":
                grp = "intended_acceptor"
            elif base in mineral:
                grp = "mineral_nitrogen"
            elif is_carbon:
                grp = "unexpected_carbon"
            elif base in ("h2", "h2s", "s", "h"):
                grp = "unexpected_energy"
            else:
                grp = "unresolved"
            leak_rows.append(dict(condition=cond, reaction=r.id, metabolite=base, flux=f, direction=("uptake" if f < 0 else "export"), group=grp))
    with (OUT / "exchange_leakage_inventory.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "metabolite", "flux", "direction", "group"], delimiter="\t")
        w.writeheader()
        for d in leak_rows:
            w.writerow(d)

    # I. H2CO3 response
    resp_rows = []
    for cond in conds:
        for h in [0.0, 0.5, 1.0, 2.0, 4.0]:
            m, sol = phenotype(model, rows, cond, h2co3=h)
            dr = donor_rxn(cond)
            resp_rows.append(dict(condition=cond, h2co3_input=h, growth=sol.objective_value, donor_uptake=sol.fluxes[dr], o2_uptake=sol.fluxes[O2], status=sol.status))
    with (OUT / "h2co3_response.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "h2co3_input", "growth", "donor_uptake", "o2_uptake", "status"], delimiter="\t")
        w.writeheader()
        for d in resp_rows:
            w.writerow(d)

    summary = dict(generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"), positive_growth={cond: float([d for d in pos_rows if d["method"] == "FBA" and d["condition"] == cond][0]["growth"]) for cond in conds}, carbon_disabled_growth={cond: float([d for d in carb_rows if d["condition"] == cond and d["h2co3"] == "disabled"][0]["growth"]) for cond in conds}, rubisco_disabled_growth={cond: float([d for d in cal_rows if d["condition"] == cond and d["knockout"] == "RUBISCO_disabled"][0]["growth"]) for cond in conds}, donor_disabled_growth={cond: float([d for d in donor_rows if d["condition"] == cond and d["scenario"] == "donor_disabled"][0]["growth"]) for cond in conds}, no_acceptor_growth={cond: float([d for d in acc_rows if d["condition"] == cond and d["scenario"] == "no_terminal_acceptor"][0]["growth"]) for cond in conds}, free_atp={cond: float([d for d in energy_rows if d["condition"] == cond][0]["free_ATP_with_all_ex_closed"]) for cond in conds})
    (OUT / "phase3_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE3 done")
    print("positive growth", summary["positive_growth"])
    print("carbon disabled growth", summary["carbon_disabled_growth"])
    print("rubisco disabled growth", summary["rubisco_disabled_growth"])
    print("donor disabled growth", summary["donor_disabled_growth"])
    print("no acceptor growth", summary["no_acceptor_growth"])
    print("free atp", summary["free_atp"])

if __name__ == "__main__":
    main()
