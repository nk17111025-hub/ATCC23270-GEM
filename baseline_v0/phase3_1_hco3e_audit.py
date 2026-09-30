#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 3.1: HCO3E bypass trace + physiology inventory cleanup."""
import csv
import datetime as dt
import json
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import pfba

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "qc" / "physiology"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
FE3 = "Ex_fe3_LSQBKT_e_RSQBKT_"
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

def read_carbon_bases():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 2")
    out = set()
    for r in range(1, s.nrows):
        abbr = s.cell_value(r, 0)
        formula = s.cell_value(r, 2)
        if abbr == "" or abbr is None:
            continue
        if "C" in str(formula):
            out.add(str(abbr).split("[")[0])
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

def coeff(r, mid):
    for k, v in r.metabolites.items():
        if k.id == mid:
            return v
    return None

def co2_reactions(model):
    prod, cons = [], []
    for r in model.reactions:
        c = coeff(r, "co2_c")
        if c is None:
            continue
        (prod if c > 0 else cons).append(r.id)
    return prod, cons

def hco3_consumers(model):
    out = []
    for r in model.reactions:
        c = coeff(r, "hco3_c")
        if c is not None and c < 0:
            out.append(r.id)
    return out

def run_pfba(model, rows, cond, hco3=2.0, hco3e_off=False, atpm=None, cap_donor=None, cap_o2=None):
    m = model.copy()
    apply_operational(m, cond_col(cond), rows)
    if hco3e_off:
        m.reactions.get_by_id("HCO3E").bounds = (0.0, 0.0)
    if hco3 is None:
        m.reactions.get_by_id(H2CO3).bounds = (0.0, 0.0)
    else:
        m.reactions.get_by_id(H2CO3).bounds = (-float(hco3), -float(hco3))
    if atpm is not None:
        m.reactions.get_by_id("ATPM").bounds = (float(atpm), float(atpm))
    dr = donor_rxn(cond)
    if cap_donor is not None:
        # cap uptake magnitude: lower bound = -cap, upper bound = +cap
        m.reactions.get_by_id(dr).bounds = (-float(cap_donor), 1000.0)
    if cap_o2 is not None:
        m.reactions.get_by_id(O2).bounds = (-float(cap_o2), 1000.0)
    set_objective(m)
    try:
        sol = pfba(m)
    except Exception:
        sol = None
    return m, sol

def total_abs_flux(sol):
    return float(sol.fluxes.abs().sum())

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = load_model()
    rows = read_table1()
    carbon_bases = read_carbon_bases()
    conds = ["FIM", "TTM", "TSM"]

    co2_prod, co2_cons = co2_reactions(model)
    hco3_cons = hco3_consumers(model)

    # 1. HCO3E bypass trace
    trace = []
    delta_rows = []
    for cond in conds:
        dr = donor_rxn(cond)
        m_n, s_n = run_pfba(model, rows, cond, hco3=2.0)
        m_off, s_off = run_pfba(model, rows, cond, hco3=2.0, hco3e_off=True)
        base = dict(condition=cond)
        rec = {}
        for key, get in [("growth", lambda s: s.fluxes[BIO]), ("HCO3tpp", lambda s: s.fluxes["HCO3tpp"]), ("HCO3E", lambda s: s.fluxes["HCO3E"]), ("RUBISCO", lambda s: s.fluxes["RUBISCO"]), ("RUBISCOX", lambda s: s.fluxes["RUBISCOX"]), ("ACCOAC", lambda s: s.fluxes["ACCOAC"]), ("BIOC2", lambda s: s.fluxes["BIOC2"]), ("CBPS", lambda s: s.fluxes["CBPS"]), ("donor", lambda s: s.fluxes[dr]), ("O2", lambda s: s.fluxes[O2]), ("ATPM", lambda s: s.fluxes["ATPM"]), ("ATPS5rpp", lambda s: s.fluxes["ATPS5rpp"]), ("total_abs_flux", total_abs_flux)]:
            rec[key + "_normal"] = get(s_n)
            rec[key + "_off"] = get(s_off)
            rec[key + "_delta"] = get(s_off) - get(s_n)
        rec["condition"] = cond
        trace.append(rec)
        # per-reaction carbon flux (normal vs off), only for carbon reactions with |delta|>1e-7
        carb_rxns = ["HCO3tpp", "HCO3E", "RUBISCO", "RUBISCOX", "ACCOAC", "BIOC2", "CBPS", "DBTS", "PPC"] + co2_prod
        for rid in sorted(set(carb_rxns)):
            if rid not in m_n.reactions:
                continue
            dn = s_n.fluxes[rid]
            do = s_off.fluxes[rid]
            if abs(do - dn) > 1e-7:
                delta_rows.append(dict(condition=cond, reaction=rid, normal_flux=dn, off_flux=do, delta=do - dn))
    with (OUT / "hco3e_flux_trace.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition"] + [k for k in trace[0].keys() if k != "condition"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in trace:
            w.writerow(d)
    with (OUT / "hco3e_flux_delta.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "normal_flux", "off_flux", "delta"], delimiter="\t")
        w.writeheader()
        for d in delta_rows:
            w.writerow(d)

    # 2. energy cost + capped tests
    cost_rows = []
    for cond in conds:
        dr = donor_rxn(cond)
        m_n, s_n = run_pfba(model, rows, cond, hco3=2.0)
        m_off, s_off = run_pfba(model, rows, cond, hco3=2.0, hco3e_off=True)
        g_n = s_n.fluxes[BIO]
        g_off = s_off.fluxes[BIO]
        donor_n = abs(s_n.fluxes[dr]); donor_off = abs(s_off.fluxes[dr])
        o2_n = abs(s_n.fluxes[O2]); o2_off = abs(s_off.fluxes[O2])
        row = dict(condition=cond, donor_per_biomass_normal=donor_n/g_n, donor_per_biomass_off=donor_off/g_off, o2_per_biomass_normal=o2_n/g_n, o2_per_biomass_off=o2_off/g_off, atps5rpp_normal=s_n.fluxes["ATPS5rpp"], atps5rpp_off=s_off.fluxes["ATPS5rpp"], total_abs_flux_normal=total_abs_flux(s_n), total_abs_flux_off=total_abs_flux(s_off))
        # capped tests
        for name, cap_d, cap_o in [("cap_donor", donor_n, None), ("cap_o2", None, o2_n), ("cap_both", donor_n, o2_n)]:
            mc, sc = run_pfba(model, rows, cond, hco3=2.0, hco3e_off=True, cap_donor=cap_d, cap_o2=cap_o)
            row[name + "_growth"] = (sc.fluxes[BIO] if sc is not None else "infeasible")
        cost_rows.append(row)
    with (OUT / "hco3e_energy_cost.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "donor_per_biomass_normal", "donor_per_biomass_off", "o2_per_biomass_normal", "o2_per_biomass_off", "atps5rpp_normal", "atps5rpp_off", "total_abs_flux_normal", "total_abs_flux_off", "cap_donor_growth", "cap_o2_growth", "cap_both_growth"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in cost_rows:
            w.writerow(d)

    # 3. zero-carbon maintenance
    maint_rows = []
    for cond in conds:
        dr = donor_rxn(cond)
        for atpm in [3.475, 0.0]:
            m, s = run_pfba(model, rows, cond, hco3=None, atpm=atpm)
            if s is None:
                maint_rows.append(dict(condition=cond, atpm=atpm, growth="infeasible", donor_uptake="", o2_uptake="", atps5rpp="", atpm_flux="", total_abs_flux=""))
            else:
                maint_rows.append(dict(condition=cond, atpm=atpm, growth=s.fluxes[BIO], donor_uptake=s.fluxes[dr], o2_uptake=s.fluxes[O2], atps5rpp=s.fluxes["ATPS5rpp"], atpm_flux=s.fluxes["ATPM"], total_abs_flux=total_abs_flux(s)))
    with (OUT / "zero_carbon_maintenance.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "atpm", "growth", "donor_uptake", "o2_uptake", "atps5rpp", "atpm_flux", "total_abs_flux"], delimiter="\t")
        w.writeheader()
        for d in maint_rows:
            w.writerow(d)

    # 4. corrected leakage classification
    biomass = model.reactions.get_by_id("Afe_biomass_mc507_WT_139p0M")
    fe2_coef = biomass.get_coefficient("fe2_c")
    fe3_coef = biomass.get_coefficient("fe3_c")
    mineral = {"k", "na1", "mg2", "ca2", "mn2", "zn2", "cu2", "pi", "so4", "so4aa", "n2", "nh4", "polypi"}
    leak = []
    for cond in conds:
        m, s = run_pfba(model, rows, cond, hco3=2.0)
        for r in m.reactions:
            if not r.id.startswith("Ex_"):
                continue
            f = float(s.fluxes[r.id])
            if abs(f) < 1e-7:
                continue
            met = next(iter(r.metabolites))
            base = met.id[:-2].replace("_DASH_", "-")
            direction = "uptake" if f < 0 else "export"
            is_carbon = base in carbon_bases
            # classify
            if direction == "export" and base == "bio":
                grp = "BIOMASS_EXPORT"
            elif direction == "export" and base == "eps_AFE":
                grp = "BIOSYNTHETIC_BYPRODUCT_EXPORT"
            elif base == "h2co3":
                grp = "INTENDED_CARBON"
            elif base == "o2":
                grp = "INTENDED_ACCEPTOR"
            elif base == "fe2" and cond == "FIM":
                grp = "INTENDED_DONOR"
            elif base == "ttton" and cond == "TTM":
                grp = "INTENDED_DONOR"
            elif base == "tsul" and cond == "TSM":
                grp = "INTENDED_DONOR"
            elif base in ("fe2", "fe3") and cond != "FIM":
                grp = "BIOMASS_MINERAL_REQUIREMENT"
            elif base == "fe3" and cond == "FIM":
                grp = "REDOX_PRODUCT_EXPORT"
            elif base == "h":
                grp = "ACID_BASE_BALANCE"
            elif base == "4hba" and direction == "export":
                grp = "BIOSYNTHETIC_BYPRODUCT_EXPORT"
            elif base in mineral:
                grp = "MINERAL_NITROGEN"
            elif base in ("so4",) and direction == "export":
                grp = "REDOX_PRODUCT_EXPORT"
            elif direction == "export":
                grp = "BIOSYNTHETIC_BYPRODUCT_EXPORT" if is_carbon else "UNRESOLVED_EXPORT"
            elif is_carbon:
                grp = "UNEXPECTED_CARBON_UPTAKE"
            else:
                grp = "UNRESOLVED"
            leak.append(dict(condition=cond, reaction=r.id, metabolite=base, flux=f, direction=direction, carbon=is_carbon, group=grp))
    with (OUT / "exchange_leakage_inventory_corrected.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "metabolite", "flux", "direction", "carbon", "group"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in leak:
            w.writerow(d)

    summary = dict(generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"), biomass_fe2_coefficient=fe2_coef, biomass_fe3_coefficient=fe3_coef, hco3_consumers=hco3_cons, co2_producers=co2_prod, co2_consumers=co2_cons, trace=trace, energy_cost=cost_rows, maintenance=maint_rows, unexpected_carbon_uptake=[d for d in leak if d["group"] == "UNEXPECTED_CARBON_UPTAKE"])
    (OUT / "phase3_1_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE3.1 done")
    print("biomass fe2/fe3 coef", fe2_coef, fe3_coef)
    print("hco3 consumers", hco3_cons)
    print("co2 consumers", co2_cons)
    for cond in conds:
        t = [x for x in trace if x["condition"] == cond][0]
        print(cond, "growth_off", round(t["growth_off"], 6), "ACCOAC_off", round(t["ACCOAC_off"], 4), "BIOC2_off", round(t["BIOC2_off"], 4), "CBPS_off", round(t["CBPS_off"], 4), "RUBISCO_off", round(t["RUBISCO_off"], 4), "total_abs_off", round(t["total_abs_flux_off"], 2))
    print("unexpected carbon uptake", summary["unexpected_carbon_uptake"])

if __name__ == "__main__":
    main()
