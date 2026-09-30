#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1.5: ATP energy-generating cycle (EGC) impact diagnostic."""
import csv
import datetime as dt
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
import libsbml
import xlrd
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "qc"
LOGS = ROOT / "baseline_v0" / "logs"

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
        rows.append(dict(reaction_id=str(rid), gra=str(cv(7)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
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
            points.append(dict(condition="FIM", donor="fe2", mu=fe2_mu, donor_uptake=fe2_up, o2=fe2_o2, co2=fe2_co2))
        if isnum(s4_mu):
            points.append(dict(condition="TTM", donor="s4o6", mu=s4_mu, donor_uptake=s4_up, o2=None, co2=s4_co2))
        if isnum(s2_mu):
            points.append(dict(condition="TSM", donor="s2o3", mu=s2_mu, donor_uptake=s2_up, o2=None, co2=s2_co2))
    book.release_resources()
    return points

def load_model():
    raw = FROZEN_SBML.read_bytes()
    doc = libsbml.readSBMLFromString(raw.decode("utf-8", "replace"))
    return _sbml_to_model(doc)

def apply_base_medium(model, cond, rows):
    cobra_ids = {r.id for r in model.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in cobra_ids:
            continue
        if not (cid.startswith("Ex_") or cid == "CYT2"):
            continue
        rx = model.reactions.get_by_id(cid)
        rx.lower_bound = float(r["lb_" + cond])
        rx.upper_bound = float(r["ub_" + cond])

def close_all_exchanges(model):
    for r in model.reactions:
        if r.id.startswith("Ex_"):
            r.lower_bound = 0.0
            r.upper_bound = 0.0

def set_objective_biomass(model):
    model.objective = model.reactions.get_by_id(BIO)

def set_max_atpm(model):
    atpm = model.reactions.get_by_id("ATPM")
    atpm.lower_bound = 0.0
    atpm.upper_bound = 1000.0
    model.objective = atpm
    return atpm

def nonzero_fluxes(model, sol, tol=1e-7):
    out = []
    for r in model.reactions:
        f = sol.fluxes[r.id]
        if abs(f) > tol:
            out.append((r.id, f))
    return out

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    rows = read_table1()
    points = read_table5()
    model = load_model()

    # ---- Section 1: minimal cycle + interruption ----
    base = model.copy()
    close_all_exchanges(base)
    set_max_atpm(base)
    sol = base.optimize()
    nz = nonzero_fluxes(base, sol)
    min_cycle = [dict(reaction_id=rid, flux=f) for rid, f in sorted(nz, key=lambda x: -abs(x[1]))]
    with (OUT / "phase1_5_minimal_egc_flux.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "flux"], delimiter="\t")
        w.writeheader()
        for d in min_cycle:
            w.writerow(d)

    interruption = []
    for block in ["Htpp", "ATPS5rpp"]:
        c = model.copy()
        close_all_exchanges(c)
        c.reactions.get_by_id(block).lower_bound = 0.0
        c.reactions.get_by_id(block).upper_bound = 0.0
        set_max_atpm(c)
        s = c.optimize()
        interruption.append(dict(blocked_reaction=block, status=s.status, atpm_max=s.objective_value))
    with (OUT / "phase1_5_egc_interruption.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["blocked_reaction", "status", "atpm_max"], delimiter="\t")
        w.writeheader()
        for d in interruption:
            w.writerow(d)

    # ---- Section 2: Htpp direction sensitivity ----
    direction = []
    for label, lb, ub in [("original", -1000.0, 1000.0), ("c_to_p_only", 0.0, 1000.0), ("p_to_c_only", -1000.0, 0.0), ("disabled", 0.0, 0.0)]:
        c = model.copy()
        close_all_exchanges(c)
        c.reactions.get_by_id("Htpp").lower_bound = lb
        c.reactions.get_by_id("Htpp").upper_bound = ub
        set_max_atpm(c)
        s = c.optimize()
        direction.append(dict(case=label, atpm_max=s.objective_value, atps5rpp_flux=s.fluxes["ATPS5rpp"], htpp_flux=s.fluxes["Htpp"], status=s.status))
    with (OUT / "phase1_5_htpp_direction.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["case", "status", "atpm_max", "atps5rpp_flux", "htpp_flux"], delimiter="\t")
        w.writeheader()
        for d in direction:
            w.writerow(d)

    # ---- Section 3: EGC usage at max growth + FVA ----
    cond_map = {"fe2": "FIM", "ttton": "TTM", "tsul": "TSM"}
    donor_rxn = {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}
    growth_rows = []
    fva_rows = []
    for cond, label in [("fe2", "FIM"), ("ttton", "TTM"), ("tsul", "TSM")]:
        c = model.copy()
        apply_base_medium(c, cond, rows)
        set_objective_biomass(c)
        s = c.optimize()
        d = dict(condition=label, biomass=s.objective_value, ATPM=s.fluxes["ATPM"], ATPS5rpp=s.fluxes["ATPS5rpp"], Htpp=s.fluxes["Htpp"])
        d["donor_uptake"] = s.fluxes[donor_rxn[label]]
        d["O2_uptake"] = s.fluxes[O2]
        d["H2CO3_uptake"] = s.fluxes[H2CO3]
        for rid in ["CYT2", "CYT1", "CYTA2", "CYTAA31", "CYTBC1", "NADHI", "SULDO", "SCCR", "TSQOC", "SQRED1", "4THASE1", "4THASE2"]:
            d[rid] = s.fluxes[rid]
        growth_rows.append(d)
        # FVA at 100% max biomass
        fva = flux_variability_analysis(c, reaction_list=["Htpp", "ATPS5rpp"], fraction_of_optimum=1.0)
        for rid in ["Htpp", "ATPS5rpp"]:
            row = fva.loc[rid]
            fva_rows.append(dict(condition=label, reaction=rid, minimum=row["minimum"], maximum=row["maximum"]))
    with (OUT / "phase1_5_baseline_growth_egc.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = list(growth_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for d in growth_rows:
            w.writerow(d)
    with (OUT / "phase1_5_egc_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "minimum", "maximum"], delimiter="\t")
        w.writeheader()
        for d in fva_rows:
            w.writerow(d)

    # ---- Section 4: biomass impact of Htpp / ATPS5rpp restriction ----
    impact = []
    for cond, label in [("fe2", "FIM"), ("ttton", "TTM"), ("tsul", "TSM")]:
        orig = model.copy()
        apply_base_medium(orig, cond, rows)
        set_objective_biomass(orig)
        so = orig.optimize()
        base_obj = so.objective_value
        for variant, lb, ub in [("Htpp_disabled", 0.0, 0.0), ("Htpp_c_to_p", 0.0, 1000.0), ("Htpp_p_to_c", -1000.0, 0.0), ("ATPS5rpp_disabled", 0.0, 0.0)]:
            c = model.copy()
            apply_base_medium(c, cond, rows)
            if variant.startswith("Htpp"):
                c.reactions.get_by_id("Htpp").lower_bound = lb
                c.reactions.get_by_id("Htpp").upper_bound = ub
            else:
                c.reactions.get_by_id("ATPS5rpp").lower_bound = 0.0
                c.reactions.get_by_id("ATPS5rpp").upper_bound = 0.0
            set_objective_biomass(c)
            sv = c.optimize()
            impact.append(dict(condition=label, variant=variant, biomass=sv.objective_value, abs_diff=sv.objective_value - base_obj, rel_diff=((sv.objective_value - base_obj) / base_obj if base_obj else 0.0), status=sv.status))
    with (OUT / "phase1_5_biomass_impact.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "variant", "status", "biomass", "abs_diff", "rel_diff"], delimiter="\t")
        w.writeheader()
        for d in impact:
            w.writerow(d)

    # ---- Section 5: Table 5 EGC pre-check ----
    t5 = []
    for p in points:
        cond = p["condition"]
        condcol = "fe2" if cond == "FIM" else ("ttton" if cond == "TTM" else "tsul")
        # original
        o = model.copy()
        apply_base_medium(o, condcol, rows)
        o.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective_biomass(o)
        so = o.optimize()
        # Htpp disabled
        h = model.copy()
        apply_base_medium(h, condcol, rows)
        h.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        h.reactions.get_by_id("Htpp").lower_bound = 0.0
        h.reactions.get_by_id("Htpp").upper_bound = 0.0
        set_objective_biomass(h)
        sh = h.optimize()
        drow = dict(condition=cond, donor=p["donor"], mu_reported=p["mu"], co2_input=p["co2"])
        drow["growth_orig"] = so.objective_value
        drow["growth_htpp_off"] = sh.objective_value
        drow["growth_diff"] = sh.objective_value - so.objective_value
        drow["donor_uptake_orig"] = so.fluxes[donor_rxn[cond]]
        drow["donor_uptake_htpp_off"] = sh.fluxes[donor_rxn[cond]]
        if p["o2"] is not None:
            drow["o2_uptake_orig"] = so.fluxes[O2]
            drow["o2_uptake_htpp_off"] = sh.fluxes[O2]
        else:
            drow["o2_uptake_orig"] = ""
            drow["o2_uptake_htpp_off"] = ""
        t5.append(drow)
    with (OUT / "phase1_5_table5_egc_precheck.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "donor", "mu_reported", "co2_input", "growth_orig", "growth_htpp_off", "growth_diff", "donor_uptake_orig", "donor_uptake_htpp_off", "o2_uptake_orig", "o2_uptake_htpp_off"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in t5:
            w.writerow(d)

    # ---- Section 6: GPR 508 vs 507 ----
    genes = Counter()
    gene_rxn = defaultdict(list)
    for r in rows:
        for g in re.findall(r"AFE_\d+", r["gra"]):
            genes[g] += 1
            gene_rxn[g].append(r["reaction_id"])
    gpr_rows = []
    for g in sorted(genes):
        gpr_rows.append(dict(gene_id=g, occurrence=genes[g], reactions=",".join(gene_rxn[g])))
    with (OUT / "phase1_5_gpr_count_check.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_id", "occurrence", "reactions"], delimiter="\t")
        w.writeheader()
        for d in gpr_rows:
            w.writerow(d)

    summary = dict(generated_at=dt.datetime.now().isoformat(timespec="seconds"), minimal_egc=min_cycle, interruption=interruption, htpp_direction=direction, baseline_growth=growth_rows, fva=fva_rows, biomass_impact=impact, table5=t5, unique_gpr_genes=len(genes))
    (OUT / "phase1_5_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE1.5 done")
    print("minimal EGC flux (nonzero reactions):", len(min_cycle), [m["reaction_id"] for m in min_cycle][:10])
    print("interruption:", interruption)
    print("htpp direction:", direction)
    print("baseline growth ATPM/ATPS5rpp/Htpp:", [(g["condition"], g["biomass"], g["ATPM"], g["ATPS5rpp"], g["Htpp"]) for g in growth_rows])
    print("fva:", fva_rows)
    print("unique GPR genes:", len(genes))

if __name__ == "__main__":
    main()
