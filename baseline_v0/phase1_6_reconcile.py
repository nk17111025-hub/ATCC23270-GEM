#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1.6: full Table 1 bound reconciliation + affected re-diagnosis."""
import csv
import datetime as dt
import json
import re
from collections import defaultdict
from pathlib import Path
import libsbml
import xlrd
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import flux_variability_analysis, find_blocked_reactions

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "qc"
SCEN = ROOT / "baseline_v0" / "scenario_definition"

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

def xml_bounds(model):
    return {r.id: (r.lower_bound, r.upper_bound) for r in model.reactions}

def apply_full_table1(model, cond, rows, skip=()):
    ids = {r.id for r in model.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids:
            continue
        if cid in skip:
            continue
        rx = model.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))

def set_objective_biomass(model):
    model.objective = model.reactions.get_by_id(BIO)

def max_sv_residual(model, sol):
    import numpy as np
    S = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    order = [r.id for r in model.reactions]
    v = [float(sol.fluxes.get(rid, 0.0)) for rid in order]
    return float(np.max(np.abs(S @ v)))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = load_model()
    rows = read_table1()
    points = read_table5()
    xb = xml_bounds(model)
    t1map = {r["reaction_id"]: r for r in rows}

    # ---- 1. A/B/C comparison ----
    conds = [("fe2", "FIM"), ("ttton", "TTM"), ("tsul", "TSM")]
    recon = []
    omitted_sets = {}
    for cond, label in conds:
        omitted = []
        for r in rows:
            cid = enc(r["reaction_id"])
            xlb, xub = xb.get(cid, (None, None))
            tlb, tub = r["lb_" + cond], r["ub_" + cond]
            differs = (xlb is None) or abs(float(tlb) - float(xlb)) > 1e-9 or abs(float(tub) - float(xub)) > 1e-9
            is_ex = cid.startswith("Ex_")
            was_applied = (is_ex or cid == "CYT2")
            if differs:
                recon.append(dict(condition=label, reaction_id=cid, table1_id=r["reaction_id"], subsystem=r["subsystem"], raw_sbml_lb=xlb, raw_sbml_ub=xub, table1_lb=tlb, table1_ub=tub, was_applied_in_phase1=was_applied, omitted=("YES" if not was_applied else "NO")))
                if not was_applied:
                    omitted.append(cid)
        omitted_sets[label] = sorted(set(omitted))
    with (OUT / "phase1_6_bound_reconciliation.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "reaction_id", "table1_id", "subsystem", "raw_sbml_lb", "raw_sbml_ub", "table1_lb", "table1_ub", "was_applied_in_phase1", "omitted"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in recon:
            w.writerow(d)

    # ---- 2. re-run with complete Table 1 bounds (except Afe_biomass = growth-relevant) ----
    rediag = []
    for cond, label in conds:
        # growth-relevant: full bounds except Afe_biomass (keep SBML 0/1000)
        g = model.copy()
        apply_full_table1(g, cond, rows, skip=("Afe_biomass_mc507_WT_139p0M",))
        set_objective_biomass(g)
        sg = g.optimize()
        fva = flux_variability_analysis(g, reaction_list=["Htpp", "ATPS5rpp"], fraction_of_optimum=1.0)
        # EGC free-energy test under published condition: close all exchanges + maximize ATPM
        eg = model.copy()
        apply_full_table1(eg, cond, rows, skip=("Afe_biomass_mc507_WT_139p0M",))
        for r in eg.reactions:
            if r.id.startswith("Ex_"):
                r.bounds = (0.0, 0.0)
        atpm = eg.reactions.get_by_id("ATPM")
        atpm.bounds = (0.0, 1000.0)
        eg.objective = atpm
        seg = eg.optimize()
        blocked = find_blocked_reactions(g, zero_cutoff=1e-7)
        rediag.append(dict(condition=label, biomass=sg.objective_value, max_sv=max_sv_residual(g, sg), Htpp_min=fva.loc["Htpp", "minimum"], Htpp_max=fva.loc["Htpp", "maximum"], ATPS5rpp_min=fva.loc["ATPS5rpp", "minimum"], ATPS5rpp_max=fva.loc["ATPS5rpp", "maximum"], EGC_atpm_max_with_Htpp_off=seg.objective_value, blocked_count=len(blocked)))
    with (OUT / "phase1_6_rediagnosis.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "biomass", "max_sv", "Htpp_min", "Htpp_max", "ATPS5rpp_min", "ATPS5rpp_max", "EGC_atpm_max_with_Htpp_off", "blocked_count"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in rediag:
            w.writerow(d)

    # literal full bounds (including Afe_biomass=0/0) -> demonstrate growth collapse
    literal = []
    for cond, label in conds:
        m = model.copy()
        apply_full_table1(m, cond, rows)
        set_objective_biomass(m)
        s = m.optimize()
        literal.append(dict(condition=label, status=s.status, biomass=s.objective_value))
    with (OUT / "phase1_6_literal_full_bounds.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "status", "biomass"], delimiter="\t")
        w.writeheader()
        for d in literal:
            w.writerow(d)

    # ---- 3. Table 5 check with complete bounds + point H2CO3 ----
    donor_rxn = {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}
    t5 = []
    for p in points:
        cond = p["condition"]
        condcol = "fe2" if cond == "FIM" else ("ttton" if cond == "TTM" else "tsul")
        m = model.copy()
        apply_full_table1(m, condcol, rows, skip=("Afe_biomass_mc507_WT_139p0M",))
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective_biomass(m)
        s = m.optimize()
        d = dict(condition=cond, donor=p["donor"], mu_reported=p["mu"], co2=p["co2"], growth=s.objective_value, donor_uptake=s.fluxes[donor_rxn[cond]])
        d["o2_uptake"] = (s.fluxes[O2] if p["o2"] is not None else "")
        d["Htpp_flux"] = s.fluxes["Htpp"]
        d["ATPS5rpp_flux"] = s.fluxes["ATPS5rpp"]
        t5.append(d)
    with (OUT / "phase1_6_table5_check.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "donor", "mu_reported", "co2", "growth", "donor_uptake", "o2_uptake", "Htpp_flux", "ATPS5rpp_flux"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in t5:
            w.writerow(d)

    # ---- 4. update scenario JSON to v3.0 ----
    old_path = SCEN / "2016_scenario_constraints.json"
    old = json.loads(old_path.read_text(encoding="utf-8"))
    overrides = {}
    for cond, label in conds:
        lst = []
        for r in rows:
            cid = enc(r["reaction_id"])
            xlb, xub = xb.get(cid, (None, None))
            tlb, tub = r["lb_" + cond], r["ub_" + cond]
            if xlb is None or abs(float(tlb) - float(xlb)) > 1e-9 or abs(float(tub) - float(xub)) > 1e-9:
                lst.append(dict(reaction_id=cid, table1_id=r["reaction_id"], raw_sbml_bounds=[xlb, xub], published_condition_bounds=[tlb, tub]))
        overrides[label] = lst
    new = dict(
        schema_version="3.0",
        title="2016 iMC507 scenario constraints (complete condition-specific bound overrides)",
        generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        note="Complete Table 1 (mmc1.xls cols 11-16) condition bounds, including internal transport/sulfur reactions. Afe_biomass published bound 0/0 is recorded but flagged as growth-blocking artifact (see special_cases).",
        source_files=old["source_files"],
        model_scope=old["model_scope"],
        raw_sbml_defaults=dict(note="mmc3.xml default bounds; see phase1_6_bound_reconciliation.tsv for every reaction where published condition bounds differ."),
        condition_bound_overrides=overrides,
        special_cases=[
            dict(reaction_id="Afe_biomass_mc507_WT_139p0M", issue="Table 1 cols 11-16 set 0/0 for all conditions, but raw SBML is 0/1000; applying 0/0 makes growth zero. Requires decision before strict reproduction."),
            dict(reaction_id="Htpp", issue="Table 1 sets 0/0 in all three conditions (published minimal media); this disables the reversible proton transport and blocks the ATP energy-generating cycle."),
        ],
        base_medium_definition=old.get("base_medium_definition"),
        experimental_validation_table5=old.get("experimental_validation_table5"),
        special_analysis_scenarios=old.get("special_analysis_scenarios"),
        energy_parameters=old.get("energy_parameters"),
    )
    json.dump(new, open(str(old_path), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    (old_path.with_suffix("").__str__())

    summary = dict(generated_at=dt.datetime.now().isoformat(timespec="seconds"), diff_counts={label: len([d for d in recon if d["condition"] == label]) for _, label in conds}, omitted=omitted_sets, rediagnosis=rediag, literal_full_bounds=literal, table5=t5)
    (OUT / "phase1_6_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE1.6 done")
    print("diff_counts", summary["diff_counts"])
    print("omitted", omitted_sets)
    print("rediagnosis", rediag)
    print("literal", literal)
    print("table5 sample", t5[:3])

if __name__ == "__main__":
    main()
