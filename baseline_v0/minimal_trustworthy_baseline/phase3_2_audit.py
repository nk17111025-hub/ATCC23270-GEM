#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 3.2 central-carbon structural sanity audit on minimal_trustworthy_baseline_v1."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import pfba, flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
V1 = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"
MMC1 = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase3_2_central_carbon_audit"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
TTTON = "Ex_ttton_LSQBKT_e_RSQBKT_"
TSUL = "Ex_tsul_LSQBKT_e_RSQBKT_"

CENTRAL = [
    # EMP / gluconeogenesis
    "PGI1", "PFK", "FBP", "FBA", "FBA3", "TPI", "GAPD1", "GAPD2", "PGK", "PGM1", "PGMT", "PGMT2", "ENO", "PYK", "PDH", "BDGK", "MTRP",
    # PPP
    "G6PDH2", "PGL", "PGDH", "RPE", "RPI", "TKT1", "TKT2", "TALA", "DDGPA",
    # CBB
    "RUBISCO", "RUBISCOX", "PRUK", "SBPASE", "GLYCH", "HCO3E", "HCO3tpp",
    # TCA
    "CS", "ACONT1", "ACONT2", "ICDHyr", "ICITRED", "SUCOAS", "SUCD", "FUM", "MDH",
    # anaplerosis
    "PPC", "PPA", "MALS",
    # pyruvate / acetyl-CoA / OAA
    "ACS", "FDH", "ACCOAC", "MCOATA", "MACPD", "ACOATA", "BIOC1", "BIOC2",
    # glyoxylate
    "GLXCL", "GLYCK", "GLXR",
    # redox
    "NADTRHD",
    # 3OAS condensation (acetyl-CoA/malonyl-CoA interface)
    "3OAS1", "3OAS27", "3OAS28",
    # PEP/PYR carboxylase family
    "PPC", "ME1", "ME2", "PPS",
]

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
        rows.append(dict(reaction_id=str(rid), name=str(cv(1)), equation=str(cv(2)), confidence=cv(3), ec=str(cv(4)), pmid=cv(5), subsystem=str(cv(6)), gra=str(cv(7)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
    book.release_resources()
    return rows

def load_v1():
    raw = V1.read_bytes()
    return _sbml_to_model(libsbml.readSBMLFromString(raw.decode("utf-8", "replace")))

def apply_cond(m, cond, rows):
    ids = {r.id for r in m.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids or cid == "Afe_biomass_mc507_WT_139p0M":
            continue
        rx = m.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))
    m.reactions.get_by_id("MACPD").bounds = (0.0, 0.0)
    m.reactions.get_by_id("ACOATA").bounds = (0.0, 1000.0)
    m.reactions.get_by_id("Htpp").bounds = (0.0, 0.0)

def cond_col(cond):
    return {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}[cond]

def donor_rxn(cond):
    return {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}[cond]

def set_objective(m):
    m.objective = m.reactions.get_by_id(BIO)

def solve(m, h2co3=2.0):
    m.reactions.get_by_id(H2CO3).bounds = (-float(h2co3), -float(h2co3)) if h2co3 is not None else (0.0, 0.0)
    set_objective(m)
    try:
        sol = pfba(m)
        return sol
    except Exception:
        return None

def present(rxns):
    return [r for r in rxns if r in globals().get("_model", None) and False] if False else [r for r in rxns]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v1_hash = hashlib.sha256(V1.read_bytes()).hexdigest()
    model = load_v1()
    rows = read_table1()
    t1 = {r["reaction_id"]: r for r in rows}

    # only keep central reactions present in model
    central = [r for r in CENTRAL if r in model.reactions]
    # add ME/PPS if present
    for extra in ["ME1", "ME2", "PPS", "ME"]:
        if extra in model.reactions and extra not in central:
            central.append(extra)
    # deduplicate preserving order
    seen = set()
    central = [r for r in central if not (r in seen or seen.add(r))]

    # inventory
    inv = []
    for rid in central:
        r = model.reactions.get_by_id(rid)
        info = t1.get(rid, {})
        inv.append(dict(reaction_id=rid, name=r.name, equation=r.build_reaction_string(), reversible=r.reversibility, lb=r.lower_bound, ub=r.upper_bound, confidence=info.get("confidence", ""), ec=info.get("ec", ""), pmid=info.get("pmid", ""), gpr=info.get("gra", ""), subsystem=info.get("subsystem", "")))
    with (OUT / "central_carbon_reaction_inventory.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "name", "equation", "reversible", "lb", "ub", "confidence", "ec", "pmid", "gpr", "subsystem"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in inv:
            w.writerow(d)

    conds = ["FIM", "TTM", "TSM"]

    # reference growth + FVA at 100/95/90
    fva_ref = []
    for cond in conds:
        m = model.copy()
        apply_cond(m, cond_col(cond), rows)
        sol = solve(m)
        ref_g = float(sol.fluxes[BIO]) if sol is not None else None
        for frac in [1.0, 0.95, 0.90]:
            fva = flux_variability_analysis(m, reaction_list=central, fraction_of_optimum=frac)
            for rid in central:
                lo = float(fva.loc[rid, "minimum"]); hi = float(fva.loc[rid, "maximum"])
                fva_ref.append(dict(condition=cond, fraction=frac, reaction=rid, minimum=lo, maximum=hi, reference_flux=(float(sol.fluxes[rid]) if sol is not None else "")))
    with (OUT / "central_carbon_fva_reference.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "fraction", "reaction", "minimum", "maximum", "reference_flux"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in fva_ref:
            w.writerow(d)

    # Test B: single-reaction KO rescue screen (FIM representative)
    ko_rxns = ["HCO3E", "PGI1", "PFK", "FBP", "GAPD1", "PYK", "PPC", "ACCOAC", "MCOATA", "G6PDH2", "TKT1", "RUBISCO", "MDH", "FUM", "CS", "PDH", "ENO", "TPI", "PGK", "RPE", "RPI", "TALA", "SBPASE", "PRUK", "ACS", "NADTRHD"]
    ko_rows = []
    for cond in ["FIM"]:
        m = model.copy()
        apply_cond(m, cond_col(cond), rows)
        base = solve(m)
        base_g = float(base.fluxes[BIO]) if base is not None else 0.0
        base_flux = {rid: float(base.fluxes[rid]) for rid in central} if base is not None else {}
        for rid in ko_rxns:
            if rid not in m.reactions:
                continue
            mm = m.copy()
            mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            s = solve(mm)
            g = float(s.fluxes[BIO]) if s is not None else "infeasible"
            # top flux changes
            diffs = []
            if s is not None:
                for r2 in central:
                    d = float(s.fluxes[r2]) - base_flux.get(r2, 0.0)
                    if abs(d) > 1e-4:
                        diffs.append((r2, round(d, 4)))
                diffs.sort(key=lambda x: -abs(x[1]))
            ko_rows.append(dict(condition=cond, knocked_out=rid, reference_growth=base_g, growth_after_ko=g, top_flux_changes=";".join([f"{a}:{b}" for a, b in diffs[:8]])))
    with (OUT / "single_reaction_rescue_screen.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "knocked_out", "reference_growth", "growth_after_ko", "top_flux_changes"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in ko_rows:
            w.writerow(d)

    # Test C: reversibility watchlist
    rev = []
    for rid in central:
        r = model.reactions.get_by_id(rid)
        if r.reversibility and r.lower_bound < 0:
            info = t1.get(rid, {})
            rev.append(dict(reaction_id=rid, equation=r.build_reaction_string(), lb=r.lower_bound, ub=r.upper_bound, confidence=info.get("confidence", ""), gpr=info.get("gra", ""), ec=info.get("ec", ""), pmid=info.get("pmid", ""), subsystem=info.get("subsystem", "")))
    with (OUT / "reversibility_watchlist.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "equation", "lb", "ub", "confidence", "gpr", "ec", "pmid", "subsystem"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in rev:
            w.writerow(d)

    # Test D: CO2/HCO3 inventory (from model, all reactions)
    co2_prod, co2_cons, hco3_prod, hco3_cons = [], [], [], []
    for r in model.reactions:
        for mid, sign in [("co2_c", 1), ("co2_c", -1), ("hco3_c", 1), ("hco3_c", -1)]:
            c = None
            for k, v in r.metabolites.items():
                if k.id == mid:
                    c = v
            if c is None:
                continue
            if mid == "co2_c" and c > 0:
                co2_prod.append(r.id)
            elif mid == "co2_c" and c < 0:
                co2_cons.append(r.id)
            elif mid == "hco3_c" and c > 0:
                hco3_prod.append(r.id)
            elif mid == "hco3_c" and c < 0:
                hco3_cons.append(r.id)
            break
    with (OUT / "co2_hco3_cycle_candidates.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["role", "reaction_ids"], delimiter="\t")
        w.writeheader()
        w.writerow(dict(role="co2_producers", reaction_ids=";".join(sorted(set(co2_prod)))))
        w.writerow(dict(role="co2_consumers", reaction_ids=";".join(sorted(set(co2_cons)))))
        w.writerow(dict(role="hco3_producers", reaction_ids=";".join(sorted(set(hco3_prod)))))
        w.writerow(dict(role="hco3_consumers", reaction_ids=";".join(sorted(set(hco3_cons)))))

    # Test E: free redox under closed exchanges (isolate internal)
    redox = {}
    for drain_met, label in [("ATPM", "free_ATP"), ("nadh_c", "free_NADH"), ("nadph_c", "free_NADPH")]:
        mm = model.copy()
        apply_cond(mm, "fe2", rows)
        for r in mm.reactions:
            if r.id.startswith("Ex_"):
                r.bounds = (0.0, 0.0)
        mm.reactions.get_by_id("MACPD").bounds = (0.0, 0.0)
        mm.reactions.get_by_id("ACOATA").bounds = (0.0, 1000.0)
        mm.reactions.get_by_id("Htpp").bounds = (0.0, 0.0)
        if drain_met == "ATPM":
            rx = mm.reactions.get_by_id("ATPM"); rx.bounds = (0.0, 1000.0); mm.objective = rx
        else:
            dr = cobra.Reaction("TMP_drain"); mm.add_reactions([dr]); dr.add_metabolites({drain_met: -1.0}); dr.bounds = (0.0, 1000.0); mm.objective = dr
        sol = mm.optimize()
        redox[label] = sol.objective_value
    with (OUT / "redox_cycle_candidates.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["test", "objective"], delimiter="\t")
        w.writeheader()
        for k, v in redox.items():
            w.writerow(dict(test=k, objective=v))

    # Test F: pyruvate/acetyl-CoA/OAA closure audit (FIM KOs)
    py_ko = ["PDH", "PYK", "PPC", "ACS", "MDH", "FUM", "CS", "ACCOAC", "MCOATA", "MALS", "MMCD", "MMM", "GLXCL", "GLYCK"]
    py_rows = []
    m = model.copy()
    apply_cond(m, "fe2", rows)
    base = solve(m)
    base_g = float(base.fluxes[BIO]) if base is not None else 0.0
    base_flux = {rid: float(base.fluxes[rid]) for rid in central} if base is not None else {}
    for rid in py_ko:
        if rid not in m.reactions:
            continue
        mm = m.copy()
        mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        s = solve(mm)
        g = float(s.fluxes[BIO]) if s is not None else "infeasible"
        diffs = []
        if s is not None:
            for r2 in central:
                d = float(s.fluxes[r2]) - base_flux.get(r2, 0.0)
                if abs(d) > 1e-4:
                    diffs.append((r2, round(d, 4)))
            diffs.sort(key=lambda x: -abs(x[1]))
        py_rows.append(dict(knocked_out=rid, reference_growth=base_g, growth_after_ko=g, top_flux_changes=";".join([f"{a}:{b}" for a, b in diffs[:10]])))
    with (OUT / "pyruvate_accoa_oaa_route_audit.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["knocked_out", "reference_growth", "growth_after_ko", "top_flux_changes"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in py_rows:
            w.writerow(d)

    # Test H: low-confidence central reactions
    low = []
    for rid in central:
        info = t1.get(rid, {})
        conf = info.get("confidence", "")
        gpr = info.get("gra", "")
        pm = info.get("pmid", "")
        if (conf == 1 or conf == "1" or conf == 1.0) or (gpr == "" and not pm):
            r = model.reactions.get_by_id(rid)
            low.append(dict(reaction_id=rid, confidence=conf, gpr=gpr, pmid=pm, ec=info.get("ec", ""), equation=r.build_reaction_string(), reversible=r.reversibility, subsystem=info.get("subsystem", "")))
    with (OUT / "low_confidence_central_reactions.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "confidence", "gpr", "pmid", "ec", "equation", "reversible", "subsystem"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in low:
            w.writerow(d)

    # FVA perturbations (targeted after selected KOs) -> summary of dormant reverse capacity
    fva_pert = []
    for cond in ["FIM"]:
        m = model.copy()
        apply_cond(m, "fe2", rows)
        for ko in ["HCO3E", "PPC", "PDH", "MDH"]:
            mm = m.copy()
            mm.reactions.get_by_id(ko).bounds = (0.0, 0.0)
            s = solve(mm)
            if s is None:
                continue
            fva = flux_variability_analysis(mm, reaction_list=central, fraction_of_optimum=1.0)
            for rid in central:
                lo = float(fva.loc[rid, "minimum"]); hi = float(fva.loc[rid, "maximum"])
                if lo < -1e-4 or hi > 1e-4:
                    fva_pert.append(dict(condition=cond, perturbation=ko, reaction=rid, minimum=lo, maximum=hi, reference_flux=float(s.fluxes[rid])))
    with (OUT / "central_carbon_fva_perturbations.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "perturbation", "reaction", "minimum", "maximum", "reference_flux"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in fva_pert:
            w.writerow(d)

    # Test A: closed-carbon internal cycle screen
    closed_carbon = []
    for cond in conds:
        m = model.copy()
        apply_cond(m, cond_col(cond), rows)
        # close all exchanges
        for r in m.reactions:
            if r.id.startswith("Ex_"):
                r.bounds = (0.0, 0.0)
        m.reactions.get_by_id("ATPM").bounds = (0.0, 0.0)  # remove maintenance so cycle detection is pure
        set_objective(m)
        try:
            sol = m.optimize()
        except Exception:
            closed_carbon.append(dict(condition=cond, biomass=0.0, central_cyclic_flux="infeasible"))
            continue
        # FVA on central reactions with biomass fixed 0
        fva = flux_variability_analysis(m, reaction_list=central, fraction_of_optimum=1.0)
        cyclic = []
        for rid in central:
            lo = float(fva.loc[rid, "minimum"]); hi = float(fva.loc[rid, "maximum"])
            if abs(lo) > 1e-6 or abs(hi) > 1e-6:
                cyclic.append(f"{rid}[{lo:.3g},{hi:.3g}]")
        closed_carbon.append(dict(condition=cond, biomass=float(sol.fluxes[BIO]), central_cyclic_flux=";".join(cyclic)))
    with (OUT / "closed_carbon_cycle_screen.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "biomass", "central_cyclic_flux"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in closed_carbon:
            w.writerow(d)

    summary = dict(generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"), v1_sha256=v1_hash, central_reactions=len(central), reference_growth={cond: None for cond in conds}, free_redox=redox, ko_rescue=ko_rows, closed_carbon=closed_carbon, low_confidence=low, reversible_watchlist=rev)
    (OUT / "phase3_2_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE3.2 audit done")
    print("v1_sha256", v1_hash)
    print("central reactions", len(central))
    print("free redox", redox)
    print("closed carbon cycles:")
    for d in closed_carbon:
        print("  ", d["condition"], "cyclic count", len(d["central_cyclic_flux"].split(";")) if d["central_cyclic_flux"] not in ("", "infeasible") else 0)
    print("low confidence", len(low))
    for d in low:
        print("  LOW", d["reaction_id"], "conf", d["confidence"], "gpr", d["gpr"], "|", d["equation"][:60])

if __name__ == "__main__":
    main()
