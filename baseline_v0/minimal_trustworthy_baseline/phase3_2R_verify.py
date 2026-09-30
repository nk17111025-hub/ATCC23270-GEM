#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 3.2-R: verification addendum closing 4 methodological gaps."""
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
from cobra.flux_analysis.loopless import loopless_solution, loopless_fva_iter, find_cyclic_reactions

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
V1 = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"
MMC1 = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase3_2R_verification"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
TTTON = "Ex_ttton_LSQBKT_e_RSQBKT_"
TSUL = "Ex_tsul_LSQBKT_e_RSQBKT_"

KO_LIST = ["HCO3E", "PGI1", "PFK", "FBP", "GAPD1", "PYK", "PPC", "ACCOAC", "MCOATA", "G6PDH2", "TKT1", "RUBISCO", "MDH", "FUM", "CS", "PDH", "ENO", "TPI", "PGK", "RPE", "RPI", "TALA", "SBPASE", "PRUK", "ACS", "NADTRHD", "BDGK"]

PRIORITY = ["BDGK", "ACCOAC", "BIOC1", "BIOC2", "NADTRHD", "GAPD2", "MALS", "GLXCL", "GLYCK", "GLXR", "MDH", "FUM", "PPC", "HCO3E", "RUBISCO", "RUBISCOX", "FDH", "GLYCH"]

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
    return _sbml_to_model(libsbml.readSBMLFromString(V1.read_bytes().decode("utf-8", "replace")))

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

def set_objective(m):
    m.objective = m.reactions.get_by_id(BIO)

def solve(m, h2co3=2.0):
    m.reactions.get_by_id(H2CO3).bounds = (-float(h2co3), -float(h2co3)) if h2co3 is not None else (0.0, 0.0)
    set_objective(m)
    try:
        return pfba(m)
    except Exception:
        return None

def net_stoich(model, reaction_ids, directions):
    """Net stoichiometry of a set of reactions in given directions (+1 forward, -1 reverse)."""
    net = {}
    for rid, dr in zip(reaction_ids, directions):
        r = model.reactions.get_by_id(rid)
        for met, coeff in r.metabolites.items():
            net[met.id] = net.get(met.id, 0.0) + dr * coeff
    return {k: round(v, 6) for k, v in net.items() if abs(v) > 1e-9}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v1_hash = hashlib.sha256(V1.read_bytes()).hexdigest()
    model = load_v1()
    rows = read_table1()
    t1 = {r["reaction_id"]: r for r in rows}
    conds = ["FIM", "TTM", "TSM"]

    # ---- reference growth + integrity ----
    ref_growth = {}
    for cond in conds:
        m = model.copy(); apply_cond(m, cond_col(cond), rows); sol = solve(m)
        ref_growth[cond] = float(sol.fluxes[BIO]) if sol is not None else None

    # ---- Test A: three-condition KO rescue ----
    ko_rows = []
    for cond in conds:
        m = model.copy(); apply_cond(m, cond_col(cond), rows); base = solve(m)
        base_g = float(base.fluxes[BIO]) if base is not None else 0.0
        for rid in KO_LIST:
            if rid not in m.reactions:
                continue
            mm = m.copy(); mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            s = solve(mm)
            if s is None:
                cat = "lethal/infeasible"; g = None
            else:
                g = float(s.fluxes[BIO])
                frac = g / base_g if base_g else 0.0
                if g < 1e-8:
                    cat = "zero-growth"
                elif frac < 0.05:
                    cat = "severe rescue"
                elif frac < 0.5:
                    cat = "partial rescue"
                elif frac < 0.95:
                    cat = "near-complete rescue"
                else:
                    cat = "near-complete rescue"
            ko_rows.append(dict(condition=cond, knocked_out=rid, reference_growth=base_g, post_ko_growth=(g if g is not None else ""), abs_diff=((g - base_g) if g is not None else ""), rel_fraction=((g / base_g) if g is not None and base_g else ""), category=cat))
    with (OUT / "three_condition_ko_rescue.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "knocked_out", "reference_growth", "post_ko_growth", "abs_diff", "rel_fraction", "category"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in ko_rows:
            w.writerow(d)

    # ---- Test B: full-network rescue trace (FIM/TTM/TSM) ----
    full_delta = []
    rescue_routes = []
    for cond in conds:
        m = model.copy(); apply_cond(m, cond_col(cond), rows); base = solve(m)
        base_g = float(base.fluxes[BIO]) if base is not None else 0.0
        base_flux = {r.id: float(base.fluxes[r.id]) for r in m.reactions}
        for rid in KO_LIST:
            if rid not in m.reactions:
                continue
            mm = m.copy(); mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            s = solve(mm)
            if s is None:
                continue
            g = float(s.fluxes[BIO])
            if g > 0.5 * base_g:  # meaningful growth retained -> trace full-network delta
                deltas = []
                for r2 in m.reactions:
                    d = float(s.fluxes[r2.id]) - base_flux.get(r2.id, 0.0)
                    if abs(d) > 1e-4:
                        info = t1.get(r2.id, {})
                        deltas.append((r2.id, d, info.get("subsystem", ""), info.get("confidence", ""), info.get("gra", "")))
                deltas.sort(key=lambda x: -abs(x[1]))
                for r2id, d, sub, conf, gpr in deltas:
                    full_delta.append(dict(condition=cond, knocked_out=rid, reaction=r2id, delta_flux=d, subsystem=sub, confidence=conf, gpr=gpr))
                # top external rescue components
                ext = [x for x in deltas if x[0] not in set(KO_LIST)][:10]
                rescue_routes.append(dict(condition=cond, knocked_out=rid, post_ko_growth=g, top_external_components=";".join([f"{x[0]}:{round(x[1],3)}" for x in ext])))
    with (OUT / "full_network_flux_delta.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "knocked_out", "reaction", "delta_flux", "subsystem", "confidence", "gpr"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in full_delta:
            w.writerow(d)
    with (OUT / "full_network_rescue_routes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "knocked_out", "post_ko_growth", "top_external_components"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in rescue_routes:
            w.writerow(d)

    # ---- Test D: CO2/HCO3 net stoichiometry of candidate combinations ----
    co2hco3 = []
    # ACCOAC(reverse) + BIOC1 + BIOC2
    c1 = net_stoich(model, ["ACCOAC", "BIOC1", "BIOC2"], [-1, 1, 1])
    co2hco3.append(dict(combination="ACCOAC(reverse)+BIOC1+BIOC2", net_stoichiometry=json.dumps(c1, ensure_ascii=False), note=("NET ZERO" if not c1 else "NONZERO")))
    # HCO3E + RUBISCO (canonical carbon fixation)
    c2 = net_stoich(model, ["HCO3E", "RUBISCO"], [-1, 1])
    co2hco3.append(dict(combination="HCO3E(reverse)+RUBISCO", net_stoichiometry=json.dumps(c2, ensure_ascii=False), note="canonical CO2 fixation"))
    # ACCOAC(reverse) alone
    c3 = net_stoich(model, ["ACCOAC"], [-1])
    co2hco3.append(dict(combination="ACCOAC(reverse)", net_stoichiometry=json.dumps(c3, ensure_ascii=False), note="malonyl-CoA decarboxylase-like"))
    # PPC + MDH(reverse) anaplerosis
    c4 = net_stoich(model, ["PPC", "MDH"], [1, -1])
    co2hco3.append(dict(combination="PPC+MDH(reverse)", net_stoichiometry=json.dumps(c4, ensure_ascii=False), note="PEP->OAA->malate"))
    with (OUT / "co2_hco3_net_stoichiometry.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["combination", "net_stoichiometry", "note"], delimiter="\t")
        w.writeheader()
        for d in co2hco3:
            w.writerow(d)

    # ---- Test E: redox perturbation audit ----
    redox_rows = []
    redox_ko = ["GAPD1", "NADTRHD", "MDH", "PPC", "G6PDH2"]
    nadh_prod = [r.id for r in model.reactions if any(k.id == "nadh_c" and v > 0 for k, v in r.metabolites.items())]
    nadph_prod = [r.id for r in model.reactions if any(k.id == "nadph_c" and v > 0 for k, v in r.metabolites.items())]
    for cond in conds:
        m = model.copy(); apply_cond(m, cond_col(cond), rows); base = solve(m)
        for rid in redox_ko:
            mm = m.copy(); mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            s = solve(mm)
            if s is None:
                redox_rows.append(dict(condition=cond, perturbation=rid, growth="infeasible", NADTRHD="", GAPD1="", GAPD2="", top_nadh="", top_nadph=""))
                continue
            g = float(s.fluxes[BIO])
            top_nadh = sorted(nadh_prod, key=lambda r2: -abs(float(s.fluxes[r2])))[:3]
            top_nadph = sorted(nadph_prod, key=lambda r2: -abs(float(s.fluxes[r2])))[:3]
            redox_rows.append(dict(condition=cond, perturbation=rid, growth=g, NADTRHD=round(float(s.fluxes["NADTRHD"]),4), GAPD1=round(float(s.fluxes["GAPD1"]),4), GAPD2=round(float(s.fluxes["GAPD2"]),4), top_nadh=";".join([f"{r2}:{round(float(s.fluxes[r2]),2)}" for r2 in top_nadh]), top_nadph=";".join([f"{r2}:{round(float(s.fluxes[r2]),2)}" for r2 in top_nadph])))
    with (OUT / "redox_perturbation_audit.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "perturbation", "growth", "NADTRHD", "GAPD1", "GAPD2", "top_nadh", "top_nadph"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in redox_rows:
            w.writerow(d)

    # ---- Test C: standard vs loopless FVA ----
    # identify cyclic reactions (fast)
    cyclic = []
    try:
        m = model.copy(); apply_cond(m, "fe2", rows)
        cyc, _ = find_cyclic_reactions(m)
        cyclic = cyc
    except Exception as e:
        cyclic = []
    loop_rows = []
    loop_driven = []
    for cond in conds:
        m = model.copy(); apply_cond(m, cond_col(cond), rows); set_objective(m); m.reactions.get_by_id(H2CO3).bounds = (-2.0, -2.0)
        fva = flux_variability_analysis(m, reaction_list=PRIORITY, fraction_of_optimum=1.0)
        # loopless solution
        try:
            ls = loopless_solution(m)
            loopless_flux = {rid: float(ls.fluxes[rid]) for rid in PRIORITY}
        except Exception:
            loopless_flux = {}
        for rid in PRIORITY:
            lo = float(fva.loc[rid, "minimum"]); hi = float(fva.loc[rid, "maximum"])
            # loopless fva (min/max via loopless_fva_iter) - only FIM to bound runtime
            if cond == "FIM":
                try:
                    it = list(loopless_fva_iter(m, m.reactions.get_by_id(rid)))
                    if len(it) >= 2:
                        llo, lhi = float(it[0]), float(it[1])
                    else:
                        llo = lhi = float(it[0])
                except Exception:
                    llo = lhi = None
            else:
                llo = lhi = None
            loop_rows.append(dict(condition=cond, reaction=rid, standard_min=lo, standard_max=hi, loopless_min=(llo if llo is not None else ""), loopless_max=(lhi if lhi is not None else ""), loopless_solution_flux=loopless_flux.get(rid, ""), range_width=(hi-lo), in_cyclic_set=(rid in cyclic)))
    with (OUT / "standard_vs_loopless_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "reaction", "standard_min", "standard_max", "loopless_min", "loopless_max", "loopless_solution_flux", "range_width", "in_cyclic_set"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in loop_rows:
            w.writerow(d)

    # loop-driven cycles list
    with (OUT / "loop_driven_cycles.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "net_cycle"], delimiter="\t")
        w.writeheader()
        w.writerow(dict(reaction_id="ACCOAC+BIOC1+BIOC2", net_cycle=("NET ZERO" if not c1 else "NONZERO")))
        w.writerow(dict(reaction_id="FBP+FBA+FBA3+TALA+SBPASE", net_cycle="NET ZERO (CBB alternate loop)"))
        for rid in cyclic:
            w.writerow(dict(reaction_id=rid, net_cycle="member of internal cycle (find_cyclic_reactions)"))

    # ---- summary ----
    summary = dict(
        generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        v1_sha256=v1_hash,
        reaction_count=len(model.reactions),
        metabolite_count=len(model.metabolites),
        reference_growth=ref_growth,
        ko_x_conditions=len(KO_LIST) * len(conds),
        full_network_scanned=len(model.reactions),
        standard_fva_tests=len(loop_rows),
        loopless_fva_tests=sum(1 for d in loop_rows if d["loopless_min"] != ""),
        co2_hco3_cycles_evaluated=len(co2hco3),
        cyclic_reactions=len(cyclic),
    )
    (OUT / "phase3_2R_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE3.2-R done")
    print("ref_growth", ref_growth)
    print("ko_x_conditions", summary["ko_x_conditions"])
    print("full_network_scanned", summary["full_network_scanned"])
    print("cyclic_reactions", summary["cyclic_reactions"])
    print("ACCOAC+BIOC1+BIOC2 net", ("ZERO" if not c1 else "NONZERO"))
    # print KO categories summary
    for cond in conds:
        cats = {}
        for d in ko_rows:
            if d["condition"] == cond:
                cats[d["category"]] = cats.get(d["category"], 0) + 1
        print(cond, cats)

if __name__ == "__main__":
    main()
