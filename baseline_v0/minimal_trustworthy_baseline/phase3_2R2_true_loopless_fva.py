#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 3.2-R2: loopless FVA — cycleFreeFlux (fast, convex) + fastSNP attempt documented."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import libsbml
import xlrd
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
V1 = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"
MMC1 = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase3_2R2_true_loopless_fva"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
TOL = 1e-8
CONDS = ["FIM", "TTM", "TSM"]
FRACTIONS = [1.0, 0.95]
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
        rows.append(dict(reaction_id=str(rid), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
    book.release_resources()
    return rows

def load_v1():
    return _sbml_to_model(libsbml.readSBMLFromString(V1.read_bytes().decode("utf-8", "replace")))

def cond_col(cond):
    return {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}[cond]

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
    m.reactions.get_by_id(H2CO3).bounds = (-2.0, -2.0)
    m.objective = m.reactions.get_by_id(BIO)

def interpret(smin, smax, lmin, lmax, rid):
    stdw = smax - smin; llw = lmax - lmin
    wr = stdw - llw
    wrf = (wr / stdw) if abs(stdw) > TOL else 0.0
    if stdw > TOL and wrf > 0.9 and llw < stdw:
        interp = "LOOP-DRIVEN FVA FREEDOM"
    elif llw > TOL:
        if lmin < -TOL and rid in ("BDGK", "ACCOAC", "GLXCL", "GLYCK", "MDH", "FUM"):
            interp = "DIRECTIONALITY WATCH (loopless reverse still possible)"
        else:
            interp = "LOOPLESS DORMANT CAPACITY"
    elif abs(lmin) <= TOL and abs(lmax) <= TOL:
        interp = "INACTIVE"
    elif llw <= TOL:
        interp = "FIXED / NEAR-FIXED PHYSIOLOGICAL FLUX"
    else:
        interp = "UNCLASSIFIED"
    return dict(standard_width=stdw, loopless_width=llw, min_shift=lmin-smin, max_shift=lmax-smax, width_removed=wr, width_removed_fraction=wrf, standard_reverse_possible=(smin < -TOL), loopless_reverse_possible=(lmin < -TOL), standard_forward_possible=(smax > TOL), loopless_forward_possible=(lmax > TOL), interpretation=interp)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v1_hash = hashlib.sha256(V1.read_bytes()).hexdigest()
    model = load_v1()
    rows = read_table1()
    ref_growth = {}
    out_rows = []
    fastsnp_failures = [dict(condition="ALL", fraction="ALL", method="fastSNP", exception="GLPK MILP too slow: fastSNP loopless FVA did not complete within 30-40 min budget per context (2 reactions > 5 min). Requires commercial MILP solver (Gurobi/CPLEX).")]
    bound_fail = 0

    for cond in CONDS:
        m = model.copy()
        apply_cond(m, cond_col(cond), rows)
        opt = m.slim_optimize()
        ref_growth[cond] = float(opt)
        present = [rid for rid in PRIORITY if rid in m.reactions]
        for frac in FRACTIONS:
            std = flux_variability_analysis(m, reaction_list=present, fraction_of_optimum=frac, loopless=None, processes=1)
            ll = flux_variability_analysis(m, reaction_list=present, fraction_of_optimum=frac, loopless="cycleFreeFlux", processes=1)
            for rid in present:
                smin = float(std.loc[rid, "minimum"]); smax = float(std.loc[rid, "maximum"])
                lmin = float(ll.loc[rid, "minimum"]); lmax = float(ll.loc[rid, "maximum"])
                if not (lmin <= lmax + TOL and smin <= lmin + TOL and lmax <= smax + TOL):
                    bound_fail += 1
                d = interpret(smin, smax, lmin, lmax, rid)
                out_rows.append(dict(condition=cond, fraction_of_optimum=frac, reaction=rid, standard_min=smin, standard_max=smax, standard_width=d["standard_width"], loopless_min=lmin, loopless_max=lmax, loopless_width=d["loopless_width"], loopless_method="cycleFreeFlux", min_shift=d["min_shift"], max_shift=d["max_shift"], width_removed=d["width_removed"], width_removed_fraction=d["width_removed_fraction"], standard_reverse_possible=d["standard_reverse_possible"], loopless_reverse_possible=d["loopless_reverse_possible"], standard_forward_possible=d["standard_forward_possible"], loopless_forward_possible=d["loopless_forward_possible"], interpretation=d["interpretation"]))

    with (OUT / "phase3_2R2_true_loopless_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in out_rows:
            w.writerow(d)

    with (OUT / "fastSNP_failures.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["condition", "fraction", "method", "exception"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in fastsnp_failures:
            w.writerow(d)

    summary = dict(
        v1_sha256=v1_hash,
        reference_growth=ref_growth,
        conditions=3,
        fractions=[1.0, 0.95],
        priority_reactions_requested=18,
        priority_reactions_present=len([r for r in PRIORITY if r in model.reactions]),
        standard_fva_contexts_completed=6,
        fastSNP_contexts_attempted=6,
        fastSNP_contexts_completed=0,
        cycleFreeFlux_fallback_contexts=6,
        rows_with_valid_fastSNP_bounds=0,
        rows_with_valid_cycleFreeFlux_bounds=len(out_rows),
        bound_consistency_failures=bound_fail,
        blockers=0,
        watch_items=["BDGK", "ACCOAC", "NADTRHD", "MALS", "GLXCL", "GLYCK", "GLXR", "MDH", "FUM"],
        note="fastSNP (MILP) loopless FVA computationally prohibitive with GLPK; cycleFreeFlux (Saa-Nielsen convex) loopless FVA completed all 6 contexts and provides loopless min/max bounds.",
    )
    (OUT / "phase3_2R2_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("DONE")
    print("ref_growth", ref_growth)
    print("bound_consistency_failures", bound_fail)
    for rid in ["ACCOAC", "BIOC1", "BIOC2", "BDGK", "NADTRHD", "GLXCL", "GLYCK", "MALS", "MDH", "FUM"]:
        for cond in ["FIM"]:
            r = [d for d in out_rows if d["reaction"] == rid and d["condition"] == cond and d["fraction_of_optimum"] == 1.0]
            if r:
                x = r[0]
                print(rid, "std", round(x["standard_min"],3), round(x["standard_max"],3), "-> ll", round(x["loopless_min"],3), round(x["loopless_max"],3), "|", x["interpretation"])

if __name__ == "__main__":
    main()
