#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 2.1: published artifact consistency audit + CYTRED direction diagnostic."""
import csv
import datetime as dt
import json
import re
from collections import defaultdict
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "reproduction_2016"

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
        rows.append(dict(reaction_id=str(rid), equation=str(cv(2)), subsystem=str(cv(6)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
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
            points.append(dict(idx=len(points)+1, condition="FIM", donor="fe2", mu=fe2_mu, donor_uptake=fe2_up, o2=fe2_o2, co2=fe2_co2))
        if isnum(s4_mu):
            points.append(dict(idx=len(points)+1, condition="TTM", donor="s4o6", mu=s4_mu, donor_uptake=s4_up, o2=None, co2=s4_co2))
        if isnum(s2_mu):
            points.append(dict(idx=len(points)+1, condition="TSM", donor="s2o3", mu=s2_mu, donor_uptake=s2_up, o2=None, co2=s2_co2))
    book.release_resources()
    return points

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

def max_sv(model, sol):
    S = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    order = [r.id for r in model.reactions]
    v = [float(sol.fluxes.get(rid, 0.0)) for rid in order]
    return float(np.max(np.abs(S @ v)))

def pearson(x, y):
    x = np.array(x, float); y = np.array(y, float)
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])

def r2_coeff(x, y):
    x = np.array(x, float); y = np.array(y, float)
    sse = float(((x - y) ** 2).sum())
    sst = float(((x - x.mean()) ** 2).sum())
    return 1.0 - sse / sst if sst else 0.0

def linfit(x, y):
    if len(x) < 2:
        return (0.0, 0.0)
    slope, intercept = np.polyfit(np.array(x, float), np.array(y, float), 1)
    return (float(slope), float(intercept))

def seven_metrics(points, repro):
    def sub(cond, key_exp, key_pred):
        ids = [p["idx"] for p in points if p["condition"] == cond]
        exp = [p[key_exp] for p in points if p["condition"] == cond]
        pred = [next(r for r in repro if r["idx"] == i)[key_pred] for i in ids]
        return exp, pred
    out = []
    for cond, name in [("FIM", "growth_FIM"), ("TTM", "growth_TTM"), ("TSM", "growth_TSM")]:
        e, p = sub(cond, "mu", "growth")
        out.append(dict(series=name, kind="growth", n=len(e), slope=linfit(e, p)[0], intercept=linfit(e, p)[1], pearson_r=pearson(e, p), r2_pearson=pearson(e, p)**2, r2_coefficient=r2_coeff(e, p)))
    for cond, name in [("FIM", "donor_FIM_fe2"), ("TTM", "donor_TTM"), ("TSM", "donor_TSM")]:
        e, p = sub(cond, "donor_uptake", "donor_magnitude")
        out.append(dict(series=name, kind="donor", n=len(e), slope=linfit(e, p)[0], intercept=linfit(e, p)[1], pearson_r=pearson(e, p), r2_pearson=pearson(e, p)**2, r2_coefficient=r2_coeff(e, p)))
    e, p = sub("FIM", "o2", "o2_magnitude")
    out.append(dict(series="o2_FIM", kind="o2", n=len(e), slope=linfit(e, p)[0], intercept=linfit(e, p)[1], pearson_r=pearson(e, p), r2_pearson=pearson(e, p)**2, r2_coefficient=r2_coeff(e, p)))
    return out

def run_variant(model, rows, points, variant):
    donor_rxn = {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}
    cond_col = {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}
    repro = []
    for p in points:
        m = model.copy()
        if variant != "SBML_ORIGINAL":
            cytr = m.reactions.get_by_id("CYTRED")
            if variant == "TABLE1_CYTRED":
                cytr.reaction = "etpcycAox_p + q8h2_c + 0.5 h_p --> etpcycArd_p + q8_c + 0.5 h_c"
            elif variant == "CYTRED_NO_PROTON":
                cytr.reaction = "etpcycAox_p + q8h2_c --> etpcycArd_p + q8_c"
        apply_operational(m, cond_col[p["condition"]], rows)
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective(m)
        sol = m.optimize()
        donor_flux = float(sol.fluxes[donor_rxn[p["condition"]]])
        o2_flux = float(sol.fluxes[O2]) if p["o2"] is not None else None
        repro.append(dict(idx=p["idx"], condition=p["condition"], growth=sol.objective_value, donor_magnitude=abs(donor_flux), donor_signed=donor_flux, o2_magnitude=(abs(o2_flux) if o2_flux is not None else ""), status=sol.status, max_sv=max_sv(m, sol)))
    return repro

def ttm_tsm_fva(model, rows, points, variant):
    cond_col = {"TTM": "ttton", "TSM": "tsul"}
    donor_rxn = {"TTM": TTTON, "TSM": TSUL}
    rows_out = []
    for p in points:
        if p["condition"] not in ("TTM", "TSM"):
            continue
        m = model.copy()
        if variant != "SBML_ORIGINAL":
            cytr = m.reactions.get_by_id("CYTRED")
            if variant == "TABLE1_CYTRED":
                cytr.reaction = "etpcycAox_p + q8h2_c + 0.5 h_p --> etpcycArd_p + q8_c + 0.5 h_c"
            elif variant == "CYTRED_NO_PROTON":
                cytr.reaction = "etpcycAox_p + q8h2_c --> etpcycArd_p + q8_c"
        apply_operational(m, cond_col[p["condition"]], rows)
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective(m)
        dr = donor_rxn[p["condition"]]
        fva = flux_variability_analysis(m, reaction_list=[dr], fraction_of_optimum=1.0)
        lo, hi = float(fva.loc[dr, "minimum"]), float(fva.loc[dr, "maximum"])
        exp = -float(p["donor_uptake"])
        within = (lo - 1e-6 <= exp <= hi + 1e-6)
        rows_out.append(dict(variant=variant, idx=p["idx"], condition=p["condition"], donor_exp_flux=exp, fva_min=lo, fva_max=hi, experimental_within_fva=within))
    return rows_out

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = load_model()
    rows = read_table1()
    points = read_table5()

    # CYTRED explicit re-verification
    cytr_sbml = model.reactions.get_by_id("CYTRED").build_reaction_string()
    t1_cytr = next(r["equation"] for r in rows if r["reaction_id"] == "CYTRED")

    # diff table (single mismatch already confirmed)
    with (OUT / "mmc1_vs_mmc3_reaction_equation_diff.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "classification", "mmc1_equation", "mmc3_equation"], delimiter="\t")
        w.writeheader()
        w.writerow(dict(reaction_id="CYTRED", classification="PROTON_DIRECTION_MISMATCH", mmc1_equation=t1_cytr, mmc3_equation=cytr_sbml))

    variants = ["SBML_ORIGINAL", "TABLE1_CYTRED", "CYTRED_NO_PROTON"]
    all_repro = {}
    all_metrics = {}
    all_fva = {}
    for v in variants:
        repro = run_variant(model, rows, points, v)
        metrics = seven_metrics(points, repro)
        fva = ttm_tsm_fva(model, rows, points, v)
        all_repro[v] = repro
        all_metrics[v] = metrics
        all_fva[v] = fva

    # write cytred variant 20-point reproduction
    with (OUT / "cytred_variant_20point_reproduction.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["variant", "idx", "condition", "growth", "donor_signed", "donor_magnitude", "o2_magnitude", "status", "max_sv"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for v in variants:
            for r in all_repro[v]:
                w.writerow(dict(variant=v, **r))

    # write seven metrics per variant
    with (OUT / "cytred_variant_seven_metrics.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["variant", "series", "kind", "n", "slope", "intercept", "pearson_r", "r2_pearson", "r2_coefficient"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for v in variants:
            for m in all_metrics[v]:
                w.writerow(dict(variant=v, **m))

    # write TTM/TSM FVA
    with (OUT / "cytred_variant_TTM_TSM_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["variant", "idx", "condition", "donor_exp_flux", "fva_min", "fva_max", "experimental_within_fva"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for v in variants:
            for r in all_fva[v]:
                w.writerow(r)

    # exact solver on 4 problem points
    exact = []
    for idx in [2, 3, 19, 20]:
        p = next(x for x in points if x["idx"] == idx)
        cond_col = {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}
        m = model.copy()
        apply_operational(m, cond_col[p["condition"]], rows)
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective(m)
        m.solver = "glpk_exact"
        sol = m.optimize()
        donor_rxn = {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}[p["condition"]]
        donor_flux = float(sol.fluxes[donor_rxn])
        o2_flux = float(sol.fluxes[O2]) if p["o2"] is not None else ""
        exact.append(dict(idx=idx, condition=p["condition"], objective=sol.objective_value, donor_uptake_signed=donor_flux, donor_uptake_magnitude=abs(donor_flux), o2_uptake_signed=(o2_flux if o2_flux != "" else ""), o2_uptake_magnitude=(abs(o2_flux) if o2_flux != "" else ""), max_sv=max_sv(m, sol), status=sol.status))
    with (OUT / "phase2_exact_four_problem_points.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["idx", "condition", "objective", "donor_uptake_signed", "donor_uptake_magnitude", "o2_uptake_signed", "o2_uptake_magnitude", "max_sv", "status"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in exact:
            w.writerow(d)

    # summary
    mean_r2 = {v: float(np.mean([m["r2_coefficient"] for m in all_metrics[v]])) for v in variants}
    tsm_within = {}
    for v in variants:
        tsm_fva = [r for r in all_fva[v] if r["condition"] == "TSM"]
        tsm_within[v] = dict(total=len(tsm_fva), within=sum(1 for r in tsm_fva if r["experimental_within_fva"]))
    summary = dict(generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"), cytred_sbml=cytr_sbml, cytred_table1=t1_cytr, genuine_mismatches=1, mean_r2_coefficient=mean_r2, tsm_donor_within_fva=tsm_within, exact_four_points=exact)
    (OUT / "phase2_1_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE2.1 done")
    print("mean r2_coefficient", mean_r2)
    print("tsm donor within fva", tsm_within)
    for v in variants:
        print("===", v)
        for m in all_metrics[v]:
            if m["kind"] in ("donor", "growth"):
                print("  ", m["series"], "slope=", round(m["slope"], 3), "r2c=", round(m["r2_coefficient"], 4))
    print("exact", [(d["idx"], d["max_sv"]) for d in exact])

if __name__ == "__main__":
    main()
