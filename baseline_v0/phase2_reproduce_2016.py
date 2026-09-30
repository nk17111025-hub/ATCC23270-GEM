#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 2: strict 2016 published-model reproduction (20 Table 5 points)."""
import csv
import datetime as dt
import json
import math
import re
from collections import defaultdict
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
OUT = ROOT / "baseline_v0" / "reproduction_2016"
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
        rows.append(dict(reaction_id=str(rid), name=str(cv(1)), equation=str(cv(2)), subsystem=str(cv(6)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
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
    doc = libsbml.readSBMLFromString(raw.decode("utf-8", "replace"))
    return _sbml_to_model(doc)

def apply_operational_bounds(model, cond, rows):
    ids = {r.id for r in model.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids:
            continue
        if cid == "Afe_biomass_mc507_WT_139p0M":
            continue  # keep SBML 0/1000 (published biomass objective)
        rx = model.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))

def set_objective_biomass(model):
    model.objective = model.reactions.get_by_id(BIO)

def max_sv(model, sol):
    S = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    order = [r.id for r in model.reactions]
    v = [float(sol.fluxes.get(rid, 0.0)) for rid in order]
    return float(np.max(np.abs(S @ v)))

def pearson(x, y):
    x = np.array(x, float)
    y = np.array(y, float)
    if x.std() == 0 or y.std() == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])

def r2_pearson(x, y):
    r = pearson(x, y)
    return r * r

def r2_coefficient(x, y):
    x = np.array(x, float)
    y = np.array(y, float)
    sse = float(((x - y) ** 2).sum())
    sst = float(((x - x.mean()) ** 2).sum())
    if sst == 0:
        return 0.0
    return 1.0 - sse / sst

def linfit(x, y):
    x = np.array(x, float)
    y = np.array(y, float)
    if len(x) < 2:
        return (0.0, 0.0)
    slope, intercept = np.polyfit(x, y, 1)
    return (float(slope), float(intercept))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    model = load_model()
    rows = read_table1()
    points = read_table5()
    t1map = {r["reaction_id"]: r for r in rows}
    donor_rxn = {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}
    cond_col = {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}

    # ---- Section 1A: Afe_biomass exception evidence ----
    ab_xml = model.reactions.get_by_id("Afe_biomass_mc507_WT_139p0M")
    ab_t1 = t1map["Afe_biomass_mc507_WT_139p0M"]

    # ---- Section 1B: count recalc ----
    recalc = {}
    for cond, label in [("fe2", "FIM"), ("ttton", "TTM"), ("tsul", "TSM")]:
        ex = tr = internal = biomass = 0
        for r in rows:
            cid = enc(r["reaction_id"])
            xb = model.reactions.get_by_id(cid)
            tlb, tub = r["lb_" + cond], r["ub_" + cond]
            if abs(float(tlb) - xb.lower_bound) > 1e-9 or abs(float(tub) - xb.upper_bound) > 1e-9:
                if cid.startswith("Ex_"):
                    ex += 1
                elif cid == "Afe_biomass_mc507_WT_139p0M":
                    biomass += 1
                else:
                    internal += 1
        recalc[label] = dict(exchange=ex, internal_other=internal, afe_biomass=biomass, total=ex + internal + biomass)

    # ---- Section 2 + 3 + 4: 20-point reproduction, FVA, pFBA ----
    repro = []
    fva_out = []
    pfba_out = []
    for p in points:
        cond = p["condition"]
        cc = cond_col[cond]
        m = model.copy()
        apply_operational_bounds(m, cc, rows)
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective_biomass(m)
        sol = m.optimize()
        donor_flux = float(sol.fluxes[donor_rxn[cond]])
        o2_flux = float(sol.fluxes[O2]) if p["o2"] is not None else None
        repro.append(dict(idx=p["idx"], condition=cond, donor=p["donor"], h2co3_input=float(p["co2"]), growth_pred=sol.objective_value, growth_exp=p["mu"], donor_pred_signed=donor_flux, donor_pred_magnitude=abs(donor_flux), donor_exp=p["donor_uptake"], o2_pred_signed=(o2_flux if o2_flux is not None else ""), o2_pred_magnitude=(abs(o2_flux) if o2_flux is not None else ""), o2_exp=(p["o2"] if p["o2"] is not None else ""), status=sol.status, max_sv=max_sv(m, sol)))

        # FVA at 100% max biomass
        fvarxns = [donor_rxn[cond], O2, "Htpp", "ATPS5rpp", "ATPM"]
        fva = flux_variability_analysis(m, reaction_list=fvarxns, fraction_of_optimum=1.0)
        for rid in fvarxns:
            lo = float(fva.loc[rid, "minimum"])
            hi = float(fva.loc[rid, "maximum"])
            fva_out.append(dict(idx=p["idx"], condition=cond, reaction=rid, minimum=lo, maximum=hi))

        # pFBA
        m2 = model.copy()
        apply_operational_bounds(m2, cc, rows)
        m2.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective_biomass(m2)
        psol = pfba(m2)
        pfba_out.append(dict(idx=p["idx"], condition=cond, donor=p["donor"], growth=psol.objective_value, donor_flux=float(psol.fluxes[donor_rxn[cond]]), donor_magnitude=abs(float(psol.fluxes[donor_rxn[cond]])), o2_flux=(float(psol.fluxes[O2]) if p["o2"] is not None else ""), o2_magnitude=(abs(float(psol.fluxes[O2])) if p["o2"] is not None else ""), Htpp=float(psol.fluxes["Htpp"]), ATPS5rpp=float(psol.fluxes["ATPS5rpp"]), ATPM=float(psol.fluxes["ATPM"])))

    with (OUT / "table5_20point_reproduction.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["idx", "condition", "donor", "h2co3_input", "growth_pred", "growth_exp", "donor_pred_signed", "donor_pred_magnitude", "donor_exp", "o2_pred_signed", "o2_pred_magnitude", "o2_exp", "status", "max_sv"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in repro:
            w.writerow(d)
    with (OUT / "table5_optimal_fva.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["idx", "condition", "reaction", "minimum", "maximum"], delimiter="\t")
        w.writeheader()
        for d in fva_out:
            w.writerow(d)
    with (OUT / "table5_pfba_diagnostic.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["idx", "condition", "donor", "growth", "donor_flux", "donor_magnitude", "o2_flux", "o2_magnitude", "Htpp", "ATPS5rpp", "ATPM"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in pfba_out:
            w.writerow(d)

    # ---- classification of validation targets ----
    classification = []
    tol = 1e-6
    for p in points:
        cond = p["condition"]
        idx = p["idx"]
        dr = donor_rxn[cond]
        # growth
        growth_unique = True
        # donor
        d_fva = [d for d in fva_out if d["idx"] == idx and d["reaction"] == dr][0]
        lo, hi = d_fva["minimum"], d_fva["maximum"]
        exp_flux = -float(p["donor_uptake"])
        if hi - lo < tol:
            dclass = "UNIQUE_OR_NEAR_UNIQUE"
        elif lo - tol <= exp_flux <= hi + tol:
            dclass = "EXPERIMENT_WITHIN_OPTIMAL_FVA"
        else:
            dclass = "EXPERIMENT_OUTSIDE_OPTIMAL_FVA"
        if p["o2"] is not None:
            o_fva = [d for d in fva_out if d["idx"] == idx and d["reaction"] == O2][0]
            olo, ohi = o_fva["minimum"], o_fva["maximum"]
            oexp = -float(p["o2"])
            if ohi - olo < tol:
                oclass = "UNIQUE_OR_NEAR_UNIQUE"
            elif olo - tol <= oexp <= ohi + tol:
                oclass = "EXPERIMENT_WITHIN_OPTIMAL_FVA"
            else:
                oclass = "EXPERIMENT_OUTSIDE_OPTIMAL_FVA"
        else:
            oclass = "N/A"
        classification.append(dict(idx=idx, condition=cond, donor_target_class=dclass, o2_target_class=oclass))

    # ---- Section 5: seven validation series ----
    series_def = []
    for cond, label in [("FIM", "growth_FIM"), ("TTM", "growth_TTM"), ("TSM", "growth_TSM")]:
        sub = [p for p in points if p["condition"] == cond]
        series_def.append(dict(series=label, kind="growth", condition=cond, n=len(sub), exp=[p["mu"] for p in sub], pred=[next(r for r in repro if r["idx"] == p["idx"])["growth_pred"] for p in sub]))
    for cond, label, dname in [("FIM", "donor_FIM_fe2", "fe2"), ("TTM", "donor_TTM_tetrathionate", "s4o6"), ("TSM", "donor_TSM_thiosulfate", "s2o3")]:
        sub = [p for p in points if p["condition"] == cond]
        series_def.append(dict(series=label, kind="donor_uptake", condition=cond, n=len(sub), exp=[p["donor_uptake"] for p in sub], pred=[next(r for r in repro if r["idx"] == p["idx"])["donor_pred_magnitude"] for p in sub]))
    sub_o2 = [p for p in points if p["condition"] == "FIM"]
    series_def.append(dict(series="o2_FIM", kind="o2_uptake", condition="FIM", n=len(sub_o2), exp=[p["o2"] for p in sub_o2], pred=[next(r for r in repro if r["idx"] == p["idx"])["o2_pred_magnitude"] for p in sub_o2]))

    series_rows = []
    metrics = []
    for sd in series_def:
        slope, intercept = linfit(sd["exp"], sd["pred"])
        r = pearson(sd["exp"], sd["pred"])
        r2p = r2_pearson(sd["exp"], sd["pred"])
        r2c = r2_coefficient(sd["exp"], sd["pred"])
        metrics.append(dict(series=sd["series"], kind=sd["kind"], n=sd["n"], slope=slope, intercept=intercept, pearson_r=r, r2_pearson=r2p, r2_coefficient_1_minus_SSE_SST=r2c))
        for i in range(sd["n"]):
            series_rows.append(dict(series=sd["series"], kind=sd["kind"], condition=sd["condition"], point=i + 1, experimental=sd["exp"][i], predicted=sd["pred"][i]))
    with (OUT / "seven_validation_series.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["series", "kind", "condition", "point", "experimental", "predicted"], delimiter="\t")
        w.writeheader()
        for d in series_rows:
            w.writerow(d)
    with (OUT / "seven_validation_metrics.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["series", "kind", "n", "slope", "intercept", "pearson_r", "r2_pearson", "r2_coefficient_1_minus_SSE_SST"], delimiter="\t")
        w.writeheader()
        for d in metrics:
            w.writerow(d)

    # ---- Section 6: published energy parameter check ----
    energy_target = {"ATPS5rpp": 5.0, "CYTAA31": 0.0, "NADHI": 5.0, "CYTBC1": 5.0, "CYTAA32": 1.3, "CYTRED": 0.5, "CYTBO3": 1.8, "CYTBD": 1.8}
    energy_rows = []
    for rid, pub in energy_target.items():
        r = model.reactions.get_by_id(rid)
        # net protons translocated across membrane = |h_p consumed - h_p produced| (or |h_c produced - h_c consumed|)
        hp_react = sum(v for k, v in r.metabolites.items() if k.id == "h_p" and v < 0)
        hp_prod = sum(v for k, v in r.metabolites.items() if k.id == "h_p" and v > 0)
        hc_react = sum(v for k, v in r.metabolites.items() if k.id == "h_c" and v < 0)
        hc_prod = sum(v for k, v in r.metabolites.items() if k.id == "h_c" and v > 0)
        net_p_to_c = abs(hp_react) - hp_prod  # protons leaving p (consumed in p) minus entering p (produced)
        magnitude = abs(net_p_to_c)
        energy_rows.append(dict(reaction_id=rid, name=r.name, published_protons=pub, encoded_protons=magnitude, h_p_reactant=abs(hp_react), h_p_product=hp_prod, h_c_reactant=abs(hc_react), h_c_product=hc_prod, match=(abs(magnitude - pub) < 1e-6)))
    # GAM / NGAM
    ab = model.reactions.get_by_id("Afe_biomass_mc507_WT_139p0M")
    atp_c = ab.get_coefficient("atp_c")
    h2o_c = ab.get_coefficient("h2o_c")
    atpm = model.reactions.get_by_id("ATPM")
    energy_rows.append(dict(reaction_id="GAM_(Afe_biomass)", name="biomass GAM", published_protons=139, encoded_protons=abs(atp_c), h_p_reactant="", h_p_product="", h_c_reactant="", h_c_product="", match=(abs(abs(atp_c) - 139) < 0.5)))
    energy_rows.append(dict(reaction_id="NGAM_(ATPM)", name="ATP maintenance", published_protons=3.475, encoded_protons=atpm.lower_bound, h_p_reactant="", h_p_product="", h_c_reactant="", h_c_product="", match=(abs(atpm.lower_bound - 3.475) < 1e-6)))
    with (OUT / "published_energy_parameter_check.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "name", "published_protons", "encoded_protons", "h_p_reactant", "h_p_product", "h_c_reactant", "h_c_product", "match"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in energy_rows:
            w.writerow(d)

    # ---- Section 7: EGC status ----
    htpp_all_zero = all(abs(float(psol["Htpp"])) < 1e-9 for psol in pfba_out) and all(abs(r["minimum"]) < 1e-9 and abs(r["maximum"]) < 1e-9 for r in fva_out if r["reaction"] == "Htpp")

    # ---- Section 8: GLPK-exact cross-check on 3 representative points ----
    exact = []
    for idx in [1, 10, 20]:
        p = next(x for x in points if x["idx"] == idx)
        cc = cond_col[p["condition"]]
        m = model.copy()
        apply_operational_bounds(m, cc, rows)
        m.reactions.get_by_id(H2CO3).bounds = (-float(p["co2"]), -float(p["co2"]))
        set_objective_biomass(m)
        sol = m.optimize()
        glpk_obj = sol.objective_value
        m.solver = "glpk_exact"
        sol_exact = m.optimize()
        exact.append(dict(idx=idx, condition=p["condition"], glpk_objective=glpk_obj, glpk_exact_objective=sol_exact.objective_value, diff=sol_exact.objective_value - glpk_obj))
    with (OUT / "table5_glpk_exact_crosscheck.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["idx", "condition", "glpk_objective", "glpk_exact_objective", "diff"], delimiter="\t")
        w.writeheader()
        for d in exact:
            w.writerow(d)

    # classification TSV
    with (OUT / "table5_target_classification.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["idx", "condition", "donor_target_class", "o2_target_class"], delimiter="\t")
        w.writeheader()
        for d in classification:
            w.writerow(d)

    summary = dict(generated_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"), afe_biomass_xml_bounds=[ab_xml.lower_bound, ab_xml.upper_bound], afe_biomass_table1_bounds=[ab_t1["lb_fe2"], ab_t1["ub_fe2"]], count_recalc=recalc, feasible_all=all(d["status"] == "optimal" for d in repro), max_sv_max=max(d["max_sv"] for d in repro), htpp_all_zero=htpp_all_zero, mean_r2_pearson=float(np.mean([m["r2_pearson"] for m in metrics])), mean_r2_coefficient=float(np.mean([m["r2_coefficient_1_minus_SSE_SST"] for m in metrics])), energy_matches=all(d["match"] for d in energy_rows), exact=exact)
    (OUT / "phase2_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("PHASE2 done")
    print("feasible_all", summary["feasible_all"], "max_sv_max", summary["max_sv_max"])
    print("mean r2_pearson", summary["mean_r2_pearson"], "mean r2_coefficient", summary["mean_r2_coefficient"])
    print("htpp_all_zero", htpp_all_zero, "energy_matches", summary["energy_matches"])
    print("metrics")
    for m in metrics:
        print("  ", m["series"], "n=", m["n"], "r=", round(m["pearson_r"], 4), "r2p=", round(m["r2_pearson"], 4), "r2c=", round(m["r2_coefficient_1_minus_SSE_SST"], 4), "slope=", round(m["slope"], 3))
    print("energy")
    for d in energy_rows:
        print("  ", d["reaction_id"], d["published_protons"], "->", d["encoded_protons"], "match=", d["match"])

if __name__ == "__main__":
    main()
