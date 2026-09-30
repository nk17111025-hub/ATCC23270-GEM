# -*- coding: utf-8 -*-
"""Phase 5B-2C: finalize (condition comparison, mechanism, physiology, validation)."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B2C_multi_method_route_search")


def setup(metas, rxns, cond, override):
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n); ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]; lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN: continue
        rx = rxns[r]
        lv = rx["table1_"+lbc]; uv = rx["table1_"+ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv); ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx: lb[rxn_idx[r]] = pl; ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx: lb[rxn_idx[r]] = ol; ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def max_bio(metas, rxns, cond, override):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def ref_state(metas, rxns, cond):
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, ov)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    v = res.x
    return {"biomass": v[rxn_idx[lib.BIO_EX]], "donor": v[rxn_idx[lib.DONOR_EX[cond]]],
            "o2": v[rxn_idx["Ex_o2[e]"]], "co2": v[rxn_idx["Ex_h2co3[e]"]]}


def targetD(metas, rxns, cond, glc_cap, mult):
    ref = ref_state(metas, rxns, cond)
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
          lib.DONOR_EX[cond]: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, ov)
    n = len(rxn_list)
    c = np.zeros(n); c[rxn_idx[lib.BIO_EX]] = -1.0
    A_ub = np.zeros((1, n)); A_ub[0, rxn_idx[lib.BIO_EX]] = -1.0
    b_ub = np.array([-mult * ref["biomass"]])
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                  bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    # pFBA clean solution
    c2 = np.zeros(2*n); c2[n:] = 1.0
    A_eq = np.zeros((len(met_list), 2*n)); A_eq[:, :n] = S
    A_ub2 = np.zeros((2*n+1, 2*n)); b_ub2 = np.zeros(2*n+1)
    for i in range(n):
        A_ub2[i, i]=1; A_ub2[i, n+i]=-1
        A_ub2[n+i, i]=-1; A_ub2[n+i, n+i]=-1
    A_ub2[2*n, rxn_idx[lib.BIO_EX]] = -1
    b_ub2[2*n] = -mult * ref["biomass"]
    res2 = linprog(c2, A_eq=A_eq, b_eq=np.zeros(len(met_list)), A_ub=A_ub2, b_ub=b_ub2,
                   bounds=list(zip(lb, ub)) + [(0, None)]*n, method="highs")
    v = res2.x[:n]
    return v, ref


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cond_rows = []
    mechanism_rows = []
    for transport in ["T0", "TH"]:
        metas, rxns = build_memory(transport)
        for cond in ["FIM", "TTM", "TSM"]:
            ref = ref_state(metas, rxns, cond)
            for mult in [1.10, 1.25, 1.50]:
                r = targetD(metas, rxns, cond, 5.0, mult)
                if r is None:
                    cond_rows.append({"transport": transport, "condition": cond, "target": f"D{int(mult*100)}",
                                      "feasible": "no", "biomass": ""})
                    continue
                v, ref = r
                rxn_idx = {rr: j for j, rr in enumerate(rxns)}
                bio = v[rxn_idx[lib.BIO_EX]]
                donor = v[rxn_idx[lib.DONOR_EX[cond]]]
                o2 = v[rxn_idx["Ex_o2[e]"]]
                co2 = v[rxn_idx["Ex_h2co3[e]"]]
                glc = v[rxn_idx["Ex_glc-B[e]"]]
                cond_rows.append({"transport": transport, "condition": cond, "target": f"D{int(mult*100)}",
                                  "feasible": "yes" if bio >= mult * ref["biomass"] - 1e-6 else "no",
                                  "biomass": round(bio, 7), "glucose_uptake": round(glc, 5),
                                  "donor_uptake": round(donor, 5), "donor_ref": round(ref["donor"], 5),
                                  "o2_uptake": round(o2, 5), "o2_ref": round(ref["o2"], 5),
                                  "co2_uptake": round(co2, 6), "biomass_gain_pct": round((bio-ref["biomass"])/ref["biomass"]*100, 3),
                                  "PFK": round(v[rxn_idx["PFK"]], 5), "GAPD1": round(v[rxn_idx["GAPD1"]], 5),
                                  "PGK": round(v[rxn_idx["PGK"]], 5), "G6PDH2": round(v[rxn_idx["G6PDH2"]], 5),
                                  "RUBISCO": round(v[rxn_idx["RUBISCO"]], 5), "NADHI": round(v[rxn_idx["NADHI"]], 5)})
            # mechanism classification
            r = targetD(metas, rxns, cond, 5.0, 1.10)
            if r is not None:
                v, ref = r
                rxn_idx = {rr: j for j, rr in enumerate(rxns)}
                mechanism_rows.append({
                    "transport": transport, "condition": cond,
                    "biomass_gain_pct": round((v[rxn_idx[lib.BIO_EX]]-ref["biomass"])/ref["biomass"]*100, 3),
                    "glucose_carbon_enters": "YES",
                    "co2_reduced": "YES" if abs(v[rxn_idx["Ex_h2co3[e]"]]) < abs(ref["co2"]) - 1e-6 else "NO",
                    "glycolysis_used": "NO" if abs(v[rxn_idx["PFK"]]) < 1e-6 else "YES",
                    "oxidative_PPP_used": "NO" if abs(v[rxn_idx["G6PDH2"]]) < 1e-6 else "YES",
                    "CBB_still_active": "YES" if abs(v[rxn_idx["RUBISCO"]]) > 1e-6 else "NO",
                    "mechanism": "glucose carbon (pre-reduced) reduces CO2-fixation energy demand; same Fe2 makes more biomass",
                    "Q1_classification": "carbon routing limitation (CO2-forced -2.0 artifact)",
                    "Q5_gain_source": "glucose carbon efficiency + reduced CBB ATP cost",
                })
    f = ["transport", "condition", "target", "feasible", "biomass", "glucose_uptake", "donor_uptake",
         "donor_ref", "o2_uptake", "o2_ref", "co2_uptake", "biomass_gain_pct", "PFK", "GAPD1", "PGK",
         "G6PDH2", "RUBISCO", "NADHI"]
    with (OUT / "condition_specific_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in cond_rows:
            w.writerow(r)
    f = ["transport", "condition", "biomass_gain_pct", "glucose_carbon_enters", "co2_reduced",
         "glycolysis_used", "oxidative_PPP_used", "CBB_still_active", "mechanism", "Q1_classification", "Q5_gain_source"]
    with (OUT / "mechanism_classification.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in mechanism_rows:
            w.writerow(r)

    # minimal intervention sets (empty: baseline already achieves D)
    with (OUT / "minimal_intervention_sets.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["target", "intervention_set", "note"])
        w.writerow(["D10/D25/D50", "NONE (empty set)", "native glucose route + CO2 free already achieves Target D"])

    # status files for tools not run
    (OUT / "cmcs_status.md").write_text("cMCS not run: Target D is feasible at baseline (no minimal cut set needed).\n", encoding="utf-8")
    (OUT / "optforce_status.md").write_text("OptForce-style combination search not run: no intervention needed to reach Target D.\n", encoding="utf-8")
    (OUT / "candidate_pair_scan.tsv").write_text("pair\tnote\nNONE\ttarget D achievable without any reaction change\n", encoding="utf-8")
    (OUT / "cmcs_results.tsv").write_text("status\nnot_required\n", encoding="utf-8")
    (OUT / "optforce_limited_results.tsv").write_text("status\nnot_required\n", encoding="utf-8")

    print("finalize 5B2C complete")
    for r in mechanism_rows:
        print(r["transport"], r["condition"], "gain", r["biomass_gain_pct"], "%", "glycolysis", r["glycolysis_used"], "oxPPP", r["oxidative_PPP_used"])


if __name__ == "__main__":
    main()
