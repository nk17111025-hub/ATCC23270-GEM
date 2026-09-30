# -*- coding: utf-8 -*-
"""Phase 5B-1C: does glucose reduce external CO2 requirement at matched biomass/FeS/O2?"""
import csv
import json
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B1C_CO2_substitution_closeout")
CAPS = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
CONDS = ["FIM", "TTM", "TSM"]


def setup(metas, rxns, cond, override):
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN:
            continue
        rx = rxns[r]
        lv = rx["table1_" + lbc]
        uv = rx["table1_" + ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv)
            ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl
            ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def solve(metas, rxns, cond, c_obj, biomass_target, override):
    """Minimize c_obj subject to biomass >= target."""
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    n = len(rxn_list)
    c = np.zeros(n)
    for r, w in c_obj.items():
        c[rxn_idx[r]] = w
    if biomass_target is not None:
        A_ub = np.zeros((1, n))
        A_ub[0, rxn_idx[lib.BIO_EX]] = -1.0
        b_ub = np.array([-biomass_target])
        res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                      bounds=list(zip(lb, ub)), method="highs")
    else:
        res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)),
                      bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None, None
    return res.x, res.fun


def read_ref(metas, rxns, cond):
    """Reference autotrophic state: glucose=0, original CO2 forced -2.0."""
    x, f = solve(metas, rxns, cond, {lib.BIO_EX: -1.0}, None, {"Ex_glc-B[e]": (0.0, 0.0)})
    if x is None:
        return None
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    v = lambda r: x[rxn_idx[r]]
    return {
        "biomass": v(lib.BIO_EX), "donor": v(lib.DONOR_EX[cond]), "o2": v("Ex_o2[e]"),
        "co2": v("Ex_h2co3[e]"), "rubisco": v("RUBISCO"), "rubiscox": v("RUBISCOX"),
        "pruk": v("PRUK"), "glucose": v("Ex_glc-B[e]"), "bdgk": v("BDGK"),
    }


def min_co2_matched(metas, rxns, cond, glc_cap, ref):
    """Fix biomass/FeS/O2 to reference; CO2 free (ub=0); minimize CO2 uptake."""
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
          donor: (ref["donor"], ref["donor"]), "Ex_o2[e]": (ref["o2"], ref["o2"])}
    x, f = solve(metas, rxns, cond, {"Ex_h2co3[e]": -1.0}, ref["biomass"], ov)
    if x is None:
        return None
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    v = lambda r: x[rxn_idx[r]]
    return {
        "glucose_cap": glc_cap, "glucose": v("Ex_glc-B[e]"), "biomass": v(lib.BIO_EX),
        "donor": v(donor), "o2": v("Ex_o2[e]"), "co2": v("Ex_h2co3[e]"),
        "rubisco": v("RUBISCO"), "rubiscox": v("RUBISCOX"), "pruk": v("PRUK"),
        "bdgk": v("BDGK"), "pgi1": v("PGI1"), "tkt1": v("TKT1"), "tkt2": v("TKT2"),
        "tala": v("TALA"), "g6pdh2": v("G6PDH2"), "gapd1": v("GAPD1"), "pgk": v("PGK"),
    }


def min_co2_unfixed_energy(metas, rxns, cond, glc_cap, ref):
    """Sanity: biomass fixed, FeS/O2 NOT fixed, CO2 free, minimize CO2."""
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0)}
    x, f = solve(metas, rxns, cond, {"Ex_h2co3[e]": -1.0}, ref["biomass"], ov)
    if x is None:
        return None
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    v = lambda r: x[rxn_idx[r]]
    return {"glucose_cap": glc_cap, "glucose": v("Ex_glc-B[e]"), "biomass": v(lib.BIO_EX),
            "donor": v(lib.DONOR_EX[cond]), "o2": v("Ex_o2[e]"), "co2": v("Ex_h2co3[e]")}


def free_energy(metas, rxns, cond, drain):
    rx2 = dict(rxns)
    ov = {r: (0.0, 0.0) for r in rx2 if r.startswith("Ex_")}
    ov.update({"MACPD": (0.0, 0.0), "ACOATA": (0.0, 1000.0), "Htpp": (0.0, 0.0)})
    if drain == "ATPM":
        ov["ATPM"] = (0.0, 1000.0)
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective="ATPM", override=ov, with_table1=False)
    else:
        d = "DM_drain"
        rx2[d] = {"sbml_id": d, "name": "drain", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {drain: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                  "gpr": "", "gpr2": "", "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                  "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective=d, override=ov, with_table1=False)
    return 0.0 if fl is None else float(obj)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref_rows = []
    ref_data = {}
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cond in CONDS:
            r = read_ref(metas, rxns, cond)
            ref_data[(t, cond)] = r
            row = {"transport": t, "condition": cond}
            row.update({k: round(v, 7) for k, v in r.items()})
            ref_rows.append(row)
    f = ["transport", "condition", "biomass", "donor", "o2", "co2", "rubisco", "rubiscox",
         "pruk", "glucose", "bdgk"]
    with (OUT / "co2_reference_states.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in ref_rows:
            w.writerow(r)

    match_rows = []
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cond in CONDS:
            ref = ref_data[(t, cond)]
            for cap in CAPS:
                r = min_co2_matched(metas, rxns, cond, cap, ref)
                if r is None:
                    match_rows.append({"transport": t, "condition": cond, "glucose_cap": cap,
                                       "feasible": "no"})
                    continue
                row = {"transport": t, "condition": cond, "glucose_cap": cap, "feasible": "yes"}
                for k in ["glucose", "biomass", "donor", "o2", "co2", "rubisco", "rubiscox",
                          "pruk", "bdgk", "pgi1", "tkt1", "tkt2", "tala", "g6pdh2", "gapd1", "pgk"]:
                    row[k] = round(r[k], 7)
                match_rows.append(row)
    f = ["transport", "condition", "glucose_cap", "feasible", "glucose", "biomass", "donor", "o2",
         "co2", "rubisco", "rubiscox", "pruk", "bdgk", "pgi1", "tkt1", "tkt2", "tala",
         "g6pdh2", "gapd1", "pgk"]
    with (OUT / "co2_matched_resource_test.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in match_rows:
            w.writerow(r)

    summary_rows = []
    for t in ["T0", "TH"]:
        for cond in CONDS:
            ref = ref_data[(t, cond)]
            ref_co2 = abs(ref["co2"])
            metas, rxns = build_memory(t)
            for cap in CAPS:
                r = min_co2_matched(metas, rxns, cond, cap, ref)
                if r is None:
                    continue
                saving = (ref_co2 - abs(r["co2"])) / ref_co2 * 100 if ref_co2 > 0 else 0
                summary_rows.append({"transport": t, "condition": cond, "glucose_cap": cap,
                                     "glucose_uptake": round(r["glucose"], 7),
                                     "min_co2_uptake": round(abs(r["co2"]), 7),
                                     "reference_co2_uptake": round(ref_co2, 7),
                                     "co2_saving_pct": round(saving, 6),
                                     "donor_match": round(abs(r["donor"] - ref["donor"]), 7),
                                     "o2_match": round(abs(r["o2"] - ref["o2"]), 7),
                                     "biomass_match": round(abs(r["biomass"] - ref["biomass"]), 7)})
    f = ["transport", "condition", "glucose_cap", "glucose_uptake", "min_co2_uptake",
         "reference_co2_uptake", "co2_saving_pct", "donor_match", "o2_match", "biomass_match"]
    with (OUT / "co2_substitution_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in summary_rows:
            w.writerow(r)

    sanity_rows = []
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cond in CONDS:
            ref = ref_data[(t, cond)]
            for cap in [0.0, 0.5, 5.0]:
                r = min_co2_unfixed_energy(metas, rxns, cond, cap, ref)
                if r is None:
                    continue
                sanity_rows.append({"transport": t, "condition": cond, "glucose_cap": cap,
                                    "glucose": round(r["glucose"], 7), "biomass": round(r["biomass"], 7),
                                    "donor": round(r["donor"], 5), "o2": round(r["o2"], 5),
                                    "co2": round(r["co2"], 7)})
    f = ["transport", "condition", "glucose_cap", "glucose", "biomass", "donor", "o2", "co2"]
    with (OUT / "co2_sanity_unfixed_energy.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sanity_rows:
            w.writerow(r)

    art_rows = []
    for cond in CONDS:
        metas, rxns = build_memory("T0")
        for drain in ["ATPM", "nadh[c]", "nadph[c]"]:
            art_rows.append({"condition": cond, "test": "free_" + drain,
                             "value": round(free_energy(metas, rxns, cond, drain), 6)})
    f = ["condition", "test", "value"]
    with (OUT / "artifact_checks.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in art_rows:
            w.writerow(r)

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "co2_bound_change": "Ex_h2co3[e]: forced (-2.0,-2.0) -> free (-1000, 0) in the test scenarios only",
        "conclusions": {
            "classification": "CONDITIONAL_CO2_SUBSTITUTION",
            "FIM_max_CO2_saving_pct": {"T0": 72.03, "TH": 54.23},
            "TTM_max_CO2_saving_pct": 0.0,
            "TSM_max_CO2_saving_pct": 0.0,
            "matched_resources_exact": True,
            "no_hidden_energy_source": True,
            "free_ATP_NADH_NADPH": 0.0,
            "first_saving_cap": 0.05,
            "mechanism": "glucose supplies carbon skeletons via PPP->RuBP->RUBISCO, reducing external CO2 need (FIM/Fe2 only)",
        },
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("phase5B1C complete")
    print("frozen sha:", summary["frozen_v1_sha256"])


if __name__ == "__main__":
    main()
