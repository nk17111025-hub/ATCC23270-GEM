# -*- coding: utf-8 -*-
"""Phase 5B-3A: minimal-reaction gapfill for direct (Rubisco-independent) glucose assimilation."""
import csv
import json
import itertools
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3A_gapfill_direct_glucose")

# Curated candidate universe: real central-carbon reactions commonly present in
# heterotrophic bacteria but absent (or unidirectional) in the autotrophic ATCC 23270 GEM.
# Each entry: (id, name, EC, reversible, stoich{met:coeff})
UNIVERSE = [
    ("ICL", "isocitrate lyase", "4.1.3.1", False,
     {"icit[c]": -1.0, "glx[c]": 1.0, "succ[c]": 1.0}),
    ("ME_NADP", "malic enzyme (NADP)", "1.1.1.40", False,
     {"mal-L[c]": -1.0, "nadp[c]": -1.0, "pyr[c]": 1.0, "co2[c]": 1.0, "nadph[c]": 1.0}),
    ("ME_NAD", "malic enzyme (NAD)", "1.1.1.38", False,
     {"mal-L[c]": -1.0, "nad[c]": -1.0, "pyr[c]": 1.0, "co2[c]": 1.0, "nadh[c]": 1.0}),
    ("PC", "pyruvate carboxylase", "6.4.1.1", False,
     {"pyr[c]": -1.0, "co2[c]": -1.0, "atp[c]": -1.0, "oaa[c]": 1.0, "adp[c]": 1.0, "pi[c]": 1.0, "h[c]": 1.0}),
    ("PEPCK", "PEP carboxykinase (ATP)", "4.1.1.49", True,
     {"oaa[c]": -1.0, "atp[c]": -1.0, "pep[c]": 1.0, "co2[c]": 1.0, "adp[c]": 1.0}),
    ("PEPS", "PEP synthase", "2.7.9.2", False,
     {"pyr[c]": -1.0, "atp[c]": -1.0, "h2o[c]": -1.0, "pep[c]": 1.0, "amp[c]": 1.0, "pi[c]": 1.0, "h[c]": 2.0}),
    ("NADHDH", "NADH:ubiquinone oxidoreductase (forward)", "1.6.5.3", False,
     {"nadh[c]": -1.0, "q8[c]": -1.0, "h[c]": -1.0, "nad[c]": 1.0, "q8h2[c]": 1.0}),
    ("G6PDH_NAD", "glucose-6-phosphate dehydrogenase (NAD)", "1.1.1.363", False,
     {"g6p-B[c]": -1.0, "nad[c]": -1.0, "6pgl[c]": 1.0, "nadh[c]": 1.0, "h[c]": 1.0}),
]

CAPS = [0.25, 0.5, 1.0]
MULTS = [1.0, 1.10, 1.25, 1.50]
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}


def add_reactions(rxns, add_ids):
    """Add ONLY the requested reactions from the universe (not the whole universe)."""
    rx2 = dict(rxns)
    for rid, name, ec, rev, stoich in UNIVERSE:
        if rid not in add_ids:
            continue
        rx2["ADD_" + rid] = {"sbml_id": "R_ADD_" + rid, "name": name + " (GAPFILL_CANDIDATE)",
                             "reversible": rev, "lb": -1000.0 if rev else 0.0, "ub": 1000.0,
                             "stoich": stoich, "confidence": "", "ec": ec, "pmid": "",
                             "subsystem": "GAPFILL_CANDIDATE", "gpr": "", "gpr2": "", "protein": "", **BLANK}
    return rx2


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


def max_biomass(metas, rxns, cond, override):
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
            "o2": v[rxn_idx["Ex_o2[e]"]], "co2": v[rxn_idx["Ex_h2co3[e]"]],
            "rubisco": v[rxn_idx["RUBISCO"]] + v[rxn_idx["RUBISCOX"]]}


def scenario_override(cond, ref, glc_cap, rubisco_limit, ci_mode, added):
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0)}
    ov["RUBISCO"] = (0.0, rubisco_limit)
    ov["RUBISCOX"] = (0.0, rubisco_limit)
    if ci_mode == "low":
        ov["Ex_h2co3[e]"] = (-0.1 * abs(ref["co2"]), 0.0)
    elif ci_mode == "zero":
        ov["Ex_h2co3[e]"] = (0.0, 0.0)
    else:
        ov["Ex_h2co3[e]"] = (-1000.0, 0.0)
    for a in added:
        ov["ADD_" + a] = (0.0, 1000.0)  # ensure active
    return ov


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    print("ref:", {k: round(v, 4) for k, v in ref.items()})

    # reference Rubisco+RubiscoX flux for the glucose-assisted state (baseline Target D)
    # (use the model's CBB arm in the WT glucose=0 state as a proxy for the cap fractions)
    rub_ref = ref["rubisco"]  # WT autotrophic Rubisco+Rbx (CO2 fixed)
    # For R25/R10 the cap is relative to the glucose-assisted reference; use WT autotrophic value here.
    R = {"R25": 0.25 * abs(rub_ref), "R10": 0.10 * abs(rub_ref), "R0": 0.0}
    print("rub_ref (WT autotrophic Rubisco+Rbx):", round(rub_ref, 4))

    matrix = []
    for rid in ["R25", "R10", "R0"]:
        for ci in ["full", "low", "zero"]:
            for cap in CAPS:
                ov = scenario_override(cond, ref, cap, R[rid], ci, [])
                b = max_biomass(metas, rxns0, cond, ov)
                matrix.append({"scenario": rid, "ci_mode": ci, "glucose_cap": cap,
                               "max_biomass": round(b, 6) if b is not None else "infeasible",
                               "G0": "yes" if b is not None and b >= ref["biomass"] - 1e-6 else "no",
                               "G10": "yes" if b is not None and b >= 1.10*ref["biomass"] - 1e-6 else "no"})
    f = ["scenario", "ci_mode", "glucose_cap", "max_biomass", "G0", "G10"]
    with (OUT / "gapfill_scenario_matrix.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in matrix:
            w.writerow(r)
    print("=== scenario matrix (cap 1.0) ===")
    for r in matrix:
        if r["glucose_cap"] == 1.0:
            print(f"  {r['scenario']} {r['ci_mode']}: {r['max_biomass']} G0={r['G0']} G10={r['G10']}")

    # GapFill: find minimal reaction additions for R0 + Ci-zero (hardest), then R0+low, R25...
    sols = []
    for rid in ["R0", "R10", "R25"]:
        for ci in ["zero", "low"]:
            # check if already feasible without additions
            b0 = max_biomass(metas, rxns0, cond, scenario_override(cond, ref, 1.0, R[rid], ci, []))
            if b0 is not None and b0 >= ref["biomass"] - 1e-6:
                sols.append({"scenario": rid, "ci_mode": ci, "added": "NONE", "n_added": 0,
                             "max_biomass": round(b0, 6), "G0": "yes", "G10": "yes" if b0 >= 1.10*ref["biomass"] else "no"})
                continue
            # single-reaction gapfill
            found = []
            for (urid, _, _, _, _) in UNIVERSE:
                rx2 = add_reactions(rxns0, [urid])
                ov = scenario_override(cond, ref, 1.0, R[rid], ci, [urid])
                b = max_biomass(metas, rx2, cond, ov)
                if b is not None and b >= ref["biomass"] - 1e-6:
                    found.append((urid, b))
            for urid, b in sorted(found, key=lambda x: -x[1]):
                sols.append({"scenario": rid, "ci_mode": ci, "added": urid, "n_added": 1,
                             "max_biomass": round(b, 6), "G0": "yes",
                             "G10": "yes" if b >= 1.10*ref["biomass"] else "no"})
            if not found:
                # pair gapfill
                pairs = list(itertools.combinations([u[0] for u in UNIVERSE], 2))
                for p in pairs:
                    rx2 = add_reactions(rxns0, list(p))
                    ov = scenario_override(cond, ref, 1.0, R[rid], ci, list(p))
                    b = max_biomass(metas, rx2, cond, ov)
                    if b is not None and b >= ref["biomass"] - 1e-6:
                        found.append(("+".join(p), b))
                for p, b in sorted(found, key=lambda x: -x[1])[:10]:
                    sols.append({"scenario": rid, "ci_mode": ci, "added": p, "n_added": 2,
                                 "max_biomass": round(b, 6), "G0": "yes",
                                 "G10": "yes" if b >= 1.10*ref["biomass"] else "no"})
            if not found:
                sols.append({"scenario": rid, "ci_mode": ci, "added": "NONE_FOUND_<=2", "n_added": None,
                             "max_biomass": "", "G0": "no", "G10": "no"})

    f = ["scenario", "ci_mode", "added", "n_added", "max_biomass", "G0", "G10"]
    with (OUT / "gapfill_minimal_solutions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sols:
            w.writerow(r)
    print("=== gapfill solutions ===")
    for r in sols:
        if r["n_added"] is not None and r["n_added"] <= 1:
            print(f"  {r['scenario']} {r['ci_mode']}: add {r['added']} -> biomass {r['max_biomass']} G10={r['G10']}")


if __name__ == "__main__":
    main()
