# -*- coding: utf-8 -*-
"""Phase 5B-5: de novo Rubisco-KO minimum-reaction rescue search (v2 model)."""
import csv
import json
import itertools
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B5")
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}

# Candidate universe: (id, name, EC, reversible, stoich, source)
CAND = [
    ("ME1", "malic enzyme NAD", "1.1.1.38", False, {"mal-L[c]": -1, "nad[c]": -1, "pyr[c]": 1, "co2[c]": 1, "nadh[c]": 1}, "BiGG"),
    ("ME2", "malic enzyme NADP", "1.1.1.40", False, {"mal-L[c]": -1, "nadp[c]": -1, "pyr[c]": 1, "co2[c]": 1, "nadph[c]": 1}, "BiGG"),
    ("PPCK", "PEP carboxykinase ATP", "4.1.1.49", True, {"oaa[c]": -1, "atp[c]": -1, "pep[c]": 1, "co2[c]": 1, "adp[c]": 1}, "BiGG"),
    ("ICL", "isocitrate lyase", "4.1.3.1", False, {"icit[c]": -1, "glx[c]": 1, "succ[c]": 1}, "BiGG"),
    ("PPS", "PEP synthase", "2.7.9.2", False, {"pyr[c]": -1, "atp[c]": -1, "h2o[c]": -1, "pep[c]": 1, "amp[c]": 1, "pi[c]": 1, "h[c]": 2}, "BiGG"),
    ("GND", "6-phosphogluconate dehydrogenase", "1.1.1.44", False, {"6pgc[c]": -1, "nadp[c]": -1, "ru5p-D[c]": 1, "co2[c]": 1, "nadph[c]": 1}, "BiGG"),
    ("NADH16", "NADH dehydrogenase (forward)", "1.6.5.3", False, {"nadh[c]": -1, "q8[c]": -1, "h[c]": -1, "nad[c]": 1, "q8h2[c]": 1}, "BiGG"),
    ("PPDK", "pyruvate phosphate dikinase", "2.7.9.1", False, {"pyr[c]": -1, "atp[c]": -1, "pi[c]": -1, "pep[c]": 1, "amp[c]": 1, "ppi[c]": 1}, "BiGG"),
    ("G6PDH_NAD", "G6PDH NAD", "1.1.1.363", False, {"g6p-B[c]": -1, "nad[c]": -1, "6pgl[c]": 1, "nadh[c]": 1, "h[c]": 1}, "curated"),
    ("PC", "pyruvate carboxylase", "6.4.1.1", False, {"pyr[c]": -1, "co2[c]": -1, "atp[c]": -1, "oaa[c]": 1, "adp[c]": 1, "pi[c]": 1, "h[c]": 1}, "curated"),
    ("PFL", "pyruvate formate lyase", "2.3.1.54", False, {"pyr[c]": -1, "coa[c]": -1, "accoa[c]": 1, "for[c]": 1}, "curated"),
    ("NADHox", "NADH oxidase", "1.6.3.4", False, {"nadh[c]": -1, "o2[c]": -1, "h[c]": -1, "nad[c]": 1, "h2o2[c]": 1}, "curated"),
    ("RPIr", "ribose-5-phosphate isomerase (reverse)", "5.3.1.6", True, {"ru5p-D[c]": -1, "r5p[c]": 1}, "curated-direction"),
    ("G6PDH2r_fwd", "G6PDH oxidative (NADP)", "1.1.1.49", False, {"g6p-B[c]": -1, "nadp[c]": -1, "6pgl[c]": 1, "nadph[c]": 1, "h[c]": 1}, "native-direction"),
]


def build_v2(t):
    metas, rxns = build_memory(t)
    rxns["NADHI"]["lb"] = -1000.0
    return metas, rxns


def add_reactions(rxns, ids):
    rx2 = dict(rxns)
    for rid, name, ec, rev, stoich, src in CAND:
        if rid not in ids:
            continue
        rx2["ADD_" + rid] = {"sbml_id": "ADD_" + rid, "name": name, "reversible": rev,
                             "lb": -1000.0 if rev else 0.0, "ub": 1000.0, "stoich": stoich,
                             "confidence": "", "ec": ec, "pmid": "", "subsystem": "candidate",
                             "gpr": "", "gpr2": "", "protein": "", **BLANK}
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


def R0_override(cond, ref):
    donor = lib.DONOR_EX[cond]
    return {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0), donor: (ref["donor"], 0.0),
            "Ex_o2[e]": (ref["o2"], 0.0), "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns0 = build_v2("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    WT = ref["biomass"]
    cand_ids = [c[0] for c in CAND]

    # candidate universe TSV
    with (OUT / "phase5B5_candidate_universe.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["candidate_id", "reaction_equation", "reversible", "source_universe", "already_in_v2", "excluded", "exclusion_reason"])
        for rid, name, ec, rev, stoich, src in CAND:
            eq = " + ".join(f"{abs(c)} {m}" for m, c in stoich.items() if c < 0) + " -> " + \
                 " + ".join(f"{abs(c)} {m}" for m, c in stoich.items() if c > 0)
            w.writerow([rid, eq, str(rev), src, "no", "no", ""])

    # baseline R0 (no additions)
    b0 = max_biomass(metas, rxns0, cond, R0_override(cond, ref))
    print("R0 baseline (no additions):", b0)

    # exhaustive search 1..6
    found = {}
    for k in [1, 2, 3, 4, 5, 6]:
        found[k] = []
        for combo in itertools.combinations(cand_ids, k):
            rx2 = add_reactions(rxns0, list(combo))
            b = max_biomass(metas, rx2, cond, R0_override(cond, ref))
            if b is not None and b > 1e-6:
                found[k].append(("+".join(sorted(combo)), b))
        if found[k]:
            print(f"=== {k}-reaction solutions: {len(found[k])} ===")
            for combo, b in sorted(found[k], key=lambda x: -x[1])[:10]:
                print(f"  {combo} -> biomass {round(b,6)}")
            break  # found minimum; stop
        else:
            print(f"no {k}-reaction solution")

    # write minimum solutions
    min_k = next((k for k in [1,2,3,4,5,6] if found[k]), None)
    with (OUT / "phase5B5_minimum_rescue_solutions.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["solution_id", "n_additions", "added_reactions", "biomass", "biomass_vs_WT"])
        if min_k is not None:
            for i, (combo, b) in enumerate(sorted(found[min_k], key=lambda x: -x[1])):
                w.writerow([f"S{i+1}", min_k, combo, round(b, 6), round(b/WT, 3)])
        else:
            w.writerow(["NONE", "", "", "", ""])
    print("minimum k:", min_k)


if __name__ == "__main__":
    main()
