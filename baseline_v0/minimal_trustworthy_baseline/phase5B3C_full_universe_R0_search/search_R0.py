# -*- coding: utf-8 -*-
"""Full-universe R0 minimal pathway search (BiGG-derived + curated candidates)."""
import csv
import json
import itertools
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3C_full_universe_R0_search")
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}

# curated extras not in BiGG (with EC)
EXTRA = {
    "G6PDH_NAD": {"name": "glucose-6-phosphate dehydrogenase (NAD)", "ec": "1.1.1.363", "rev": False,
                  "stoich": {"g6p-B[c]": -1, "nad[c]": -1, "6pgl[c]": 1, "nadh[c]": 1, "h[c]": 1}},
    "PC": {"name": "pyruvate carboxylase", "ec": "6.4.1.1", "rev": False,
           "stoich": {"pyr[c]": -1, "co2[c]": -1, "atp[c]": -1, "oaa[c]": 1, "adp[c]": 1, "pi[c]": 1, "h[c]": 1}},
    "PFL": {"name": "pyruvate formate lyase", "ec": "2.3.1.54", "rev": False,
            "stoich": {"pyr[c]": -1, "coa[c]": -1, "accoa[c]": 1, "for[c]": 1}},
    "NADHox": {"name": "NADH oxidase", "ec": "1.6.3.4", "rev": False,
               "stoich": {"nadh[c]": -1, "o2[c]": -1, "h[c]": -1, "nad[c]": 1, "h2o2[c]": 1}},
}


def load_bigg():
    d = json.load(open("_bigg_reactions.json", encoding="utf-8"))
    out = {}
    for r in d:
        if r.get("mappable") == "yes":
            out[r["bigg_id"]] = {"name": r["name"], "ec": "", "rev": r.get("reversible") == "True",
                                 "stoich": json.loads(r["stoich"])}
    return out


def build_universe(rxns0):
    metas, _ = build_memory("T0")
    model_mets = set(metas.keys())
    bigg = load_bigg()
    # native reactions by normalized stoichiometry
    def norm(stoich):
        return frozenset((m, round(c, 6)) for m, c in stoich.items() if abs(c) > 1e-9)
    native_norms = set()
    for r, rx in rxns0.items():
        native_norms.add(norm(rx["stoich"]))
        native_norms.add(norm({m: -c for m, c in rx["stoich"].items()}))  # reverse
    universe = {}
    for rid, spec in {**bigg, **EXTRA}.items():
        if rid in ("HEX1", "PGI", "PFK", "FUM", "PPC", "MALS", "PGL", "EDD", "EDA", "TKT1", "TKT2",
                   "TALA", "RPI", "RPE", "GLUSy", "ASPTA", "ALATA_L", "ACONT", "F6PA", "PPA", "PPKr",
                   "CYTBO3_4pp", "SUCDi", "G6PDH2r"):
            continue  # already native-equivalent or not relevant
        stoich = spec["stoich"]
        if not all(m in model_mets for m in stoich):
            continue  # uses a metabolite not in the model
        if norm(stoich) in native_norms:
            continue  # native-equivalent (skip)
        # keep only reactions using model metabolites
        universe[rid] = spec
    return universe


def add_reactions(rxns, ids, universe):
    rx2 = dict(rxns)
    for rid in ids:
        spec = universe[rid]
        rev = spec.get("rev", False)
        rx2["ADD_" + rid] = {"sbml_id": "R_ADD_" + rid, "name": spec["name"] + " (NON_NATIVE)",
                             "reversible": rev, "lb": -1000.0 if rev else 0.0, "ub": 1000.0,
                             "stoich": spec["stoich"], "confidence": "", "ec": spec.get("ec", ""),
                             "pmid": "", "subsystem": "NON_NATIVE", "gpr": "", "gpr2": "", "protein": "", **BLANK}
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


def R0_override(cond, ref, glc_cap, ci_zero):
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
          "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}
    ov["Ex_h2co3[e]"] = (0.0, 0.0) if ci_zero else (-0.1 * abs(ref["co2"]), 0.0)
    return ov


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    universe = build_universe(rxns0)
    ids = list(universe.keys())
    print("non-native universe size:", len(ids))
    print("candidates:", ids)

    # baseline R0 (no additions)
    for ci in [True, False]:
        b = max_biomass(metas, rxns0, cond, R0_override(cond, ref, 1.0, ci))
        print(f"R0 baseline ci_zero={ci}: {b}")

    # search 1..4 reactions
    found = {}
    for k in [1, 2, 3, 4, 5, 6]:
        found[k] = []
        for combo in itertools.combinations(ids, k):
            rx2 = add_reactions(rxns0, list(combo), universe)
            b_ci0 = max_biomass(metas, rx2, cond, R0_override(cond, ref, 1.0, True))
            b_cilow = max_biomass(metas, rx2, cond, R0_override(cond, ref, 1.0, False))
            ok0 = b_ci0 is not None and b_ci0 >= ref["biomass"] - 1e-6
            oklow = b_cilow is not None and b_cilow >= ref["biomass"] - 1e-6
            if ok0 or oklow:
                found[k].append({"reactions": "+".join(sorted(combo)), "n": k,
                                 "ci_zero_biomass": round(b_ci0, 6) if b_ci0 is not None else "infeasible",
                                 "ci_low_biomass": round(b_cilow, 6) if b_cilow is not None else "infeasible"})
        if found[k]:
            print(f"=== {k}-reaction solutions found: {len(found[k])} ===")
            for r in sorted(found[k], key=lambda x: -max(x["ci_zero_biomass"] if isinstance(x["ci_zero_biomass"], float) else 0, x["ci_low_biomass"] if isinstance(x["ci_low_biomass"], float) else 0))[:15]:
                print("  ", r)
        else:
            print(f"no {k}-reaction solution")

    # write minimal_reaction_sets
    with (OUT / "minimal_reaction_sets.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["n", "reactions", "ci_zero_biomass", "ci_low_biomass"])
        for k in [1, 2, 3, 4, 5, 6]:
            for r in found[k]:
                w.writerow([r["n"], r["reactions"], r["ci_zero_biomass"], r["ci_low_biomass"]])

    # save universe for other scripts
    with open("_universe.json", "w", encoding="utf-8") as f:
        json.dump(universe, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
