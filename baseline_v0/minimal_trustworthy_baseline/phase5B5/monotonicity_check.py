# -*- coding: utf-8 -*-
"""Phase 5B-5 monotonicity sanity check: ALL 14 candidates + R0."""
import json
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory
from search_rescue import CAND, add_reactions, setup, ref_state, R0_override, BLANK

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B5")


def build_v2(t):
    metas, rxns = build_memory(t)
    rxns["NADHI"]["lb"] = -1000.0
    return metas, rxns


def max_biomass(metas, rxns, cond, override):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def demand(metas, rxns, cond, met, override):
    rx2 = dict(rxns)
    rx2["DM"] = {"sbml_id": "DM", "name": "demand", "reversible": False, "lb": 0.0, "ub": 1000.0,
                 "stoich": {met: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                 "gpr": "", "gpr2": "", "protein": "", **BLANK}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rx2, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx["DM"]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx["DM"]]


def free_energy(metas, rxns, cond, drain):
    rx2 = dict(rxns)
    ov = {r: (0.0, 0.0) for r in rx2 if r.startswith("Ex_")}
    ov.update({"MACPD": (0.0, 0.0), "ACOATA": (0.0, 1000.0), "Htpp": (0.0, 0.0),
               "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)})
    if drain == "ATPM":
        ov["ATPM"] = (0.0, 1000.0)
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective="ATPM", override=ov, with_table1=False)
    else:
        d = "DM_drain"
        rx2[d] = {"sbml_id": d, "name": "drain", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {drain: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                  "gpr": "", "gpr2": "", "protein": "", **BLANK}
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective=d, override=ov, with_table1=False)
    return 0.0 if fl is None else float(obj)


def main():
    metas, rxns0 = build_v2("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    rx2 = add_reactions(rxns0, [c[0] for c in CAND])
    ov = R0_override(cond, ref)

    bio = max_biomass(metas, rx2, cond, ov)
    print("ALL-14 biomass:", round(bio, 7) if bio is not None else "INFEASIBLE")

    targets = ["3pg[c]", "r5p[c]", "e4p[c]", "accoa[c]", "oaa[c]", "akg[c]"]
    reach = {}
    for m in targets:
        v = demand(metas, rx2, cond, m, ov)
        reach[m] = (round(v, 6) if v is not None else "infeasible")
        print(f"  {m}: {reach[m]}")

    fe = {
        "free_ATP": free_energy(metas, rx2, cond, "ATPM"),
        "free_NADH": free_energy(metas, rx2, cond, "nadh[c]"),
        "free_NADPH": free_energy(metas, rx2, cond, "nadph[c]"),
    }
    print("artifact tests:", fe)

    # NADH/NADPH balance: net production under R0 max-biomass pFBA
    # (report via a pFBA-like solution is complex; report free-energy + reachability as the evidence)
    result = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "all_14_biomass": round(bio, 7) if bio is not None else "infeasible",
        "reachability": reach,
        "artifact_tests": fe,
    }
    decision = "CANDIDATE_UNIVERSE_PROVEN_INSUFFICIENT" if (bio is None or bio <= 1e-6) else "CANDIDATE_UNIVERSE_SUFFICIENT_BUT_MINIMUM_IS_GREATER_THAN_6"
    result["decision"] = decision
    (OUT / "phase5B5_monotonicity_check.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("DECISION:", decision)


if __name__ == "__main__":
    main()
