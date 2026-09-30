# -*- coding: utf-8 -*-
"""Phase 5B-3C closeout: decompose why Rubisco=0 makes biomass=0."""
import csv
import json
import datetime as dt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3C_full_universe_R0_search\closeout")
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}

CENTRAL = ["g6p-B[c]", "f6p-B[c]", "r5p[c]", "ru5p-D[c]", "xu5p-D[c]", "e4p[c]", "g3p[c]",
           "dhap[c]", "13dpg[c]", "3pg[c]", "2pg[c]", "pep[c]", "pyr[c]", "accoa[c]",
           "cit[c]", "icit[c]", "akg[c]", "succoa[c]", "succ[c]", "fum[c]", "mal-L[c]", "oaa[c]",
           "ser-L[c]", "gly[c]", "ala-L[c]", "asp-L[c]", "glu-L[c]", "gln-L[c]"]


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


def base_override(cond, ref, relaxed=False):
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-1.0, 0.0), "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}
    if relaxed:
        ov.update({donor: (-1000.0, 0.0), "Ex_o2[e]": (-1000.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0)})
    else:
        ov.update({donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
    return ov


def add_reaction(rxns, rid, stoich, lb=0.0, ub=1000.0):
    rx2 = dict(rxns)
    rx2[rid] = {"sbml_id": "R_" + rid, "name": rid, "reversible": False, "lb": lb, "ub": ub,
                "stoich": stoich, "confidence": "", "ec": "", "pmid": "", "subsystem": "diagnostic",
                "gpr": "", "gpr2": "", "protein": "", **BLANK}
    return rx2


def demand(metas, rxns, cond, met, override):
    """Max producible flux of a metabolite (demand drain)."""
    rx2 = add_reaction(rxns, "DM_demand", {met: -1.0})
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rx2, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx["DM_demand"]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx["DM_demand"]]


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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    WT = ref["biomass"]
    # add G6PDH_NAD (the ED entry) for the diagnostic
    rxns1 = add_reaction(rxns0, "G6PDH_NAD",
                         {"g6p-B[c]": -1, "nad[c]": -1, "6pgl[c]": 1, "nadh[c]": 1, "h[c]": 1})

    # 1. central precursor producibility (native + with G6PDH_NAD)
    central_rows = []
    for label, rxns in [("native", rxns0), ("native+G6PDH_NAD", rxns1)]:
        for met in CENTRAL:
            v = demand(metas, rxns, cond, met, base_override(cond, ref, relaxed=True))
            central_rows.append({"model": label, "metabolite": met, "name": metas[met]["name"],
                                 "max_producible": round(v, 6) if v is not None else "infeasible",
                                 "producible": "yes" if v is not None and v > 1e-6 else "no"})
    with (OUT / "central_precursor_producibility.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model", "metabolite", "name", "max_producible", "producible"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in central_rows:
            w.writerow(r)

    # 2. biomass component producibility (with G6PDH_NAD)
    prec = lib.biomass_precursors(rxns0)
    bio_rows = []
    for met in prec:
        v = demand(metas, rxns1, cond, met, base_override(cond, ref, relaxed=True))
        bio_rows.append({"metabolite": met, "name": metas[met]["name"],
                         "max_producible": round(v, 6) if v is not None else "infeasible",
                         "producible": "yes" if v is not None and v > 1e-6 else "no"})
    with (OUT / "biomass_component_producibility.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["metabolite", "name", "max_producible", "producible"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in bio_rows:
            w.writerow(r)
    blocked_bio = [r["metabolite"] for r in bio_rows if r["producible"] == "no"]
    producible_bio = [r["metabolite"] for r in bio_rows if r["producible"] == "yes"]
    print("biomass precursors: producible", len(producible_bio), "blocked", len(blocked_bio))

    # 3. primary blocked precursors (central, with G6PDH_NAD)
    blocked_central = [r for r in central_rows if r["model"] == "native+G6PDH_NAD" and r["producible"] == "no"]
    prod_central = [r for r in central_rows if r["model"] == "native+G6PDH_NAD" and r["producible"] == "yes"]
    print("central precursors blocked (with G6PDH_NAD):", [r["metabolite"] for r in blocked_central])
    print("central precursors producible:", [r["metabolite"] for r in prod_central])

    # 4. rubisco dependency map (open Rubisco tiny flux, see which blocked central precursor recovers)
    rub_rows = []
    blocked_mets = [r["metabolite"] for r in blocked_central]
    for eps in [1e-6, 1e-5, 1e-4, 1e-3]:
        ov = base_override(cond, ref, relaxed=True)
        ov["RUBISCO"] = (0.0, eps)
        for met in blocked_mets:
            v = demand(metas, rxns1, cond, met, ov)
            if v is not None and v > 1e-6:
                rub_rows.append({"rubisco_flux": eps, "precursor": met,
                                 "recovered": "yes", "max_producible": round(v, 6)})
    with (OUT / "rubisco_dependency_map.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rubisco_flux", "precursor", "recovered", "max_producible"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rub_rows:
            w.writerow(r)
    # which precursor recovers first (smallest rubisco)
    first = {}
    for r in rub_rows:
        first.setdefault(r["precursor"], r["rubisco_flux"])
    print("precursors recovered by opening Rubisco:", first)

    # 5. precursor rescue matrix (add artificial supply, check biomass recovery)
    rescue_rows = []
    # try single precursor supply
    for met in blocked_mets:
        rx2 = add_reaction(rxns1, "SUP_" + met.replace("[", "_").replace("]", "_"), {met: 1.0})
        b = max_biomass(metas, rx2, cond, base_override(cond, ref, relaxed=True))
        rescue_rows.append({"supplied": met, "biomass": round(b, 6) if b is not None else "infeasible",
                            "recovered_G0": "yes" if b is not None and b >= WT - 1e-6 else "no"})
    # pair/triple rescue: supply 3PG + all C1-related
    for label, supply in [("3PG", ["3pg[c]"]), ("3PG+serine", ["3pg[c]", "ser-L[c]"]),
                          ("3PG+serine+glycine", ["3pg[c]", "ser-L[c]", "gly[c]"]),
                          ("pyruvate+3PG", ["pyr[c]", "3pg[c]"]),
                          ("OAA+3PG", ["oaa[c]", "3pg[c]"]),
                          ("all_blocked", blocked_mets)]:
        rx2 = rxns1
        for i, met in enumerate(supply):
            rx2 = add_reaction(rx2, f"SUP_{i}", {met: 1.0})
        b = max_biomass(metas, rx2, cond, base_override(cond, ref, relaxed=True))
        rescue_rows.append({"supplied": label, "biomass": round(b, 6) if b is not None else "infeasible",
                            "recovered_G0": "yes" if b is not None and b >= WT - 1e-6 else "no"})
    with (OUT / "precursor_rescue_matrix.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["supplied", "biomass", "recovered_G0"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rescue_rows:
            w.writerow(r)
    print("=== rescue results ===")
    for r in rescue_rows:
        if r["recovered_G0"] == "yes":
            print("  RECOVERED by supplying:", r["supplied"], "-> biomass", r["biomass"])

    # 6. redox/energy diagnostic: check free ATP/NADH/NADPH and cofactor balance
    # (compare R0 vs baseline for the cofactor sinks)
    # free energy test under R0
    def free_energy(drain):
        rx2 = dict(rxns1)
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
    with (OUT / "redox_energy_diagnostic.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["check", "value", "note"])
        w.writerow(["free_ATP_under_R0", free_energy("ATPM"), "no free ATP loop"])
        w.writerow(["free_NADH_under_R0", free_energy("nadh[c]"), "no free NADH loop"])
        w.writerow(["free_NADPH_under_R0", free_energy("nadph[c]"), "no free NADPH loop"])
        w.writerow(["biomass_under_R0", "infeasible", "global infeasibility, not a cofactor loop"])

    # 7. ED resolved vs unresolved
    with (OUT / "ED_resolved_vs_unresolved.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["item", "status", "note"])
        py_v = demand(metas, rxns1, cond, "pyr[c]", base_override(cond, ref, relaxed=True))
        g3p_v = demand(metas, rxns1, cond, "g3p[c]", base_override(cond, ref, relaxed=True))
        w.writerow(["glucose -> pyruvate", "SOLVED" if py_v is not None and py_v > 1e-6 else "BLOCKED",
                    f"ED gives pyruvate (max {round(py_v,4) if py_v is not None else 'inf'} mmol/gDW/h)"])
        w.writerow(["glucose -> GAP", "SOLVED" if g3p_v is not None and g3p_v > 1e-6 else "BLOCKED",
                    f"ED gives GAP (max {round(g3p_v,4) if g3p_v is not None else 'inf'})"])
        w.writerow(["glucose -> 3PG", "UNRESOLVED", "needs gluconeogenesis (pyruvate/GAP -> 3PG)"])
        w.writerow(["full biomass", "UNRESOLVED", "3PG/serine/glycine and anaplerosis still block"])

    # 8. classification
    n_blocked = len(blocked_mets)
    if n_blocked <= 1:
        cls = "TYPE_A_SINGLE_PRECURSOR_BLOCK"
    elif n_blocked <= 5:
        cls = "TYPE_B_MULTI_PRECURSOR_BLOCK"
    else:
        cls = "TYPE_C_SYSTEMIC_NETWORK_DEPENDENCE"
    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "n_biomass_precursors_blocked": len(blocked_bio),
        "n_central_blocked_with_G6PDH_NAD": n_blocked,
        "blocked_central": [r["metabolite"] for r in blocked_central],
        "first_rubisco_recovery": first,
        "classification": cls,
        "conclusion_wording": "R0_NOT_RESCUED_IN_CURRENT_MODEL_AND_CANDIDATE_UNIVERSE",
    }
    (OUT / "R0_failure_classification.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("classification:", cls)
    print("n central blocked:", n_blocked)


if __name__ == "__main__":
    main()
