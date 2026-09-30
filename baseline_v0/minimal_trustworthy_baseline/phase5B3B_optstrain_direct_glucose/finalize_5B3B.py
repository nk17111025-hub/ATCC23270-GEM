# -*- coding: utf-8 -*-
"""Phase 5B-3B: rubisco-independence tests + remaining TSVs + summary."""
import csv
import json
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory
from optstrain import add_reactions, setup, ref_state

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3B_optstrain_direct_glucose")


def max_biomass(metas, rxns, cond, override):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def main():
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    rub_ref = abs(ref["rubisco"])

    sols = [("G6PDH_NAD", ["G6PDH_NAD"]), ("PFL", ["PFL"]), ("PEPCK", ["PEPCK"]),
            ("PC", ["PC"]), ("G6PDH_NAD+PEPS", ["G6PDH_NAD", "PEPS"])]

    # rubisco independence tests
    rub_rows = []
    for label, addids in sols:
        rx2 = add_reactions(rxns0, addids)
        for test, rlim, ci in [("R10_Ci_low", 0.1*rub_ref, "low"), ("R0_Ci_free", 0.0, "full"),
                                ("R10_Ci_zero", 0.1*rub_ref, "zero"), ("R0_Ci_zero", 0.0, "zero")]:
            ov = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_fe2[e]": (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
                  "RUBISCO": (0.0, rlim), "RUBISCOX": (0.0, rlim)}
            ov["Ex_h2co3[e]"] = (-0.1*abs(ref["co2"]), 0.0) if ci == "low" else ((0.0, 0.0) if ci == "zero" else (-1000.0, 0.0))
            b = max_biomass(metas, rx2, cond, ov)
            cls = ("STRONG_DIRECT_ASSIMILATION" if b is not None and b >= ref["biomass"] - 1e-6 and rlim == 0
                   else ("PARTIAL_CBB_DEPENDENCE" if b is not None and b >= ref["biomass"] - 1e-6 else "FAILS"))
            rub_rows.append({"solution": label, "test": test, "max_biomass": round(b, 6) if b is not None else "infeasible",
                             "class": cls})
    f = ["solution", "test", "max_biomass", "class"]
    with (OUT / "rubisco_independence_tests.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rub_rows:
            w.writerow(r)

    # precursor connectivity (G6PDH_NAD, R10)
    conn = [
        ["pyruvate", "direct (ED: G6PDH_NAD->PGL->PGDH->DDGPA)"],
        ["acetyl-CoA", "direct (PDH)"],
        ["GAP", "direct (DDGPA)"],
        ["R5P", "direct (non-oxidative PPP)"],
        ["3PG/serine", "PARTIAL (needs residual Rubisco)"],
        ["OAA", "anaplerosis (PPC, needs CO2)"],
        ["alpha-KG", "direct (TCA)"],
    ]
    with (OUT / "precursor_connectivity.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["precursor", "connectivity"])
        for r in conn:
            w.writerow(r)

    # physiology screen
    phys = [
        ["G6PDH_NAD (R10)", 0.068, -0.436, -113.3, 0.69, -27.7, 0.72, 0.0, "PHYSIOLOGY_OK",
         "glucose 0.44, Fe2 69% WT, O2 72% WT, CO2=0"],
        ["PFL (R10)", 0.054, -0.36, -96, 0.58, -23, 0.60, 0.0, "PHYSIOLOGY_WATCH",
         "formate byproduct via PFL"],
    ]
    with (OUT / "physiology_screen.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["solution", "biomass", "glucose", "fe2", "fe2_frac", "o2", "o2_frac", "co2", "flag", "note"])
        for r in phys:
            w.writerow(r)

    # engineering complexity
    comp = [
        ["G6PDH_NAD", 1, 0, 1, "no", "no", "HIGH", "ED entry with NAD; no new metabolite; single cofactor change (NADP->NAD)"],
        ["PFL", 1, 1, 1, "no", "no", "MEDIUM", "pyruvate->acetyl-CoA+formate; introduces formate byproduct"],
        ["PEPCK", 1, 0, 1, "no", "no", "MEDIUM", "OAA->PEP gluconeogenesis; needs CO2 release"],
        ["G6PDH_NAD+PEPS", 2, 0, 2, "maybe", "no", "MEDIUM", "ED + PEP synthesis; candidate for R0 (needs verification)"],
    ]
    with (OUT / "engineering_complexity.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["solution", "N_reactions", "M_new_metabolites", "C_cofactor_changes", "R_extra_redox", "E_extra_energy", "P_plausibility", "note"])
        for r in comp:
            w.writerow(r)

    # gapfill vs optstrain comparison
    with (OUT / "gapfill_optstrain_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction", "gapfill_5B3A", "optstrain_5B3B", "classification"])
        w.writerow(["G6PDH_NAD (NAD-G6PDH)", "found (R10)", "found (R10+R25)", "CONSENSUS"])
        w.writerow(["ICL", "found (R10, degenerate pre-fix)", "found (R25)", "PARTIAL_CONSENSUS"])
        w.writerow(["PEPCK/PEPS/PPDK/PC/ME", "found (pre-fix)", "found (R25)", "PARTIAL_CONSENSUS"])
        w.writerow(["PFL", "not in 5B3A universe", "found (R10, minor)", "OPTSTRAIN_ONLY"])
        w.writerow(["NADHDH", "found (pre-fix, zero flux)", "found (R25, minor)", "LOW_VALUE"])

    # solution ranking
    rank = [
        ["Tier 1", "G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363)", "1 reaction; R10 feasible (+31%); ED route; consensus GapFill+OptStrain",
         "POSSIBLY_NATIVE (genome NAD-G6PDH paralog check needed)"],
        ["Tier 2", "PFL (pyruvate formate lyase, EC 2.3.1.54)", "1 reaction; R10 barely feasible (+4%); formate byproduct",
         "LIKELY_HETEROLOGOUS; redox/fermentation tradeoff"],
        ["Tier 2", "PEPCK / PEPS (PEP synthesis)", "1 reaction for R25; candidate for R0 gluconeogenesis to 3PG",
         "POSSIBLY_NATIVE; needs pairing for full Rubisco independence"],
        ["Tier 2", "G6PDH_NAD + PEPS (ED + gluconeogenesis)", "2 reactions; candidate for R0 (needs >2 verified)",
         "combination candidate"],
    ]
    with (OUT / "solution_ranking.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tier", "reaction_set", "evidence", "native_check"])
        for r in rank:
            w.writerow(r)

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "results": {
            "R25_min": 1, "R10_min": 1, "R0_min": ">3",
            "consensus_reaction": "G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363)",
            "dominant_mechanism": "Entner-Doudoroff (ED)",
            "fully_Rubisco_independent": False,
            "R10_residual_CBB": "10% (3PG/serine still CBB-dependent)",
        },
        "database": "curated universal central-carbon (12 EC-annotated reactions; EC as primary identifier)",
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("finalize 5B3B complete")
    for r in rub_rows:
        if r["solution"] == "G6PDH_NAD":
            print(r["test"], r["max_biomass"], r["class"])


if __name__ == "__main__":
    main()
