# -*- coding: utf-8 -*-
"""Phase 5B-3A: write remaining TSVs + validation summary."""
import csv
import json
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory
from gapfill import add_reactions, setup, ref_state, scenario_override, UNIVERSE

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3A_gapfill_direct_glucose")


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

    # reaction universe summary
    with (OUT / "reaction_universe_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "name", "EC", "reversible", "equation", "source"])
        for rid, name, ec, rev, stoich in UNIVERSE:
            eq = " + ".join(f"{abs(c)} {m}" for m, c in stoich.items() if c < 0) + " -> " + \
                 " + ".join(f"{abs(c)} {m}" for m, c in stoich.items() if c > 0)
            w.writerow([rid, name, ec, str(rev), eq, "curated central-carbon (EC-annotated)"])

    # rubisco dependency tests for the G6PDH_NAD solution (R10)
    rx2 = add_reactions(rxns0, ["G6PDH_NAD"])
    dep = []
    tests = [
        ("R10_G6PDH_NAD", {"RUBISCO": (0.0, 0.1*rub_ref), "RUBISCOX": (0.0, 0.1*rub_ref), "Ex_h2co3[e]": (-1000.0, 0.0)}),
        ("RUBISCO0_Ci_free", {"RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0)}),
        ("R10_Ci_zero", {"RUBISCO": (0.0, 0.1*rub_ref), "RUBISCOX": (0.0, 0.1*rub_ref), "Ex_h2co3[e]": (0.0, 0.0)}),
        ("RUBISCO0_Ci_zero", {"RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)}),
    ]
    for label, extra in tests:
        ov = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_fe2[e]": (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0)}
        ov.update(extra)
        b = max_biomass(metas, rx2, cond, ov)
        dep.append({"solution": "G6PDH_NAD", "test": label,
                    "max_biomass": round(b, 6) if b is not None else "infeasible",
                    "biomass_vs_WT": (round(b/ref["biomass"], 3) if b is not None else "n/a"),
                    "class": "RUBISCO_DEPENDENT" if (b is None or b < ref["biomass"]) and "RUBISCO0" in label else ("DIRECT_ASSIMILATION" if b is not None and b >= ref["biomass"] else "n/a")})
    f = ["solution", "test", "max_biomass", "biomass_vs_WT", "class"]
    with (OUT / "rubisco_dependency_tests.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in dep:
            w.writerow(r)

    # direct precursor accessibility (from the R10 G6PDH_NAD pFBA solution)
    # glucose reaches: pyruvate (ED, direct), acetyl-CoA (PDH), 3PG (partly via 10% Rubisco)
    acc = [
        ["pyruvate", "YES (ED: G6PDH_NAD->PGL->PGDH->DDGPA)", "direct"],
        ["acetyl-CoA", "YES (PDH)", "direct"],
        ["GAP", "YES (DDGPA produces GAP)", "direct"],
        ["3PG", "PARTIAL (residual 10% Rubisco required)", "CBB-dependent"],
        ["serine/glycine", "PARTIAL (needs 3PG -> 10% Rubisco)", "CBB-dependent"],
        ["R5P", "YES (non-oxidative PPP)", "direct"],
        ["E4P", "YES (non-oxidative PPP)", "direct"],
        ["OAA", "PARTIAL (PPC needs CO2)", "anaplerosis"],
        ["alpha-KG", "YES (TCA via acetyl-CoA)", "direct"],
    ]
    with (OUT / "direct_precursor_accessibility.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["precursor", "accessibility", "route_class"])
        for r in acc:
            w.writerow(r)

    # existing model overlap
    overlap = [
        ["G6PDH_NAD", "TRUE_NETWORK_ADDITION_CANDIDATE",
         "model has only NADP-dependent G6PDH2 (EC 1.1.1.49); NAD-dependent G6PDH (EC 1.1.1.363) absent"],
        ["ICL", "TRUE_NETWORK_ADDITION_CANDIDATE", "MALS present but isocitrate lyase absent (glyoxylate shunt incomplete)"],
        ["ME_NADP/ME_NAD", "TRUE_NETWORK_ADDITION_CANDIDATE", "malic enzyme absent"],
        ["PC", "TRUE_NETWORK_ADDITION_CANDIDATE", "pyruvate carboxylase absent (only ACCOAC)"],
        ["PEPCK/PEPS", "TRUE_NETWORK_ADDITION_CANDIDATE", "PEP synthesis from pyruvate/OAA absent (only PPC PEP->OAA)"],
        ["NADHDH", "POTENTIAL_MODEL_OMISSION", "NADHI exists but only reverse-ETC direction (NAD+->NADH); forward NADH oxidation unidirectional"],
    ]
    with (OUT / "existing_model_overlap.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["candidate", "classification", "note"])
        for r in overlap:
            w.writerow(r)

    # solution ranking
    rank = [
        ["Tier 1", "G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363)", "1 reaction enables R10 (Rubisco<=10%), biomass +31%",
         "enables Entner-Doudoroff pathway with NAD; moderate glucose; PHYSIOLOGY_OK"],
        ["Tier 2", "NADHDH (NADH:Q oxidoreductase forward)", "alone insufficient (zero flux under R10); needs combination",
         "redox-direction candidate; needs pairing with a carbon route"],
        ["Tier 2", "PEPCK / PEPS (PEP synthesis)", "candidate for R0 gluconeogenesis to 3PG",
         "needs >2 reactions for full Rubisco independence"],
    ]
    with (OUT / "solution_ranking.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tier", "reaction_set", "evidence", "note"])
        for r in rank:
            w.writerow(r)

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "method": "minimal-reaction gapfill (curated EC-annotated central-carbon universe, 8 candidates)",
        "results": {
            "R25": "feasible with 0 additions",
            "R10": "feasible with 1 addition: G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363)",
            "R0": "needs >2 additions (not found in universe)",
            "min_reactions": {"R25": 0, "R10": 1, "R0": ">2"},
            "key_addition": "G6PDH_NAD -> enables Entner-Doudoroff (glucose -> pyruvate + GAP) with NAD",
            "fully_Rubisco_independent": False,
            "R10_residual_Rubisco": "10% (3PG/serine still CBB-dependent)",
        },
        "note": "candidate universe is curated (8 EC-annotated reactions); not the full BioCyc PGDB (1313 reactions) because compound-ID mapping was not automated",
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("finalize 5B3A complete")
    for r in dep:
        print(r["test"], r["max_biomass"], r["class"])


if __name__ == "__main__":
    main()
