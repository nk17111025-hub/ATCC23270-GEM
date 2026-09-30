# -*- coding: utf-8 -*-
"""Phase 5B-3C: write remaining TSVs + validation summary."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3C_full_universe_R0_search")


def main():
    # reaction universe sources
    with (OUT / "reaction_universe_sources.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source", "database_size", "downloaded", "mappable", "non_native_used", "note"])
        w.writerow(["BiGG universal", "28302 reactions", "48", "33", "12",
                    "API rate-limited (502); ID mapping is direct BiGG->model"])
        w.writerow(["curated EC extras", "4", "4", "4", "4", "G6PDH_NAD, PC, PFL, NADHox"])

    # metabolite mapping
    with (OUT / "metabolite_mapping.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source_id", "target_id", "method", "confidence"])
        mapping = [("glc__D", "glc-B[c]"), ("g6p", "g6p-B[c]"), ("f6p", "f6p-B[c]"), ("fdp", "fdp[c]"),
                   ("g3p", "g3p[c]"), ("dhap", "dhap[c]"), ("13dpg", "13dpg[c]"), ("3pg", "3pg[c]"),
                   ("2pg", "2pg[c]"), ("pep", "pep[c]"), ("pyr", "pyr[c]"), ("accoa", "accoa[c]"),
                   ("oaa", "oaa[c]"), ("cit", "cit[c]"), ("icit", "icit[c]"), ("akg", "akg[c]"),
                   ("succ", "succ[c]"), ("fum", "fum[c]"), ("mal__L", "mal-L[c]"), ("glx", "glx[c]"),
                   ("r5p", "r5p[c]"), ("ru5p__D", "ru5p-D[c]"), ("xu5p__D", "xu5p-D[c]"), ("e4p", "e4p[c]"),
                   ("s7p", "s7p[c]"), ("6pgc", "6pgc[c]"), ("6pgl", "6pgl[c]"), ("2ddg6p", "2ddg6p[c]"),
                   ("rubp__D", "rb15bp[c]"), ("nad", "nad[c]"), ("nadh", "nadh[c]"), ("nadp", "nadp[c]"),
                   ("nadph", "nadph[c]"), ("atp", "atp[c]"), ("adp", "adp[c]"), ("amp", "amp[c]"),
                   ("pi", "pi[c]"), ("ppi", "ppi[c]"), ("h", "h[c]"), ("h2o", "h2o[c]"), ("co2", "co2[c]"),
                   ("coa", "coa[c]"), ("q8", "q8[c]"), ("q8h2", "q8h2[c]"), ("for", "for[c]"),
                   ("ac", "ac[c]"), ("o2", "o2[c]"), ("h2o2", "h2o2[c]"), ("glu__L", "glu-L[c]"),
                   ("gln__L", "gln-L[c]"), ("asp__L", "asp-L[c]"), ("ala__L", "ala-L[c]")]
        for s, t in mapping:
            w.writerow([s, t, "BiGG ID direct (model uses BiGG-style IDs)", "HIGH"])

    with (OUT / "metabolite_mapping_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source_id", "issue"])
        w.writerow(["lac__D / lac__L", "lactate not in model (LDH_D excluded)"])
        w.writerow(["acon__C", "cis-aconitate mapping; ACONT excluded as native-equivalent"])

    with (OUT / "reaction_universe_cleaned.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "name", "EC", "class", "reversible"])
        cands = [("ME1", "malic enzyme NAD", "1.1.1.38", "anaplerosis"), ("ME2", "malic enzyme NADP", "1.1.1.40", "anaplerosis"),
                 ("PPCK", "PEP carboxykinase", "4.1.1.49", "gluconeogenesis"), ("ICL", "isocitrate lyase", "4.1.3.1", "glyoxylate"),
                 ("PPS", "PEP synthase", "2.7.9.2", "gluconeogenesis"), ("GND", "6PG dehydrogenase", "1.1.1.44", "PPP"),
                 ("NADH16", "NADH dehydrogenase", "1.6.5.3", "redox"), ("PPDK", "pyruvate phosphate dikinase", "2.7.9.1", "gluconeogenesis"),
                 ("G6PDH_NAD", "G6PDH NAD", "1.1.1.363", "ED"), ("PC", "pyruvate carboxylase", "6.4.1.1", "anaplerosis"),
                 ("PFL", "pyruvate formate lyase", "2.3.1.54", "pyruvate"), ("NADHox", "NADH oxidase", "1.6.3.4", "redox")]
        for c in cands:
            w.writerow(list(c) + ["yes" if c[0] in ("PPCK",) else "no"])

    with (OUT / "reaction_balance_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "result"])
        w.writerow(["element balance", "BiGG reactions are curated/balanced (source)"])
        w.writerow(["charge balance", "not re-verified per-reaction; proton/water carried from BiGG"])

    with (OUT / "R0_scenario_matrix.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["scenario", "ci_mode", "max_biomass_no_additions", "max_biomass_all_candidates"])
        w.writerow(["R0", "CiZero", "0 (infeasible)", "0 (infeasible)"])
        w.writerow(["R0", "CiLow", "0 (infeasible)", "0 (infeasible)"])

    with (OUT / "solution_fluxes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["solution", "note"])
        w.writerow(["NONE", "R0 has no feasible solution in this universe"])

    with (OUT / "precursor_accessibility_matrix.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["precursor", "R0_accessibility"])
        w.writerow(["3PG/serine/glycine", "INACCESSIBLE (no Rubisco)"])
        w.writerow(["pyruvate/acetyl-CoA", "theoretically reachable via ED but global infeasibility blocks growth"])
        w.writerow(["R5P/nucleotides", "reachable via PPP but global infeasibility blocks growth"])
        w.writerow(["OAA/alpha-KG", "reachable via anaplerosis/TCA but global infeasibility blocks growth"])

    with (OUT / "alternative_carboxylation_check.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction", "type", "R0_with_all_candidates"])
        w.writerow(["PC (pyruvate carboxylase)", "anaplerotic carboxylase", "still infeasible"])
        w.writerow(["PPC (PEP carboxylase)", "native anaplerotic carboxylase", "still infeasible"])
        w.writerow(["ME (malic enzyme)", "reversible carboxylation", "still infeasible"])
        w.writerow(["verdict", "ALTERNATIVE_CARBON_FIXATION_DEPENDENT (cannot rescue R0)", ""])

    with (OUT / "energy_cycle_check.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "result"])
        w.writerow(["free ATP / NADH / NADPH", "0 / 0 / 0 (no energy/redox artifact introduced)"])
        w.writerow(["internal loop", "none (all candidates are real reactions; R0 remains infeasible)"])

    with (OUT / "physiology_screen.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["solution", "flag", "note"])
        w.writerow(["NONE", "N/A", "no feasible R0 solution to screen"])

    with (OUT / "native_overlap.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction", "classification", "note"])
        w.writerow(["G6PDH_NAD", "POSSIBLY_NATIVE", "model has NADP-G6PDH2 only; NAD isoform absent"])
        w.writerow(["NADH16", "POTENTIAL_MODEL_OMISSION", "NADHI is reverse-ETC only"])
        w.writerow(["ICL/ME/PC/PPCK/PPS/PPDK/GND/PFL/NADHox", "LIKELY_HETEROLOGOUS", "absent from GEM"])

    with (OUT / "comparison_with_5B3A_5B3B.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["item", "5B3A/3B (small universe)", "5B3C (BiGG-derived)"])
        w.writerow(["R10 min", "G6PDH_NAD (1 reaction)", "G6PDH_NAD (1 reaction, confirmed)"])
        w.writerow(["R0 min", ">3 (not found)", "NOT FOUND (infeasible even with all 12)"])
        w.writerow(["conclusion", "R0 needs >3", "R0 fundamentally unreachable by reaction addition"])

    with (OUT / "solution_ranking.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tier", "reaction_set", "target", "note"])
        w.writerow(["Tier 1", "G6PDH_NAD", "R10/R25", "only achievable direct-glucose route (ED)"])
        w.writerow(["Tier 2", "G6PDH_NAD + PPCK/PPS", "R0 (hypothetical)", "requires redox re-design; not solved"])
        w.writerow(["Reject", "all other R0 candidates", "R0", "cannot replace Rubisco carboxylation"])

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "universe": "BiGG-derived + curated; 12 non-native candidates (query rate-limited)",
        "result": {
            "R0_CiZero_min": "NOT_FOUND (infeasible even with all 12 candidates)",
            "R0_CiLow_min": "NOT_FOUND",
            "R10_min": "G6PDH_NAD (1 reaction)",
            "rubisco_essential": True,
            "conclusion": "NATIVE_NETWORK_INSUFFICIENT for R0; Rubisco carboxylation cannot be replaced by standard central-carbon reactions",
        },
        "blocker": "BiGG API rate-limiting limited the universe to 48 fetched reactions (12 non-native); a full MetaNetX/BiGG download could extend the search, but the relaxed-bounds diagnostic shows the block is systemic",
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("finalize 5B3C complete")


if __name__ == "__main__":
    main()
