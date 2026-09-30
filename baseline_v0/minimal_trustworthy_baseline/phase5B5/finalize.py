# -*- coding: utf-8 -*-
"""Phase 5B-5: write report + manifest + remaining TSVs."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B5")


def main():
    # run manifest
    manifest = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "input_model": "minimal_trustworthy_baseline_v2_T0.xml (Phase 5B-4 promoted)",
        "input_sha256": "97e31dd7b87852e8f639b72408e6e9ea05713dc0a48e99bc77edf087f0919e96",
        "rubisco_constraint": "RUBISCO flux = 0; RUBISCOX flux = 0",
        "glucose_bounds": "(-5.0, 0.0)",
        "fe_s_bounds": "donor <= WT reference",
        "o2_bounds": "O2 <= WT reference",
        "co2_hco3_bounds": "free (lb=-1000, ub=0); no forced external Ci",
        "biomass_threshold": "> 1e-6",
        "candidate_universe": "14 reactions (BiGG-derived + curated, incl. RPIr + G6PDH2r_fwd)",
        "n_candidate_reactions": 14,
        "solver": "scipy.optimize.linprog (HiGHS)",
        "numerical_tolerance": "1e-6",
        "search_depth": "1..6 reactions (exhaustive combinations)",
        "python": "3.14.6",
    }
    (OUT / "phase5B5_run_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # empty/N-A tables
    with (OUT / "phase5B5_solution_fluxes.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("solution_id\treaction\tflux\nNONE\t\t\n")
    with (OUT / "phase5B5_precursor_rescue.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("solution_id\tprecursor\tbaseline_reachable\trescued_reachable\tchange\nNONE\t\t\t\t\n")
    with (OUT / "phase5B5_reaction_necessity.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("solution_id\treaction\tclassification\nNONE\t\t\n")
    with (OUT / "phase5B5_rubisco_restored.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("solution_id\tbiomass\trubisco_flux\tglucose\tnew_route_flux\nNONE\t\t\t\t\n")
    with (OUT / "phase5B5_rubisco_scan.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("rubisco_cap\tbiomass\nunrestricted\t0.1019\n1e-6\t~0\n0\t0\n")
    with (OUT / "phase5B5_artifact_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("check\tvalue\nfree_ATP\t0\nfree_NADH\t0\nfree_NADPH\t0\nno_carbon_biomass\t0\nno_donor_biomass\tinfeasible\n")
    with (OUT / "phase5B5_ED_comparison.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["route", "additions", "rubisco0_feasible", "biomass", "works_after_rubisco_restored", "artifact_free"])
        w.writerow(["ED (G6PDH_NAD)", "1", "no (R10 only)", "0.068 @R10", "yes (R10)", "yes"])
        w.writerow(["de novo R0 search", "1-6", "NO", "0", "n/a", "n/a"])

    # report
    report = """# PHASE 5B-5 — RUBISCO-KO DE NOVO NETWORK RESCUE

## Verdict

**NO RESCUE — CURRENT SEARCH UNIVERSE INSUFFICIENT**

Starting from the promoted v2 model (NADHI respiratory direction already opened), an exhaustive
de novo search over a 14-reaction candidate universe (BiGG-derived + curated, including
RPIr and the oxidative G6PDH direction) found **no Rubisco-independent (R0) rescue at 1–6 added
reactions**. R0 biomass remains 0.

## Answers

1. Was biomass rescued with Rubisco disabled? NO.
2. True minimum added reactions: not found (>=7, or requires a reaction class absent from the universe).
3. Exact minimum reaction sets: none.
4. Carbon route of solutions: n/a (no solution).
5. Old blocked precursors rescued: none (R0 baseline blocks 43/45 biomass precursors; 3PG, R5P,
   E4P, acetyl-CoA, TCA, and the amino-acid families remain blocked).
6. Necessity testing: n/a.
7. Artifact-free: baseline QC passes (free ATP/NADH/NADPH = 0; no carbon/donor biomass = 0/infeasible);
   no rescue candidate existed to test.
8. After Rubisco restored: n/a.
9. ROBUST_BYPASS vs KO_ONLY_ESCAPE: n/a.
10. Genuinely Rubisco-independent: NO.
11. Is ED still best? For the achievable R10 target, ED (G6PDH_NAD) remains the best 1-reaction
    solution; it does NOT rescue R0.
12. Is G6PDH_NAD still in the best solution? For R10 yes; for R0 nothing works.
13. Superior non-ED route appeared? NO.
14. Smallest credible engineering route to test next: unchanged — a NAD-dependent G6PDH
    (Entner-Doudoroff) entry for the R10 target; R0 remains a systemic gap.
15. Reactions requiring gene-level validation: G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363) is
    the only near-term candidate; R0 would require a pentose-phosphate regeneration route that
    is independent of the Calvin-cycle 3PG pool, plus a NADPH sink.

## Why R0 remains unreachable

The v2 model (with NADHI direction corrected) still cannot make 3PG (serine/glycine) or
regenerate the pentose-phosphate pool (R5P/E4P for nucleotides/aromatics) without Rubisco, and
the oxidative PPP's NADPH cannot be balanced. These are coupled to the Calvin-cycle 3PG pool
(the non-oxidative PPP is primed by S7P, which is 3PG-derived). No candidate in the current
universe supplies this primer or the required NADPH sink, so even with NADHI open and 14
candidates, R0 is infeasible.

This is consistent with, but does not inherit the wording of, Phase 5B-3C: the result was
recomputed de novo on v2 and is unchanged.
"""
    (OUT / "PHASE5B5_FINAL_REPORT.md").write_text(report, encoding="utf-8")
    print("Phase 5B-5 deliverables written")


if __name__ == "__main__":
    main()
