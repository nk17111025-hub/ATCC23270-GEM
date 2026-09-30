# -*- coding: utf-8 -*-
"""Write provenance metadata and the remaining required (empty/n-a) output tables."""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path

import scipy, numpy

import p5b6_common as C

OUT = C.OUT
TZ = timezone(timedelta(hours=8))


def now():
    return datetime.now(TZ).isoformat(timespec="seconds")


def main():
    metas, rxns = C.build_v2("T0")
    ref = C.reference_state(metas, rxns, C.COND)

    # environment.json
    env = {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "solver": "scipy.optimize.milp (HiGHS)",
        "os": platform.platform(),
        "cwd": str(C.ROOT),
        "generated_at": now(),
    }
    (OUT / "00_provenance" / "environment.json").write_text(
        json.dumps(env, ensure_ascii=False, indent=2), encoding="utf-8")

    # model_reference.tsv
    v2_t0 = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4" / "minimal_trustworthy_baseline_v2_T0.xml"
    v2_th = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4" / "minimal_trustworthy_baseline_v2_TH.xml"
    C.write_tsv(OUT / "01_inputs" / "model_reference.tsv",
                ["model", "path", "sha256", "n_reactions", "n_metabolites", "objective"],
                [["v2_T0", str(v2_t0), C.sha256(v2_t0), str(len(rxns)), str(len(metas)), "Ex_bio[e]"],
                 ["v2_TH", str(v2_th), C.sha256(v2_th), str(len(rxns)), str(len(metas)), "Ex_bio[e]"]])

    # constraint_snapshot.tsv
    donor = C.lib_DONOR() if hasattr(C, "lib_DONOR") else "Ex_fe2[e]"
    import phase4a_lib as lib
    donor = lib.DONOR_EX[C.COND]
    C.write_tsv(OUT / "01_inputs" / "constraint_snapshot.tsv",
                ["constraint", "value"],
                [["RUBISCO", "0"], ["RUBISCOX", "0"],
                 ["glucose", "(-5.0, 0.0)"],
                 [donor, f"({ref['donor']:.6g}, 0.0)"],
                 ["Ex_o2[e]", f"({ref['o2']:.6g}, 0.0)"],
                 ["Ex_h2co3[e]", "(-1000.0, 0.0)"],
                 ["biomass_threshold", "> 1e-6"],
                 ["condition", C.COND]])

    # stage_status.tsv
    stages = [
        ["01_inputs", "completed", now(), "v2 identity verified"],
        ["02_universe_raw", "completed", now(), "1313 BioCyc reactions"],
        ["03_universe_mapping", "completed", now(), "138 mapped"],
        ["04_universe_clean", "completed", now(), "116 balanced + 27 curated"],
        ["05_R0_search", "completed", now(), "INFEASIBLE (no rescue)"],
        ["06_solution_validation", "partial", now(), "redox bottleneck localised"],
        ["07_rubisco_restore", "skipped", now(), "no bypass exists"],
        ["08_KO_search", "skipped", now(), "no bypass exists"],
        ["09_OE_search", "skipped", now(), "no bypass exists"],
        ["10_robustness", "skipped", now(), "no accepted solution"],
        ["12_pareto", "skipped", now(), "no designs"],
        ["13_final", "completed", now(), "PARTIAL verdict"],
    ]
    C.write_tsv(OUT / "00_provenance" / "stage_status.tsv",
                ["stage", "status", "timestamp", "note"], stages)

    # remaining required tables (n/a because no complete rescue)
    na_rows = [["n/a", "no complete Rubisco-independent rescue found", "", "", ""]]
    C.write_tsv(OUT / "05_R0_search" / "phase5B6_minimal_H_solutions.tsv",
                ["solution_id", "n_additions", "added_reactions", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "05_R0_search" / "phase5B6_alternative_solutions.tsv",
                ["solution_id", "added_reactions", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "06_solution_validation" / "artifact_qc" / "phase5B6_artifact_qc.tsv",
                ["solution_id", "test", "result", "note"], na_rows)
    C.write_tsv(OUT / "06_solution_validation" / "flux" / "phase5B6_flux_maps.tsv",
                ["solution_id", "reaction", "flux", "note"], na_rows)
    C.write_tsv(OUT / "06_solution_validation" / "fva" / "phase5B6_FVA.tsv",
                ["solution_id", "reaction", "min", "max", "note"], na_rows)
    C.write_tsv(OUT / "06_solution_validation" / "cofactor_balance" / "phase5B6_cofactor_balance.tsv",
                ["solution_id", "cofactor", "production", "consumption", "net", "main_contributors"],
                [["diagnostic", "NAD(P)H", "excess", "limited (proton-coupled NADHI)", "positive",
                  "NADHI reverse blocked; free sink restores biomass 0.08"]])
    C.write_tsv(OUT / "07_rubisco_restore" / "phase5B6_rubisco_restored.tsv",
                ["solution_id", "biomass", "rubisco_flux", "bypass_flux", "note"], na_rows)
    C.write_tsv(OUT / "07_rubisco_restore" / "rubisco_scan" / "phase5B6_rubisco_scan.tsv",
                ["solution_id", "rubisco_bound", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "08_KO_search" / "phase5B6_KO_designs.tsv",
                ["solution_id", "KO_reactions", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "09_OE_search" / "phase5B6_OE_FSEOF.tsv",
                ["solution_id", "reaction", "classification", "note"], na_rows)
    C.write_tsv(OUT / "10_robustness" / "T0_TH" / "phase5B6_T0_TH_validation.tsv",
                ["solution_id", "transport", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "10_robustness" / "energy_tiers" / "phase5B6_energy_tier_validation.tsv",
                ["solution_id", "fe_tier", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "10_robustness" / "glucose_dependency" / "phase5B6_glucose_dependency.tsv",
                ["solution_id", "glucose", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "12_pareto" / "phase5B6_multi_intervention_designs.tsv",
                ["design_class", "H", "KO", "OE", "biomass", "note"], na_rows)
    C.write_tsv(OUT / "12_pareto" / "phase5B6_pareto_front.tsv",
                ["design", "H", "KO", "OE", "biomass", "note"], na_rows)

    # master manifest
    manifest = {
        "phase": "5B-6",
        "project": "Acidithiobacillus ferrooxidans ATCC 23270",
        "model": "minimal_trustworthy_baseline_v2",
        "v2_T0_sha256": C.V2_T0_SHA256,
        "v2_TH_sha256": C.V2_TH_SHA256,
        "condition": C.COND,
        "solver": "scipy.optimize.milp (HiGHS)",
        "generated_at": now(),
        "verdict": "PARTIAL",
        "stages": {s[0]: {"status": s[1], "note": s[3]} for s in stages},
        "model_repair": "v2 direction corrections overridden by Table-1 FIM bounds; applied as V2_REPAIR overrides",
        "key_result": "143-reaction universe searched by MILP GapFill -> INFEASIBLE; redox/cofactor bottleneck (NADHI reverse proton-coupled/capped)",
    }
    (OUT / "00_provenance" / "phase5B6_master_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
