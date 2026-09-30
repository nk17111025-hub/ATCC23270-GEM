# -*- coding: utf-8 -*-
"""Phase 5B-4: write route-search results + comparison + report."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B4")


def main():
    # GapFill results (unchanged from Phase 5B-3A, since v2 phenotype is unchanged)
    with (OUT / "phase5B4_gapfill_results.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["scenario", "min_additions", "added_reactions", "biomass", "note"])
        w.writerow(["R25", "0", "NONE", "0.0994 (G10)", "native already feasible"])
        w.writerow(["R10", "1", "G6PDH_NAD (EC 1.1.1.363)", "0.068 (+31%)", "ED entry; unchanged from 5B-3A"])
        w.writerow(["R0", "NOT_FOUND", ">3 (infeasible even with all candidates)", "0", "NADHI opening does not rescue"])
    with (OUT / "phase5B4_optstrain_results.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["scenario", "min_interventions", "best_sets", "note"])
        w.writerow(["R25", "1", "G6PDH_NAD / PEPCK / PC / ICL / PFL (1-reaction)", "unchanged"])
        w.writerow(["R10", "1", "G6PDH_NAD (primary), PFL (weak)", "unchanged"])
        w.writerow(["R0", "NOT_FOUND", "none", "infeasible with <=6 reactions; NADHI open does not help"])
    with (OUT / "phase5B4_full_universe_R0_results.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["scenario", "result", "note"])
        w.writerow(["R0-CiZero", "infeasible", "even with all 12 candidates + NADHI open"])
        w.writerow(["R0-CiLow", "infeasible", "unchanged from Phase 5B-3C"])
    with (OUT / "phase5B4_route_classification.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction", "class", "note"])
        w.writerow(["NADHI", "N1 (native, direction corrected)", "respiratory NADH oxidation enabled; safe but no phenotype change"])
        w.writerow(["G6PDH_NAD", "H (heterologous candidate)", "still the R10 minimum; not native"])
        w.writerow(["PPC (GPR AFE_1883)", "N1 (native, GPR corrected)", "metadata only; not an engineering intervention"])
    with (OUT / "phase5B4_old_vs_new.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["question", "phase5B-3C old model", "phase5B-4 v2", "changed"])
        rows = [
            ["native glucose route", "CBB (glucose->PPP->RuBP->Rubisco->3PG)", "CBB (same)", "NO"],
            ["maximum biomass (donor/O2<=WT)", "0.1019", "0.1019", "NO"],
            ["direct lower glycolysis", "not used (PFK/GAPD=0)", "not used (same)", "NO"],
            ["oxidative PPP", "not used (G6PDH2=0)", "not used (same)", "NO"],
            ["NADH disposal", "limited (NADHI reverse-ETC only)", "enabled (NADHI reversible)", "YES (capability)"],
            ["first blocked precursor", "3PG", "3PG", "NO"],
            ["Rubisco dependence", "essential", "essential", "NO"],
            ["R25 min intervention", "0", "0", "NO"],
            ["R10 min intervention", "G6PDH_NAD (1)", "G6PDH_NAD (1)", "NO"],
            ["R0 feasible", "no", "no", "NO"],
            ["best added reaction(s)", "G6PDH_NAD", "G6PDH_NAD", "NO"],
            ["G6PDH_NAD still useful", "yes (R10)", "yes (R10)", "NO"],
            ["ED still preferred", "yes (R10)", "yes (R10)", "NO"],
            ["smallest heterologous gap", "1 reaction (R10)", "1 reaction (R10)", "NO"],
        ]
        for r in rows:
            w.writerow(r)
    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "verdict": "v2 promoted (minimal_trustworthy_baseline_v2)",
        "NADHI_QC": "PASS (no free-energy artifact)",
        "key_result": "NADHI direction correction is safe but does NOT change glucose phenotype, Rubisco dependence, or R0 feasibility",
        "phase5B3C_overturned": [],
        "phase5B3C_upheld": ["R0 infeasible", "3PG first block", "CBB route dominant", "G6PDH_NAD = R10 minimum"],
    }
    (OUT / "phase5B4_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("Phase 5B-4 finalize complete")


if __name__ == "__main__":
    main()
