# -*- coding: utf-8 -*-
"""Write Phase 5B-3E deliverables."""
import csv
import json
import datetime as dt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3E")

ROWS = [
    ["G6PDH2 / zwf", "G6PDH2", "AFE_2025", "RU820_RS09360", "WP_012536916.1",
     "G510 UniProt + Phase4.1G", "B7J4N6", "6pgl + H + NADPH -> g6p-B + NADP (lb 0)",
     "lb -1000 (opened)", "EXISTING_CORRECTION_CONFIRMED", "retain Phase 4.1G", "EC 1.1.1.49, NADP; oxidative direction opened in Phase 4.1G"],
    ["GAPD1/GAPD2", "GAPD1/GAPD2", "AFE_3251", "RU820_RS15125", "WP_009567621.1",
     "G510 UniProt + Phase4.1G", "B7JBC9", "13dpg + H + NADH -> g3p + NAD + Pi (lb 0)",
     "lb -1000 (opened)", "EXISTING_CORRECTION_CONFIRMED", "retain Phase 4.1G", "type-I GAPDH reversible; glycolytic direction opened"],
    ["PGK", "PGK", "AFE_3250", "RU820_RS15120", "WP_012537677.1",
     "G510 UniProt + Phase4.1G", "B7JBC8", "3pg + ATP -> 13dpg + ADP (lb 0)",
     "lb -1000 (opened)", "EXISTING_CORRECTION_CONFIRMED", "retain Phase 4.1G", "reviewed reversible PGK (GO glycolysis+gluconeogenesis)"],
    ["RpiA", "RPI", "AFE_0629", "RU820_RS03035", "WP_009560936.1",
     "G52 g52_gene_decisions_draft.json", "row 586", "r5p -> ru5p-D (lb 0)",
     "unchanged", "HOLD_INSUFFICIENT_EVIDENCE", "no edit", "G52: maintain 2016 GPR, no formula change, direction unconfirmed (Rhea RHEA:14657 reversible but not independently confirmed)"],
    ["PPC", "PPC", "AFE_1810", "RU820_RS08340", "WP (2016)", "G56 模型动作清单/未决点", "row 1581", "co2 + h2o + pep -> h + oaa + pi (GPR AFE_1810)",
     "unchanged", "SUMMARY_WAS_WRONG_NO_EDIT", "no edit", "G56: keep AFE_1810; AFE_1883 insufficient evidence (don't write OR)"],
    ["PPC candidate", "PPC", "AFE_1883", "RU820_RS08690", "WP_012536819.1", "G56 未决点", "row 1647", "candidate PEPC",
     "not added", "HOLD_INSUFFICIENT_EVIDENCE", "no edit", "G56: AFE_1883 isozyme needs independent evidence"],
    ["GHMT3 / GlyA", "GHMT3", "AFE_0295", "RU820_RS01430", "WP (2016)", "G51/G54 跨基因核查", "GHMT3 row 280", "nad + thf + gly -> co2 + nadh + mlthf + nh4 (GPR AFE_0295 GlyA)",
     "GPR AFE_0295 removed (metadata)", "VERIFIED_GPR_ERROR_FIXED", "remove AFE_0295 from GHMT3 GPR", "GlyA is SHMT (GHMT2), not the glycine-cleavage complex"],
    ["GMHEPAT/GMHEPK", "GMHEPAT/GMHEPK", "AFE_1406", "?", "?", "G55 (MISSING)", "?", "?", "unchanged", "UNRESOLVED", "no edit", "G55 source not found in workspace"],
    ["ARGDC", "ARGDC", "AFE_1417/AFE_1471", "?", "?", "G55 (MISSING)", "?", "?", "unchanged", "UNRESOLVED", "no edit", "G55 source not found"],
    ["DB4PS", "DB4PS", "AFE_1494/AFE_0299", "?", "?", "G55 (MISSING)", "?", "?", "unchanged", "UNRESOLVED", "no edit", "G55 source not found; DB4PS also appears in G51 (needs verification)"],
    ["NADHI / Complex I", "NADHI", "multiple Nuo", "?", "?", "G58 (audit)", "?", "nad + q8h2 -> nadh + q8 (reverse ETC)", "unchanged", "UNRESOLVED", "audit only", "Stage 5: forward NADH oxidation not promoted without evidence"],
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["item", "reaction_id", "old_AFE", "current_RU820", "WP", "G5X_source_file", "G5X_row_or_record",
            "2016_equation", "current_corrected", "verification_status", "proposed_action", "reason"]
    with (OUT / "phase5B3E_source_verification.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in ROWS:
            w.writerow(r)

    with (OUT / "phase5B3E_applied_changes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction", "field_changed", "old_value", "new_value", "source", "verification_status", "biological_reason"])
        w.writerow(["GHMT3", "GPR", "AFE_0295 (GlyA)", "(removed)", "G51/G54", "VERIFIED_GPR_ERROR_FIXED",
                    "GlyA is serine hydroxymethyltransferase (GHMT2), not the glycine-cleavage complex (GHMT3); GPR is metadata-only (not in SBML)"])

    held = [r for r in ROWS if r[9] in ("HOLD_INSUFFICIENT_EVIDENCE", "UNRESOLVED", "SUMMARY_WAS_WRONG_NO_EDIT")]
    with (OUT / "phase5B3E_hold_list.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["item", "reaction_id", "verification_status", "reason"])
        for r in held:
            w.writerow([r[0], r[1], r[9], r[11]])

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "verdict": "CONDITIONAL_PASS",
        "applied_sbml_changes": [],
        "applied_metadata_changes": ["GHMT3 GPR: remove AFE_0295 (metadata-only)"],
        "new_model_sha256": {
            "T0": "3b7d6fda5f5b1b30369f65db44a2036c919401470a43b2be35792207afc50345",
            "TH": "ea404ac16cdeb42b90ea2581f6592d8e7aa28dc4b9553707be9a770f88428933",
        },
        "new_model_identical_to_phase4_1G": True,
        "R0_diagnosis_unchanged": True,
        "key_findings": {
            "verified_model_error_fixed": ["G6PDH2", "GAPD1", "GAPD2", "PGK"],  # already fixed in Phase 4.1G
            "verified_GPR_error_fixed": ["GHMT3 (AFE_0295 GlyA removal)"],
            "summary_was_wrong_no_edit": ["PPC (AFE_1810/AFE_1883)"],
            "hold_insufficient_evidence": ["RPI", "GMHEPAT/GMHEPK", "ARGDC", "DB4PS"],
            "unresolved": ["NADHI"],
        },
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("Phase 5B-3E deliverables written")


if __name__ == "__main__":
    main()
