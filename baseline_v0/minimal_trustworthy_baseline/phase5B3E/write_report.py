# -*- coding: utf-8 -*-
"""Write the Phase 5B-3E verification report."""
from pathlib import Path

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3E")

REPORT = """# PHASE 5B-3E — G5X-VERIFIED CENTRAL-CARBON STRUCTURAL CORRECTION

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c` (UNCHANGED).

## Verdict

**CONDITIONAL PASS** — only ONE correction is independently verified (GHMT3/GlyA GPR), and it is
a metadata (Table-1 GPR) change that does not alter the SBML. The two headline structural
hypotheses in the prompt (RPI directionality, PPC GPR) are **not supported** by the G5X sources.
The new versioned models are therefore byte-identical to the Phase 4.1G corrected models.

## Source re-verification results

| item | reaction | verification | action |
|---|---|---|---|
| G6PDH2 / zwf (AFE_2025) | G6PDH2 | EXISTING_CORRECTION_CONFIRMED | retain Phase 4.1G (lb -1000) |
| GAPD1/GAPD2 (AFE_3251) | GAPD1/GAPD2 | EXISTING_CORRECTION_CONFIRMED | retain Phase 4.1G (lb -1000) |
| PGK (AFE_3250) | PGK | EXISTING_CORRECTION_CONFIRMED | retain Phase 4.1G (lb -1000) |
| RpiA (AFE_0629) | RPI | HOLD_INSUFFICIENT_EVIDENCE | no edit (direction unconfirmed) |
| PPC (AFE_1810 / AFE_1883) | PPC | SUMMARY_WAS_WRONG_NO_EDIT | no edit (keep AFE_1810) |
| GHMT3 / GlyA (AFE_0295) | GHMT3 | VERIFIED_GPR_ERROR_FIXED | remove AFE_0295 from GHMT3 GPR |
| GMHEPAT / GMHEPK (AFE_1406) | GMHEPAT/GMHEPK | UNRESOLVED | G55 source missing |
| ARGDC (AFE_1417 / AFE_1471) | ARGDC | UNRESOLVED | G55 source missing |
| DB4PS (AFE_1494 / AFE_0299) | DB4PS | UNRESOLVED | G55 source missing |
| NADHI / Complex I | NADHI | UNRESOLVED | audit only (Stage 5) |

## Key evidence findings

### PPC (contradicts the prompt)

G56 (`G56_模型动作清单.tsv`, `G56_未决点与验证建议.tsv`) states: AFE_1810 (RU820_RS08340) is
"保留现状" (old model already associates PPC with AFE_1810; database annotation alone cannot
rewrite the reaction or GPR). AFE_1883 (RU820_RS08690) is "证据不足暂缓" (PPC exists but the
2016 GPR is AFE_1810; AFE_1883 isozyme capability needs independent evidence, "暂不写OR").
The prompt's hypothesis (remove AFE_1810 as AckA-like, add AFE_1883) is rejected.

### GHMT3 / GlyA (confirms the prompt)

G51/G54 (`G54_跨基因核查记录.md`, `G54_阶段验收报告.md`, `G51_组内验收说明.md`) state: 2016
GHMT3 = glycine-cleavage system (`nad + thf + gly -> co2 + nadh + mlthf + nh4`, EC 2.1.2.1)
with GPR AFE_0295 / GlyA. AFE_0295 = RU820_RS01430, annotated serine hydroxymethyltransferase
(SHMT), which is GHMT2 (serine<->glycine), NOT the glycine-cleavage complex. Conclusion: the
GHMT3 GPR is a mismatch; the normal SHMT reaction (GHMT2) remains correctly associated with
AFE_0295.

### RPI (does not support the prompt)

G52 (`g52_gene_decisions_draft.json`) states: AFE_0629 = RpiA (RU820_RS03035, WP_009560936.1);
"维持 2016 历史 GPR 候选…本轮不提升分数、不改式"; "方向缺少独立实验证实". The Rhea chemistry
(RHEA:14657) is reversible, but the G5X reviewer explicitly holds the direction, so no RPI bound
change is applied.

## Why the new models are byte-identical to Phase 4.1G

The only verified correction (GHMT3 GPR removal) lives in the mmc1.xls Table 1 annotation, which
the SBML XML does not carry. RPI and PPC (the only SBML-visible candidates) were not supported.
Therefore no LOWER_BOUND / UPPER_BOUND / stoichiometry edit is warranted, and the two new
versioned models have the same SHA256 as the Phase 4.1G corrected models.

## Stage 9 — R0 re-diagnosis

Because the SBML is unchanged, the Phase 5B-3C R0 decomposition is reproduced exactly: R0
(Rubisco=0) remains infeasible; 3PG is still the first blocked precursor; R5P/E4P remain
independently blocked; GAP->1,3-BPG remains blocked by NADH disposal; acetyl-CoA/TCA remain
blocked. The Phase 5B-3C conclusion (TYPE_C_SYSTEMIC_NETWORK_DEPENDENCE,
R0_NOT_RESCUED_IN_CURRENT_MODEL_AND_CANDIDATE_UNIVERSE) is unchanged.

## Answers to the required questions

Q1. Genuine 2016 model-structure problems: G6PDH2 (reductive-only), GAPD1/2 (gluconeogenic-only),
    PGK (gluconeogenic-only) — all already corrected in Phase 4.1G. RPI and PPC are NOT verified
    errors; NADHI is a direction concern but unresolved.
Q2. R0 remains Rubisco-dependent: YES (unchanged).
Q3. 3PG still the first blocked precursor: YES.
Q4. R5P/E4P still independently blocked: YES.
Q5. GAP->1,3-BPG still blocked by NADH disposal: YES.
Q6. Acetyl-CoA still blocked: YES.
Q7. Phase 5B-3C conclusion materially changes: NO.
Q8. True heterologous ED completion still required: YES.
Q9. Remaining engineering gaps: 3PG/serine supply + NADH sink, R5P/E4P pentose balance,
    acetyl-CoA/TCA anaplerosis (unchanged from Phase 5B-3C).

## Files

PHASE5B3E_G5X_VERIFICATION_REPORT.md, phase5B3E_source_verification.tsv,
phase5B3E_applied_changes.tsv, phase5B3E_hold_list.tsv,
minimal_trustworthy_baseline_v1_glucose_G5Xverified_T0.xml / _TH.xml,
validation_summary.json, read_g5x.py, finalize.py.
"""


def main():
    (OUT / "PHASE5B3E_G5X_VERIFICATION_REPORT.md").write_text(REPORT, encoding="utf-8")
    print("report written")


if __name__ == "__main__":
    main()
