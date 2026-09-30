# updatedv3.1 — final G51-G56 integration report

## Verdict

`PASS_G51_G56_INTEGRATION`

This PASS is limited to the accepted G51-G56 evidence scope. G57-G510 are not integrated.

## Parent and output

- parent: `minimal_trustworthy_baseline_v3_T0.xml`
- parent SHA256: `5ba88e6dd31d375ea0acbd94442e190b0de9f67d2e61169f8b73380d57becb77`
- output SHA256: `026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3`
- reactions: 620 -> 620
- metabolites: 575 -> 575
- actual stoichiometry/bound diff vs T0: **NONE**

## What was integrated

- all G51-G56 final gene decisions: **1758/1758** stored in `GENE_REVIEW_RECORDS`
- all parsed existing-reaction review records: **410** stored in `ALL_REVIEW_RECORDS`
- unique existing T0 reactions reviewed by G51-G56: **328**
- reactions with numeric 2026 confidence evidence: **240**
- deterministic field-level actions: **35** total: **34 metadata/GPR/EC actions applied**, **1 structural candidate held**

The XML is SBML Level 2, so GPR/confidence/identifier updates are encoded in reaction `<notes>` and mirrored in the master workbook. No SBML Level conversion was performed.

## Important corrections implemented

- FCLT: current GPR `RU820_RS00855`; current EC **4.98.1.1** (legacy 4.99.1.1 not retained as current EC).
- PGM1: `RU820_RS02555` branch retained; `RU820_RS02995` branch remains HOLD; current branch EC **5.4.2.12**; full replacement GPR not fabricated.
- ADNSE: current EC nomenclature **3.13.2.1** recorded; chemistry unchanged.
- NADS1: `RU820_RS02075` withdrawn; NADS2 retains that current locus as a candidate.
- G55 corrections implemented: GMHEPAT/GMHEPK, ARGDC, DB4PS, including GMHEPK EC 2.7.1.167.
- Other deterministic G51/G52/G54 GPR decisions are listed line-by-line in `EXECUTABLE_ACTIONS` / `updatedv3.1_action_audit.tsv`.

## Structural HOLD policy

No G51-G56 candidate was forced into reaction stoichiometry or bounds without a final, network-safe instruction.

`TRDR` is the key case: G52 supplies a balanced corrected candidate equation and score 4 for the catalytic chemistry, but explicitly requires network validation. Applying that candidate alone collapsed the validated FIM autotrophic state, so it is recorded as **HOLD_NOT_APPLIED**, not silently forced into v3.1.

`R_NDH2_NATIVE` and `R_NOX_NATIVE` are **retained**. G56 classifies their exact-reaction assignment as insufficient/HOLD; that does not authorize deletion.

## Confidence conflicts

Six reactions have cross-record score conflicts and are deliberately stored as conflicts rather than guessed:

- CU2tpp [1,2]
- ZN2tpp [1,2]
- S7PI [1,2]
- PGM1 [1,2]
- FDH [1,2]
- Pitpp [1,2]

## Independent QC

- XML parse: PASS
- unique reaction/species IDs: PASS
- all species references valid: PASS
- all LB <= UB: PASS
- unexplained structural diffs: **0**
- G57-G510 data used: **0**
- 12 paired regression checks: all T0 and updatedv3.1 statuses/objectives identical
- glucose-closed benchmark: FIM `0.0520763871`, TTM `0.0520763871`, TSM `0.0520763871`, identical to T0

The raw WT rows in the regression helper leave glucose open and therefore are not the canonical autotrophic benchmark; the glucose-closed rows reproduce the validated v3 autotrophic state.

## Files

- `updatedv3.1.xml` — executable SBML plus G51-G56 reaction metadata
- `updatedv3.1_master.xlsx` — complete 620-reaction state plus all G51-G56 gene/reaction review records and validation
- `updatedv3.1_action_audit.tsv` — exact deterministic field actions and HOLD status
- `regression_summary.tsv/json` — independent T0-vs-v3.1 regression comparison
