# updatedv3.1 Build Report

## Base

- Source: `baseline_v0/minimal_trustworthy_baseline/phase5B6/00_provenance/v3/minimal_trustworthy_baseline_v3_T0.xml`
- SHA256: `5ba88e6dd31d375ea0acbd94442e190b0de9f67d2e61169f8b73380d57becb77` (verified)
- 620 reactions / 575 metabolites (unchanged: no reaction added or removed)

## Confidence standard (2016 iMC507, Methods 2.1)

- 4: biochemically characterized enzyme
- 3: genetic knockout or physiological evidence
- 2: indirect evidence / sequence homology
- 1: gap-filling / network-function reaction with no other evidence

## Applied changes

### 1. Confidence update (comprehensive)

- 330 reactions carry `CONFIDENCE_2026` in reaction notes + master.
- Distribution: score 2 = 281, score 1 = 24, score 3 = 19, score 4 = 6.
- R_NDH2_NATIVE and R_NOX_NATIVE are KEPT; G56 classified both loci as 证据不足暂缓 with no
  exact reaction, so both are set to confidence 1 (gap-filling, no direct evidence) with `G5X_SOURCE=G56`.

### 2. GPR (current RU820 loci) + EC

| reaction_id | GENE_ASSOCIATION_CURRENT | LEGACY_AFE_ASSOCIATION | EC |
|---|---|---|---|
| R_GMHEPAT | NONE (AFE_1406 removed) | AFE_1406 |  |
| R_GMHEPK | RU820_RS06520 | AFE_1406 | 2.7.1.167 |
| R_ARGDC | RU820_RS06805 | AFE_1471 |  |
| R_DB4PS | RU820_RS01450 | AFE_0299 |  |

Equations and bounds for these four are unchanged.

## Counts

- IMPLEMENTED_ACTIONS = 335 (330 confidence + 4 GPR + 1 EC)
- ALREADY_SATISFIED = 0
- UNRESOLVED_NONEXECUTABLE = 1244 (G51-G56 records held/insufficient/non-executable)
- UNEXPLAINED_DIFFS = 0
- QC_ERRORS = 0

## Anti-fabrication QC

- T0 -> updatedv3.1 diff: 330 reaction `<notes>` added, 0 structural changes.
- No reaction added or removed; no equation/bound/stoichiometry change.
- Confidence values come from the G5X 2026 reaction score / 2016 Confidence Level columns.
- GPR values are valid `RU820_RS...` loci; legacy AFE kept only as annotation.
- Only G51-G56 conclusions used; G57-G510 excluded.

## Regression note

No structural change to the network; the only differences are metadata (confidence/GPR/EC).
FBA phenotypes are therefore unchanged from T0. A full FBA re-run was not executed
(no local solver/test runner); the structural diff confirms no phenotype-affecting change.
