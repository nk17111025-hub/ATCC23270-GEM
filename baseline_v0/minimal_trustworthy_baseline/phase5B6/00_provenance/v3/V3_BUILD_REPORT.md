# V3 BUILD REPORT — minimal_trustworthy_baseline_v3

## Verdict

`V3_PROMOTED — EXECUTION BUG REMOVED`

## Parentage

v3 is built only from the promoted Phase 5B-4 `minimal_trustworthy_baseline_v2`:

| file | v2 SHA256 (parent) | verified |
|---|---|---|
| T0 | `97e31dd7b87852e8f639b72408e6e9ea05713dc0a48e99bc77edf087f0919e96` | yes |
| TH | `2a9c86d309d3f1e53e461f5a5df9cad65de974718db95ddc1bbe0945d917be39` | yes |

The v3 XML content equals the v2 XML content (the model reactions and base bounds are
unchanged); the SHA256 of the v3 XML is therefore identical to v2. The v3 change is a
**bound-execution architecture correction**, not a reaction-content change.

## Defect removed

The historical Table-1 FIM/TTM/TSM reaction-specific bounds could overwrite approved
internal-reaction corrections (GAPD1, GAPD2, PGK, G6PDH2, NADHI were re-closed to lb=0, and
BDGK was re-opened to lb=-1000). v3 makes internal-reaction base bounds authoritative.

## New execution rule

`load v3 -> apply environmental condition (exchange-only) -> apply experiment/search override`

Historical Table-1 internal-reaction bounds are never applied. Any attempt to apply one over
v3 raises `BOUND_PRECEDENCE_VIOLATION` (see `scripts/v3_runtime.py`).

## Condition definitions

`v3_condition_definition.tsv` classifies every Table-1 entry as either:

- `ENVIRONMENTAL_CONSTRAINT` — the 28 exchange reactions (Fe2/O2/CO2/NH4/SO4/…).
- `MODEL_INTERNAL_BOUND` — 586 internal reactions, marked `legacy_reference_only`, NOT-APPLIED.

## Effective bounds (all PASS)

| reaction | v3 base | FIM | TTM | TSM | status |
|---|---|---|---|---|---|
| BDGK | (0, 1000) | (0, 1000) | (0, 1000) | (0, 1000) | internal_base_authoritative |
| G6PDH2 | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | internal_base_authoritative |
| GAPD1 | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | internal_base_authoritative |
| GAPD2 | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | internal_base_authoritative |
| PGK | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | internal_base_authoritative |
| NADHI | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | (-1000, 1000) | internal_base_authoritative |

## Regression QC (all PASS)

| test | v3 | expected | result |
|---|---|---|---|
| FIM WT biomass | 0.0520764 | 0.052076 | PASS |
| TTM WT biomass | 0.0520764 | 0.052076 | PASS |
| TSM WT biomass | 0.0520764 | 0.052076 | PASS |
| no-donor (FIM/TTM/TSM) | 0 | 0 | PASS |
| no-carbon | 0 | 0 | PASS |
| glucose-closed | 0.0520764 | 0.052076 | PASS |
| free ATP | 0 | 0 | PASS |
| free NADH | 0 | 0 | PASS |
| free NADPH | 0 | 0 | PASS |

No collateral phenotype change. The WT autotrophic biomass is unchanged (the approved
direction corrections do not affect autotrophic CO2 fixation), as required.

## Native R0 (zero heterologous additions) — new starting point

Rubisco = RubiscoX = 0, glucose open, Fe2 ≤ WT, O2 ≤ WT, no artificial carbon/ATP/redox:

**biomass = 0 (infeasible)**

Under the corrected v3 execution the central-carbon EMP core is flux-capable, but native R0
remains infeasible for the reasons diagnosed in Phase 5B-6 (amino-acid/cofactor biosynthesis
and the NADH/NADPH redox imbalance). This is now the authoritative R0 starting point.

## Promotion

v3 promoted as the sole executable baseline for subsequent Phase 5B-6/5B-7 work.
v1, Phase 4.1G and v2 remain unchanged and are preserved for traceability.
