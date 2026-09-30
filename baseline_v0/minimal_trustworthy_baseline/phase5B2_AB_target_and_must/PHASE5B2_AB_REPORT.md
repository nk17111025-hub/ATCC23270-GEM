# PHASE 5B-2A / 5B-2B — TARGET PHENOTYPE SPACE + MUST REACTION DISCOVERY

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). No reactions added, no bounds changed in the model definition; all target
constraints are applied as scenario overrides only.

## Q1 — Which targets are feasible?

| Target | FIM | TTM | TSM |
|---|---|---|---|
| A (conservative) | feasible | feasible | feasible |
| B (medium) | feasible | feasible | feasible |
| C (strong) | feasible | feasible | feasible |
| D (+10/25/50% biomass at fixed Fe/S+O2) | **infeasible** | **infeasible** | **infeasible** |

Targets A/B/C become feasible at low glucose caps (A=0.1, B=0.25, C=0.5 mmol/gDW/h).
Target D is infeasible in every condition: glucose does not raise biomass when both Fe/S
and O2 are pinned to the autotrophic reference values.

## Q2 — FIM: MUST reactions to strengthen glucose → CO2 substitution

- **MUST-ON** (off in WT, required in target): `Ex_glc-B[e]`, `GLCtex`, `GLCtpp`
  (glucose uptake + transport), and `BDGK` (glucose → G6P; WT can also use it via maltose,
  so it is a strong range shift rather than a strict MUST-ON).
- **MUST-UP** (direction flip / level): `PGI1` (F6P → G6P reverse in WT flips to G6P → F6P
  forward in the target) and `Ex_h2co3[e]` (CO2 uptake −2.0 → −1.5…0).
- The CBB arm (`RUBISCO`, `PRUK`, `TKT/TALA`) and respiratory chain show downward range
  shifts, but these are confounded by WT FVA degeneracy (internal cycles), so they are
  recorded as `STRONG_RANGE_SHIFT`, not strict MUST.

## Q3 — TTM/TSM: reactions that make glucose uptake possible

The same core changes are required in TTM and TSM as in FIM: glucose uptake/transport
`MUST-ON`, `PGI1` forward `MUST-UP`, CO2 `MUST-UP`. No sulfur-specific extra MUST reaction
appears at the single-reaction level.

## Q4 — Cross-condition core MUST

Yes. Three universal reactions/features appear in FIM, TTM, and TSM:
1. glucose uptake + transport (MUST-ON),
2. `PGI1` forward (G6P → F6P) direction (MUST-UP),
3. external CO2 decrease (MUST-UP).

## Q5 — Opposite-direction (FIM↑ vs TTM/TSM↓) reactions

None found at the MUST level. The core MUST set is identical across conditions. The only
quantitative difference is the extent of CO2 substitution (FIM can reach CO2 ≈ 0; TTM/TSM
reach ~91–93% saving in pFBA, leaving residual CO2 ≈ −0.14…−0.18).

## Q6 — MUST reactions that depend on high flux / low confidence / abnormal behavior

- The core MUST reactions use only well-annotated, GPR-bearing reactions (`PGI1` AFE_2924,
  `BDGK`; the glucose transporter is still GPR-unresolved by design).
- Glucose uptake required is moderate (~0.33–0.45 mmol/gDW/h), not extreme.
- The donor is reduced to the target's minimum (90/70/50% WT), not bound-saturated.
- No MUST reaction is a confidence-1 or no-GPR reaction; the transport reactions are
  explicitly `NATIVE_GLUCOSE_UPTAKE_GPR_UNRESOLVED`.
- Physiology flags: A/B in FIM are `WATCH` only because CO2 approaches 0; nothing is
  `IMPLAUSIBLE`.

## Q7 — Reaction-priority tiers

- **Tier 1 (strong MUST + physiologically reasonable)**: glucose uptake + transport,
  `PGI1` forward direction (G6P → F6P), external CO2 reduction.
- **Tier 2 (strong MUST + needs thermodynamic/physiological verification)**: the Fe/S
  oxidation-capacity reduction implied by "donor ≥ 90/70/50% WT", and the non-oxidative-PPP
  directionality (TKT/TALA/RPI/RPE range shifts).
- **Tier 3 (mathematically valid but suspicious)**: the large `STRONG_RANGE_SHIFT` set on
  RUBISCO/PRUK/GAPDH/respiration, which is confounded by WT FVA degeneracy and should not be
  read as MUST until a loopless/parsimonious FVA is run.

## Caveats / blockers

- FVA ranges are degenerate (internal loops widen WT ranges); strict MUST-UP/DOWN on the CBB
  and respiratory reactions is therefore uncertain. The MUST-ON and MUST-UP (direction-flip)
  results above are robust because they are tight.
- Target D infeasibility should be re-checked after a loopless-FVA pass; it is consistent
  with Phase 5B-1 (no biomass gain at matched energy).
- No gene-level mapping was performed; reaction-level MUST ≠ gene KO/OE.

## Files

target_phenotype_spec.md, reference_states.tsv, target_feasibility_matrix.tsv,
target_flux_summary.tsv, wt_fva.tsv, target_fva_<cond>_<target>.tsv, must_single_reactions.tsv,
strong_range_shifts.tsv, must_condition_comparison.tsv, must_pair_candidates.tsv,
pfba_target_sanity.tsv, physiology_flags.tsv, validation_summary.json + scripts.
