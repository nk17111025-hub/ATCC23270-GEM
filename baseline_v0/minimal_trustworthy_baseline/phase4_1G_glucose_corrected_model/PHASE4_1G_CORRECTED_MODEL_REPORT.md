# PHASE 4.1G CORRECTED GLUCOSE MODEL — REPORT

Frozen baseline: `minimal_trustworthy_baseline_v1.xml`
(SHA256 `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c` — UNCHANGED).

## 1. Glucokinase locus resolution

**GLUCOKINASE_CANDIDATE_ONLY.** AFE_2841 → RU820_RS13140 → WP_012537424.1 is the current
ROK-family/KEGG K25026 candidate, but it was added as BDGK GPR only in the 2024 model
(the 2016 model had no GPR). No catalytic evidence. BDGK GPR left unresolved.
See `glucokinase_locus_resolution.md`.

## 2. Exact model corrections made (evidence-corrected candidate, v1 untouched)

Bound corrections only (no new GPR, no heterologous genes, no stoichiometry changes):

| reaction | change | reason |
|---|---|---|
| BDGK | lb -1000 → 0 (forward-only) | block reverse ATP-generating glucose formation |
| G6PDH2 | lb 0 → -1000 (open oxidative direction) | reaction is G6PDH EC 1.1.1.49 (oxidative PPP); reductive-only was a constraint |
| GAPD1 | lb 0 → -1000 | type-I GAPDH reversible; open G3P→1,3-BPG |
| GAPD2 | lb 0 → -1000 | same gene (AFE_3251), NADP-dependent |
| PGK | lb 0 → -1000 | reviewed bidirectional PGK; open 1,3-BPG→3PG |

PFK already carried the correct GPR (AFE_1807 / PfkB) — no change.

Added (C1): extracellular/periplasmic beta-D-glucose + provisional native uptake
(exchange + outer porin + inner transporter), with **no transporter GPR**
(`NATIVE_GLUCOSE_UPTAKE_GPR_UNRESOLVED`). Two energetics variants:

- **T0**: facilitated (non-energized) glucose transport.
- **TH**: H+-coupled glucose symport.

## 3. Baseline regression (glucose uptake closed)

FIM / TTM / TSM all return **0.052076** for both T0 and TH — identical to the frozen v1
reference. Free ATP / NADH / NADPH = 0 (no new energy/redox cycle). HCO3E KO stays
infeasible. Fe/S oxidation unchanged. The bound corrections did not perturb autotrophic growth.

## 4. Growth vs glucose uptake (FIM; TTM/TSM qualitatively identical)

| glucose cap | growth | glucose uptake | Fe2 uptake | O2 uptake | CO2 uptake |
|---|---|---|---|---|---|
| 0 | 0.052 | 0 | -164.5 | -38.7 | -2.0 |
| 0.05 | 0.060 | -0.050 | -175.1 | | -2.0 |
| 0.1 | 0.068 | -0.100 | -185.7 | | -2.0 |
| 0.25 | 0.091 | -0.250 | -217.6 | | -2.0 |
| 0.5 | 0.130 | -0.500 | -270.7 | | -2.0 |
| 1 | 0.208 | -1.000 | -376.9 | | -2.0 |
| 2 | 0.365 | -2.000 | -589.3 | | -2.0 |
| 5 | 0.667 | -3.93 | -1000 (bound) | | -2.0 |

Biomass increases monotonically with glucose. Fe2 oxidation increases (not decreases) —
glucose supplies carbon while Fe/S oxidation supplies energy/ATP.

## 5. T0 vs TH comparison

Quantitatively near-identical (e.g. FIM cap 5: T0 = 0.667 vs TH = 0.664). The transport
energetics choice does not change the qualitative conclusion.

## 6. Fe/S oxidation retained — YES

Fe2 / tetrathionate / thiosulfate uptake remains nonzero (and increases) in every glucose
condition.

## 7. Direct lower glycolysis active — NO

GAPD1 and PGK remain in the gluconeogenic direction (positive flux) even with reversible
bounds. The optimal solution does not route glucose through G3P→1,3-BPG→3PG.

## 8. Oxidative PPP active — NO

G6PDH2 flux is ~0; glucose does not enter the oxidative PPP / Entner-Doudoroff branch.

## 9. RUBISCO dependence — unchanged

Glucose-assisted growth is entirely RUBISCO-dependent. RUBISCO KO is **infeasible** even at
glucose cap 5. The earlier glucose → PPP → RuBP → RUBISCO → 3PG detour **does not disappear**
after the correction: it is the model's optimal route, not an artifact of the one-way bounds.
Inorganic-carbon demand does not decrease (CO2 uptake stays at the fixed -2.0 bound).

## 10. NADTRHD dependence

NADTRHD (confidence 1, no GPR) carries growing flux with glucose (0.95 → 9.31 mmol/gDW/h at
cap 5), i.e. substantial redox-balancing load at high glucose uptake. Flagged (see below).

## 11. Artifact flags

- NADTRHD large flux (up to 9.3) at high glucose caps — redox-balancing, confidence 1 / no GPR.
- Glucose "growth" is RUBISCO/CO2-fixation dependent, not heterotrophic; the CO2-closed test
  still grows (0.71) because RUBISCO recycles internally produced CO2 — not biological heterotrophy.
- Fe2 uptake saturates its -1000 bound at high glucose (energy-limited solution).
- No new free ATP/NADH/NADPH cycle (all 0). No glyoxylate-shunt activation (MALS=0,
  GLXCL/GLYCK/GLXR ≈ 0).

## 12. Final status

**PHASE4_1G_PASS_WITH_TRANSPORT_UNCERTAINTY**

Qualitative phenotype reproduced: glucose availability → increased biomass while Fe/S
oxidation remains active. Transport mechanism/GPR remains unresolved (T0 vs TH sensitivity).

## 13. Files produced

- glucokinase_locus_resolution.md
- glucose_model_corrections.tsv
- minimal_trustworthy_baseline_v1_glucose_corrected_T0.xml / _TH.xml
- glucose_baseline_regression.tsv
- glucose_uptake_sensitivity.tsv
- glucose_flux_route_analysis.tsv
- glucose_artifact_tests.tsv
- glucose_watch_reaction_fluxes.tsv
- phase4_1G_validation_summary.json
- build_glucose_model.py / run_glucose.py / phase4a_lib.py (reproducible)

## Key result for CENTRAL REVIEW

Opening the gluconeogenic GAPDH/PGK and reductive G6PDH bounds did **not** switch glucose
to glycolysis or the oxidative PPP. The model remains an obligate-autotroph-style network
that uses glucose as a CO2-acceptor (RuBP) for RUBISCO, with Fe/S oxidation supplying energy.
Whether this reflects real biology (obligate autotrophy with glucose supplementation) or a
remaining model limitation (over-reversible RUBISCO / incomplete heterotrophic glucose
catabolism) requires CENTRAL REVIEW adjudication. No transporter GPR was assigned and no
strain was designed.
