# PHASE 4.1G-4A — WHY GLUCOSE IS STILL RUBISCO-DEPENDENT

Diagnosis only. The T0/TH candidate models were inspected without any modification
(no reactions added, no bounds changed, no repair).

## 1. Actual dominant glucose route (FIM; TTM/TSM consistent)

glucose → G6P (BDGK) → F6P (PGI1) → **non-oxidative PPP** (TKT1/TKT2/TALA) →
Ru5P → RuBP (PRUK) → **RUBISCO (carboxylase, CO2 fixation)** → 3PG → 2PG (PGM1) →
PEP (ENO) → pyruvate (PYK) → acetyl-CoA (PDH) → TCA / biomass.

Lower EMP (GAPD1/PGK glycolytic direction) and the oxidative PPP (G6PDH2) carry
**zero flux** in every glucose condition. GAPD1/PGK stay in the gluconeogenic
direction (positive flux), i.e. the model is still running gluconeogenesis, not glycolysis.

## 2. First real lower-EMP bottleneck

Step-by-step forced-glycolytic-flux audit (FIM, glucose cap 1):

| step | forcing glycolytic direction | result |
|---|---|---|
| F6P→FDP (PFK) | force forward | feasible (growth 0.208) |
| FDP→GAP+DHAP (FBA) | force forward | feasible |
| GAP→1,3-BPG (GAPD1) | force reverse | feasible at 1e-6 / 1e-4 |
| **1,3-BPG→3PG (PGK)** | force reverse | **INFEASIBLE at all magnitudes** |

The first true bottleneck is **PGK (1,3-BPG → 3PG)**. Forcing even 1e-6 glycolytic
flux makes the model infeasible, and simple cofactor relaxations (free NADH / ATP /
NADPH drain, open proton exchange, free ATPM) do NOT restore it. The block is not a
bound artifact — GAPD1/PGK were already opened to reversible in the candidate model.

## 3. Is NADTRHD essential for the glucose benefit? — NO

NADTRHD KO leaves glucose growth unchanged (0.208). When NADTRHD is removed, the model
simply shifts redox balancing to more gluconeogenic GAPD1 (2.97 → 17.5) and more RUBISCO
(2.66 → 11.39). NADTRHD is a convenient redox shuttle, not a required route.

## 4. Complete ED / gluconate route? — No usable route

- Entner-Doudoroff: reaction-wise **COMPLETE** (G6PDH2, PGL, PGDH, DDGPA), but the
  oxidative entry (G6PDH2) carries zero flux — present but never used.
- oxidative PPP (G6P→6PG→Ru5P+CO2): **INCOMPLETE** — there is no decarboxylating
  6-phosphogluconate dehydrogenase (PGDH is the ED dehydratase).
- non-oxidative PPP: COMPLETE and actively used.
- gluconate pathway: ABSENT. Glyoxylate shunt: PARTIAL (MALS only, no ICL).

## 5. Why RUBISCO KO is infeasible (three layers)

1. **Artifact layer**: RUBISCO KO with forced CO2 uptake (Ex_h2co3 = −2.0) is infeasible
   simply because CO2 has no sink. With CO2 free, RUBISCO KO alone is feasible because the
   model switches to **RUBISCOX** (the RuBP oxygenase) + the glycolate salvage pathway
   (PGLYCP → GLYCH → GLXCL → GLYCK) to keep making 3PG from RuBP.
2. **Structural layer**: the model requires the **RuBP carboxylase/oxygenase system**.
   Closing BOTH RUBISCO and RUBISCOX makes glucose completely unusable (growth 0,
   glucose uptake 0).
3. **Redox layer (the real reason)**: the model is an Fe2 oxidizer whose redox economy is
   Fe2 oxidation (reducing-power SOURCE) balanced by CO2 fixation (reducing-power SINK).
   Lower-EMP glycolysis and the oxidative PPP are reducing-power SOURCES, so the model
   never uses them — they would add NADH on top of the already-excess Fe2-derived NADH.
   Glucose can therefore only be assimilated through the CO2-fixation machinery (as the
   RuBP acceptor), not by direct catabolism. Fe2 oxidation is also irreplaceable
   (Fe2 closed + glucose → infeasible).

## 6. Low-confidence-feature dependence

The phenotype does not depend on a low-confidence reaction. The dominant route uses
well-established enzymes (BDGK, PGI1, TKT/TALA, PRUK, RUBISCO/RUBISCOX, PGM/ENO/PYK/PDH).
NADTRHD (confidence 1, no GPR) carries large flux but is dispensable. The dependence is on
the autotrophic redox architecture, not on any single low-confidence reaction.

## 7. Unresolved model issues (report to CENTRAL REVIEW, not fixed)

- The model cannot generate energy from glucose (Fe2 closed → infeasible); whether this is
  true obligate-autotrophy or a missing heterotrophic energy route is unresolved.
- The forced CO2 uptake (−2.0) creates a spurious "RUBISCO-KO infeasible" result.
- The redox-based refusal of glycolysis is a model property that may or may not reflect
  real ATCC 23270 physiology (13C-glucose labelling would be needed to adjudicate).

## Files

glucose_full_flux_trace.tsv, lower_EMP_block_diagnosis.tsv, glucose_dependency_tests.tsv,
alternative_route_inventory.tsv, rubisco_dependency_mechanism.tsv, validation_summary.json,
diagnose.py (reproducible).
