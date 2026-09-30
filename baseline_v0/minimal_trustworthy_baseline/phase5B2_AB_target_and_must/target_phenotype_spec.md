# Target Phenotype Specification (Phase 5B-2A)

WT reference = autotrophic optimum (glucose = 0, CO2 forced −2.0): biomass 0.052076,
Fe2 −164.509 / tetrathionate −5.216 / thiosulfate −9.128, O2 −38.663 / −15.792 / −15.792.
CO2 in target scenarios is freed (allow uptake, no forced value).

| Target | biomass | Fe/S donor | O2 | CO2 | glucose |
|---|---|---|---|---|---|
| A (conservative mixotrophy) | ≥ WT | ≥ 90% WT | ≤ 120% WT | ≤ 75% WT (−25%) | > 0 |
| B (medium mixotrophy) | ≥ WT | ≥ 70% WT | ≤ 120% WT | ≤ 50% WT (−50%) | > 0 |
| C (strong mixotrophy) | ≥ WT | ≥ 50% WT | unbounded | → 0 | > 0 |
| D10/D25/D50 (fixed-energy gain) | > WT (+10/25/50%) | = WT | = WT | free | > 0 |

Glucose caps scanned: 0.05, 0.1, 0.25, 0.5, 1, 2, 5 mmol/gDW/h.

## Feasibility

| Target | FIM | TTM | TSM | first feasible glucose cap |
|---|---|---|---|---|
| A | yes | yes | yes | 0.1 |
| B | yes | yes | yes | 0.25 |
| C | yes | yes | yes | 0.5 |
| D10/D25/D50 | **no** | **no** | **no** | — |

Target D (biomass gain at exactly matched Fe/S + O2) is infeasible: glucose does not raise
biomass when both energy resources are pinned to the autotrophic values.
