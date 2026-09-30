# E2 KI-only StrainDesign search

- Authoritative host: `D:\嗜酸氧化亚铁硫杆菌\before final\updatedv3.1_verified\updatedv3.1.xml`
- Host SHA256: `026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3`
- E1 candidates loaded: 2424
- Search candidate direction variants: 2427
- Original FIM+WT baseline: `{"status": "optimal", "biomass": 0.673529948005056, "glucose_uptake": -3.977833303964595, "rubisco_flux": 4.143052406221354, "co2_flux": -2.0, "fe2_flux": -1000.0, "o2_flux": -242.00093172816875}`
- Relaxed CO2 search baseline: `{"status": "optimal", "biomass": 0.7302586898902531, "glucose_uptake": -4.674279012799309, "rubisco_flux": 2.32355316518396, "co2_flux": 0.0, "fe2_flux": -1000.0, "o2_flux": -243.495654128335}`

StrainDesign was called through `SDModule('protect')` and `compute_strain_designs` with SCIP, explicit `ko_cost={}`, unit `ki_cost` for each optional source-supported candidate reaction, and fixed seed `1`. Opposite one-way donor directions are separate cost-1 variants; a reversible source supports one reversible KI. Candidate sources and bounds are reconstructed from the three local donor GEMs. The search caps glucose uptake at the relaxed host baseline and disallows O2 secretion. A small positive Fe2 uptake guard preserves the Fe oxidation condition without setting a high Fe2 fraction. The model-default GLPK precheck was bypassed because it reported infeasible on the large network; StrainDesign rechecks feasibility with selected SCIP.

## Search results

| Scenario | Max cost | Status | Runtime (s) | Solutions | Rubisco cap | Biomass | Glucose | Rubisco | CO2 | Fe2 | O2 | Full-pool glucose-off max | Glucose threshold | Proof | KI-set glucose-off max | On-off delta |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| reference | 0 | optimal | 0.0 | 1 | None | 0.7302586898902741 | -4.674279012799262 | 2.3235531651840797 | 0.0 | -1000.0 | -243.4956541283347 | 0.05374048269487021 | 0.5842069519122025 | True | 0.05207638710273988 | 0.6781823027875342 |
| rubisco_75 | 1 | infeasible | 4.183287620544434 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_75 | 2 | infeasible | 4.081525087356567 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_75 | 3 | infeasible | 4.1593451499938965 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_75 | 4 | infeasible | 4.132860422134399 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_75 | 5 | infeasible | 3.9113759994506836 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_75 | 6 | infeasible | 4.299245357513428 | 0 | 1.74266487388797 | None | None | None | None | None | None | 0.05374048269501412 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 1 | infeasible | 3.990018129348755 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 2 | infeasible | 4.095384120941162 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 3 | infeasible | 4.083797216415405 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 4 | infeasible | 3.9171335697174072 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 5 | infeasible | 4.913416862487793 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_50 | 6 | infeasible | 9.08707308769226 | 0 | 1.16177658259198 | None | None | None | None | None | None | 0.05374048269501471 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 1 | infeasible | 12.486059665679932 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 2 | infeasible | 12.04314136505127 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 3 | infeasible | 11.649184942245483 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 4 | infeasible | 10.411291360855103 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 5 | infeasible | 4.090828895568848 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
| rubisco_25 | 6 | infeasible | 4.1244049072265625 | 0 | 0.58088829129599 | None | None | None | None | None | None | 0.05374048269501886 | 0.5842069519122025 | True | None | None |
- `rubisco_75`: SCIP proved no viable KI set with at most 6 additions.
- `rubisco_50`: SCIP proved no viable KI set with at most 6 additions.
- `rubisco_25`: SCIP proved no viable KI set with at most 6 additions.

All three full-pool preflight LPs were feasible (biomass approximately 0.8073 with every candidate available). Thus these results give a lower bound of **more than 6 KIs** for this candidate pool and phenotype, without establishing a finite minimum. No constrained scenario has a selected KI combination or a solution phenotype. All 18 bounded searches completed as `infeasible`; none timed out.

## Minimum KI sets

- **reference (0 KI):**  —  (status `optimal`).
  Phenotype: biomass=0.7302586898902741, glucose=-4.674279012799262, Rubisco=2.3235531651840797, CO2=0.0, Fe2=-1000.0, O2=-243.4956541283347; glucose-off max biomass=0.05207638710273988, biomass on-off delta=0.6781823027875342.

Scenario checkpoints under `checkpoints/` are authoritative for interrupted or unresolved runs. `time_limit` and `error` statuses are retained and must not be interpreted as proven infeasibility.
