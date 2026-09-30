# MINIMAL_TRUSTWORTHY_BASELINE_V1_CHANGELOG

generated_at=2026-09-26T08:10:20+08:00
source_model=baseline_v0\original\mmc3.xml
source_sha256=0dbaacfc07aafb892db0731a517160cf46006bc502b0dedd56f950c14b840246
v1_model=baseline_v0\minimal_trustworthy_baseline\freeze\minimal_trustworthy_baseline_v1.xml
v1_sha256=bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c

## Centrally approved structural changes (3 reactions)

### 1. MACPD (R_MACPD) - disable
- bounds: [0.0, 1000.0] -> [0.0, 0.0]
- reason: malonyl-ACP decarboxylase; Confidence=1; no GPR/PMID; EC 2.3.1.38 inconsistent; closes artificial HCO3E bypass
- rollback: restore UPPER_BOUND=1000.0

### 2. ACOATA (R_ACOATA) - disable reverse
- bounds: [-1000.0, 1000.0] -> [0.0, 1000.0]
- reason: reverse direction closes artificial HCO3E bypass; forward retained as conservative assumption
- rollback: restore LOWER_BOUND=-1000.0

### 3. Htpp (R_Htpp) - disable by default
- bounds: [-1000.0, 1000.0] -> [0.0, 0.0]
- reason: unconstrained reversible Htpp creates free-ATP EGC (Htpp+ATPS5rpp+ATPM); published Table 1 already sets 0/0; no GPR/PMID for alternative direction
- status: DISABLED_BY_DEFAULT_DUE_TO_UNCONSTRAINED_EGC; BIOLOGICAL_PROTON_LEAK_DIRECTION_AND_CAPACITY_UNRESOLVED
- rollback: restore LOWER_BOUND=-1000.0, UPPER_BOUND=1000.0

## Deferred (not repaired in this freeze)
- CYTRED mmc1/mmc3 proton-direction discrepancy
- CA2tpp imbalance
- remaining charge imbalances
- TSM quantitative thiosulfate-uptake limitation
- 508 vs 507 GPR count
- uncertain ETC proton stoichiometries
- uncertain biomass composition
- unfinished genome-wide G5X review