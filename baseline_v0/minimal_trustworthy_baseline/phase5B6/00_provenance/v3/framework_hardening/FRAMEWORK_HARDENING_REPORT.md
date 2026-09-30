# FRAMEWORK HARDENING REPORT — Phase 5B-6R.3

## Verdict

FRAMEWORK_HARDENING_PASS

## Results

- authoritative model identity: PASS
- XML round-trip: PASS
- complete 2016 regression: 4/4
- artifact QC: 7/7
- MILP unit tests: 6/6
- runtime loads actual v3 XML: yes
- new intermediate metabolites supported: yes (unit test 2)

## R0_STRICT definition

Ex_h2co3[e] = (0, 1000) (no external inorganic-carbon uptake, secretion allowed)
Ex_glc-B[e] = (-5.0, -0.0001) (glucose uptake required)
RUBISCO = RUBISCOX = (0, 0)
donor <= WT, O2 <= WT (deterministic pFBA reference)

## Donor/O2 reference

WT biomass 0.0520764. Donor/O2 FVA ranges are NOT unique; deterministic pFBA reference is used.
See wt_donor_o2_reference.tsv.

## R0 sanity (condition semantics only, not biological conclusion)
- OLD_BUGGED_R0: OPTIMAL, biomass -0, glucose -0
- R0_STRICT: OPTIMAL, biomass -0, glucose -0.149014
- R0_CO2_BIDIRECTIONAL_SENSITIVITY: OPTIMAL, biomass -0, glucose -0.149014

## Next

The real universal-network search is deferred until after this clean PASS.
