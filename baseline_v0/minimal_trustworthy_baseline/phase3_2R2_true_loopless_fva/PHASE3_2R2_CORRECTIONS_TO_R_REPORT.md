# PHASE 3.2-R2 CORRECTIONS TO PHASE 3.2-R REPORT

## Correction 1 — loopless FVA completion
- Previous wording: "standard vs loopless FVA completed" (Phase 3.2-R).
- Correct wording: Phase 3.2-R produced standard FVA + representative loopless solutions only; it did NOT produce true reaction-wise loopless min/max bounds (loopless_fva_tests = 0). Phase 3.2-R2 performs the bounds calculation.

## Correction 2 — NADTRHD KO compensation
- Previous wording: "NADTRHD KO is compensated by ETC/proton exchange".
- Correct wording: NADTRHD KO induces distributed whole-network rerouting including ETC, folate metabolism, nucleotide metabolism, amino-acid metabolism (glycine/serine/glutamate) and central-carbon reactions; no single artificial rescue route was identified.

## Correction 3 — BDGK dormant capacity
- Previous wording: "BDGK 0..55.7 is true dormant capacity" (based on standard FVA).
- Correct wording: standard FVA gave BDGK [0, 55.7]; true loopless FVA (cycleFreeFlux, fastSNP unavailable with GLPK) confirms BDGK loopless [0, 55.7] — forward-only, no reverse. So "0..55.7 is true loopless forward dormant capacity" is now CONFIRMED by cycleFreeFlux; the reverse (G6P -> glucose + ATP) is NOT loopless-feasible.

## Correction 4 — fastSNP method status
- fastSNP (MILP-based loopless FVA) did not complete with the available GLPK solver (single-reaction loopless FVA > minutes; 6 contexts did not finish within 40 min). cycleFreeFlux (Saa-Nielsen convex loopless formulation) completed all 6 contexts and provides loopless min/max bounds. fastSNP requires a commercial MILP solver (Gurobi/CPLEX, as used in the original 2016 paper).
