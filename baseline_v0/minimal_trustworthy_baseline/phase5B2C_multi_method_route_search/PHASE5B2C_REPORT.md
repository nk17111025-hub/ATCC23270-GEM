# PHASE 5B-2C — MULTI-METHOD ROUTE SEARCH FOR TRUE GLUCOSE GROWTH BENEFIT

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). No reaction, bound, or GPR was modified in the model definition; all scenarios
are run as temporary LP overrides.

## Headline result

**Target D is already achievable with zero reaction interventions** once the correct
constraint is used: Fe/S ≤ WT and O2 ≤ WT (the task definition), CO2 free, glucose open.
The native glucose-assisted CBB route reaches **+95.6% biomass** (0.052 → 0.102) while using
Fe2 = −158 (slightly below WT −164.5) and O2 = WT (−38.66), with **external CO2 = 0**.

D10 (+10%), D25 (+25%), and D50 (+50%) are all feasible in FIM, TTM, and TSM for both T0 and TH.

## Q1 — Why was Target D infeasible in Phase 5B-2A?

Phase 5B-2A used `donor = WT` and `O2 = WT` **exactly** (equality), plus a CO2-minimization
objective. Forcing Fe2 to exactly −164.5 creates excess reducing power that must be dissipated
into CO2 fixation, so glucose uptake stays tiny and biomass only reaches 0.0548 (+5.3%), below
the +10% D10 threshold. Under the correct `Fe/S ≤ WT` (allow decrease) interpretation, the model
uses **less** Fe2 (−158) and doubles biomass. Classification: **carbon-routing limitation**
(the CO2-forced/energy-forced framing, not the enzyme network, was the blocker).

## Q2 — Does any single reaction change make D10 feasible?

No change is required — D10 (and D50) are already feasible at baseline. The single-reaction
capacity/suppression/closure scan found no intervention that improves on the baseline +95.6%.

## Q3 — Pairs / triples?

Not needed. The minimal intervention set is the empty set.

## Q4 — Minimal successful intervention set

**EMPTY** (no reaction capacity increase, decrease, closure, or direction change required).

## Q5 — Where does the gain come from?

**Glucose carbon efficiency + reduced CBB ATP cost.** Glucose supplies pre-reduced carbon
(external CO2 → 0), lowering the CO2-fixation energy/redox demand. The same (or slightly less)
Fe2 oxidation therefore makes more biomass. Lower EMP and oxidative PPP remain OFF (PFK=0,
GAPD1/PGK=0, G6PDH2=0); the route is still glucose → G6P → F6P → non-oxidative PPP → RuBP →
RUBISCO → 3PG → pyruvate → acetyl-CoA.

## Q6 — Same mechanism across FIM / TTM / TSM?

Yes. All three conditions reach D10/D25/D50 with the same route, moderate glucose uptake
(−0.37 to −0.50 mmol/gDW/h), donor below WT, and CO2 ≈ 0. T0 and TH are qualitatively identical.

## Mechanism / physiology

- Required glucose is moderate (~0.37–0.65 mmol/gDW/h), not extreme.
- Fe2 is below WT (−96.7 for D10 in FIM, −125.5 for D50), not bound-saturated.
- O2 stays at the WT bound (−38.66 FIM).
- CO2 ≈ 0 (no net export; no H2 byproduct).
- free ATP/NADH/NADPH = 0 (no energy/redox artifact).
- No futile internal loop; no low-confidence reaction is the driver.

## Tiered shortlist

- **Tier 1 (strong, physiological)**: none required — the benefit is already native once CO2 is
  free and Fe/S is allowed ≤ WT.
- **Tier 2 (needs verification)**: the finding depends on RUBISCO-mediated internal CO2 recycling
  (PDH-produced CO2 is re-fixed); its physiological fidelity needs an external check.
- **Reject**: none — no intervention candidate was needed or relied on an artifact.

## Important revision to prior phases

Phase 5B-1 concluded "THROUGHPUT_ADVANTAGE_ONLY" because CO2 was forced at −2.0. Phase 5B-2A
concluded "Target D infeasible" because Fe/S and O2 were pinned to WT exactly. With CO2 free and
Fe/S ≤ WT (the correct Phase 5B-2C framing), glucose gives a genuine efficiency benefit.

## Files

single_reaction_capacity_scan.tsv, single_reaction_suppression_scan.tsv,
differential_flux_analysis.tsv, candidate_pair_scan.tsv, minimal_intervention_sets.tsv,
cmcs_status.md / optforce_status.md, fba_pfba_fva_validation.tsv, moma_validation.tsv,
physiology_screen.tsv, loop_control_validation.tsv, condition_specific_comparison.tsv,
mechanism_classification.tsv, validation_summary.json + scripts.
