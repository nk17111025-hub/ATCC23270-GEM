# PHASE 5B-1 — QUANTIFY THE BENEFIT OF THE GLUCOSE-ASSISTED CBB ROUTE

Models: Phase 4.1G glucose candidate (T0 primary; TH qualitatively identical).
Frozen v1 baseline unchanged.

## Headline result

The glucose-assisted CBB route is **robust and real**, but its biomass benefit is a
**throughput** effect, not an efficiency effect. When the Fe/S donor (and O2) is held at the
autotrophic value, glucose is not taken up at all and biomass does not increase.

## 1. Route robustness (FBA / pFBA / FVA)

Across FVA at 30–100% of maximum biomass (FIM, glucose cap 5), the identified route is
mandatory in every optimal solution:

- BDGK, PRUK, RUBISCO always positive (RUBISCO min 2.64 at 30% → 4.12 at 100%);
- G6PDH2 always 0 (oxidative PPP never used);
- GAPD1/PGK stay gluconeogenic (never glycolytic).

Verdict: **ROBUSTLY_ACTIVE**, not a single-solution artifact.

## 2. Unconstrained growth

FIM: biomass 0.052 → 0.667 at glucose cap 5, while Fe2 uptake rises from −164.5 to −1000
(bound-saturated) and O2 from −38.7 to −242. The growth increase is accompanied by a ~6×
energy-throughput increase.

## 3. Donor-matched (the decisive test)

Fe/S donor fixed at the autotrophic optimum (−164.5 FIM, −5.22 TTM, −9.13 TSM):

- glucose uptake = 0;
- biomass unchanged (0.052076);
- relative gain = 0% at every glucose cap, for FIM/TTM/TSM.

## 4. Oxygen-matched

O2 fixed at autotrophic value (−38.66 FIM): glucose uptake = 0, biomass unchanged.

## 5. Donor + O2 matched

Both fixed: glucose uptake = 0, biomass unchanged.

## 6. Fixed-biomass resource saving

At 30/50/70/90% of autotrophic biomass, the minimized Fe/S, O2, and CO2 uptakes are
**identical** with and without glucose (e.g. 90% biomass: Fe2 −159.75, O2 −37.46, CO2 −2.0 in
both cases). Glucose does not reduce any inorganic resource at fixed biomass.

## 7. Carbon-fixation burden

Per-biomass metrics fall with glucose (CO2/biomass 38.4 → 3.0; RUBISCO/biomass 41.6 → 6.2;
PRUK/biomass 41.6 → 6.2), but this is a consequence of higher biomass. Absolute RUBISCO flux
rises (2.17 → 4.12) and CO2 uptake stays fixed at −2.0. Classification: **MIXED** — per-biomass
CBB burden falls, but total CBB throughput rises; at fixed biomass there is no CO2 saving.

## 8. Energy and redox mechanism

ATP synthase 20.8 → 157.7, respiration (CYTAA31/O2) and Fe2 rise in proportion to biomass.
The benefit is carbon-driven, not energy-driven: the autotrophic model is CO2-carbon-limited,
so extra Fe2 energy cannot raise biomass. Glucose supplies carbon skeletons (via the
PPP→RuBP→RUBISCO route), lifting the carbon limitation and letting the extra Fe/S energy be
converted into biomass.

## 9. Production envelopes

At Fe2 = autotrophic (−164.5), the glc0 and glc5 envelopes are identical (0.052). At higher Fe2
(−246, −329), glucose shifts the frontier (0.113, 0.173 vs 0.052). Glucose shifts the feasible
frontier only when Fe/S energy is allowed to rise.

## 10. Artifacts

Free ATP/NADH/NADPH = 0; growth without electron donor = infeasible; glucose-only with Fe/S
closed = infeasible. Fe2 saturation at high glucose is flagged (bound-saturated growth is not
an efficiency advantage). Conclusions are identical for T0 and TH.

## Overall interpretation

**THROUGHPUT_ADVANTAGE_ONLY.** Under matched Fe/S (or Fe/S+O2) resources glucose gives zero
biomass benefit; the apparent advantage is due entirely to higher Fe/S + O2 energy throughput,
which glucose enables by supplying carbon skeletons that lift the autotrophic CO2 limitation.

## Files

unconstrained_growth_comparison.tsv, donor_matched_comparison.tsv,
oxygen_matched_comparison.tsv, donor_oxygen_matched_comparison.tsv,
fixed_biomass_resource_saving.tsv, carbon_fixation_burden.tsv,
energy_redox_comparison.tsv, production_envelopes.tsv, route_flux_robustness.tsv,
artifact_robustness_tests.tsv, validation_summary.json, quantify_benefit.py.
