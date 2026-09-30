# PHASE 5A — TCA COMPLETION HYPOTHESIS TEST

Hypothesis: closing the canonical 2-oxoglutarate → succinyl-CoA break in the TCA cycle
lets glucose contribute more effectively to biomass.

H0 = Phase 4.1G glucose candidate model.
H1 = H0 + one hypothetical OGDH reaction
`2-oxoglutarate + CoA + NAD+ → succinyl-CoA + CO2 + NADH`
(`akg[c] + coa[c] + nad[c] → succoa[c] + co2[c] + nadh[c]`, lb 0/1000, no GPR,
labelled HYPOTHETICAL_ENGINEERING_REACTION).

No other reaction, bound, or GPR was changed; the frozen v1 baseline was not touched.

## Result: OGDH IS NEVER USED

In every condition (FIM/TTM/TSM; glucose caps 0–5), the hypothetical OGDH flux is
**exactly 0.0**. H1 growth is numerically identical to H0 in every case. The added reaction
does not change biomass, glucose uptake, Fe/S uptake, O2/CO2 exchange, Rubisco/RubiscoX,
NADTRHD, PDH, or any TCA flux.

| condition | glucose cap | H0 growth | H1 growth | OGDH flux |
|---|---|---|---|---|
| FIM | 0 | 0.052076 | 0.052076 | 0 |
| FIM | 1 | 0.208306 | 0.208306 | 0 |
| FIM | 5 | 0.666610 | 0.666610 | 0 |
| TTM/TSM | 5 | 0.833222 | 0.833222 | 0 |

## Why (mechanism)

The model is an Fe2-oxidizing chemolithoautotroph with **excess reducing power**
(Fe2 oxidation drives reverse electron transport → NADH/NADPH). The oxidative TCA
(including OGDH) is a **NADH source**, so the model never runs it — it would add reducing
power on top of the already-excess Fe2-derived NADH. The model instead uses the TCA
fragments biosynthetically:

- citrate → isocitrate → 2-oxoglutarate (ICDHyr) → **glutamate / amino acids**;
- succinyl-CoA (SUCOAS) → lysine / heme;
- fumarate → malate → OAA (MDH) → aspartate / citrate synthase.

The 2-OG → succinyl-CoA gap is therefore **not a bottleneck**: the cell does not need a
complete oxidative TCA cycle, so a hypothetical OGDH is simply unused.

## Task-by-task conclusions

1. **Baseline regression** — PASS. H1 without glucose matches H0 (0.052076), free
   ATP/NADH/NADPH = 0 (no new energy/redox cycle).
2. **Glucose sensitivity** — H1 identical to H0 at every cap and condition.
3. **TCA closure** — NO_MEANINGFUL_TCA_FLUX for OGDH (flux 0). The existing TCA halves
   remain separate biosynthetic branches.
4. **Rubisco dependency** — unchanged. RUBISCO KO → 0.648, RUBISCOX KO → 0.742,
   double KO → 0, identical for H0 and H1.
5. **Route shift** — NO_ROUTE_CHANGE. Lower EMP / oxidative PPP / TCA fluxes unchanged.
6. **Biomass benefit** — none (Δbiomass = 0 at every point; H1 = H0 exactly).
7. **Redox** — unchanged. OGDH (a NADH source) is not used, so it neither helps nor
   worsens the reducing-power state.
8. **Artifact tests** — PASS. Free ATP/NADH/NADPH = 0; no new loop; glucose-only growth
   with Fe/S closed remains infeasible (Fe/S still required).

## Overall interpretation

**TCA_COMPLETION_NO_MAJOR_EFFECT**

Adding OGDH does not close a functional TCA cycle in the FBA solution and does not change
glucose-assisted growth. The limiting factor is not the missing OGDH reaction but the
autotrophic redox economy (Fe2 oxidation = reducing-power source, CO2 fixation = sink),
which makes all reducing-power-generating catabolic steps (glycolysis, oxidative PPP, OGDH)
unused.

## Files

H0_H1_growth_comparison.tsv, TCA_flux_comparison.tsv, rubisco_dependency_comparison.tsv,
glucose_route_shift.tsv, redox_comparison.tsv, artifact_tests.tsv, validation_summary.json,
build_h1.py + test_ogdh.py (+ phase4a_lib.py, build_glucose_model.py).
