# PHASE 5B-3A — MINIMAL GLUCOSE ASSIMILATION GAPFILL SEARCH

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). No reaction is permanently added to the model; candidates are applied as
temporary LP scenarios only.

## Method

Minimal-reaction gapfill over a curated, EC-annotated central-carbon universe (8 candidates:
ICL, ME_NADP, ME_NAD, PC, PEPCK, PEPS, NADHDH, G6PDH_NAD). Targets: biomass ≥ WT (G0) to
≥1.5×WT (G50), glucose uptake 0.25–1.0 mmol/gDW/h, Fe2 ≤ WT, O2 ≤ WT, and three Rubisco caps
(R25 = ≤25%, R10 = ≤10%, R0 = 0) × three external-Ci modes (full / low / zero).

## Q1 — R25 / R10 / R0 feasibility

| scenario | no-addition max biomass | gapfill result |
|---|---|---|
| R25 (Rubisco ≤ 25%) | 0.0994 (G10 yes) | **feasible with 0 additions** |
| R10 (Rubisco ≤ 10%) | 0.0519 (below WT) | **feasible with 1 addition** |
| R0 (Rubisco = 0) | 0 (no growth) | **needs >2 additions (not found in universe)** |

## Q2 — Minimum reaction count

0 for R25, 1 for R10, >2 for R0.

## Q3 — Most common / key addition

**G6PDH_NAD — NAD-dependent glucose-6-phosphate dehydrogenase (EC 1.1.1.363).**
`g6p-B + NAD → 6pgl + NADH + H`. This is the only single reaction that turns R10 feasible
(biomass 0.052 → 0.068, +31%).

## Q4 — Mainly lower glycolysis? — NO

The gap is not lower EMP. The model already contains GAPD1/GAPD2/PGK (reversible in the
glucose candidate), but they are not the limiting feature. The effective addition is the
oxidative-PPP / Entner-Doudoroff entry with NAD: G6PDH_NAD → PGL → PGDH → DDGPA gives
`glucose → pyruvate + GAP` directly (DDGPA flux +0.338, PYK +0.095, PDH +0.214).

## Q5 — Fully Rubisco-independent solution? — NO

The G6PDH_NAD solution reaches pyruvate and acetyl-CoA directly (via ED), but 3PG → serine /
glycine still requires the residual 10% Rubisco. Closing Rubisco to 0 makes the solution
infeasible (`RUBISCO_DEPENDENT`). R0 needs >2 reactions (e.g. a pyruvate→PEP / OAA→PEP step
plus anaplerosis), which the 8-reaction universe did not resolve.

## Q6 — Which solutions still depend on internal CO2 recycling?

The R10 solution keeps Rubisco at its 10% cap (RUBISCO flux = 0.2166) for 3PG regeneration;
it is therefore `CBB_RECYCLING_DEPENDENT` for 3PG/serine, but glucose → pyruvate is direct.

## Q7 — Which additions are likely current-GEM omissions?

- **G6PDH_NAD (EC 1.1.1.363)** — `TRUE_NETWORK_ADDITION_CANDIDATE` (model has only
  NADP-dependent G6PDH2).
- ICL (EC 4.1.3.1) — glyoxylate shunt incomplete (MALS only).
- ME, PC, PEPCK/PEPS — missing anaplerotic / gluconeogenic nodes.
- NADHDH — `POTENTIAL_MODEL_OMISSION` (NADHI present but only reverse-ETC direction).

## Q8 — Best candidate for downstream gene/literature audit

**G6PDH_NAD (NAD-dependent glucose-6-phosphate dehydrogenase).** It is a single, well-defined
reaction (EC 1.1.1.363) that enables a direct glucose → pyruvate route (Entner-Doudoroff)
under 90% Rubisco reduction, with moderate glucose (−0.44), Fe2 (−113 vs WT −164.5), O2
(−27.7 vs WT −38.7), and no CO2 uptake. It should be checked against the ATCC 23270 genome
for a NAD-dependent G6PDH locus.

## Tiered ranking

- **Tier 1**: G6PDH_NAD (1 reaction → R10, +31%, direct pyruvate).
- **Tier 2**: NADHDH (redox direction, needs a carbon-route partner); PEPCK/PEPS (candidates
  for R0 gluconeogenesis to 3PG, needs combination).
- **Reject**: none (no candidate relied on a free-energy loop; a bug was found and fixed —
  the first pass accidentally activated the whole universe, which was corrected).

## Caveats

- The reaction universe is a curated 8-reaction central-carbon set, not the full BioCyc PGDB
  (1313 reactions), because BioCyc compound-frameid → model-metabolite mapping was not
  automated. A full-universe gapfill may find smaller/different solutions for R0.
- Carbon-atom attribution is network-route inference, not proof; 13C labelling would be
  required to confirm direct glucose-carbon entry.

## Files

reaction_universe_summary.tsv, gapfill_scenario_matrix.tsv, gapfill_minimal_solutions.tsv,
rubisco_dependency_tests.tsv, direct_precursor_accessibility.tsv, existing_model_overlap.tsv,
solution_ranking.tsv, validation_summary.json + scripts.
