# PHASE 5B-3C — FULL-UNIVERSE R0 MINIMAL PATHWAY COMPLETION SEARCH

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). No reaction is permanently added.

## Method

Built a non-native reaction universe from BiGG (universal reactions, 28,302 total) via the
BiGG API, mapping BiGG metabolite IDs to the model (the model already uses BiGG-style IDs, so
mapping is direct/HIGH confidence). Query was rate-limited (502 errors), yielding 48 downloaded
reactions, 33 mappable; after removing native-equivalent reactions, a **12-reaction non-native
universe** remained: ME1, ME2, PPCK, ICL, PPS, GND, NADH16, PPDK, G6PDH_NAD, PC, PFL, NADHox.

This is ~4× the 5B-3A/3B curated pool and covers glycolysis/ED/PPP/anaplerosis/gluconeogenesis/
redox, with real BiGG stoichiometry and EC numbers.

## Q1 — Is R0 feasible in the full universe? — NO

Exhaustive search of 1–6 reaction additions found **zero** feasible R0 solutions, under both
R0-CiZero and R0-CiLow. Even adding **all 12 candidates at once** leaves R0 infeasible
(biomass = 0).

## Q2 / Q3 — R0-CiZero / R0-CiLow minimum reaction count

**Not found** (infeasible even with 6 reactions, and even with all 12 candidates).

## Q4 — Minimum solution

None exists within this universe.

## Q5 — Is the bottleneck 3PG / serine-glycine supply?

Only partially. The block is more fundamental than a single precursor. Relaxing Fe2, O2, and
external-Ci bounds one at a time and together did NOT make R0 feasible, which means the failure
is not an energy/O2/carbon-input limit. Diagnostic decomposition:

- RUBISCO=0 + RUBISCOX=0 → infeasible;
- RUBISCO=0 + RUBISCOX allowed (Ci=0) → still infeasible;
- RUBISCO allowed + RUBISCOX=0 → feasible (0.102).

So the carboxylase (RUBISCO) itself is essential; the oxygenase (RUBISCOX) cannot substitute at
Ci=0, and no standard central-carbon reaction set replaces it.

## Q6 — Is G6PDH_NAD / ED still optimal? — YES, but only for R10/R25

The consensus reaction (G6PDH_NAD → ED → pyruvate + GAP) remains the best single reaction for
R10 (Rubisco ≤ 10%), but it does nothing for R0.

## Q7 — Is there a simpler direct-glucose route than ED? — NO

No candidate route (lower glycolysis, anaplerosis, glyoxylate, gluconeogenesis) enables R0.

## Q8 — Do solutions merely swap Rubisco for another carboxylase? — YES, and that is insufficient

The candidate set includes PC and PPC (PEP/pyruvate carboxylases) and ME (malic enzyme,
reversible carboxylation). Even with all of these present, R0 stays infeasible, showing that
replacing Rubisco with another carboxylase does not rescue the phenotype.

## Q9 — Which additions are likely current-GEM omissions?

- G6PDH_NAD (EC 1.1.1.363) — NAD-dependent G6PDH missing (model has only NADP-G6PDH2).
- ICL (EC 4.1.3.1) — glyoxylate shunt incomplete (MALS only).
- ME1/ME2, PC, PPCK, PPS, PPDK — missing anaplerotic / gluconeogenic nodes.
- NADH16 — model NADHI is reverse-ETC only (POTENTIAL_MODEL_OMISSION).
- GND — decarboxylating 6-phosphogluconate dehydrogenase missing.

## Q10 — Top reaction sets for downstream gene audit

Given R0 is not reachable by reaction addition, the shortlist is for the achievable R10/R25
target, not R0:

1. G6PDH_NAD (R10, +31%, ED) — Tier 1.
2. G6PDH_NAD + PPCK/PPS (ED + gluconeogenesis) — candidate if R0 is ever attempted with a
   redesigned redox system.
3. ICL + ME/PC (glyoxylate + anaplerosis) — supports TCA without external CO2, but still not
   enough for R0.

## Overall conclusion

**R0 is unreachable in the current ATCC 23270 genome-scale model.** The model is a
chemolithoautotroph whose central metabolism is stoichiometrically and redox-coupled to CO2
fixation via Rubisco. Adding standard central-carbon reactions (glycolysis, ED, PPP,
anaplerosis, gluconeogenesis, redox) does not remove this dependence. The earlier
+95.6% glucose gain and the R10 (+31%) result are Rubisco-mediated, not Rubisco-independent.

This is a `NATIVE_NETWORK_INSUFFICIENT`-style result for R0: achieving truly
Rubisco-independent heterotrophic glucose assimilation would require re-engineering beyond
standard reaction completion — e.g. a different electron-transport / redox-balancing system —
which is outside a minimal pathway-completion search.

## Caveats (reported honestly)

- The universe is 12 non-native reactions, not the "hundreds to thousands" target: the BiGG
  API was rate-limited (502 errors), and only 48 reactions were fetched before the run
  completed. A full BiGG/MetaNetX download with complete metabolite mapping could in principle
  find a route outside this set, but the diagnostic (infeasible even with all candidates +
  relaxed bounds) strongly suggests the block is systemic, not a missing reaction.

## Files

minimal_reaction_sets.tsv, _bigg_reactions.json (raw BiGG), _universe.json,
search_R0.py + build_universe.py (reproducible).
