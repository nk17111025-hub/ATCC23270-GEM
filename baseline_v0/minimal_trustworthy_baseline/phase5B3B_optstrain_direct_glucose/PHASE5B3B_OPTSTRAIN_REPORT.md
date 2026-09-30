# PHASE 5B-3B — OPTSTRAIN SEARCH FOR DIRECT GLUCOSE ASSIMILATION

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). No reaction is permanently added; all candidates are temporary LP scenarios.

## Method

OptStrain-style minimal-heterologous-reaction search over a curated 12-reaction universal
central-carbon database (EC-numbered; ED / PPP / glyoxylate / anaplerosis / gluconeogenesis /
redox / pyruvate-formate classes). Targets: biomass ≥ WT to ≥1.5×WT, glucose 0.25–1.0
mmol/gDW/h, Fe2 ≤ WT, O2 ≤ WT, Rubisco ≤ 25%/10%/0% (R25/R10/R0) × external-Ci full/low/zero.

## 1. Minimum heterologous reaction count

| scenario | minimum non-native reactions | best reaction set |
|---|---|---|
| R25 (Rubisco ≤ 25%) | 1 | G6PDH_NAD (or PEPCK/PC/PFL/ICL/PEPS/PPDK/NADHDH) |
| R10 (Rubisco ≤ 10%) | 1 | **G6PDH_NAD** (PFL is a weak second) |
| R0 (Rubisco = 0) | **>3** | not found in universe |

## 2. Fully Rubisco-independent solution? — NO

G6PDH_NAD reaches pyruvate/acetyl-CoA directly (Entner-Doudoroff), but 3PG → serine/glycine
still needs the residual Rubisco. At Rubisco = 0 it fails (`FAILS`), so it is
`PARTIAL_CBB_DEPENDENCE`. R0 requires >3 reactions (a pyruvate/OAA → PEP gluconeogenic step
plus anaplerosis/redox support) and was not resolved.

## 3. Most common external reaction

**G6PDH_NAD — NAD-dependent glucose-6-phosphate dehydrogenase (EC 1.1.1.363).** This is the
Consensus hit: found by both the Phase 5B-3A GapFill and this OptStrain search.

## 4. Dominant mechanism

**Entner-Doudoroff (Route B).** G6PDH_NAD → PGL → PGDH → DDGPA converts glucose → pyruvate + GAP
directly. Secondary mechanisms: gluconeogenesis (PEPCK/PEPS/PPDK), anaplerosis (PC/ME),
glyoxylate (ICL), pyruvate-formate (PFL). Classical lower glycolysis (Route A) is not the
limiting feature — GAPD1/GAPD2/PGK already exist (reversible).

## 5. Solutions that reduce CO2 / Fe2 / O2 together

G6PDH_NAD (R10): CO2 → 0, Fe2 69% of WT, O2 72% of WT, glucose −0.44 mmol/gDW/h (within the
experiment-supported 0.37–0.65 range).

## 6. Closest to experimental glucose uptake

G6PDH_NAD (glucose −0.44) and PFL (−0.36) are within the physiological range.

## 7. GapFill vs OptStrain convergence

**Consensus**: G6PDH_NAD (NAD-G6PDH). Partial consensus: ICL, PEPCK/PEPS/PPDK/PC/ME (found in
both, but only under R25 or as degenerate hits). OptStrain-only: PFL (not in the 5B-3A
universe). NADHDH appeared in both but carries little/zero flux — low value.

## 8. Possibly-native vs likely-heterologous

- **G6PDH_NAD** — `POSSIBLY_NATIVE` (model has only NADP-G6PDH2; a genome NAD-G6PDH paralog
  should be checked in GCF_049532655.1).
- **NADHDH** — `POTENTIAL_MODEL_OMISSION` (NADHI is reverse-ETC only).
- **PFL, PEPCK, PEPS, PPDK, PC, ICL, ME** — `POSSIBLY_NATIVE` / `LIKELY_HETEROLOGOUS`
  (absent from the GEM; genome paralog/promiscuity needs checking).

## 9. Top reaction sets for downstream gene-level audit

1. **G6PDH_NAD** — Tier 1: 1 reaction, R10 feasible (+31%), ED route, consensus.
2. **PFL** — Tier 2: 1 reaction, R10 marginal (+4%), formate byproduct.
3. **PEPCK / PEPS** — Tier 2: 1 reaction for R25; candidate for R0 gluconeogenesis to 3PG.
4. **G6PDH_NAD + PEPS** — Tier 2: 2 reactions; candidate for fuller Rubisco independence
   (needs further verification).

## Caveats

- The database is a curated 12-reaction universal set (EC-annotated); a full MetaNetX/BiGG
  universe was not downloaded (no automated BiGG/MetaNetX compound-ID mapping), so R0 may have
  smaller solutions outside this set.
- Carbon-atom attribution is network-route inference, not proof; 13C-MFA is required to
  confirm direct glucose-carbon entry.

## Files

reaction_database_summary.tsv, native_vs_non_native_mapping.tsv, optstrain_scenario_matrix.tsv,
optstrain_minimal_solutions.tsv, rubisco_independence_tests.tsv, precursor_connectivity.tsv,
physiology_screen.tsv, engineering_complexity.tsv, gapfill_optstrain_comparison.tsv,
solution_ranking.tsv, validation_summary.json + scripts.
