# PHASE 5B-3C CLOSEOUT — R0 ESSENTIALITY DECOMPOSITION

Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
(UNCHANGED). Diagnostic supply/demand reactions are temporary only.

## Result wording

**R0_NOT_RESCUED_IN_CURRENT_MODEL_AND_CANDIDATE_UNIVERSE** (not "R0 impossible", not
"Rubisco biologically essential").

## Q1 — First blocked precursor

**3PG (3-phosphoglycerate).** With Rubisco = 0 + glucose + G6PDH_NAD, the reachable pool is
G6P, F6P, Ru5P, Xu5P, GAP, DHAP, pyruvate, alanine. The first blocked node is **3PG** (together
with 13dpg / 2PG / PEP), which recovers at the smallest tested Rubisco flux (1e-6).

## Q2 — Is 3PG / serine / glycine the main bottleneck?

It is the **primary** bottleneck (serine and glycine are downstream of 3PG and stay blocked),
but it is not the only one — see Q3.

## Q3 — Other independent precursor blocks

- **R5P** (ribose-5-phosphate): Ru5P is reachable, but RPI is irreversible in the
  Ru5P→R5P direction; R5P is needed for nucleotides.
- **E4P** (erythrose-4-phosphate): the non-oxidative PPP cannot supply it without 3PG-derived
  carbon (S7P pool).
- **acetyl-CoA**: pyruvate is reachable (via ED), but PDH / CoA / redox coupling is blocked.
- **TCA / anaplerotic nodes**: OAA, citrate, isocitrate, α-KG, succinate, fumarate, malate.
- **Amino acids**: aspartate, glutamate, glutamine (downstream of OAA / α-KG), plus
  serine/glycine (downstream of 3PG).

## Q4 — Does supplying one precursor restore biomass? — NO

No single precursor supply rescued biomass (G0). Pairs (3PG+serine, 3PG+serine+glycine,
pyruvate+3PG, OAA+3PG) also did not.

## Q5 — Minimum precursor pools to restore biomass

Not achievable by a small precursor set: **43 of 45 biomass constituents are blocked**, and
even supplying all 20 blocked central metabolites did not restore biomass because the remaining
blocked biomass precursors (nucleotides, lipids, cofactors) are also downstream of the same
3PG / PPP / acetyl-CoA network.

## Q6 — Failure class

**TYPE_C — SYSTEMIC_NETWORK_DEPENDENCE.** Rubisco in the current GEM is not merely a 3PG
source; it supplies the 3PG/pentose/TCA/redox architecture that all 43 biomass precursors
depend on. The failure is a **multi-factor coupling** of carbon topology (3PG → lower
glycolysis/PPP/TCA), redox (glycolytic NADH with no sink), and precursor supply — not a single
metabolite or a pure ATP/redox loop.

## Q7 — What ED (G6PDH_NAD) already solves

ED solves **glucose → pyruvate + GAP** (and pyruvate → alanine). It does NOT solve:
glucose → 3PG/2PG/PEP (needs gluconeogenesis with a NADH sink), glucose → R5P/E4P (needs the
3PG-fed pentose-phosphate balance), or glucose → acetyl-CoA/TCA/OAA/αKG.

In one sentence: **ED solves glucose → pyruvate, but R0 still fails because the model cannot
make 3PG (serine/glycine), R5P/E4P (nucleotides/aromatics), and acetyl-CoA/TCA (lipids and
the aspartate/glutamate families) without Rubisco-mediated CO2 fixation.**

## Q8 — Priority for any future design

1. A **3PG supply route** (gluconeogenic pyruvate/GAP → 3PG) together with a **NADH sink**.
2. A **pentose-phosphate balance** independent of the 3PG/Calvin pool (R5P and E4P).
3. An **acetyl-CoA / TCA / anaplerotic** route that does not depend on Rubisco-derived carbon.

Because these are coupled (not a single gap), the realistic next step is a network-level
re-design (different electron transport / redox balancing), not a small reaction addition.

## Files

biomass_component_producibility.tsv, central_precursor_producibility.tsv,
primary_blocked_precursors.tsv, rubisco_dependency_map.tsv, precursor_rescue_matrix.tsv,
redox_energy_diagnostic.tsv, ED_resolved_vs_unresolved.tsv, R0_failure_classification.json.
