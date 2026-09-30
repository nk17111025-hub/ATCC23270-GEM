# PHASE 5B-5 — RUBISCO-KO DE NOVO NETWORK RESCUE

## Verdict

**NO RESCUE — CURRENT SEARCH UNIVERSE INSUFFICIENT**

Starting from the promoted v2 model (NADHI respiratory direction already opened), an exhaustive
de novo search over a 14-reaction candidate universe (BiGG-derived + curated, including
RPIr and the oxidative G6PDH direction) found **no Rubisco-independent (R0) rescue at 1–6 added
reactions**. R0 biomass remains 0.

## Answers

1. Was biomass rescued with Rubisco disabled? NO.
2. True minimum added reactions: not found (>=7, or requires a reaction class absent from the universe).
3. Exact minimum reaction sets: none.
4. Carbon route of solutions: n/a (no solution).
5. Old blocked precursors rescued: none (R0 baseline blocks 43/45 biomass precursors; 3PG, R5P,
   E4P, acetyl-CoA, TCA, and the amino-acid families remain blocked).
6. Necessity testing: n/a.
7. Artifact-free: baseline QC passes (free ATP/NADH/NADPH = 0; no carbon/donor biomass = 0/infeasible);
   no rescue candidate existed to test.
8. After Rubisco restored: n/a.
9. ROBUST_BYPASS vs KO_ONLY_ESCAPE: n/a.
10. Genuinely Rubisco-independent: NO.
11. Is ED still best? For the achievable R10 target, ED (G6PDH_NAD) remains the best 1-reaction
    solution; it does NOT rescue R0.
12. Is G6PDH_NAD still in the best solution? For R10 yes; for R0 nothing works.
13. Superior non-ED route appeared? NO.
14. Smallest credible engineering route to test next: unchanged — a NAD-dependent G6PDH
    (Entner-Doudoroff) entry for the R10 target; R0 remains a systemic gap.
15. Reactions requiring gene-level validation: G6PDH_NAD (NAD-dependent G6PDH, EC 1.1.1.363) is
    the only near-term candidate; R0 would require a pentose-phosphate regeneration route that
    is independent of the Calvin-cycle 3PG pool, plus a NADPH sink.

## Why R0 remains unreachable

The v2 model (with NADHI direction corrected) still cannot make 3PG (serine/glycine) or
regenerate the pentose-phosphate pool (R5P/E4P for nucleotides/aromatics) without Rubisco, and
the oxidative PPP's NADPH cannot be balanced. These are coupled to the Calvin-cycle 3PG pool
(the non-oxidative PPP is primed by S7P, which is 3PG-derived). No candidate in the current
universe supplies this primer or the required NADPH sink, so even with NADHI open and 14
candidates, R0 is infeasible.

This is consistent with, but does not inherit the wording of, Phase 5B-3C: the result was
recomputed de novo on v2 and is unchanged.
