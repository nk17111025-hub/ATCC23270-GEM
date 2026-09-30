# PHASE 5B-4 — EVIDENCE-CORRECTED MODEL REBUILD & DE NOVO ROUTE REDISCOVERY

## Verdict

**v2 PROMOTED** (`minimal_trustworthy_baseline_v2`). NADHI direction correction passes QC with no
free-energy artifact, but it does **not** change the glucose phenotype, Rubisco dependence, or
any route-discovery result. Phase 5B-3C conclusions are UPHELD, not overturned.

## Model provenance and SHA256

- frozen v1: `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c` (unchanged)
- Phase4.1G T0 source: `3b7d6fda…c50345` (matched)
- Phase4.1G TH source: `ea404ac1…8428933` (matched)
- final v2 T0: `97e31dd7…0919e96`
- final v2 TH: `2a9c86d3…917be39`

## Applied corrections

| reaction | change | evidence |
|---|---|---|
| NADHI | LOWER_BOUND 0 -> -1000 (open respiratory NADH-oxidizing direction) | CURRENT_GENOME_DATABASE_SUPPORTED |
| PPC | GPR AFE_1810 -> AFE_1883 / RU820_RS08690 (metadata only) | CURRENT_GENOME_DATABASE_SUPPORTED |
| GHMT3 | GPR: remove AFE_0295 / GlyA (metadata only) | VERIFIED_GPR_ERROR_FIXED |

Only NADHI is an SBML-visible edit; PPC and GHMT3 GPR corrections are Table-1 metadata (the
SBML does not encode GPR).

## QC (all PASS)

free ATP / free NADH / free NADPH = 0 / 0 / 0; baseline (glucose closed) FIM/TTM/TSM = 0.052076;
no-donor biomass = infeasible; no-carbon biomass = 0. NADHI reversibility did not create an
ATP-synthase / quinone / PMF loop.

## Recomputed native phenotype (v2, FIM, donor<=WT, O2<=WT, CO2 free, glucose)

Maximum biomass = **0.1019** (vs WT 0.052076) — identical to Phase 5B-2C. Route is unchanged:
glucose -> G6P (BDGK) -> F6P (PGI1) -> non-oxidative PPP (TKT/TALA/RPI/RPE) -> Ru5P -> RuBP
(PRUK) -> RUBISCO -> 3PG -> 2PG (PGM1) -> PEP (ENO) -> pyruvate (PYK) -> acetyl-CoA (PDH).
PFK = 0, GAPD1/GAPD2 = 0, G6PDH2 = 0: lower EMP and oxidative PPP remain unused even with NADHI
open. NADHI = +2.36 (still reverse-ETC direction).

## Precursor reachability (Rubisco=0, glucose)

Producible: G6P, F6P, Ru5P, GAP, DHAP. Blocked: 3PG (first), R5P, E4P, 13dpg, 2PG, PEP, pyruvate,
acetyl-CoA, TCA (cit/icit/akg/succoa/succ/fum/mal/oaa), serine, glycine, alanine, aspartate,
glutamate, glutamine. First blocked precursor remains **3PG**.

## Rubisco dependence

unrestricted -> 0.1019; 1e-6 -> ~0; 1e-7/1e-8/0 -> 0. Rubisco remains quantitatively essential
for glucose-assisted growth.

## Recomputed route discovery

- R25: 0 additions (native feasible).
- R10: G6PDH_NAD (EC 1.1.1.363), 1 reaction, +31% — unchanged.
- R0: infeasible even with all 12 candidates + NADHI open.

G6PDH_NAD remains the R10 minimum; ED remains the preferred engineering route; no new route
became superior.

## Answers to the required questions

1. Optimal native glucose route in v2: CBB-assisted (glucose -> PPP -> RuBP -> Rubisco -> 3PG), unchanged.
2. Glucose through lower glycolysis: NO (PFK/GAPD = 0).
3. NADH disposal still limiting: YES (the NADHI direction is now available, but the redox/energy
   economy still routes glucose through CBB; NADH disposal alone does not unlock glycolysis).
4. First blocked precursor: 3PG (unchanged).
5. Rubisco still required: YES.
6. R0 feasible natively: NO.
7. R25 minimum: 0.
8. R10 minimum: G6PDH_NAD (1).
9. R0 minimum: not found (infeasible).
10. G6PDH_NAD still selected: YES (R10).
11. ED still preferred: YES.
12. Different route superior: NO.
13. Truly heterologous selected reactions: G6PDH_NAD (NAD-dependent G6PDH); NADHI is a native
    direction correction (N1), not an engineering addition.
14. Smallest genuine engineering gap: 1 reaction (R10); R0 remains a systemic gap.
15. Predicted biomass improvement: +96% (0.052 -> 0.102) under accepted constraints (unchanged).
16. Phase 5B-3C conclusions overturned: NONE.

## Project-level verdict

Correcting the NADHI respiratory direction (and the PPC/GHMT3 GPR metadata) is safe but
phenotypically neutral: the corrected model still routes glucose through the Rubisco-dependent
CBB cycle, still cannot use lower glycolysis, and still cannot achieve Rubisco-independent (R0)
growth. The genuine engineering gap therefore remains what Phase 5B-3C found — a systemic
3PG/pentose/redox dependency on CO2 fixation — and the minimal near-term intervention remains a
NAD-dependent G6PDH (Entner-Doudoroff) entry for the R10 target, not a native-network fix.
