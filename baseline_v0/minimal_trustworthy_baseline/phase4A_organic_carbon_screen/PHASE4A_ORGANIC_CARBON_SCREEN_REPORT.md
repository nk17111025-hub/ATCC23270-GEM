# PHASE 4A — LOW-COST ORGANIC-CARBON NATIVE-CAPABILITY SCREEN

Model: `baseline_v0/minimal_trustworthy_baseline/freeze/minimal_trustworthy_baseline_v1.xml`

## 1. Frozen-model integrity

- SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c` (MATCH, byte-identical).
- 615 reactions, 573 metabolites, 3 compartments (c/p/e).
- Approved v1 patches reproduced: MACPD=0/0, ACOATA=0/1000, Htpp=0/0.
- Reference growth reproduced with a scipy/HiGHS FBA engine:
  FIM=0.052076, TTM=0.052076, TSM=0.052076.
- Free-energy test (all exchanges closed): free ATP = 0.0, free NADH = 0, free NADPH = 0.
  No free-energy / free-redox cycle detected.

## 2. Candidate substrates found

Intracellular forms present in the model:

- glucose (`glc-A[c]`, `glc-B[c]`), glycerol (`glyc[c]`), acetate (`ac[c]`),
  pyruvate (`pyr[c]`), formate (`for[c]`).
- 4-hydroxybenzyl alcohol (`4hba[c/p/e]`) — the only organic compound with a native
  exchange + transport path (`Ex_4hba[e]`, `4HBAtex`, `4HBAtpp`).

Absent from the model entirely (no metabolite, no enzyme):

- ethanol, acetaldehyde, lactate, methanol, free fructose, sucrose.

The only organic-carbon exchange reactions in the model are `Ex_4hba[e]` and `Ex_eps_AFE[e]`
(EPS is a secreted product, not a substrate). There is NO extracellular exchange for
glucose, glycerol, acetate, pyruvate, or formate.

## 3. Substrates with complete native paths

None usable. `4hba` has complete native transport (exchange + outer/inner-membrane transport),
but it is a thiamine-biosynthesis intermediate with NO catabolic route to central carbon.
The native feasibility test confirms `4hba` uptake is ~0 and growth is unchanged.

## 4. Substrates one functional step away

- acetate — add one transporter; native `ACS` (acetyl-CoA synthetase, conf 2, AFE_1969)
  converts acetate to acetyl-CoA.
- pyruvate — add one transporter; pyruvate is already a central metabolite.
- glucose — add one transporter; native `BDGK` links glucose to G6P. But the lower EMP
  (GAPD1/GAPD2/PGK) is written gluconeogenic and irreversible, so glucose reaches
  pyruvate/acetyl-CoA only via the PPP→RuBP→RUBISCO→3PG route (see artifact flags).

## 5. Substrates requiring >= 2 new functions

- glycerol — transporter + reverse G3PD2 (glycerol-3-phosphate oxidation to DHAP).
- fructose — transporter + fructokinase.
- sucrose — transporter + invertase / sucrose phosphorylase.
- lactate — transporter + lactate dehydrogenase.

## 6. Substrates requiring major pathway construction

- ethanol — transporter + ethanol dehydrogenase + acetaldehyde dehydrogenase (no ethanol/acetaldehyde present).
- methanol — transporter + methanol dehydrogenase + full C1 assimilation (absent).
- formate — transporter + reverse FTHFL + full C1 assimilation (absent). Formate currently only
  oxidizes to CO2 (FDH) or donates formyl to purines (GART).

## 7. Native feasibility results (Task 4, no new reactions)

Only `4hba` has a native exchange/transport path. Under FIM/TTM/TSM and uptake caps
0.1–10 mmol/gDW/h, `4hba` uptake remains ~0 and biomass stays 0.052076 (identical to baseline).
`4hba` is therefore NOT usable as an organic-carbon source by the frozen model.

## 8. Biomass changes (hypothetical single-transporter probes, mixotrophic, CO2 intact)

| substrate | baseline | +substrate (FIM) | +substrate (TTM/TSM) |
|---|---|---|---|
| glucose | 0.052076 | 0.667 | 1.614 |
| acetate | 0.052076 | 0.077 | 0.077 |
| pyruvate | 0.052076 | 0.195 | 0.195 |
| glycerol | 0.052076 | 0.0525 | 0.0525 |
| formate | 0.052076 | 0.312 | 0.312 |
| 4hba | 0.052076 | 0.052076 | 0.052076 |

These are model predictions with an added hypothetical cytosolic exchange, NOT experimental proof.

## 9. Fe/S oxidation retained or lost

In every mixotrophic probe the Fe/S donor uptake remains active (Fe2 for FIM, tetrathionate for TTM,
thiosulfate for TSM) and oxygen uptake remains active. Organic carbon availability does not shut off
Fe/S oxidation in these solutions.

## 10. Organic carbon actually entering central metabolism

- acetate — YES (ACS -> acetyl-CoA -> TCA/lipids).
- pyruvate — YES (already central; PDH -> acetyl-CoA).
- glucose — PARTIAL. Reaches G6P/F6P and pentose-phosphate (nucleotides/glycogen), but pyruvate
  is reached only through the PPP→RuBP→RUBISCO→3PG route, not direct glycolysis (glycolysis is blocked).
- glycerol — essentially NO (only the phospholipid backbone; G3PD2 oxidation blocked).
- formate — NO as organic carbon; formate carbon enters only as CO2 after FDH oxidation, then RUBISCO fixation.
- 4hba — NO.

## 11. Artifact flags

- formate "growth" is an artifact of formate -> CO2 (FDH) -> RUBISCO fixation; there is no C1 assimilation route.
- glucose growth is partly sustained by RUBISCO running on internally recycled CO2 (RUBISCO is reversible);
  the CO2-closed "heterotrophic" result is therefore not a clean heterotrophic prediction.
- acetate utilization activates low-confidence glyoxylate-cycle alternatives (MALS, reverse GLXCL/GLYCK).
- NADTRHD (confidence 1, no GPR) carries large flux in the glucose and formate probes (redox balancing).
- All CO2-closed heterotrophic tests are affected by internal CO2 recycling; treat them as sensitivity
  points, not heterotrophy claims.

## 12. WATCH reactions activated

- BDGK: forward glucose -> G6P in the glucose probe (expected direction).
- NADTRHD: large flux in glucose (+9.3) and formate (+4.4) probes.
- MALS: active in the acetate probe (+0.23).
- GLXCL / GLYCK: reverse flux in the acetate probe (-0.12).
- ACCOAC: active (forward) in all probes.
- MDH/FUM: normal TCA-cycle flux (FUM negative = fumarate -> malate).

None of these WATCH reactions were modified; they are only observed to carry flux.

## 13. Minimal intervention counts

| substrate | minimum functional interventions |
|---|---|
| acetate | 1 (transporter) |
| pyruvate | 1 (transporter) |
| glucose | 1 (transporter); +2 direction reversals (GAPD1, PGK) for direct glycolysis |
| glycerol | 2 (transporter + reverse G3PD2) |
| fructose / sucrose / lactate | 2 (transporter + one enzyme) |
| ethanol | 3 (transporter + 2 enzymes) |
| methanol / formate | >=3 incl. full C1 assimilation (major) |

## 14. Evidence gaps requiring CENTRAL REVIEW

- BDGK (glucokinase-like) confidence 1 / no GPR — glucose phosphorylation remains weak.
- G3PD2 direction — glycerol-3-phosphate oxidation is not represented.
- GAPD1/GAPD2/PGK and G6PDH2 are written gluconeogenic/reductive and irreversible;
  catabolic directionality (glycolysis, oxidative PPP) is absent.
- NADTRHD confidence 1 / no GPR — large redox-balancing flux.
- MALS / GLXCL / GLYCK / GLXR — low-confidence glyoxylate-cycle alternatives.
- No transporter for any candidate substrate (genome-level transporter annotation needed).
- No C1 assimilation module (serine cycle / RuMP) for formate or methanol.
- Acetate acid-stress toxicity is adjudicated externally (not assessed here).

## 15. Files produced

- PHASE4A_ORGANIC_CARBON_SCREEN_REPORT.md (this file)
- organic_carbon_entry_inventory.tsv
- substrate_pathway_gap_map.tsv
- organic_carbon_engineering_distance.tsv
- native_organic_carbon_feasibility.tsv
- organic_carbon_contribution_analysis.tsv
- minimal_completion_candidates.tsv
- phase4A_validation_summary.json
- phase4a_lib.py, task1_inventory.py, task2_3_paths.py, task4_feasibility.py,
  task5_6_probe.py, task6_minimal.py, finalize.py (reproducible pipeline)

No final substrate is selected here. CENTRAL REVIEW combines these model results with
substrate price, toxicity, experimental evidence, and industrial practicality.
