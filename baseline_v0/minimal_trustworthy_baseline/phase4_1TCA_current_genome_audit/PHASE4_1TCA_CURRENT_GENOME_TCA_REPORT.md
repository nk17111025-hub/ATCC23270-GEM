# PHASE 4.1TCA — Current-genome TCA completeness audit

## Executive verdict

**Canonical oxidative TCA is not supported as complete in GCF_049532655.1.** No canonical `sucA`/OGDH E1 (K00164, EC 1.2.4.2) or `sucB`/OGDH E2 (K00658, EC 2.3.1.61) is present in the current RefSeq annotation. The generic 2-oxoacid-dehydrogenase candidates resolve to pyruvate/other 2-oxoacid complexes by KEGG assignment, current gene symbol, domain class and genomic context. Shared E3-like proteins do not establish 2-oxoglutarate substrate specificity.

**No complete alternative 2-OG closure is supported.** `RU820_RS07175 / WP_012536595.1` is a genuine NAD-dependent succinate-semialdehyde dehydrogenase and maps to historical `AFE_1546 / K00135`; it is also detected in type-strain proteomics. However, the required upstream 2-oxoglutarate decarboxylase is absent; the GABA-shunt glutamate decarboxylase and GABA transaminase are not supported; and no `korA/korB` or OOR complex is present. Therefore SSADH alone is not a bypass.

A second issue emerged: current RefSeq does not support a canonical `sdhABCD` set. Historical `AFE_2550` was sometimes labeled `sdhC`, but the 2009 sulfur-oxidation study placed `AFE_2550–2555` in the sulfur-induced HdrABC locus, and the current ortholog `RU820_RS11745 / WP_012537254.1` is annotated as an HdrB-family Fe-S protein. No SdhA/B catalytic core was identified by the current annotation audit. This makes succinate→fumarate unsupported at the current-gene level unless a highly divergent complex is later demonstrated.

Overall biological status: **TCA_TRULY_INCOMPLETE_SUPPORTED**. Model status: **MODEL_TCA_INCOMPLETE_MATCHES_BIOLOGY**, with an important local warning that `R_SUCD` is not genomically supported as currently represented.

## Current TCA inventory and route

| Step | Current evidence | Status |
|---|---|---|
| OAA + acetyl-CoA → citrate | RU820_RS14240 / WP_012537576.1 citrate synthase (legacy AFE_3065) | PRESENT |
| citrate ↔ isocitrate | RU820_RS02050 / WP_012536129.1 aconitate hydratase (AFE_0423) | PRESENT |
| isocitrate → 2-OG | RU820_RS02055 / WP_012606555.1 `icd` (AFE_0424) | PRESENT |
| 2-OG → succinyl-CoA (canonical OGDH) | no sucA/K00164; no sucB/K00658 | ABSENT |
| 2-OG → succinate by complete alternative bypass | SSADH present, upstream closure genes absent | ABSENT |
| succinyl-CoA ↔ succinate | RU820_RS02060 `sucC` + RU820_RS02070 `sucD` | PRESENT |
| succinate → fumarate | no supported canonical SdhABCD; old AFE_2550/sdhC alias is HdrB | ABSENT/UNRESOLVED DIVERGENT-COMPLEX POSSIBILITY |
| fumarate ↔ malate | RU820_RS12345 `fumC` | PRESENT |
| malate ↔ OAA | RU820_RS13930 / WP_012537534.1 MDH | PRESENT |

## Canonical OGDH candidate audit

The two strongest false-positive risks are the `AFE_1811–1813` and `AFE_3068–3070` complexes. Current loci are `RU820_RS08345–08355` and `RU820_RS14250–14260`. KEGG assigns AFE_1813 and AFE_3070 to K00161 (pyruvate dehydrogenase E1 alpha), and AFE_1812/AFE_3069 to K00162 (PDH E1 beta). `RU820_RS14260` is explicitly `pdhA`, EC 1.2.4.1. AFE_3068 is K00382/E3-like but is embedded in the PDH locus. These proteins cannot be promoted to OGDH without substrate-specific evidence.

`RU820_RS08740 / WP_012536825.1` is a true dihydrolipoyl dehydrogenase (EC 1.8.1.4), and `RU820_RS04600` is another E3-family protein. Both are shared-enzyme candidates only. Neither has neighboring `sucA/sucB`, so neither rescues OGDH.

## KEGG reconciliation

The present KEGG organism entry `afr` is **not based on GCF_049532655.1**; it still states GenBank assembly `GCA_000021485.1`, chromosome `CP001219`, created in 2008. The KEGG organism page listing `M00009 Citrate cycle` is therefore not evidence that the current genome has a complete module, and in any case a module title being listed on the organism pathway-module page is not the same as a satisfied completeness call.

For the disputed module segment, M00009 requires a valid 2-OG oxidation route: canonical K00164 + K00658 plus K00382, or accepted alternatives such as K01616 with K00382, or OOR K00174 + K00175. `afr` contains AFE_3068/K00382, but that is only a shared E3-like component and is in a pyruvate-dehydrogenase context. The required OGDH E1/E2 or OOR pair is absent. Therefore there is **no exact KEGG gene/KO that legitimately closes 2-OG→succinyl-CoA** in ATCC 23270.

## Historical versus current interpretation

Valdés et al. (2008) explicitly concluded that all TCA enzymes except the E1–E3 alpha-ketoglutarate dehydrogenase subunits were predicted, making the TCA incomplete. Osorio et al. (2013) still referred to the incomplete TCA cycle. The 2016 iMC507 reconstruction inherited that biological picture. Current GCF_049532655.1 does **not** reverse the OGDH conclusion: `sucA/sucB` remain absent.

The current reannotation actually weakens one other old TCA assignment. `AFE_2550` appeared in older work with an `sdhC`-like description, but sulfur-expression/context evidence treated it as part of HdrABC, and current RefSeq calls the ortholog an HdrB-family protein. Thus the old claim “everything except OGDH is present” should not be repeated without qualification; succinate→fumarate now also needs explicit re-verification.

## Expression/activity evidence

No direct ATCC 23270 OGDH enzyme-activity measurement demonstrating activity was found. Proteomics instead supports many flanking TCA components: aconitase, ICDH, SucC/D, FumC and citrate synthase are detected in the type-strain-equivalent DSM 14882T dataset. SSADH AFE_1546 is also detected and enriched on pyrite relative to Fe(II), but there is no gene-supported upstream route that converts 2-OG to succinate semialdehyde/GABA. Expression therefore does not convert SSADH into a complete bypass.

No glucose-supplementation-specific induction/activity evidence for a cryptic OGDH or complete bypass was located. This remains an experimental question rather than a model-repair justification.

## v1 GEM reconciliation

The frozen `minimal_trustworthy_baseline_v1.xml` contains `R_CS`, `R_ACONT1/2`, `R_ICDHyr`, `R_SUCOAS`, `R_SUCD`, `R_FUM`, and `R_MDH`, but no OGDH/AKGDH reaction. The SBML does not encode gene associations for these reactions, so all GPRs require external current-genome reconciliation.

The absence of OGDH in the model matches current biology and is **not** presently a model-repair target. `R_SUCD`, however, is problematic: it is encoded as `H+ + NADH + fumarate <-> NAD+ + succinate`, not as a canonical quinone-linked succinate dehydrogenase reaction, and no current canonical SdhABCD set supports it. This reaction should be flagged for Central Review, not edited during this phase. `R_ICDHyr` is fully reversible in v1 despite the current enzyme annotation supporting oxidative NADP-dependent ICDH; its directionality is likewise only flagged for review.

## Final classification

- Canonical OGDH: **NO**
- Complete alternative 2-OG closure: **NO**
- Biological classification: **TCA_TRULY_INCOMPLETE_SUPPORTED**
- Model classification: **MODEL_TCA_INCOMPLETE_MATCHES_BIOLOGY**
- Minimum missing biochemical function: **a substrate-supported 2-OG oxidation/closure function** (canonical OGDH E1/E2 with an E3, or a complete alternative pathway). In addition, **succinate→fumarate lacks a currently supported canonical Sdh complex and remains a second genomic/model-review issue**.

No GEM modification, FBA, gap filling or strain design was performed.
