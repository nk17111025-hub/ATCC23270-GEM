# PHASE 4.1G-2A — RAPID GLUCOSE-ROUTE EVIDENCE RETRIEVAL

Scope: publication-backed evidence for four glucose-interpretation questions only.
No model edits, no FBA, no GPR assignment, no strain design.

Genome note: the AFE_ locus tags are from the original ATCC 23270 genome
(Valdés et al. 2008, GCF_000021485). The RU820_RS tags belong to the newer
GCF_049532655.1 assembly. Protein-level annotations were retrieved from UniProt
keyed to the AFE_ names / WP_ accessions (strain ATCC 23270).

## Q1 — GLUCOSE UPTAKE

### Per-candidate evidence

| gene | UniProt annotation | glucose-specific evidence |
|---|---|---|
| RU820_RS11605 / AFE_2522 / WP_012537238.1 | Carbohydrate-selective porin, OprB family (TrEMBL) | NO (annotation only) |
| RU820_RS10385 / AFE_2250 / WP_012537080.1 | Carbohydrate-selective porin, OprB family (TrEMBL) | NO (annotation only) |
| RU820_RS10670 / AFE_2312 / WP_012537116.1 | Sugar transporter family protein (MFS, TC 2.A.1.1) | NO (annotation only) |
| RU820_RS09115 / AFE_1971 / WP_012536876.1 | Transporter, putative (MFS, PF12832) | NO (annotation only) |

No paper reports a direct glucose-transport assay, a glucose-specific porin
substrate test, glucose-induced expression of a named transporter, or a
transporter knockout/complementation in ATCC 23270.

The single strongest physiological statement is Wang et al. 2012
(PMID 22210219), which reports that the pfkB mutant "did not completely prevent
A. ferrooxidans from assimilating exogenous glucose" — i.e. glucose can enter and
be used, but the responsible transporter is NOT identified.

The PTS-like cluster (AFE_3018–AFE_3023 / RU820_RS14005–14030) was not found in any
publication, and its completeness as a glucose PTS (EI/HPr/EIIA/EIIB/EIIC) is unverified.

**Conclusion: NO DIRECT GLUCOSE-SPECIFIC TRANSPORTER IDENTIFIED.**
The OprB porins (AFE_2522, AFE_2250) and the sugar MFS (AFE_2312) are the
strongest annotation-only candidates, but "carbohydrate porin" / "sugar transporter"
must not be read as glucose specificity without an experiment.

## Q2 — GLUCOKINASE / BDGK (AFE_2841 / RU820_RS13140 / WP_012537424.1)

- UniProt: "ROK family protein" (TrEMBL, unreviewed). No EC number, no
  "glucokinase" / "glucose kinase" annotation, no catalytic-activity line.
- Literature: zero paper hits for AFE_2841.
- Wang et al. 2012 performed glucose supplementation and RT-qPCR of
  "central carbohydrate metabolism" genes, but the published abstract does not
  name glucokinase / AFE_2841; AFE_1807 (pfkB) is the enzyme functionally validated
  in that study, and it is phosphofructokinase, NOT a glucose kinase.

ROK-family membership is compatible with glucokinase (e.g. E. coli Glk is ROK-family),
but the family also includes N-acetylglucosamine kinases and ROK repressors.

**Classification: ANNOTATION_ONLY.**
**Sufficient to assign BDGK GPR? NO** (the current annotation is only "ROK family protein";
glucokinase activity is not demonstrated for AFE_2841).

## Q3 — OXIDATIVE PPP / ZWF (AFE_2025 / RU820_RS09360 / WP_012536916.1)

- UniProt: "Glucose-6-phosphate 1-dehydrogenase (G6PD) EC 1.1.1.49", gene zwf,
  with both G6PDH domains (NAD-binding PF00479 + C-terminal PF02781) and GO terms
  "glucose metabolic process" + "pentose-phosphate shunt, oxidative branch".
- No ATCC 23270-specific biochemical G6PDH assay was found.
- Expression/proteomics support is indirect: zwf/AFE_2025 appears in the proteome
  (Bellenberg 2019, PMID 30984136) and in the QS transcriptome (Mamani 2016,
  PMID 27683573), but neither shows glucose -> G6P -> oxidative-PPP flux.

**Is oxidative-PPP participation in glucose metabolism supported? PARTIAL.**
The enzyme identity is solid (EC 1.1.1.49 with complete domain architecture, HAMAP-based),
but there is no direct functional or flux evidence that glucose feeds the oxidative PPP
in ATCC 23270. Note also that the 2016 GEM represents the G6PDH step in the reductive
direction; enzyme annotation does not resolve model directionality.

## Q4 — LOWER EMP DIRECTIONALITY (GAP / PGK)

- AFE_3251 / RU820_RS15125 / WP_009567621.1 — "Glyceraldehyde-3-phosphate dehydrogenase,
  type I (EC 1.2.1.-)", NAD(P)-binding domain, GO "glucose metabolic process".
  Type-I GAPDH is the canonical reversible phosphorylating enzyme.
- AFE_3250 / RU820_RS15120 / WP_012537677.1 — "Phosphoglycerate kinase (EC 2.7.2.3)",
  REVIEWED (Swiss-Prot), catalytic activity "3-phosphoglycerate + ATP =
  3-phospho-glyceroyl phosphate + ADP", GO terms for BOTH "glycolytic process" AND
  "gluconeogenesis".

The reviewed PGK record explicitly assigns the enzyme to both directions with a reversible
reaction, which contradicts a biological one-way gluconeogenic constraint. Wang et al. 2012
also demonstrates that exogenous glucose is assimilated (requires the G3P -> 1,3-BPG -> 3PG
direction to proceed). Esparza et al. 2010 and Osorio et al. 2013 show ATCC 23270 is an
obligate autotroph running the CBB cycle, which explains why a reconstruction would write
central carbon in the gluconeogenic direction — as a reflection of autotrophic carbon flow,
not enzyme irreversibility.

**Does the 2016 irreversible gluconeogenic direction appear biologically supported?**
**LIKELY_MODEL_CONSTRAINT.** The strongest evidence is the reviewed bidirectional PGK
annotation (both GO directions + reversible reaction) together with the demonstrated
glucose assimilation; no evidence supports one-way enzyme directionality.

## Source priority applied

Primary weight was given to: Wang et al. 2012 (direct ATCC 23270 glucose experiment),
Valdés et al. 2008 (ATCC 23270 genome), Esparza et al. 2010 / 2019 (carbon fixation),
and Osorio et al. 2013 (central-carbon vs CBB regulation). Closely-related-species papers
(Acidithiobacillus Ameehan, Thiomonas, Leptospirillum) were used only as background.

## Files

- phase4_1G_2A_rapid_evidence_report.md (this file)
- phase4_1G_2A_paper_evidence.tsv (9 paper records with extraction fields)
- phase4_1G_2A_unresolved_questions.tsv (5 open items)
