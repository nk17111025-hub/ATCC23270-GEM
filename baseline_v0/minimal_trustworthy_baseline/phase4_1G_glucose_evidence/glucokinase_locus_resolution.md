# GLUCOKINASE LOCUS RESOLUTION

## Question

Which ATCC 23270 locus was called "glucokinase" in the genome / central-carbon
reconstruction literature, and is AFE_2841 truly the historical glucokinase candidate?

## Finding

- Legacy locus: **AFE_2841** (old RefSeq locus **AFE_RS13045**).
- Current mapping: **RU820_RS13140** → protein **WP_012537424.1**.
- The 2016 GEM (Campodonico) BDGK reaction carried **no GPR** (empty).
- The 2024 GEM (Khaleque) assigned **AFE_2841** as the BDGK GPR (candidate E4-1-BDGK-GPR).
- That assignment rests on cross-database annotation only:
  - KEGG `afr:AFE_2841` = **K25026** glucokinase, EC **2.7.1.2** (RHEA:17825, via KEGG R00299 — a reaction link, not a gene experiment).
  - BioCyc PGDB (GCF_000021485) `AFE_RS13045-MONOMER`.
  - UniProt "ROK family protein" (TrEMBL, unreviewed, no EC line).
- No Valdés-2008 explicit "glucokinase" gene call was found in the project copies; the
  "glucokinase" label is a later inference from the ROK-family / KEGG K25026 annotation.

## Catalytic evidence

None. For AFE_2841 there is no purified-enzyme assay, no substrate (beta-D-glucose) or
phosphate-donor (ATP) specificity test, no Km/Vmax/kcat, no G6P anomer test, no
knockout/complementation, and no glucose-induced expression measurement. The Wang 2012
glucose experiment (PMID 22210219) studied **AFE_1807 / pfkB** (phosphofructokinase), not
AFE_2841, and must not be transferred.

## Conclusion

**GLUCOKINASE_CANDIDATE_ONLY**

AFE_2841 / RU820_RS13140 / WP_012537424.1 is the current glucokinase candidate by
ROK-family/KEGG annotation, but it is not experimentally validated. BDGK GPR remains
unresolved (GLUCOKINASE_GPR_UNRESOLVED); no catalytic direction can be asserted from this locus.
