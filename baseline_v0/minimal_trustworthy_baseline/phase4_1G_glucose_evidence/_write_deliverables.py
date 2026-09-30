# -*- coding: utf-8 -*-
"""Write the two TSV deliverables."""
import csv
from pathlib import Path

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4_1G_glucose_evidence")

PAPERS = [
    ["Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications",
     "2008", "10.1186/1471-2164-9-597", "19077236",
     "A. ferrooxidans ATCC 23270 (type strain)", "whole genome (AFE_ loci)", "genome sequencing + annotation",
     "Central-carbon-metabolism model; carbohydrate-porin and transporter genes annotated only",
     "NO", "NO", "NO",
     "First complete ATCC 23270 genome; presents initial in-silico central-carbon model but assigns only predicted (bioinformatic) functions.",
     "No enzymatic or transport assay; all glucose-related assignments are prediction."],

    ["Development of a markerless gene replacement system for Acidithiobacillus ferrooxidans and construction of a pfkB mutant",
     "2012", "10.1128/aem.07230-11", "22210219",
     "A. ferrooxidans ATCC 23270", "AFE_1807 (pfkB)", "gene knockout + heterologous expression + RT-qPCR",
     "AFE_1807 encodes a functional phosphofructokinase (PFK-B); ATCC 23270 can assimilate exogenous glucose",
     "YES", "NO", "NO",
     "pfkB (AFE_1807) complemented PFK-deficient E. coli (weak, ~800-fold lower than E. coli pfkB); pfkB mutant grew slower on S(0) but still assimilated exogenous glucose; RT-qPCR of central-carbohydrate genes with/without glucose.",
     "Identifies PFK (fructose-6-P -> FBP), NOT a glucose kinase or glucose transporter; does not name the uptake system."],

    ["Genes and pathways for CO2 fixation in the obligate, chemolithoautotrophic acidophile, Acidithiobacillus ferrooxidans",
     "2010", "10.1186/1471-2180-10-229", "20799944",
     "A. ferrooxidans ATCC 23270", "cbb operons (RubisCO, PRK)", "operon mapping + RT-PCR + EMSA",
     "ATCC 23270 fixes CO2 via the Calvin-Benson-Bassham (CBB) cycle (form I RubisCO + carboxysomes); obligate chemolithoautotroph",
     "NO", "NO", "NO",
     "Four cbb operons encode CBB enzymes; CbbR regulates cbb1-3; confirms the organism's carbon flows through autotrophic CO2 fixation.",
     "Carbon-flow context only; does not test glycolysis direction or glucose use."],

    ["Effect of CO2 Concentration on Uptake and Assimilation of Inorganic Carbon in the Extreme Acidophile Acidithiobacillus ferrooxidans",
     "2019", "10.3389/fmicb.2019.00603", "31019493",
     "A. ferrooxidans ATCC 23270", "cbb operons, can-sulP", "RT-qPCR + Western blot",
     "Inorganic-carbon (CO2/HCO3-) uptake and CBB assimilation; no classical bicarbonate transporter (sulP proposed instead)",
     "NO", "NO", "NO",
     "Transcript abundance of five cbb operons and can-sulP tracked vs CO2; conceptual model of CO2 fixation; reinforces autotrophic (CO2) carbon supply.",
     "Inorganic carbon only; does not address organic-carbon (glucose) uptake or glycolysis."],

    ["Anaerobic sulfur metabolism coupled to dissimilatory iron reduction in the extremophile Acidithiobacillus ferrooxidans",
     "2013", "10.1128/aem.03057-12", "23354702",
     "A. ferrooxidans ATCC 23270", "central-carbon + CBB genes", "microarray transcriptomics + proteomics",
     "Central carbon (glycolytic) pathways up under aerobic growth; CBB components up under anaerobic growth",
     "NO", "NO", "NO",
     "Aerobic growth correlates with upregulated central-carbon pathways and higher growth; anaerobic with CBB — showing EMP and CBB are conditionally regulated, not one-way.",
     "Does not measure individual GAPDH/PGK direction; uses sulfur-grown cells, not glucose."],

    ["Proteomics Reveal Enhanced Oxidative Stress Responses and Metabolic Adaptation in Acidithiobacillus ferrooxidans Biofilm Cells on Pyrite",
     "2019", "10.3389/fmicb.2019.00592", "30984136",
     "A. ferrooxidans DSM 14882 (= ATCC 23270)", "zwf (AFE_2025) among detected proteins", "shotgun proteomics",
     "G6PDH (zwf) protein detected; carbon-fixation and oxidative-phosphorylation adaptations reported",
     "NO", "NO", "NO",
     "1157 proteins quantified; zwf/AFE_2025 appears in the proteome, consistent with a functional oxidative-PPP gene; no activity assay.",
     "Presence/abundance only; does not demonstrate glucose -> G6P -> oxidative PPP flux."],

    ["Insights into the Quorum Sensing Regulon of the Acidophilic Acidithiobacillus ferrooxidans Revealed by Transcriptomic in the Presence of an AHL Superagonist Analog",
     "2016", "10.3389/fmicb.2016.01365", "27683573",
     "A. ferrooxidans ATCC 23270", "QS regulon incl. zwf (AFE_2025)", "DNA microarray",
     "zwf (AFE_2025) among QS-regulated genes (biofilm context)",
     "NO", "NO", "NO",
     "141 genes (4.5% of genome) QS-regulated; zwf/G6PDH implicated in QS/biofilm response, not glucose catabolism.",
     "Regulation only; no directionality or kinetic evidence for G6PDH."],

    ["Characterize the Growth and Metabolism of Acidithiobacillus ferrooxidans under Electroautotrophic and Chemoautotrophic Conditions",
     "2024", "10.3390/microorganisms12030590", "38543641",
     "A. ferrooxidans ATCC 23270", "porin / galactose-metabolism genes", "RNA-seq + metabolomics",
     "Galactose metabolism and porin/pilin expression enhanced during electroautotrophic growth (EPS production)",
     "NO", "NO", "NO",
     "493 DEGs; increased porin/transmembrane-protein expression under electroautotrophy; galactose pathway up (EPS), not glucose uptake.",
     "Not a glucose-transport test; porin regulation is context-dependent (electrode adhesion)."],

    ["Genome sequencing and metabolic network reconstruction of a novel sulfur-oxidizing bacterium Acidithiobacillus Ameehan",
     "2023", "10.3389/fmicb.2023.1277847", "38053556",
     "Acidithiobacillus Ameehan (NOT ATCC 23270)", "central-carbon genes", "genome + GEM + Biolog",
     "Ameehan cannot grow on glucose or yeast extract as carbon source; incomplete CBB cycle",
     "YES", "NO", "NO",
     "Related Acidithiobacillus sp. unable to use glucose; 88.7% Biolog agreement; supports that Acidithiobacillus spp. are largely obligate autotrophs.",
     "Different species; used only as genus-level context, not ATCC 23270 evidence."],
]

UNRESOLVED = [
    ["glucose transporter identity",
     "Which of AFE_2522/AFE_2250 (OprB porins), AFE_2312 (MFS sugar porter), AFE_1971 (MFS) actually transports glucose; no direct glucose-transport assay exists",
     "carbohydrate/sugar transport assay, porin substrate-specificity, or glucose-induced expression of a specific transporter in ATCC 23270"],
    ["PTS-like cluster completeness (AFE_3018-AFE_3023 / RU820_RS14005-14030)",
     "Whether the PTS-like cluster encodes a complete glucose PTS (EI/HPr/EIIA/EIIB/EIIC) or is degenerate",
     "gene-cluster inspection + demonstration of a membrane IIC/IIB glucose component"],
    ["AFE_2841 glucokinase activity",
     "ROK-family protein annotation only; no EC 2.7.1.2 assignment and no glucose-phosphorylation assay",
     "biochemical glucose kinase assay, glucose-induced expression, or knockout phenotype of AFE_2841"],
    ["Zwf (AFE_2025) glucose oxidation",
     "Strong EC 1.1.1.49 annotation but no ATCC 23270-specific demonstration that glucose -> G6P -> oxidative PPP flux occurs",
     "G6PDH enzyme assay or 13C-glucose labelling to confirm oxidative-PPP flux"],
    ["GAPDH (AFE_3251) directionality",
     "Type-I GAPDH annotation implies reversibility, but no ATCC 23270 assay of G3P -> 1,3-BPG",
     "enzyme assay or glucose-fed flux/labelling study in ATCC 23270"],
]


def main():
    pf = ["title", "year", "doi", "pmid", "organism_strain", "target_gene_locus",
          "experiment_type", "exact_supported_function", "glucose_directly_tested",
          "transport_directly_tested", "enzyme_direction_tested", "key_result_1_3_sentences", "limitation"]
    with (OUT / "phase4_1G_2A_paper_evidence.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(pf)
        for row in PAPERS:
            w.writerow(row)

    uf = ["question", "issue", "evidence_needed"]
    with (OUT / "phase4_1G_2A_unresolved_questions.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(uf)
        for row in UNRESOLVED:
            w.writerow(row)

    print("wrote", len(PAPERS), "paper rows;", len(UNRESOLVED), "unresolved rows")


if __name__ == "__main__":
    main()
