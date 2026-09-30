"""Create Alpha-style source register for the G51 working review."""
import concurrent.futures
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
curated = json.loads((HERE / "curated_findings.json").read_text(encoding="utf-8"))
DOIS = sorted({doi.strip() for item in curated.values() for doi in item.get("evidence_doi", "").split(";") if doi.strip()})
CACHE = HERE / "crossref_sources_cache.json"
old = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}


def fetch(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "G51-evidence-review/1.0 (research data)"})
            with urllib.request.urlopen(req, timeout=30) as response:
                return doi, json.load(response)["message"]
        except Exception as exc:
            if attempt == 3:
                return doi, {"error": str(exc)}
            time.sleep(2 ** attempt)


with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for doi, data in pool.map(fetch, [d for d in DOIS if d not in old]):
        old[doi] = data
CACHE.write_text(json.dumps(old, ensure_ascii=False, indent=2), encoding="utf-8")

def item(source_id, citation, year, doi, source_type, strain, experiment, genes, reactions, measured, score, limit, url):
    return {"Source ID": source_id, "Full citation": citation, "Year": year, "DOI": doi, "Source type": source_type,
            "Strain": strain, "Experimental type": experiment, "Genes / proteins studied": genes,
            "Related Reaction IDs / candidate IDs": reactions, "What was directly measured": measured,
            "Applicable confidence level": score, "Evidence limitation": limit, "URL": url}

sources = [
    item("G51-S001", "NCBI RefSeq GCF_049532655.1, genomic.gff, annotation RS_2025_04_13", 2025, "", "Official annotation", "ATCC 23270", "Genome annotation", "RU820_RS00005–RU820_RS01505", "G51 identity", "Chromosome positions, locus, WP and product annotations", "Not a reaction score", "Protein function labels are computational annotations.", "https://drive.google.com/file/d/19ijnng4VELvWQH8QID6cII9TqYMUhes4/view"),
    item("G51-S002", "B1 old AFE to current RU820 mapping, project verified table", 2026, "", "Sequence/identity mapping", "ATCC 23270", "Protein ID exact matches; ambiguous mappings retained", "G51 old AFE", "Legacy tracking", "Old/current protein ID relationships", "Not a reaction score", "Protein identity does not prove catalysis; nonunique WP mapping is unresolved.", "B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv"),
    item("G51-S003", "Campodonico et al. 2016 iMC507 original Supplementary Table 1", 2016, "10.1016/j.meteno.2016.03.003", "Model supplementary table", "ATCC 23270", "Model reconstruction", "Original AFE GPRs", "55 distinct G51-associated reactions", "Reaction formula, score, EC, PMID, GPR and three media bounds", "Original scores only", "Old model score does not establish newly proposed gene links or revised chemistry.", "https://drive.google.com/file/d/1u23dFDLoG27LJvEi6BjYKIK3Tr-zR2-f/view"),
    item("G51-S004", "UniProtKB ATCC 23270 taxon 243159 protein records, accessed 2026-09-25", 2026, "", "Protein function database", "ATCC 23270", "Reviewed/unreviewed annotation and cross-references", "WP matched to UniProt accession", "All genes with UniProt hit", "Protein names, evidence codes, EC, InterPro/Pfam links", "2 for sequence-based precise reaction only", "Reviewed entries can still use sequence rules; paper evidence must be checked per assertion.", "https://rest.uniprot.org/uniprotkb/search"),
    item("G51-S005", "InterPro/Pfam family records from UniProtKB cross-references, accessed 2026-09-25", 2026, "", "Protein domain database", "Multiple organisms", "Domain/family classification", "G51 UniProt-linked proteins", "Domain candidates", "Family/domain accessions", "2 at most for precise reaction", "Family identity alone does not fix substrate, location or direction.", "https://www.ebi.ac.uk/interpro/"),
    item("G51-S006", "KEGG afr historical gene/KO/EC/reaction exports, project frozen 2026-09-22", 2026, "", "Metabolic annotation database", "ATCC 23270 old assembly", "Computational pathway mapping", "Old AFE matched by B1", "Candidate reaction IDs", "Gene-to-KO/EC/reaction links", "2 at most", "Old locus/version and inherited annotations require independent checking.", "https://www.kegg.jp/entry/afr"),
    item("G51-S007", "BioCyc AFER243159 PGDB 30.0 project export", 2026, "", "Metabolic annotation database", "ATCC 23270 old assembly", "PathoLogic computational predictions", "Old AFE matched by B1", "Candidate reaction IDs", "Old gene/protein/reaction associations", "2 at most", "BioCyc/MetaCyc Pathway Tools links share provenance and are not independent experiments.", "https://biocyc.org/AFER243159/"),
    item("G51-S008", "Rhea reaction definitions and ChEBI cross-references, project frozen 2026-09-23", 2026, "", "Reaction chemistry database", "General enzyme definitions", "Curated reaction stoichiometry", "EC-linked candidate proteins", "RHEA accessions", "Equations, ChEBI IDs and EC links", "No gene-specific score", "Exact chemical definition does not prove this locus catalyzes it.", "https://www.rhea-db.org/"),
    item("G51-S009", "BRENDA organism/EC project exports, accessed 2026-09-23", 2026, "", "Enzyme database", "Organism-dependent", "EC and reaction records", "G51 EC candidates", "Candidate enzyme chemistry", "EC-level data", "Depends on underlying paper", "Organism-level records cannot be assigned to an exact locus without gene/protein evidence.", "https://www.brenda-enzymes.org/"),
    item("G51-S010", "Europe PMC per-gene ID and functional phrase searches, accessed 2026-09-25", 2026, "", "Literature index", "Mixed strains", "Search index only", "293 G51 loci", "Paper leads", "DOI/title/abstract candidates", "No reaction score", "Zero hits do not prove no study; hits require strain, gene and assay verification.", "https://www.ebi.ac.uk/europepmc/webservices/rest/search"),
]

context = {
 "10.1016/j.jbiotec.2007.08.030": ("ATCC 23270", "Native enzyme purification, N-terminal protein identification", "AFE_0029/TetH", "4THASE1;4THASE2", "Tetrathionate-dependent activity and purified protein N-terminus", "4 for measured hydrolysis; exact alternative formula ratios require check", "No proof for all legacy product ratios."),
 "10.1128/JB.01472-13": ("ATCC 23270", "tetH knockout, complementation, overexpression and qPCR", "AFE_0029;AFE_0048", "4THASE1;4THASE2;TSQOC", "Tetrathionate-growth phenotype and doxD transcript changes", "3 for TetH physiological hydrolysis", "Transcription of doxD does not prove its precise reaction or complex composition."),
 "10.1128/AEM.02251-12": ("ATCC 23270", "Purified/recombinant thiosulfate dehydrogenase assay", "AFE_0042", "TSQOC", "Thiosulfate-to-tetrathionate assay with ferricyanide; quinone negative", "4 for the measured artificial-acceptor assay", "Does not validate q8-dependent 2016 TSQOC equation."),
 "10.1007/s00253-014-5830-4": ("ATCC 23270; A. caldus MTH-04", "Recombinant enzyme assay and ATCC 23270 knockout/overexpression", "AFE_0269", "RHEA:12981;SULDO", "GSH-dependent SDO activity, sulfur-growth/activity effects", "4 for measured assay; lower for exact GSSH model equation", "GSSH in vivo and compartment were not established."),
 "10.1074/jbc.M700374200": ("Other bacterial PlsY", "Membrane topology and active-site biochemical study", "PlsY homolog", "G3PAT1 candidate", "Acylphosphate:glycerol-3-phosphate acyltransferase specificity", "2 for ATCC 23270 inference", "Other organism; no ATCC 23270 substrate assay."),
 "10.1021/bi00136a016": ("Escherichia coli", "Purified PurE/PurK/PurC biochemistry", "PurE;PurK homologs", "AIRC;RHEA:13193;RHEA:19317", "Two-step ATP-dependent AIR-to-CAIR chemistry", "2 for ATCC 23270 inference", "Other organism; no direct ATCC 23270 assay."),
 "10.3389/fmicb.2019.00603": ("ATCC 23270", "Whole-cell CO2 uptake and assimilation", "AFE_0287 prediction", "HCO3E", "CO2 uptake/assimilation under different conditions", "2 for HCO3E annotation", "Whole-cell rates do not demonstrate the AFE_0287 catalytic step."),
 "10.1074/jbc.M113.480368": ("Escherichia coli", "UbiI gene perturbation and ubiquinone intermediate analysis", "UbiI/UbiB homolog context", "OPHHX", "UbiI hydroxylation role and correction of old UbiB assignment", "2 for ATCC 23270 functional inference", "Other organism; does not identify the ATCC 23270 hydroxylase."),
}

for n, doi in enumerate(DOIS, 1):
    meta = old[doi]
    authors = ", ".join(" ".join(filter(None, [a.get("given", ""), a.get("family", "")])) for a in meta.get("author", []))
    title = (meta.get("title") or [doi])[0]
    journal = (meta.get("container-title") or [""])[0]
    year = (meta.get("published", {}).get("date-parts") or [[""]])[0][0]
    citation = f"{authors} ({year}). {title}. {journal}. https://doi.org/{doi}"
    strain, exp, genes, reactions, measured, score, limit = context[doi]
    sources.append(item(f"G51-P{n:03d}", citation, year, doi, "Original research", strain, exp, genes, reactions, measured, score, limit, "https://doi.org/" + doi))

(HERE / "G51_sources.json").write_text(json.dumps(sources, ensure_ascii=False, indent=2), encoding="utf-8")
(HERE / "G51_doi_source_ids.json").write_text(json.dumps({d: f"G51-P{i:03d}" for i, d in enumerate(DOIS, 1)}, ensure_ascii=False, indent=2), encoding="utf-8")
print("sources", len(sources), "DOIs", len(DOIS), "errors", [d for d in DOIS if "error" in old[d]])
