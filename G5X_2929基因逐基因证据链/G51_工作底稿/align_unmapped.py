"""Compare the current G51 proteins lacking exact-WP AFE mapping to frozen 2008 proteins."""
import json
from collections import Counter, defaultdict
from pathlib import Path
from Bio.Align import PairwiseAligner

here = Path(__file__).resolve().parent
root = Path(__file__).resolve().parents[2] / "01_原始数据" / "01_基因组与官方注释"
old_dir = root / "GCF_000021485.1_20260922"
new_dir = root / "GCF_049532655.1_20260922"
inventory = json.loads((here / "G51_evidence_inventory.json").read_text(encoding="utf-8"))

def fasta(path):
    header = ""
    chunks = []
    for line in path.open(encoding="utf-8"):
        if line.startswith(">"):
            if header:
                yield header.split()[0][1:], "".join(chunks)
            header, chunks = line.strip(), []
        else:
            chunks.append(line.strip())
    if header:
        yield header.split()[0][1:], "".join(chunks)

old_protein = dict(fasta(next(old_dir.rglob("protein.faa"))))
new_protein = dict(fasta(next(new_dir.rglob("protein.faa"))))
old_tags = {}
for line in next(old_dir.rglob("genomic.gff")).open(encoding="utf-8"):
    if line.startswith("#"):
        continue
    parts = line.strip().split("\t")
    if len(parts) != 9 or parts[2] != "CDS":
        continue
    attrs = dict(pair.split("=", 1) for pair in parts[8].split(";") if "=" in pair)
    protein_id = attrs.get("protein_id")
    if protein_id:
        old_tags[protein_id] = {"old_refseq_locus": attrs.get("locus_tag", ""), "old_product": attrs.get("product", ""),
                                "old_start": int(parts[3]), "old_end": int(parts[4])}
gene_to_afe = {}
for line in next(old_dir.rglob("genomic.gff")).open(encoding="utf-8"):
    if line.startswith("#"):
        continue
    parts = line.strip().split("\t")
    if len(parts) != 9 or parts[2] != "gene":
        continue
    attrs = dict(pair.split("=", 1) for pair in parts[8].split(";") if "=" in pair)
    gene_to_afe[attrs.get("locus_tag")] = attrs.get("old_locus_tag", "")

kmers = defaultdict(set)
for old_id, seq in old_protein.items():
    for i in range(len(seq) - 5):
        kmers[seq[i:i+6]].add(old_id)

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

out = []
for gene in inventory:
    if gene["old_afe"]:
        continue
    seq = new_protein.get(gene["wp"])
    if not seq:
        out.append({"locus": gene["locus"], "error": "current WP absent from frozen protein.faa"})
        continue
    counts = Counter()
    for i in range(len(seq) - 5):
        for old_id in kmers.get(seq[i:i+6], ()):
            counts[old_id] += 1
    candidates = []
    for old_id, shared in counts.most_common(15):
        old_seq = old_protein[old_id]
        aln = aligner.align(seq, old_seq)[0]
        stats = aln.counts()
        aligned = stats.identities + stats.mismatches
        meta = old_tags.get(old_id, {})
        candidates.append({"old_wp": old_id, "old_afe": gene_to_afe.get(meta.get("old_refseq_locus", ""), ""),
                           "old_refseq_locus": meta.get("old_refseq_locus", ""), "old_product": meta.get("old_product", ""),
                           "old_start": meta.get("old_start"), "old_end": meta.get("old_end"),
                           "old_length": len(old_seq), "current_length": len(seq), "identity_fraction": round(stats.identities / aligned, 5) if aligned else 0,
                           "current_coverage": round(aligned / len(seq), 5), "old_coverage": round(aligned / len(old_seq), 5),
                           "aligned_aa": aligned, "matching_aa": stats.identities, "shared_6mer_count": shared})
    candidates.sort(key=lambda x: (x["identity_fraction"]*min(x["current_coverage"],x["old_coverage"]), x["aligned_aa"]), reverse=True)
    out.append({"locus": gene["locus"], "current_wp": gene["wp"], "product": gene["product"], "candidates": candidates[:5]})
(here / "G51_unmapped_sequence_alignment.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
for row in out:
    print(row["locus"], row.get("candidates", [])[:2])
