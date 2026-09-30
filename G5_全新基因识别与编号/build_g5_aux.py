#!/usr/bin/env python3
"""Rebuild G5 auxiliary old-protein alignment and synteny tables.

Run from the project root with:
    python G5_全新基因识别与编号/build_g5_aux.py

Requires Biopython. This script only writes the two G5_辅助*.tsv files.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


ROOT = Path.cwd()
G5_DIR = next(p for p in ROOT.glob("G5_*") if p.is_dir())
B1_DIR = next(p for p in ROOT.glob("B1_*") if p.is_dir())
DATA_DIR = B1_DIR / "标准化数据"
NEW_DIR = DATA_DIR / "GCF_049532655.1"
OLD_DIR = DATA_DIR / "GCF_000021485.1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def parse_gene_gff(path: Path) -> list[dict[str, object]]:
    genes = []
    for line in path.open(encoding="utf-8"):
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
        genes.append({
            "seq": fields[0], "start": int(fields[3]), "end": int(fields[4]),
            "strand": fields[6], "id": attrs.get("ID", ""),
            "locus": attrs.get("locus_tag", ""),
            "old": attrs.get("old_locus_tag", ""),
            "gene": attrs.get("gene", ""),
        })
    return genes


def protein_parent_map(path: Path) -> dict[str, list[str]]:
    parents: dict[str, list[str]] = defaultdict(list)
    for line in path.open(encoding="utf-8"):
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
        protein, parent = attrs.get("protein_id"), attrs.get("Parent")
        if protein and parent:
            parents[protein].append(parent)
    return parents


def build_alignment_rows(audit_rows: list[dict[str, str]]) -> list[list[str]]:
    new_faa = next(NEW_DIR.glob("protein.faa"))
    old_faa = next(OLD_DIR.glob("protein.faa"))
    new_proteins = {record.id.split("|")[0]: str(record.seq)
                    for record in SeqIO.parse(new_faa, "fasta")}
    old_records = list(SeqIO.parse(old_faa, "fasta"))
    old_proteins = [str(record.seq) for record in old_records]
    old_ids = [record.id.split("|")[0] for record in old_records]
    old_gene_by_id = {str(gene["id"]): gene for gene in parse_gene_gff(OLD_DIR / "genomic.gff")}
    old_parents = protein_parent_map(OLD_DIR / "genomic.gff")

    # Five-residue exact seeds narrow local alignment to a small candidate set.
    seed_index: dict[str, list[int]] = defaultdict(list)
    for index, sequence in enumerate(old_proteins):
        for seed in {sequence[i:i + 5] for i in range(max(0, len(sequence) - 4))}:
            seed_index[seed].append(index)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    output = []
    for row in audit_rows:
        locus, protein_id = row["当前locus"], row["当前protein ID"]
        query = new_proteins.get(protein_id)
        if not query:
            continue
        seed_counts: dict[int, int] = defaultdict(int)
        for i in range(max(0, len(query) - 4)):
            for index in seed_index.get(query[i:i + 5], ()):
                seed_counts[index] += 1

        best = None
        for index in sorted(seed_counts, key=seed_counts.get, reverse=True)[:30]:
            subject = old_proteins[index]
            alignment = aligner.align(query, subject)[0]
            identical = aligned = 0
            for (q1, q2), (s1, s2) in zip(*alignment.aligned):
                length = min(q2 - q1, s2 - s1)
                aligned += length
                identical += sum(query[q1 + offset] == subject[s1 + offset]
                                 for offset in range(length))
            if not aligned:
                continue
            identity = identical / aligned
            min_coverage = aligned / min(len(query), len(subject))
            current_coverage = aligned / len(query)
            old_coverage = aligned / len(subject)
            candidate = (alignment.score, identity, min_coverage, aligned,
                         old_ids[index], current_coverage, old_coverage)
            if best is None or candidate[0] > best[0]:
                best = candidate
        if best:
            old_loci = sorted({str(old_gene_by_id[parent]["locus"])
                               for parent in old_parents.get(best[4], ())
                               if parent in old_gene_by_id})
            output.append([locus, protein_id, best[4], ";".join(old_loci), f"{best[1]:.4f}",
                           f"{best[2]:.4f}", str(best[3]), f"{best[5]:.4f}",
                           f"{best[6]:.4f}", f"{best[0]:.1f}"])
    return output


def build_neighborhood_rows(audit_rows: list[dict[str, str]],
                            alignment_rows: list[list[str]]) -> list[list[str]]:
    current_genes = parse_gene_gff(NEW_DIR / "genomic.gff")
    old_genes = parse_gene_gff(OLD_DIR / "genomic.gff")
    current_by_locus = {str(g["locus"]): g for g in current_genes}
    old_by_locus = {str(g["locus"]): g for g in old_genes}
    old_by_id = {str(g["id"]): g for g in old_genes}
    old_parent_map = protein_parent_map(OLD_DIR / "genomic.gff")

    # Only unambiguous B1 exact mappings are anchors.
    b1_rows = read_tsv(next((B1_DIR / "标准化数据").glob("旧新基因映射.tsv")))
    anchors: dict[str, str] = {}
    for row in b1_rows:
        old_locus, current_locus = row.get("旧RefSeq位点", ""), row.get("当前GCF049位点", "")
        if (row.get("映射状态") == "完全匹配" and old_locus.startswith("AFE_RS")
                and "|" not in old_locus + current_locus and current_locus in current_by_locus):
            anchors[current_locus] = old_locus
    old_anchor = {current: old_by_locus.get(old) for current, old in anchors.items()}

    current_order: dict[str, list[dict[str, object]]] = defaultdict(list)
    old_order: dict[str, list[dict[str, object]]] = defaultdict(list)
    for gene in current_genes:
        current_order[str(gene["seq"])].append(gene)
    for gene in old_genes:
        old_order[str(gene["seq"])].append(gene)
    for genes in (*current_order.values(), *old_order.values()):
        genes.sort(key=lambda gene: int(gene["start"]))
    current_positions = {str(gene["locus"]): (genes, index)
                         for genes in current_order.values()
                         for index, gene in enumerate(genes)}
    alignment_by_locus = {row[0]: row for row in alignment_rows}

    results = []
    for audit in audit_rows:
        locus = audit["当前locus"]
        genes, position = current_positions[locus]
        left = right = None
        for step in range(1, len(genes) + 1):
            candidate = genes[(position - step) % len(genes)]
            if candidate["locus"] in anchors and old_anchor.get(str(candidate["locus"])):
                left = candidate
                break
        for step in range(1, len(genes) + 1):
            candidate = genes[(position + step) % len(genes)]
            if candidate["locus"] in anchors and old_anchor.get(str(candidate["locus"])):
                right = candidate
                break

        alignment = alignment_by_locus.get(locus)
        top_protein = alignment[2] if alignment else ""
        candidates = []
        for parent in old_parent_map.get(top_protein, ()):
            gene = old_by_id.get(parent)
            if gene and gene not in candidates:
                candidates.append(gene)

        inside = []
        interval = ""
        interval_available = False
        old_left = old_anchor.get(str(left["locus"])) if left else None
        old_right = old_anchor.get(str(right["locus"])) if right else None
        if old_left and old_right and old_left["seq"] == old_right["seq"]:
            ordered = old_order[str(old_left["seq"])]
            left_index = next((i for i, gene in enumerate(ordered)
                               if gene["locus"] == old_left["locus"]), None)
            right_index = next((i for i, gene in enumerate(ordered)
                                if gene["locus"] == old_right["locus"]), None)
            if left_index is not None and right_index is not None:
                between = set()
                index = (left_index + 1) % len(ordered)
                while index != right_index:
                    between.add(str(ordered[index]["locus"]))
                    index = (index + 1) % len(ordered)
                    if len(between) > len(ordered):
                        break
                inside = [gene for gene in candidates if str(gene["locus"]) in between]
                interval = f"{old_left['locus']}..{old_right['locus']} ({len(between)}旧基因)"
                interval_available = True

        def describe(gene: dict[str, object]) -> str:
            return f"{gene['locus']}[old={gene['old'] or '空'}]"

        candidate_text = ";".join(describe(gene) for gene in candidates) or "无旧CDS候选"
        inside_text = ";".join(describe(gene) for gene in inside) or "无"
        left_text = f"{left['locus']}->{anchors[str(left['locus'])]}" if left else "无"
        right_text = f"{right['locus']}->{anchors[str(right['locus'])]}" if right else "无"
        identity = float(alignment[4]) if alignment else 0.0
        coverage = float(alignment[5]) if alignment else 0.0

        if len(inside) == 1:
            consistent = "是：top候选在锚定区间"
            advice = ("建议排除：序列top hit唯一落在双侧锚区间，旧AFE对应"
                      if str(inside[0]["old"]).startswith("AFE_") else
                      "建议排除：top hit唯一落在双侧锚区间，旧old_locus_tag为空/非AFE")
        elif len(inside) > 1:
            consistent, advice = "部分一致：多个候选落在锚定区间", "保留人工复核：锚区间内有多个同WP旧locus"
        elif candidates and interval_available:
            consistent, advice = "否：top候选均在区间外", "待复核：序列候选均不在锚定区间，可能旁系同源/结构变化"
        elif interval_available:
            consistent, advice = "无序列候选", "初步支持新locus：旧蛋白top hit无旧CDS定位"
        else:
            consistent, advice = "无法判定", "待复核：两侧锚点无法建立同一旧contig区间"

        results.append([locus, audit["当前protein ID"], left_text, right_text, interval,
                        candidate_text, inside_text, consistent, advice,
                        f"{identity:.4f}", f"{coverage:.4f}"])
    return results


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    d3_dir = next(p for p in ROOT.glob("D3_*") if p.is_dir())
    d3_table = next(d3_dir.rglob("D3-1-R_全基因跨数据库交叉引用.tsv"))
    audit_rows = [row for row in read_tsv(d3_table)
                  if row["基因类型"] == "protein_coding"
                  and row["B1映射状态"] != "完全匹配"]
    assert len(audit_rows) == 324
    alignment_rows = build_alignment_rows(audit_rows)
    alignment_header = ["当前locus", "当前protein ID", "top旧protein ID", "旧AFE候选(旧GFF locus_tag)",
                        "identity", "短序列覆盖(按min长度)", "比对氨基酸数", "当前序列覆盖", "旧序列覆盖",
                        "BLOSUM62局部比对分数"]
    write_tsv(G5_DIR / "G5_辅助序列比对.tsv", alignment_header, alignment_rows)
    neighborhood_rows = build_neighborhood_rows(audit_rows, alignment_rows)
    neighborhood_header = ["当前locus", "当前protein ID", "左侧最近B1锚点(current→old AFE_RS)",
                           "右侧最近B1锚点(current→old AFE_RS)", "旧基因组锚点区间",
                           "序列top旧locus候选及old_locus_tag", "落入锚区间候选", "邻域是否一致",
                           "可判分类建议", "identity", "短序列覆盖"]
    write_tsv(G5_DIR / "G5_辅助邻域核验.tsv", neighborhood_header, neighborhood_rows)

    high = [row for row in neighborhood_rows
            if float(row[9]) >= 0.90 and float(row[10]) >= 0.90]
    concordant = [row for row in high if row[7].startswith("是：")]
    afe_tag = [row for row in concordant if "[old=AFE_" in row[6]]
    no_afe_tag = len(concordant) - len(afe_tag)
    print(f"alignment rows: {len(alignment_rows)}")
    print(f"neighborhood rows: {len(neighborhood_rows)}")
    print(f"high similarity (identity/short coverage >= 0.90): {len(high)}")
    print(f"high similarity and neighborhood concordant: {len(concordant)}")
    print(f"  explicit old_locus_tag AFE_xxxx: {len(afe_tag)}")
    print(f"  old_locus_tag empty/non-AFE: {no_afe_tag}")
    print(f"high-similarity neighborhood conflicts/ambiguous: {len(high) - len(concordant)}")


if __name__ == "__main__":
    main()
