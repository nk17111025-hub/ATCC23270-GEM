#!/usr/bin/env python3
"""Audit 324 current loci against the original GCA AFE_xxxx/ACK annotation.

Run from project root:
    python G5_*/build_g5_raw_afe.py

Requires Biopython. Writes only G5_\u8f85\u52a9\u539f\u59cbAFE\u6838\u9a8c.tsv.
The old GCA and GCF chromosome FASTAs are asserted to be the same 2,982,397 bp
sequence before projecting GCF coordinates to original GCA CDS features.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

ROOT = Path.cwd()
G5_DIR = next(p for p in ROOT.glob("G5_*") if p.is_dir())
B1_DIR = next(p for p in ROOT.glob("B1_*") if p.is_dir())
STD_DIR = next(p for p in B1_DIR.iterdir()
               if p.is_dir() and (p / "GCF_049532655.1").is_dir())
CUR_DIR = STD_DIR / "GCF_049532655.1"
OLD_REF_DIR = STD_DIR / "GCF_000021485.1"
ORIGINAL_DATA_DIR = next(p for p in ROOT.glob("01_*") if p.is_dir())
OLD_GCA_DIR = next(p for p in ORIGINAL_DATA_DIR.rglob("GCA_000021485.1")
                   if p.is_dir() and (p / "genomic.gff").exists()
                   and (p / "protein.faa").exists())
OLD_GCF_FASTA = next(ORIGINAL_DATA_DIR.rglob("GCF_000021485.1_ASM2148v1_genomic.fna"))


def read_rows(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.reader(f, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    return {r.id.split("|")[0]: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def attrs(s: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in s.split(";") if "=" in part)


def cds_rows(path: Path) -> list[dict[str, object]]:
    rows = []
    for line in path.open(encoding="utf-8"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        a = attrs(f[8])
        protein = a.get("protein_id") or a.get("Name", "")
        locus = a.get("locus_tag", "")
        if protein and locus:
            rows.append({"seq": f[0], "start": int(f[3]), "end": int(f[4]),
                         "strand": f[6], "protein": protein, "locus": locus,
                         "parent": a.get("Parent", "")})
    return rows


def chromosome_sequence(gff: Path, fasta_candidates: list[Path]) -> tuple[str, str]:
    seqids = set()
    for line in gff.open(encoding="utf-8"):
        if not line.startswith("#"):
            seqids.add(line.split("\t", 1)[0])
    seqids.discard("")
    if len(seqids) != 1:
        raise ValueError(f"Expected one chromosome in {gff}; saw {seqids}")
    fasta = next((p for p in fasta_candidates if p.exists()), None)
    if fasta is None:
        raise FileNotFoundError(f"No chromosome FASTA found among {fasta_candidates}")
    records = list(SeqIO.parse(fasta, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {fasta}; got {len(records)}")
    return next(iter(seqids)), str(records[0].seq).upper()


def alignment(query: str, subject: str, aligner: PairwiseAligner):
    aln = aligner.align(query, subject)[0]
    aligned = identical = 0
    for (q1, q2), (s1, s2) in zip(*aln.aligned):
        n = min(q2 - q1, s2 - s1)
        aligned += n
        identical += sum(query[q1 + i] == subject[s1 + i] for i in range(n))
    if not aligned:
        return aln.score, 0, 0.0, 0.0, 0.0
    return aln.score, aligned, identical / aligned, aligned / len(query), aligned / len(subject)


def main() -> None:
    # Verify coordinate transfer is valid despite GenBank/RefSeq contig IDs differing.
    old_ref_seqid, old_ref_seq = chromosome_sequence(
        OLD_REF_DIR / "genomic.gff",
        [OLD_GCF_FASTA])
    old_gca_seqid, old_gca_seq = chromosome_sequence(
        OLD_GCA_DIR / "genomic.gff",
        [OLD_GCA_DIR / "GCA_000021485.1_ASM2148v1_genomic.fna"])
    if len(old_ref_seq) != 2_982_397 or len(old_gca_seq) != 2_982_397 or old_ref_seq != old_gca_seq:
        raise ValueError("Old GCA/GCF chromosome sequences differ; coordinate projection stopped")

    cur_proteins = read_fasta(CUR_DIR / "protein.faa")
    raw_proteins = read_fasta(OLD_GCA_DIR / "protein.faa")
    ref_cds = cds_rows(OLD_REF_DIR / "genomic.gff")
    raw_cds = [x for x in cds_rows(OLD_GCA_DIR / "genomic.gff")
               if str(x["locus"]).startswith("AFE_")]

    ref_by_wp: dict[str, list[dict[str, object]]] = defaultdict(list)
    for feature in ref_cds:
        ref_by_wp[str(feature["protein"])].append(feature)
    ref_by_locus = {str(feature["locus"]): feature for feature in ref_cds}

    # Transfer RefSeq old-coordinate CDSs to GCA raw AFE CDSs. Sequence equality
    # above validates the coordinate system; seqid is deliberately normalized.
    raw_by_strand: dict[str, list[dict[str, object]]] = defaultdict(list)
    for raw in raw_cds:
        raw_by_strand[str(raw["strand"])].append(raw)
    for group in raw_by_strand.values():
        group.sort(key=lambda x: int(x["start"]))
    raw_projection: dict[str, list[tuple[dict[str, object], int]]] = defaultdict(list)
    reverse_projection: dict[str, set[str]] = defaultdict(set)
    for ref in ref_cds:
        for raw in raw_by_strand[str(ref["strand"])]:
            if int(raw["start"]) > int(ref["end"]):
                break
            if int(raw["end"]) < int(ref["start"]):
                continue
            bp = min(int(raw["end"]), int(ref["end"])) - max(int(raw["start"]), int(ref["start"])) + 1
            if bp > 0:
                raw_projection[str(ref["locus"])].append((raw, bp))
                reverse_projection[str(raw["locus"])].add(str(ref["locus"]))

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    raw_items = list(raw_proteins.items())
    raw_seed_index: dict[str, list[int]] = defaultdict(list)
    for i, (_, seq) in enumerate(raw_items):
        for seed in {seq[j:j + 5] for j in range(max(0, len(seq) - 4))}:
            raw_seed_index[seed].append(i)
    raw_id_indices: dict[str, list[int]] = defaultdict(list)
    for i, (protein_id, _) in enumerate(raw_items):
        raw_id_indices[protein_id].append(i)

    d3_dir = next(p for p in ROOT.glob("D3_*") if p.is_dir())
    d3_table = next(d3_dir.rglob("D3-1-R_全基因跨数据库交叉引用.tsv"))
    d3_rows = read_rows(d3_table)
    d3_header = d3_rows[0]
    audit_rows = [
        [row[d3_header.index("当前locus")], row[d3_header.index("当前protein ID")]]
        for row in d3_rows[1:]
        if row[d3_header.index("基因类型")] == "protein_coding"
        and row[d3_header.index("B1映射状态")] != "完全匹配"
    ]
    assert len(audit_rows) == 324
    aln_rows = read_rows(next(G5_DIR.glob("G5_\u8f85\u52a9\u5e8f\u5217\u6bd4\u5bf9.tsv")))[1:]
    neigh_rows = read_rows(next(G5_DIR.glob("G5_\u8f85\u52a9\u90bb\u57df\u6838\u9a8c.tsv")))[1:]
    aln_by_locus = {r[0]: r for r in aln_rows}
    neigh_by_locus = {r[0]: r for r in neigh_rows}
    output = []

    for audit in audit_rows:
        locus, current_wp = audit[0], audit[1]
        query = cur_proteins.get(current_wp, "")
        if not query:
            continue
        old_wp = aln_by_locus.get(locus, ["", "", ""])[2]
        ref_candidates = list(ref_by_wp.get(old_wp, ()))
        inside = set(re.findall(r"AFE_RS\d+", neigh_by_locus.get(locus, [""] * 7)[6]))
        selected_ref = ([x for x in ref_candidates if str(x["locus"]) in inside]
                        if inside else ref_candidates)

        # Search original AFE CDSs in the interval bounded by current locus's
        # nearest exact B1 anchors after mapping those anchors to old GCF coords.
        neigh = neigh_by_locus.get(locus, [""] * 4)
        left_ids = re.findall(r"AFE_RS\d+", neigh[2]) if len(neigh) > 2 else []
        right_ids = re.findall(r"AFE_RS\d+", neigh[3]) if len(neigh) > 3 else []
        left_anchor = ref_by_locus.get(left_ids[-1]) if left_ids else None
        right_anchor = ref_by_locus.get(right_ids[-1]) if right_ids else None
        interval_cds = []
        interval_label = "no_two_sided_anchor_interval"
        if left_anchor and right_anchor:
            left_edge, right_edge = int(left_anchor["end"]), int(right_anchor["start"])
            interval_label = (f"{left_anchor['locus']}:{left_anchor['start']}..{left_anchor['end']}"
                              f"->{right_anchor['locus']}:{right_anchor['start']}..{right_anchor['end']}")
            for raw in raw_cds:
                midpoint = (int(raw["start"]) + int(raw["end"])) / 2
                in_interval = (left_edge < midpoint < right_edge if left_edge < right_edge
                               else midpoint > left_edge or midpoint < right_edge)
                if in_interval:
                    interval_cds.append(raw)
        interval_hits = []
        for raw in interval_cds:
            ack = str(raw["protein"])
            subject = raw_proteins.get(ack, "")
            if not subject:
                continue
            score, aa, ident, current_cov, old_cov = alignment(query, subject, aligner)
            interval_hits.append((score, aa, ident, current_cov, old_cov, raw))
        interval_hits.sort(key=lambda x: x[0], reverse=True)
        interval_best = interval_hits[0] if interval_hits else None
        interval_hit_text = ";".join(
            f"{raw['locus']}({raw['protein']})[aa={aa};id={ident:.4f};curCov={qcov:.4f};oldCov={scov:.4f}]"
            for _, aa, ident, qcov, scov, raw in interval_hits if aa >= 15 and ident >= .30
        ) or "no_interval_hit_above_minimum"
        interval_best_text = (f"{interval_best[5]['locus']}({interval_best[5]['protein']})"
                              f"[id={interval_best[2]:.4f};curCov={interval_best[3]:.4f};oldCov={interval_best[4]:.4f}]"
                              if interval_best else "none")

        projected = []
        for ref in selected_ref:
            for raw, bp in raw_projection.get(str(ref["locus"]), ()):
                projected.append((ref, raw, bp))

        seed_counts: dict[int, int] = defaultdict(int)
        for j in range(max(0, len(query) - 4)):
            for i in raw_seed_index.get(query[j:j + 5], ()):
                seed_counts[i] += 1
        candidate_indices = set(sorted(seed_counts, key=seed_counts.get, reverse=True)[:40])
        for _, raw, _ in projected:
            candidate_indices.update(raw_id_indices.get(str(raw["protein"]), ()))
        ranked = []
        for i in candidate_indices:
            ack, subject = raw_items[i]
            ranked.append((*alignment(query, subject, aligner), ack))
        ranked.sort(reverse=True)
        best = ranked[0] if ranked else (0.0, 0, 0.0, 0.0, 0.0, "")
        top_ack = best[5]
        top_loci = sorted({str(raw["locus"]) for raw in raw_cds
                           if str(raw["protein"]) == top_ack})

        projection_text = []
        split_flags = set()
        substantial_raw_loci = set()
        for ref, raw, bp in projected:
            ack = str(raw["protein"])
            subject = raw_proteins.get(ack, "")
            if not subject:
                continue
            _, aa, ident, current_cov, old_cov = alignment(query, subject, aligner)
            projection_text.append(
                f"{ref['locus']}({ref['protein']})=>{raw['locus']}({ack})"
                f"[bp={bp};aa={aa};id={ident:.4f};curCov={current_cov:.4f};oldCov={old_cov:.4f}]"
            )
            mapped_ref_loci = reverse_projection[str(raw["locus"])]
            substantial_ref_loci = {str(ref_feature["locus"])
                                    for ref_feature, overlap_bp in raw_projection[str(ref["locus"])]
                                    if str(ref_feature["locus"]) == str(raw["locus"]) and overlap_bp >= 30}
            if len(mapped_ref_loci) > 1 and int(bp) >= 90:
                split_flags.add(
                    f"old annotation topology: {raw['locus']} overlaps {len(mapped_ref_loci)} old-RefSeq CDS (not evidence of current split)"
                )
            if int(bp) >= 90 and aa >= 30 and ident >= .30:
                substantial_raw_loci.add(str(raw["locus"]))
        if len(substantial_raw_loci) > 1:
            split_flags.add(
                f"possible current/old annotation merge: substantial sequence support at {len(substantial_raw_loci)} original AFE loci"
            )

        if projected and split_flags:
            conclusion = "原始AFE坐标对应；旧注释拓扑差异/可能merge，按覆盖与多拷贝复核"
        elif projected:
            conclusion = "旧RefSeq坐标投射命中原始AFE CDS；查看identity及双向coverage"
        elif best[2] >= .9 and best[3] >= .9 and best[4] >= .9:
            conclusion = "高相似原始AFE蛋白命中但未由旧RefSeq坐标投射"
        else:
            conclusion = "无原始AFE坐标投射命中；结合最佳ACK序列覆盖复核"

        output.append([locus, current_wp, old_wp or "无",
                       ";".join(str(x["locus"]) for x in selected_ref) or "无",
                       interval_label, interval_hit_text, interval_best_text,
                       ";".join(projection_text) or "无", ";".join(top_loci) or "无",
                       top_ack or "无", str(best[1]), f"{best[2]:.4f}",
                       f"{best[3]:.4f}", f"{best[4]:.4f}",
                       ";".join(sorted(split_flags)) or "无", conclusion])

    outpath = G5_DIR / "G5_\u8f85\u52a9\u539f\u59cbAFE\u6838\u9a8c.tsv"
    header = ["\u5f53\u524dlocus", "\u5f53\u524dprotein ID", "\u65e7RefSeq top WP",
              "\u90bb\u57df\u7b5b\u9009\u540e\u65e7RefSeq locus_tag",
              "B1\u53cc\u4fa7\u951a\u70b9\u65e7GCF\u5750\u6807\u533a\u95f4",
              "\u533a\u95f4\u539f\u59cbAFE CDS\u5019\u9009\u5e8f\u5217\u6bd4\u5bf9",
              "\u533a\u95f4\u6700\u4f73AFE\u5e8f\u5217\u547d\u4e2d",
              "\u65e7GCF\u5750\u6807\u6295\u5c04\u5230\u539f\u59cbGCA AFE CDS\u53ca\u5e8f\u5217identity/\u53cc\u5411coverage",
              "\u5168\u5e93\u6700\u4f73\u539f\u59cbAFE_xxxx", "\u6700\u4f73ACK protein ID",
              "\u6bd4\u5bf9\u6c28\u57fa\u9178\u6570", "identity", "current coverage", "old coverage",
              "annotation/topology candidates (not confirmed events)", "\u6838\u9a8c\u7ed3\u8bba"]
    with outpath.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(output)

    projected = [r for r in output if r[7] != "\u65e0"]
    split_rows = [r for r in output if r[14] != "\u65e0"]
    print(f"old chromosome IDs: RefSeq={old_ref_seqid}, GenBank={old_gca_seqid}; length={len(old_ref_seq)}; identical=True")
    print(f"rows={len(output)}; old-coordinate projection to original AFE CDS={len(projected)}")
    print(f"annotation/topology candidates (not confirmed events)={len(split_rows)}")
    print(f"output={outpath}")


if __name__ == "__main__":
    main()
