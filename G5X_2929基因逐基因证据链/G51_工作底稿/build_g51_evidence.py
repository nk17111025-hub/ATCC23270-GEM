"""Build a traceable first-pass G51 evidence inventory from frozen source files.

The inventory records observations; it does not infer biochemical reactions from
product names or transfer evidence from homologs to this strain.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
GFF = ROOT / "B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff"
MAP = ROOT / "B1_最新基因组与旧AFE编号映射/标准化数据/旧新基因映射.tsv"
MODEL = ROOT / "A0_复现2016与2024模型/A0_主任务/脚本与环境/中间/mmc1.xlsx"
KEGG = ROOT / "01_原始数据/04_数据库原始导出/KEGG"
D3 = ROOT / "D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv"
RHEA_DIR = ROOT / "01_原始数据/04_数据库原始导出/Rhea"


def attrs(text):
    return dict(part.split("=", 1) for part in text.split(";") if "=" in part)


def read_gff():
    genes = []
    cds = {}
    for line in GFF.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        a = attrs(cols[8])
        locus = a.get("locus_tag")
        if not locus:
            continue
        if cols[2] == "CDS":
            cds[locus] = {"wp": a.get("protein_id", ""), "cds_product": urllib.parse.unquote(a.get("product", "")), "cds_start": int(cols[3]), "cds_end": int(cols[4]), "cds_strand": cols[6], "cds_partial": a.get("partial", "")}
        elif cols[2] == "gene" and a.get("gene_biotype") == "protein_coding":
            genes.append({"locus": locus, "chromosome": cols[0], "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6], "gene_attrs": a})
    genes.sort(key=lambda g: (g["chromosome"], g["start"]))
    for i, g in enumerate(genes, 1):
        g["ordinal"] = i
        g.update(cds.get(g["locus"], {}))
    return genes[:293]


def load_mapping():
    out = defaultdict(list)
    with MAP.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["当前GCF049位点"]:
                out[row["当前GCF049位点"]].append(row)
    return out


def load_old_reactions():
    df = pd.read_excel(MODEL, sheet_name="Table 1", header=1).fillna("")
    out = defaultdict(list)
    for _, row in df.iterrows():
        assoc = str(row["Gene-Reaction Association"])
        for afe in set(re.findall(r"AFE_\d{4,5}", assoc)):
            out[afe].append({str(k): str(v) for k, v in row.items()})
    by_ec = defaultdict(list)
    for _, row in df.iterrows():
        ec = str(row.get("EC Number", "")).strip()
        if ec:
            by_ec[ec].append({str(k): str(v) for k, v in row.items()})
    return out, by_ec


def load_kegg():
    relations = defaultdict(lambda: defaultdict(list))
    for filename, kind in [("KEGG_link_ko_afr.txt", "ko"), ("KEGG_link_enzyme_afr.txt", "ec"), ("KEGG_link_reaction_afr.txt", "reaction")]:
        with (KEGG / filename).open(encoding="utf-8") as f:
            for line in f:
                left, right = line.strip().split("\t")[:2]
                relations[left.replace("afr:", "")][kind].append(right)
    return relations


def load_rhea():
    out = {}
    for path in RHEA_DIR.glob("*.tsv"):
        with path.open(encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                accession = row.get("Reaction identifier")
                if accession and accession.startswith("RHEA:") and row.get("Equation"):
                    out[accession] = {**row, "source_file": str(path.relative_to(ROOT))}
    return out


def uniprot_results(wps):
    cache = OUT / "uniprot_atcc23270_cache.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    out = []
    for offset in range(0, len(wps), 15):
        batch = wps[offset:offset+15]
        query = "(" + " OR ".join("xref:RefSeq-" + wp for wp in batch) + ") AND organism_id:243159"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": query, "format": "json", "size": 100})
        for retry in range(4):
            try:
                with urllib.request.urlopen(url, timeout=45) as r:
                    data = json.load(r)
                out.extend(data.get("results", []))
                break
            except Exception as e:
                if retry == 3:
                    print("UniProt batch failed", offset, str(e), flush=True)
                time.sleep(2**retry)
        print("UniProt", min(offset+15, len(wps)), "/", len(wps), flush=True)
        time.sleep(0.2)
    cache.write_text(json.dumps(out, ensure_ascii=False), encoding="utf-8")
    return out


def main():
    genes = read_gff()
    assert len(genes) == 293 and genes[0]["locus"] == "RU820_RS00005" and genes[-1]["locus"] == "RU820_RS01505"
    mapping = load_mapping()
    old_reactions, old_by_ec = load_old_reactions()
    kegg = load_kegg()
    rhea = load_rhea()
    with D3.open(encoding="utf-8-sig", newline="") as f:
        d3 = {r["当前locus"]: r for r in csv.DictReader(f, delimiter="\t")}
    uni = defaultdict(list)
    for entry in uniprot_results([g["wp"] for g in genes if g.get("wp")]):
        refs = entry.get("uniProtKBCrossReferences", [])
        for wp in {r.get("id") for r in refs if r.get("database") == "RefSeq"}:
            uni[wp].append(entry)
    results = []
    for g in genes:
        maps = mapping.get(g["locus"], [])
        exact_maps = [m for m in maps if m["当前蛋白ID"] == g.get("wp") and m["映射方法"] == "RefSeq protein_id exact" and m["映射状态"] == "完全匹配"]
        afe = sorted({m["旧AFE位点"] for m in exact_maps})
        model_rows = [r for a in afe for r in old_reactions.get(a, [])]
        kegg_obs = {a: dict(kegg.get(a, {})) for a in afe}
        uni_obs = []
        for u in uni.get(g.get("wp"), []):
            refs = u.get("uniProtKBCrossReferences", [])
            uni_obs.append({
                "accession": u.get("primaryAccession"),
                "entry_type": u.get("entryType"),
                "annotation_score": u.get("annotationScore"),
                "protein_existence": u.get("proteinExistence"),
                "name": u.get("proteinDescription", {}).get("recommendedName") or u.get("proteinDescription", {}).get("submissionNames"),
                "ec": sorted({e.get("value", "") for item in u.get("proteinDescription", {}).values() if isinstance(item, dict) for e in item.get("ecNumbers", []) if isinstance(e, dict)}),
                "crossrefs": [{"db": r["database"], "id": r["id"]} for r in refs if r["database"] in {"InterPro", "Pfam", "Rhea", "KEGG", "TCDB", "RefSeq"}],
                "url": "https://www.uniprot.org/uniprotkb/" + u.get("primaryAccession", ""),
                "entry_version": u.get("entryAudit", {}).get("entryVersion"),
            })
        status = "逐基因初筛完成；原始论文和精确反应待核"
        if not afe:
            status += "；旧 AFE 未确认"
        rhea_ids = re.findall(r"RHEA:\d+", d3.get(g["locus"], {}).get("Rhea精确交叉引用", ""))
        rhea_defs = [rhea[r] for r in rhea_ids if r in rhea]
        ec_candidates = {r["EC number"].replace("EC:", "") for r in rhea_defs if r.get("EC number")}
        ec_candidates.update(re.findall(r"\d+\.\d+\.\d+\.\d+", d3.get(g["locus"], {}).get("KEGG EC", "")))
        same_ec_old = {r["Reaction ID"]: r for ec in ec_candidates for r in old_by_ec.get(ec, [])}
        results.append({
            "ordinal": g["ordinal"], "locus": g["locus"], "wp": g.get("wp", ""),
            "chromosome": g["chromosome"], "start": g["start"], "end": g["end"], "strand": g["strand"],
            "product": g.get("cds_product", ""), "cds_match_gene": g.get("cds_start") == g["start"] and g.get("cds_end") == g["end"] and g.get("cds_strand") == g["strand"],
            "cds_partial": g.get("cds_partial", ""), "old_afe": afe, "map_rows": maps,
            "uniprot": uni_obs, "kegg": kegg_obs, "d3_crossref": d3.get(g["locus"], {}), "rhea_definitions": rhea_defs, "model_2016": model_rows, "model_2016_same_ec": list(same_ec_old.values()),
            "review_status": status, "exact_reaction_verified": False, "action": "证据不足暂缓", "reaction_score": "不评分",
        })
    (OUT / "G51_evidence_inventory.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print("genes", len(results), "unique_loci", len({g['locus'] for g in results}), "exact_afe", sum(bool(g['old_afe']) for g in results), "old_model_genes", sum(bool(g['model_2016']) for g in results), "uniprot_genes", sum(bool(g['uniprot']) for g in results), flush=True)


if __name__ == "__main__":
    main()
