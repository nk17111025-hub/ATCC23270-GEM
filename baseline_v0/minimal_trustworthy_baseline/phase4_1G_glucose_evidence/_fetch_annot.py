# -*- coding: utf-8 -*-
"""Fetch UniProt annotation for candidate glucose genes (ATCC 23270 / taxid 243159)."""
import json
import time
import urllib.parse
import urllib.request

GENES = [
    ("AFE_2522", "RU820_RS11605", "WP_012537238.1", "carbohydrate porin"),
    ("AFE_2250", "RU820_RS10385", "WP_012537080.1", "carbohydrate porin"),
    ("AFE_2312", "RU820_RS10670", "WP_012537116.1", "sugar porter MFS"),
    ("AFE_1971", "RU820_RS09115", "WP_012536876.1", "MFS transporter"),
    ("AFE_2841", "RU820_RS13140", "WP_012537424.1", "ROK family protein"),
    ("AFE_2025", "RU820_RS09360", "WP_012536916.1", "G6PDH zwf"),
    ("AFE_3251", "RU820_RS15125", "WP_009567621.1", "gap"),
    ("AFE_3250", "RU820_RS15120", "WP_012537677.1", "pgk"),
]


def uniprot_search(query, fields):
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "tsv", "fields": fields, "size": "5"}
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=40) as r:
        return r.read().decode("utf-8")


def main():
    fields = "accession,reviewed,protein_name,gene_names,organism_name,length,ec,cc_catalytic_activity,cc_function,cc_caution,ft_domain,protein_families,go_p,xref_refseq"
    out = []
    for afe, ru, wp, annot in GENES:
        row = {"AFE": afe, "RU820": ru, "WP": wp, "given_annotation": annot}
        # search by gene_exact (AFE name)
        for q in [f"gene_exact:{afe} AND taxonomy_id:243159",
                  f"xref:refseq:{wp}",
                  f"{afe}"]:
            try:
                txt = uniprot_search(q, fields)
                lines = txt.strip().split("\n")
                if len(lines) > 1:
                    row["uniprot_query"] = q
                    row["uniprot_hits"] = "\n".join(lines[1:6])
                    break
            except Exception as e:
                row.setdefault("errs", []).append(str(e))
            time.sleep(0.4)
        out.append(row)
        print("=== ", afe, "|", wp, "|", annot)
        print(row.get("uniprot_hits", "NO HIT"))
        print()
        time.sleep(0.4)
    with open("_annot.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
