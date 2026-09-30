# -*- coding: utf-8 -*-
"""Targeted follow-up: full metadata + abstracts for key papers, plus extra searches."""
import json
import time
import urllib.parse
import urllib.request

QUERIES = [
    'TITLE:"Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications"',
    'TITLE:"Development of a markerless gene replacement system for Acidithiobacillus ferrooxidans"',
    'TITLE:"Effect of CO2 Concentration on Uptake and Assimilation of Inorganic Carbon" Acidithiobacillus',
    'TITLE:"Characterize the Growth and Metabolism of Acidithiobacillus ferrooxidans"',
    'Acidithiobacillus ferrooxidans glucose metabolism',
    'Acidithiobacillus ferrooxidans glucose oxidation',
    'Acidithiobacillus ferrooxidans fructose',
    'Acidithiobacillus ferrooxidans phosphotransferase PTS',
    'Acidithiobacillus ferrooxidans glucose supplementation',
    'Acidithiobacillus ferrooxidans mixotrophic',
    'Acidithiobacillus ferrooxidans RuBisCO carbon dioxide fixation',
    'Acidithiobacillus ferrooxidans Embden Meyerhof Parnas',
    'Acidithiobacillus ferrooxidans glucokinase OR "glucose kinase" OR "ROK"',
    'Acidithiobacillus ferrooxidans glucose-6-phosphate isomerase',
    'Acidithiobacillus ferrooxidans fructose-bisphosphate',
]


def search(q, n=20):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": q, "format": "json", "pageSize": str(n), "resultType": "core",
              "fields": "title,authorString,year,doi,pmid,pubType,journalTitle,abstractText"}
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def main():
    out = []
    for q in QUERIES:
        try:
            d = search(q)
            n = d.get("hitCount", 0)
            print(f"\n===== {q} -> {n}")
            for r in d.get("resultList", {}).get("result", [])[:6]:
                rec = {
                    "query": q, "title": r.get("title", ""), "year": r.get("year", ""),
                    "doi": r.get("doi", ""), "pmid": r.get("pmid", ""),
                    "journal": r.get("journalTitle", ""), "authors": r.get("authorString", ""),
                    "abstract": (r.get("abstractText") or ""),
                }
                out.append(rec)
                print(f"  [{r.get('year','')}] {r.get('title','')[:120]} DOI={r.get('doi','')}")
        except Exception as e:
            print("ERR", q, e)
        time.sleep(0.3)
    with open("_lit_hits2.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\nTOTAL", len(out))


if __name__ == "__main__":
    main()
