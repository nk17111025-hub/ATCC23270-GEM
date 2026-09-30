# -*- coding: utf-8 -*-
"""Search Europe PMC for ATCC 23270 glucose-related evidence."""
import json
import time
import urllib.parse
import urllib.request

QUERIES = [
    "Acidithiobacillus ferrooxidans glucose",
    "Acidithiobacillus ferrooxidans glucose transport",
    "Acidithiobacillus ferrooxidans glucose uptake",
    "Acidithiobacillus ferrooxidans glucokinase",
    "Acidithiobacillus ferrooxidans glucose-6-phosphate dehydrogenase",
    "Acidithiobacillus ferrooxidans pentose phosphate pathway",
    "Acidithiobacillus ferrooxidans glycolysis",
    "Acidithiobacillus ferrooxidans gluconeogenesis",
    "Acidithiobacillus ferrooxidans Calvin cycle",
    "Acidithiobacillus ferrooxidans carbohydrate porin",
    "Acidithiobacillus ferrooxidans OprB",
    "Acidithiobacillus ferrooxidans sugar transporter",
    "Acidithiobacillus ferrooxidans Embden-Meyerhof",
    "Acidithiobacillus ferrooxidans EMP pathway",
    "Acidithiobacillus ferrooxidans phosphoglycerate kinase",
    "Acidithiobacillus ferrooxidans glyceraldehyde-3-phosphate dehydrogenase",
    "Acidithiobacillus ferrooxidans heterotrophic growth",
    "Acidithiobacillus ferrooxidans organic carbon",
    "AFE_2841",
    "AFE_2025",
    "AFE_3251",
    "AFE_2522",
]


def search(q, n=25):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": q, "format": "json", "pageSize": str(n),
              "resultType": "core", "fields": "title,authorString,year,doi,pmid,pubType,journalTitle,abstractText"}
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def main():
    allrows = []
    for q in QUERIES:
        try:
            d = search(q)
            n = d.get("hitCount", 0)
            print(f"=== {q} -> {n} hits")
            for r in d.get("resultList", {}).get("result", [])[:8]:
                t = r.get("title", "")[:110]
                ab = (r.get("abstractText") or "")[:160]
                allrows.append({
                    "query": q, "title": r.get("title", ""),
                    "year": r.get("year", ""), "doi": r.get("doi", ""),
                    "pmid": r.get("pmid", ""), "journal": r.get("journalTitle", ""),
                    "authors": r.get("authorString", ""), "abstract": ab,
                })
                print(f"   [{r.get('year','')}] {t}")
        except Exception as e:
            print("ERR", q, e)
        time.sleep(0.3)
    with open("_lit_hits.json", "w", encoding="utf-8") as f:
        json.dump(allrows, f, ensure_ascii=False, indent=2)
    print("TOTAL rows", len(allrows))


if __name__ == "__main__":
    main()
