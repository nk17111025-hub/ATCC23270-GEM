# -*- coding: utf-8 -*-
"""Fetch full metadata + abstracts for key papers by DOI."""
import json
import time
import urllib.parse
import urllib.request

DOIS = [
    "10.1186/1471-2164-9-597",        # Valdes 2008 genome
    "10.1128/aem.07230-11",           # pfkB mutant
    "10.3389/fmicb.2019.00603",       # CO2 concentration / carbon assimilation
    "10.3390/microorganisms12030590", # electroautotrophic/chemoautotrophic
    "10.3389/fmicb.2019.00592",       # proteomics oxidative stress
    "10.1128/aem.01424-20",           # mixotrophic Thiomonas
    "10.3389/fmicb.2023.1277847",     # novel Acidithiobacillus genome
    "10.1186/1471-2164-11-404",       # biofilm/planktonic transcriptome
    "10.3389/fmicb.2024.1374800",     # acidophilic heterotrophs
    "10.3390/microorganisms8071076",  # arsenic stress
]


def fetch_doi(doi):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": f'DOI:"{doi}"', "format": "json", "pageSize": "1", "resultType": "core",
              "fields": "title,authorString,year,doi,pmid,pubType,journalTitle,abstractText,affiliation"}
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=40) as r:
        d = json.load(r)
    res = d.get("resultList", {}).get("result", [])
    return res[0] if res else None


def main():
    out = []
    for doi in DOIS:
        try:
            r = fetch_doi(doi)
            if r:
                out.append(r)
                print("=" * 100)
                print("TITLE:", r.get("title"))
                print("YEAR:", r.get("year"), "| DOI:", r.get("doi"), "| PMID:", r.get("pmid"), "| J:", r.get("journalTitle"))
                print("AUTHORS:", r.get("authorString"))
                print("ABSTRACT:", (r.get("abstractText") or "")[:2500])
        except Exception as e:
            print("ERR", doi, e)
        time.sleep(0.4)
    with open("_papers.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\nDONE", len(out))


if __name__ == "__main__":
    main()
