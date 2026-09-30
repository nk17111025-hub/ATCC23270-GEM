# -*- coding: utf-8 -*-
"""Deeper targeted searches + full-text of Wang 2012."""
import json
import re
import time
import urllib.parse
import urllib.request


def epmc_search(q, n=15):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {"query": q, "format": "json", "pageSize": str(n), "resultType": "core",
              "fields": "title,authorString,year,doi,pmid,pmcid,journalTitle,abstractText"}
    req = urllib.request.Request(url + "?" + urllib.parse.urlencode(params))
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def main():
    queries = [
        'Acidithiobacillus ferrooxidans "glucose transporter"',
        'Acidithiobacillus ferrooxidans "glucose" AND (porin OR MFS OR permease OR transporter)',
        'Acidithiobacillus ferrooxidans iMC507 OR "metabolic model"',
        'Acidithiobacillus ferrooxidans "quorum sensing" transcriptomic',
        'Acidithiobacillus ferrooxidans "central carbon metabolism"',
        'Acidithiobacillus ferrooxidans "glucose-6-phosphate"',
        'Acidithiobacillus ferrooxidans "phosphofructokinase" OR pfkB',
        'Acidithiobacillus ferrooxidans heterotroph glucose grow',
        'Acidithiobacillus ferrooxidans "EMP" OR "glycolytic" OR "glycolysis" glucose',
    ]
    out = []
    for q in queries:
        try:
            d = epmc_search(q)
            print(f"\n===== {q} -> {d.get('hitCount')}")
            for r in d.get("resultList", {}).get("result", [])[:8]:
                out.append({"query": q, "title": r.get("title"), "year": r.get("year"),
                            "doi": r.get("doi"), "pmid": r.get("pmid"), "pmcid": r.get("pmcid"),
                            "journal": r.get("journalTitle"), "abstract": r.get("abstractText") or ""})
                print(f"  [{r.get('year')}] {r.get('title','')[:100]} DOI={r.get('doi')} PMC={r.get('pmcid')}")
        except Exception as e:
            print("ERR", q, e)
        time.sleep(0.3)

    # full text of Wang 2012
    try:
        d = epmc_search('DOI:"10.1128/aem.07230-11"')
        r = d.get("resultList", {}).get("result", [{}])[0]
        pmcid = r.get("pmcid", "")
        print("\nWang 2012 PMCID:", pmcid)
        if pmcid:
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "codex-evidence-retrieval")
            with urllib.request.urlopen(req, timeout=60) as rr:
                txt = rr.read().decode("utf-8", "replace")
            with open("_wang2012_fulltext.xml", "w", encoding="utf-8") as f:
                f.write(txt)
            print("full text len", len(txt))
            for term in ["glucokinase", "glucose kinase", "AFE_2841", "AFE_2025", "zwf", "AFE_3251", "gap",
                         "glyceraldehyde", "phosphoglycerate kinase", "AFE_1807", "pfkB", "glucose-6-phosphate",
                         "AFE_2522", "AFE_2250", "AFE_2312", "AFE_1971", "porin", "transporter"]:
                idx = [m.start() for m in re.finditer(re.escape(term), txt, re.IGNORECASE)]
                if idx:
                    print(f"  [{term}] {len(idx)} hits")
    except Exception as e:
        print("full text ERR", e)

    with open("_lit_hits3.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
