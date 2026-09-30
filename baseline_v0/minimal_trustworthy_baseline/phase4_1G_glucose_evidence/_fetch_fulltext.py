# -*- coding: utf-8 -*-
"""Fetch Wang 2012 full text (PMC) + two more key abstracts."""
import re
import json
import time
import urllib.request


def get(url):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "codex-evidence-retrieval")
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def main():
    # Wang 2012 full text via NCBI eutils efetch (PMC)
    txt = ""
    for base in [
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3298148&rettype=full&retmode=xml",
    ]:
        try:
            txt = get(base)
            if len(txt) > 5000:
                break
        except Exception as e:
            print("efetch ERR", e)
        time.sleep(1)
    if txt:
        open("_wang2012_fulltext.xml", "w", encoding="utf-8").write(txt)
        print("Wang full text len", len(txt))
        plain = re.sub(r"<[^>]+>", " ", txt)
        for term in ["glucokinase", "glucose kinase", "AFE_2841", "AFE_2025", "zwf", "AFE_3251",
                     "glyceraldehyde-3-phosphate dehydrogenase", "gap", "phosphoglycerate kinase", "pgk",
                     "AFE_1807", "pfkB", "glucose-6-phosphate", "AFE_2522", "AFE_2250", "AFE_2312",
                     "AFE_1971", "porin", "glucose transporter", "ROK"]:
            idx = [m.start() for m in re.finditer(re.escape(term), plain, re.IGNORECASE)]
            if idx:
                # print first context
                s = idx[0]
                ctx = plain[max(0, s - 120): s + 220]
                print(f"\n### [{term}] {len(idx)} hits")
                print(ctx.strip()[:340])

    # fetch abstracts for two more papers
    for doi in ["10.1186/1471-2180-10-229", "10.3389/fmicb.2016.01365", "10.1128/aem.03057-12"]:
        try:
            url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" + \
                  f'DOI:"{doi}"&format=json&pageSize=1&resultType=core&fields=title,year,doi,pmid,abstractText,authorString,journalTitle'
            d = json.loads(get(url))
            r = d.get("resultList", {}).get("result", [{}])[0]
            print("\n" + "=" * 100)
            print("TITLE:", r.get("title"))
            print("DOI:", r.get("doi"), "PMID:", r.get("pmid"), "AUTHORS:", r.get("authorString"))
            print("ABS:", (r.get("abstractText") or "")[:2200])
        except Exception as e:
            print("ERR", doi, e)
        time.sleep(0.4)


if __name__ == "__main__":
    main()
