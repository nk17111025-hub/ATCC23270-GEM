# -*- coding: utf-8 -*-
import json
import time
import urllib.parse
import urllib.request


def epmc(query):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = f'query={urllib.parse.quote(query)}&format=json&pageSize=1&resultType=core&fields=title,year,doi,pmid,abstractText,authorString,journalTitle'
    req = urllib.request.Request(url + "?" + params)
    req.add_header("User-Agent", "codex")
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def main():
    for doi in ["10.1186/1471-2180-10-229", "10.3389/fmicb.2016.01365", "10.1128/aem.03057-12"]:
        try:
            d = epmc('DOI:"' + doi + '"')
            r = d.get("resultList", {}).get("result", [{}])[0]
            print("=" * 100)
            print("TITLE:", r.get("title"))
            print("DOI:", r.get("doi"), "PMID:", r.get("pmid"), "AUTH:", r.get("authorString"))
            print("ABS:", (r.get("abstractText") or "")[:2200])
        except Exception as e:
            print("ERR", doi, e)
        time.sleep(0.4)


if __name__ == "__main__":
    main()
