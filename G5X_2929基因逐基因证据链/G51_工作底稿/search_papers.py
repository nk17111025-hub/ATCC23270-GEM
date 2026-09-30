"""Search Europe PMC per G51 legacy locus; retain raw candidate hits for review."""
import concurrent.futures
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROWS = json.loads((HERE / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
DEST = HERE / "europe_pmc_gene_search.json"


def search(row):
    terms = row["old_afe"] or [row["locus"], row["wp"]]
    hits = []
    for term in terms:
        query = '"' + term + '"'
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": query, "format": "json", "resultType": "core", "pageSize": 25})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=35) as response:
                    data = json.load(response)
                hits.append({"term": term, "url": url, "hit_count": data.get("hitCount", 0), "results": [{k: r.get(k) for k in ("id", "source", "title", "authorString", "pubYear", "doi", "pmid", "pmcid", "journalTitle", "pubTypeList", "abstractText")} for r in data.get("resultList", {}).get("result", [])]})
                break
            except Exception as exc:
                if attempt == 3:
                    hits.append({"term": term, "url": url, "error": str(exc), "results": []})
                time.sleep(2 ** attempt)
    return row["locus"], hits


if __name__ == "__main__":
    if DEST.exists():
        old = json.loads(DEST.read_text(encoding="utf-8"))
    else:
        old = {}
    todo = [r for r in ROWS if r["locus"] not in old]
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        futures = [pool.submit(search, r) for r in todo]
        for n, future in enumerate(concurrent.futures.as_completed(futures), 1):
            locus, hits = future.result()
            old[locus] = hits
            if n % 20 == 0 or n == len(todo):
                DEST.write_text(json.dumps(old, ensure_ascii=False), encoding="utf-8")
                print("searched", n, "/", len(todo), flush=True)
    DEST.write_text(json.dumps(old, ensure_ascii=False, indent=2), encoding="utf-8")
    print("total_loci", len(old), "hits", sum(sum(h.get("hit_count", 0) for h in v) for v in old.values()), flush=True)
