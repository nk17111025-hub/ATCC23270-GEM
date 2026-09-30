"""Check both current RU820 and WP identifiers in Europe PMC for all 293 G51 loci."""
import concurrent.futures
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

here = Path(__file__).resolve().parent
rows = json.loads((here / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
dest = here / "europe_pmc_current_id_search.json"
cache = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {}

def search(term):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({
        "query": '"' + term + '"', "format": "json", "resultType": "core", "pageSize": 25})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=35) as response:
                data = json.load(response)
            return term, {"url": url, "hit_count": data.get("hitCount", 0), "results": [
                {k: result.get(k) for k in ("id", "source", "title", "authorString", "pubYear", "doi", "pmid", "pmcid", "journalTitle", "abstractText")}
                for result in data.get("resultList", {}).get("result", [])]}
        except Exception as exc:
            if attempt == 3:
                return term, {"url": url, "error": str(exc), "results": []}
            time.sleep(2 ** attempt)

terms = [term for row in rows for term in (row["locus"], row["wp"]) if term not in cache]
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
    for n, pair in enumerate(pool.map(search, terms), 1):
        term, item = pair
        cache[term] = item
        if n % 40 == 0 or n == len(terms):
            dest.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
            print("searched", n, "/", len(terms), flush=True)
dest.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")
print("total", len(cache), "hits", sum(v.get("hit_count", 0) for v in cache.values()), "errors", sum("error" in v for v in cache.values()))
