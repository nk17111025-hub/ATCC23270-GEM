"""Fetch only Rhea ECs absent from the prior D3-1 download, preserving raw responses."""
import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = HERE / "Rhea新增原始响应"
RAW.mkdir(exist_ok=True)


def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


current = tsv(HERE / "D3-1-R_全基因跨数据库交叉引用.tsv")
previous = tsv(ROOT / "06_Rhea/D3-1/Rhea_modelBioCyc_EC_crossrefs.tsv")
wanted = {ec for row in current for field in ("KEGG EC", "BioCyc EC") for ec in row[field].split(";") if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec)}
covered = {ec.removeprefix("EC:") for row in previous for ec in row["EC number"].split(";") if ec}
missing = sorted(wanted - covered)
columns = "rhea-id,equation,ec,chebi-id,reaction-xref(KEGG),reaction-xref(MetaCyc)"
manifest = []
response_rows = []
for batch_number, start in enumerate(range(0, len(missing), 12), 1):
    batch = missing[start:start + 12]
    name = f"Rhea_new_EC_batch_{batch_number:03d}.tsv"
    path = RAW / name
    url = "https://www.rhea-db.org/rhea/?" + urllib.parse.urlencode({"query": " OR ".join("ec:" + ec for ec in batch), "columns": columns, "format": "tsv", "limit": 1000})
    status = "cached"
    if not path.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "D3-1-R ATCC23270 full-genome audit"}), timeout=45) as response:
                    data = response.read()
                    status = str(response.status)
                if not data.startswith(b"Reaction identifier\t"):
                    raise ValueError("unexpected Rhea response header")
                path.write_bytes(data)
                break
            except Exception as exc:
                status = f"ERROR {type(exc).__name__}: {exc}"
                if attempt == 2:
                    break
                time.sleep(2 * (attempt + 1))
    if path.exists():
        response_rows.extend(tsv(path))
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        size = path.stat().st_size
    else:
        digest, size = "", 0
    manifest.append({"文件": path.relative_to(ROOT).as_posix(), "查询EC": ";".join(batch), "URL/API": url, "获取时间": datetime.now().astimezone().isoformat(), "数据库版本": "Rhea REST；响应未给release字段", "HTTP状态": status, "字节数": size, "SHA256": digest})
    if batch_number % 10 == 0:
        print(f"Rhea gap {batch_number}/{(len(missing)+11)//12}", flush=True)

with (HERE / "Rhea_新增查询清单.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(manifest)
by_id = {row["Reaction identifier"]: row for row in [*previous, *response_rows] if row.get("Reaction identifier")}
with (HERE / "Rhea_全基因EC查询汇总.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(previous[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(by_id.values())
print(json.dumps({"missing_EC_requested": len(missing), "batches": len(manifest), "failed_batches": sum(not (ROOT / row["文件"]).exists() for row in manifest), "total_Rhea_reactions": len(by_id)}, ensure_ascii=False), flush=True)
