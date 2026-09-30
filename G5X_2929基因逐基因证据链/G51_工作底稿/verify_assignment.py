import json
from pathlib import Path

here = Path(__file__).resolve().parent
rows = json.loads((here / "G51_evidence_inventory.json").read_text(encoding="utf-8"))
assignment = (here / "G51_drive_assignment.txt").read_text(encoding="utf-8").splitlines()
errors = []
for i, (row, line) in enumerate(zip(rows, assignment), 1):
    p = line.split("｜")
    expected = [f"{i:04d}", row["locus"], row["wp"], f"{row['chromosome']}:{row['start']}-{row['end']}({row['strand']})", row["product"]]
    if p != expected:
        errors.append((i, p, expected))
print("assignment_rows", len(assignment), "inventory_rows", len(rows), "mismatches", len(errors))
for error in errors[:10]:
    print(error)
