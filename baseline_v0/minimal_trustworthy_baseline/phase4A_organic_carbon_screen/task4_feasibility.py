# -*- coding: utf-8 -*-
"""Task 4: native organic-carbon feasibility (no new reactions added)."""
import csv
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

CAPS = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

CENTRAL_RXNS = ["RUBISCO", "PPC", "CS", "MDH", "FUM", "PDH", "PYK", "PFK", "FBA",
                "GAPD1", "PGK", "PGI1", "BDGK", "GLYK", "ACS", "FDH", "GART", "FTHFL"]


def record(metas, rxns, cond, label, override=None):
    fluxes, bounds, opt = lib.run_pfba(metas, rxns, cond, override=override)
    if fluxes is None:
        return None
    donor = lib.DONOR_EX[cond]
    row = {
        "condition": cond,
        "case": label,
        "growth": round(opt, 9),
        "organic_uptake": round(fluxes.get("Ex_4hba[e]", 0.0), 6),
        "donor_uptake": round(fluxes.get(donor, 0.0), 6),
        "co2_hco3_uptake": round(fluxes.get("Ex_h2co3[e]", 0.0), 6),
        "oxygen_uptake": round(fluxes.get("Ex_o2[e]", 0.0), 6),
        "ATPM": round(fluxes.get("ATPM", 0.0), 6),
    }
    for r in CENTRAL_RXNS:
        row[r] = round(fluxes.get(r, 0.0), 6) if r in fluxes else ""
    return row


def main():
    metas, rxns = lib.parse_sbml_meta()
    rows = []
    for cond in ["FIM", "TTM", "TSM"]:
        base = record(metas, rxns, cond, "baseline_inorganic")
        if base:
            rows.append(base)
        for cap in CAPS:
            r = record(metas, rxns, cond, f"4hba_avail_cap{cap}", override={"Ex_4hba[e]": (-cap, 1000.0)})
            if r:
                rows.append(r)

    fields = ["condition", "case", "growth", "organic_uptake", "donor_uptake",
              "co2_hco3_uptake", "oxygen_uptake", "ATPM"] + CENTRAL_RXNS
    with (OUT / "native_organic_carbon_feasibility.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("wrote native feasibility rows", len(rows))
    for r in rows:
        print(r["condition"], r["case"], "growth", r["growth"], "org_up", r["organic_uptake"],
              "donor", r["donor_uptake"], "co2", r["co2_hco3_uptake"], "o2", r["oxygen_uptake"])


if __name__ == "__main__":
    main()
