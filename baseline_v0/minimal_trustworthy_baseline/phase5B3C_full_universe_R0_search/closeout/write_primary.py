# -*- coding: utf-8 -*-
"""Write primary_blocked_precursors.tsv (upstream-trace of first blocks)."""
import csv
from pathlib import Path

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3C_full_universe_R0_search\closeout")

ROWS = [
    ["3pg[c]", "g3p[c] (reachable)", "GAPD1/GAPD2 glycolytic (g3p->13dpg)", "GAPD1/GAPD2",
     "reverse (glycolytic)", "NAD+", "glycolytic GAP->13dpg produces NADH with no sink in the autotrophic redox economy"],
    ["13dpg[c]", "g3p[c]", "g3p -> 13dpg", "GAPD1/GAPD2", "reverse", "NAD+", "same redox block as 3PG"],
    ["2pg[c]/pep[c]", "3pg[c]", "3pg -> 2pg -> pep", "PGM1/ENO", "forward", "none", "downstream of 3PG (PGM/ENO reversible but no 3PG supply)"],
    ["r5p[c]", "ru5p-D[c] (reachable)", "ru5p -> r5p", "RPI", "reverse (blocked, lb=0)", "none", "RPI written r5p->ru5p irreversible; R5P for nucleotides not supplied"],
    ["e4p[c]", "g3p[c]+s7p[c]", "TKT/TALA (S7P requires 3PG-derived carbon)", "TKT1/TKT2/TALA", "forward", "none", "non-oxidative PPP cannot regenerate E4P without the 3PG/Calvin pool"],
    ["accoa[c]", "pyr[c] (reachable)", "pyr -> accoa + CO2", "PDH", "forward", "CoA, NAD", "PDH blocked by CoA/redox coupling without Rubisco-supported TCA"],
    ["oaa[c]/akg[c]", "accoa[c]", "TCA + anaplerosis", "CS/ICDHyr/PPC", "forward", "CO2 (PPC)", "TCA/anaplerosis need acetyl-CoA and carboxylation"],
    ["ser-L[c]/gly[c]", "3pg[c]", "3pg -> ser -> gly", "PGCD/...", "forward", "NAD/NADP", "downstream of 3PG"],
    ["asp-L[c]/glu-L[c]/gln-L[c]", "oaa[c]/akg[c]", "transamination", "ASPTA/GLUSy", "forward", "none", "downstream of OAA/alpha-KG"],
]


def main():
    with (OUT / "primary_blocked_precursors.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["precursor", "upstream_reachable", "first_blocked_conversion", "reaction_id",
                    "reaction_direction", "cofactor_requirement", "reason"])
        for r in ROWS:
            w.writerow(r)
    print("primary_blocked_precursors written")


if __name__ == "__main__":
    main()
