# -*- coding: utf-8 -*-
"""Tasks 5 and 6: hypothetical-import probes (carbon contribution) + minimal completion."""
import csv
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

# substrates with an existing cytosolic form (probed by adding a hypothetical cytosolic exchange)
PROBE = [
    ("glucose", "glc-B[c]"),
    ("glycerol", "glyc[c]"),
    ("acetate", "ac[c]"),
    ("pyruvate", "pyr[c]"),
    ("formate", "for[c]"),
    ("4hba", "4hba[c]"),
]


def add_exchange(reactions, met_key, rxn_id):
    rx2 = dict(reactions)
    rx2[rxn_id] = {
        "sbml_id": rxn_id, "name": "hypothetical exchange", "reversible": True,
        "lb": -1000.0, "ub": 1000.0, "stoich": {met_key: -1.0},
        "confidence": "", "ec": "", "pmid": "", "subsystem": "hypothetical",
        "gpr": "", "gpr2": "", "protein": "",
        "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
        "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": "",
    }
    return rx2


def probe(metas, rxns, cond, met_key, sub, cap=10.0, close_co2=False):
    rxn_id = "EX_hypo_" + sub
    rx2 = add_exchange(rxns, met_key, rxn_id)
    ov = {rxn_id: (-cap, 1000.0)}
    if close_co2:
        ov["Ex_h2co3[e]"] = (0.0, 0.0)
    fluxes, bounds, opt = lib.run_pfba(metas, rx2, cond, override=ov)
    if fluxes is None:
        return None
    return {
        "growth": opt,
        "organic_uptake": fluxes.get(rxn_id, 0.0),
        "donor_uptake": fluxes.get(lib.DONOR_EX[cond], 0.0),
        "co2_uptake": fluxes.get("Ex_h2co3[e]", 0.0),
        "oxygen_uptake": fluxes.get("Ex_o2[e]", 0.0),
        "ATPM": fluxes.get("ATPM", 0.0),
        "RUBISCO": fluxes.get("RUBISCO", 0.0),
        "PDH": fluxes.get("PDH", 0.0),
        "ACS": fluxes.get("ACS", 0.0),
        "GLYK": fluxes.get("GLYK", 0.0),
        "BDGK": fluxes.get("BDGK", 0.0),
        "PYK": fluxes.get("PYK", 0.0),
        "CS": fluxes.get("CS", 0.0),
        "FDH": fluxes.get("FDH", 0.0),
    }


def main():
    metas, rxns = lib.parse_sbml_meta()
    contrib_rows = []
    for cond in ["FIM", "TTM", "TSM"]:
        for sub, met in PROBE:
            base = probe(metas, rxns, cond, met, sub, cap=10.0, close_co2=False)
            het = probe(metas, rxns, cond, met, sub, cap=10.0, close_co2=True)
            if base is None:
                continue
            uptake = base["organic_uptake"]
            used = uptake < -1e-6
            het_g = het["growth"] if het else 0.0
            donor_nonzero = base["donor_uptake"] < -1e-6
            rubisco_active = abs(base["RUBISCO"]) > 1e-6
            notes = ""
            if not used:
                cls = "NO_UTILIZATION"
            elif sub == "formate":
                cls = "ARTIFACT_SUSPECTED"
                notes = "formate is oxidized to CO2 (FDH) and the CO2 is re-fixed by RUBISCO; no C1 assimilation route"
            elif sub == "glycerol":
                cls = "UPTAKE_WITHOUT_BIOMASS_CONTRIBUTION"
                notes = "glycerol->glyc3p feeds only the lipid backbone (G3PD2 oxidation blocked); negligible growth change"
            elif donor_nonzero and rubisco_active:
                cls = "MIXOTROPHY_LIKE"
            elif donor_nonzero:
                cls = "UPTAKE_WITHOUT_BIOMASS_CONTRIBUTION"
            elif rubisco_active:
                cls = "LITHOHETEROTROPHY_LIKE"
            else:
                cls = "POTENTIAL_HETEROTROPHY"
            if sub == "glucose":
                notes = "glucose reaches pyruvate/acetyl-CoA only via PPP->RuBP->RUBISCO->3PG route (glycolysis and oxidative-PPP blocked); CO2 fixation stays active"
            contrib_rows.append({
                "condition": cond,
                "substrate": sub,
                "probe": "mixotrophic_co2_open",
                "growth": round(base["growth"], 9),
                "organic_uptake": round(uptake, 6),
                "donor_uptake": round(base["donor_uptake"], 6),
                "co2_uptake": round(base["co2_uptake"], 6),
                "oxygen_uptake": round(base["oxygen_uptake"], 6),
                "RUBISCO": round(base["RUBISCO"], 6),
                "PDH": round(base["PDH"], 6),
                "ACS": round(base["ACS"], 6),
                "GLYK": round(base["GLYK"], 6),
                "BDGK": round(base["BDGK"], 6),
                "classification": cls,
                "notes": notes,
            })
            contrib_rows.append({
                "condition": cond,
                "substrate": sub,
                "probe": "heterotrophic_co2_closed",
                "growth": round(het_g, 9) if het else "",
                "organic_uptake": round(het["organic_uptake"], 6) if het else "",
                "donor_uptake": round(het["donor_uptake"], 6) if het else "",
                "co2_uptake": round(het["co2_uptake"], 6) if het else "",
                "oxygen_uptake": round(het["oxygen_uptake"], 6) if het else "",
                "RUBISCO": round(het["RUBISCO"], 6) if het else "",
                "PDH": round(het["PDH"], 6) if het else "",
                "ACS": round(het["ACS"], 6) if het else "",
                "GLYK": round(het["GLYK"], 6) if het else "",
                "BDGK": round(het["BDGK"], 6) if het else "",
                "classification": "",
                "notes": "CO2/HCO3 exchange closed; tests sole-carbon-source feasibility",
            })

    fields = ["condition", "substrate", "probe", "growth", "organic_uptake", "donor_uptake",
              "co2_uptake", "oxygen_uptake", "RUBISCO", "PDH", "ACS", "GLYK", "BDGK",
              "classification", "notes"]
    with (OUT / "organic_carbon_contribution_analysis.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in contrib_rows:
            w.writerow(r)

    print("wrote contribution rows", len(contrib_rows))
    for r in contrib_rows:
        if r["probe"] == "mixotrophic_co2_open":
            print(r["condition"], r["substrate"], "growth", r["growth"], "up", r["organic_uptake"],
                  "co2", r["co2_uptake"], "het_g", (r["classification"]))


if __name__ == "__main__":
    main()
