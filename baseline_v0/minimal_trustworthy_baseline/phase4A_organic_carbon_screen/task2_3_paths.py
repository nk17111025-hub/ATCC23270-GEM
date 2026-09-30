# -*- coding: utf-8 -*-
"""Tasks 2 and 3: canonical substrate-to-central-carbon pathway gap map + engineering distance."""
import csv
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

# Canonical catabolic paths. Each step: (from, to, rxn_id_or_None, kind, needs_reverse)
# kind in {exchange, transport, conversion}
PATHWAYS = {
    "glucose": [
        ("glc-B[e]", "glc-B[e]", None, "exchange", False),
        ("glc-B[e]", "glc-B[c]", None, "transport", False),
        ("glc-B[c]", "g6p-B[c]", "BDGK", "conversion", False),
        ("g6p-B[c]", "f6p-B[c]", "PGI1", "conversion", False),
        ("f6p-B[c]", "fdp[c]", "PFK", "conversion", False),
        ("fdp[c]", "g3p[c]", "FBA", "conversion", False),
        ("g3p[c]", "13dpg[c]", "GAPD1", "conversion", True),
        ("13dpg[c]", "3pg[c]", "PGK", "conversion", True),
        ("3pg[c]", "2pg[c]", "PGM1", "conversion", False),
        ("2pg[c]", "pep[c]", "ENO", "conversion", False),
        ("pep[c]", "pyr[c]", "PYK", "conversion", False),
        ("pyr[c]", "accoa[c]", "PDH", "conversion", False),
    ],
    "glycerol": [
        ("glyc[e]", "glyc[e]", None, "exchange", False),
        ("glyc[e]", "glyc[c]", None, "transport", False),
        ("glyc[c]", "glyc3p[c]", "GLYK", "conversion", False),
        ("glyc3p[c]", "dhap[c]", "G3PD2", "conversion", True),
        ("dhap[c]", "g3p[c]", "TPI", "conversion", False),
        ("g3p[c]", "13dpg[c]", "GAPD1", "conversion", True),
        ("13dpg[c]", "3pg[c]", "PGK", "conversion", True),
        ("3pg[c]", "2pg[c]", "PGM1", "conversion", False),
        ("2pg[c]", "pep[c]", "ENO", "conversion", False),
        ("pep[c]", "pyr[c]", "PYK", "conversion", False),
        ("pyr[c]", "accoa[c]", "PDH", "conversion", False),
    ],
    "ethanol": [
        ("etoh[e]", "etoh[e]", None, "exchange", False),
        ("etoh[e]", "etoh[c]", None, "transport", False),
        ("etoh[c]", "acald[c]", None, "conversion", False),
        ("acald[c]", "accoa[c]", None, "conversion", False),
    ],
    "acetate": [
        ("ac[e]", "ac[e]", None, "exchange", False),
        ("ac[e]", "ac[c]", None, "transport", False),
        ("ac[c]", "accoa[c]", "ACS", "conversion", False),
    ],
    "fructose": [
        ("fru[e]", "fru[e]", None, "exchange", False),
        ("fru[e]", "fru[c]", None, "transport", False),
        ("fru[c]", "f6p-B[c]", None, "conversion", False),
    ],
    "sucrose": [
        ("sucr[e]", "sucr[e]", None, "exchange", False),
        ("sucr[e]", "sucr[c]", None, "transport", False),
        ("sucr[c]", "glc-B[c]", None, "conversion", False),
    ],
    "lactate": [
        ("lac-L[e]", "lac-L[e]", None, "exchange", False),
        ("lac-L[e]", "lac-L[c]", None, "transport", False),
        ("lac-L[c]", "pyr[c]", None, "conversion", False),
    ],
    "pyruvate": [
        ("pyr[e]", "pyr[e]", None, "exchange", False),
        ("pyr[e]", "pyr[c]", None, "transport", False),
        ("pyr[c]", "accoa[c]", "PDH", "conversion", False),
    ],
    "formate": [
        ("for[e]", "for[e]", None, "exchange", False),
        ("for[e]", "for[c]", None, "transport", False),
        ("for[c]", "10fthf[c]", "FTHFL", "conversion", True),
        ("10fthf[c]", "central-carbon", None, "conversion", False),
    ],
    "methanol": [
        ("meoh[e]", "meoh[e]", None, "exchange", False),
        ("meoh[e]", "meoh[c]", None, "transport", False),
        ("meoh[c]", "for[c]", None, "conversion", False),
        ("for[c]", "central-carbon", None, "conversion", False),
    ],
    "4hba": [
        ("4hba[e]", "4hba[e]", "Ex_4hba[e]", "exchange", False),
        ("4hba[e]", "4hba[p]", "4HBAtex", "transport", False),
        ("4hba[p]", "4hba[c]", "4HBAtpp", "transport", False),
        ("4hba[c]", "central-carbon", None, "conversion", False),
    ],
}

# Curated engineering-distance (canonical, as-is allowed directions).
# Numeric = number of conversion reactions; "" = blocked/missing in needed direction.
DIST = {
    "glucose": ("1", "", "", "4",
                "reaches G6P/F6P via BDGK+PGI1; lower EMP blocked (GAPD1+PGK written gluconeogenic/irreversible); "
                "reaches biomass via upper-EMP + PPP (nucleotides) and glycogen"),
    "glycerol": ("", "", "", "2",
                 "GLYK->glyc3p (1) feeds lipid backbone (G3PAT1/PGSA); G3PD2 written DHAP->glyc3p (irreversible) so "
                 "glycerol oxidation to DHAP needs direction reversal; no route to pyruvate/acetyl-CoA as-is"),
    "ethanol": ("", "", "", "",
                "no ethanol/acetaldehyde metabolite or enzyme in model"),
    "acetate": ("", "", "1", "1",
                "ACS (acetate->acetyl-CoA) is present, active, irreversible forward; no route back to G6P (gluconeogenesis blocked)"),
    "fructose": ("", "", "", "",
                 "no free-fructose metabolite; only f6p-B/fdp present"),
    "sucrose": ("", "", "", "",
                "no sucrose metabolite or sucrose-cleavage enzyme"),
    "lactate": ("", "", "", "",
                "no lactate metabolite or lactate dehydrogenase"),
    "pyruvate": ("", "0", "1", "1",
                 "already a central metabolite; PDH -> acetyl-CoA present"),
    "formate": ("", "", "", "1",
                "formate->10fthf (FTHFL) is written in deformylation direction and irreversible (blocked); "
                "only oxidation to CO2 (FDH) or formyl donation to purines (GART); no C1 assimilation module"),
    "methanol": ("", "", "", "",
                 "no methanol metabolite or methanol dehydrogenase; no C1 assimilation module"),
    "4hba": ("", "", "", "",
             "complete native transport but no catabolic route from 4hba to central carbon (thiamine intermediate)"),
}


def step_status(rxns, rxn, kind, needs_reverse):
    if rxn is None:
        if kind == "exchange":
            return "MISSING_EXCHANGE"
        if kind == "transport":
            return "MISSING_TRANSPORT"
        return "MISSING_CONVERSION"
    if rxn not in rxns:
        return "MISSING_CONVERSION"
    rx = rxns[rxn]
    conf = rx["confidence"]
    try:
        confv = float(conf) if conf != "" else 99
    except Exception:
        confv = 99
    if needs_reverse and rx["lb"] >= 0:
        return "PRESENT_BUT_CLOSED"
    if confv <= 1:
        return "PRESENT_LOW_CONFIDENCE"
    if rx["gpr"] == "":
        return "PRESENT_NO_GPR"
    return "PRESENT_ACTIVE"


def required_action(st):
    return {
        "PRESENT_ACTIVE": "none",
        "PRESENT_BUT_CLOSED": "ENABLE_NATIVE_REACTION_OR_REVERSE_DIRECTION",
        "PRESENT_LOW_CONFIDENCE": "VERIFY_GENE",
        "PRESENT_NO_GPR": "VERIFY_GENE",
        "MISSING_TRANSPORT": "ADD_TRANSPORTER",
        "MISSING_CONVERSION": "ADD_ENZYME",
        "MISSING_EXCHANGE": "ADD_EXCHANGE",
    }.get(st, "UNKNOWN")


def main():
    metas, rxns = lib.parse_sbml_meta()

    gap_rows = []
    for spec in lib.SUBSTRATE_SPECS:
        sub = spec["substrate"]
        steps = PATHWAYS.get(sub, [])
        for i, (frm, to, rxn, kind, rev) in enumerate(steps):
            st = step_status(rxns, rxn, kind, rev)
            gap_rows.append({
                "substrate": sub,
                "step_index": i + 1,
                "from_metabolite": frm,
                "to_metabolite": to,
                "reaction_id": rxn if rxn else "",
                "status": st,
                "confidence": rxns[rxn]["confidence"] if (rxn and rxn in rxns) else "",
                "GPR": rxns[rxn]["gpr"] if (rxn and rxn in rxns) else "",
                "required_action": required_action(st),
                "comments": step_comment(sub, rxn),
            })

    dist_rows = []
    for spec in lib.SUBSTRATE_SPECS:
        sub = spec["substrate"]
        d = DIST.get(sub, ("", "", "", "", ""))
        n_missing = sum(1 for s in PATHWAYS.get(sub, []) if s[2] is None and s[3] == "conversion")
        n_missing_tr = sum(1 for s in PATHWAYS.get(sub, []) if s[2] is None and s[3] == "transport")
        n_missing_ex = sum(1 for s in PATHWAYS.get(sub, []) if s[2] is None and s[3] == "exchange")
        n_lowconf = 0
        n_reverse = 0
        n_open = 0
        for (frm, to, rxn, kind, rev) in PATHWAYS.get(sub, []):
            if rxn and rxn in rxns:
                rx = rxns[rxn]
                try:
                    cv = float(rx["confidence"]) if rx["confidence"] != "" else 99
                except Exception:
                    cv = 99
                if rx["gpr"] == "" or cv <= 1:
                    n_lowconf += 1
                if rev:
                    n_reverse += 1
                    if rx["lb"] >= 0:
                        n_open += 1
        transport_complete = spec["ex_met"] is not None and len(spec["transport_rxns"]) >= 1 and spec["c_met"] is not None
        # central connection: canonical path reaches G6P/F6P, pyruvate, or acetyl-CoA (a central-carbon hub)
        central = d[0] != "" or d[1] != "" or d[2] != "" or sub == "pyruvate"
        dist_rows.append({
            "substrate": sub,
            "A_missing_reactions": n_missing,
            "B_missing_transport": n_missing_tr,
            "C_missing_exchange": n_missing_ex,
            "D_low_confidence_or_noGPR_required": n_lowconf,
            "E_direction_reversal_required": n_reverse,
            "F_closed_reactions_to_open": n_open,
            "G_dist_to_G6P_F6P": d[0],
            "G_dist_to_pyruvate": d[1],
            "G_dist_to_acetyl_CoA": d[2],
            "G_dist_to_biomass_precursor": d[3],
            "native_path_complete": "TRUE" if transport_complete and (d[1] != "" or d[2] != "" or d[3] != "") else "FALSE",
            "transport_complete": "TRUE" if transport_complete else "FALSE",
            "central_connection_complete": "TRUE" if central else "FALSE",
            "distance_notes": d[4],
        })

    gf = ["substrate", "step_index", "from_metabolite", "to_metabolite", "reaction_id",
          "status", "confidence", "GPR", "required_action", "comments"]
    with (OUT / "substrate_pathway_gap_map.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gf, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in gap_rows:
            w.writerow(r)

    df = ["substrate", "A_missing_reactions", "B_missing_transport", "C_missing_exchange",
          "D_low_confidence_or_noGPR_required", "E_direction_reversal_required",
          "F_closed_reactions_to_open", "G_dist_to_G6P_F6P", "G_dist_to_pyruvate",
          "G_dist_to_acetyl_CoA", "G_dist_to_biomass_precursor", "native_path_complete",
          "transport_complete", "central_connection_complete", "distance_notes"]
    with (OUT / "organic_carbon_engineering_distance.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=df, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in dist_rows:
            w.writerow(r)

    print("gap_map rows", len(gap_rows), "distance rows", len(dist_rows))
    for r in dist_rows:
        print(r["substrate"], "G6P", r["G_dist_to_G6P_F6P"], "pyr", r["G_dist_to_pyruvate"],
              "accoa", r["G_dist_to_acetyl_CoA"], "bio", r["G_dist_to_biomass_precursor"],
              "transport", r["transport_complete"], "central", r["central_connection_complete"])


def step_comment(sub, rxn):
    c = {
        ("glycerol", "G3PD2"): "G3PD2 written DHAP->glyc3p (reduction) and irreversible; glycerol oxidation needs reversal",
        ("formate", "FTHFL"): "FTHFL written 10fthf->formate (deformylation) and irreversible; formate->10fthf blocked",
        ("glucose", "GAPD1"): "GAPD1/2 written 13dpg->g3p (gluconeogenic) and irreversible; glycolysis needs reversal",
        ("glucose", "PGK"): "PGK written 3pg->13dpg (gluconeogenic) and irreversible; glycolysis needs reversal",
        ("glucose", "BDGK"): "BDGK glucokinase-like, confidence 1, no GPR (WATCH item); forward glucose->G6P supported",
    }
    return c.get((sub, rxn), "")


if __name__ == "__main__":
    main()
