# -*- coding: utf-8 -*-
"""Download BiGG central-carbon reactions and map to the model metabolite space."""
import json
import time
import urllib.request

import phase4a_lib as lib
from build_glucose_model import build_memory

# BiGG reaction IDs relevant to central carbon / glucose assimilation / anaplerosis / redox
BIGG_RXNS = [
    "HEX1", "GLK", "PGI", "PFK", "FBP", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK",
    "PDH", "PFL", "CS", "ACONTa", "ACONTb", "ICDHyr", "AKGDH", "SUCOAS", "SUCD1", "SUCDi",
    "FUM", "MDH", "ME1", "ME2", "PC", "PEPCK", "PPC", "PPCK", "ICL", "MALS", "PPS",
    "G6PDH2r", "G6PDHr", "PGL", "GND", "EDD", "EDA", "TKT1", "TKT2", "TALA", "RPI", "RPE",
    "NADH16", "NADH17pp", "NADH12pp", "CYOOm2", "CYTBO3_4pp", "GLUDy", "GLNS", "GLUSy",
    "ASPTA", "ALATA_L", "ACONT", "F6PA", "F6PP", "G1PACT", "ACCoAC", "ALCD2x", "ALDD2x",
    "ACKr", "PTAr", "ACACT1r", "LDH_D", "D_LACD2", "G3PD1ir", "GLYCK", "GLYCDx", "GTHOr",
    "PPA", "PPKr", "PPDK", "EX_glc__D_e",
]

MET_MAP = {
    "glc__D": "glc-B[c]", "g6p": "g6p-B[c]", "g6p__B": "g6p-B[c]", "f6p": "f6p-B[c]",
    "fdp": "fdp[c]", "g3p": "g3p[c]", "dhap": "dhap[c]", "13dpg": "13dpg[c]", "3pg": "3pg[c]",
    "2pg": "2pg[c]", "pep": "pep[c]", "pyr": "pyr[c]", "accoa": "accoa[c]", "oaa": "oaa[c]",
    "cit": "cit[c]", "icit": "icit[c]", "akg": "akg[c]", "succoa": "succoa[c]", "succ": "succ[c]",
    "fum": "fum[c]", "mal__L": "mal-L[c]", "glx": "glx[c]", "glyc": "glyc[c]", "glyc3p": "glyc3p[c]",
    "r5p": "r5p[c]", "ru5p__D": "ru5p-D[c]", "xu5p__D": "xu5p-D[c]", "e4p": "e4p[c]", "s7p": "s7p[c]",
    "6pgc": "6pgc[c]", "6pgl": "6pgl[c]", "2ddg6p": "2ddg6p[c]", "rubp__D": "rb15bp[c]",
    "g1p": "g1p[c]", "g1p__B": "g1p-B[c]", "6pgl": "6pgl[c]", "acon__C": "acon-C[c]",
    "nad": "nad[c]", "nadh": "nadh[c]", "nadp": "nadp[c]", "nadph": "nadph[c]",
    "atp": "atp[c]", "adp": "adp[c]", "amp": "amp[c]", "pi": "pi[c]", "ppi": "ppi[c]",
    "h": "h[c]", "h2o": "h2o[c]", "co2": "co2[c]", "coa": "coa[c]", "q8": "q8[c]", "q8h2": "q8h2[c]",
    "for": "for[c]", "ac": "ac[c]", "o2": "o2[c]", "h2o2": "h2o2[c]", "glu__L": "glu-L[c]",
    "gln__L": "gln-L[c]", "asp__L": "asp-L[c]", "ala__L": "ala-L[c]", "h2": "h2[c]",
    "lac__D": "lac-D[c]", "lac__L": "lac-L[c]", "glyc__R": "glyc-R[c]", "dha": "dhap[c]",
}


def fetch_reaction(bigg_id):
    url = f"http://bigg.ucsd.edu/api/v2/universal/reactions/{bigg_id}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "codex"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)
    except Exception as e:
        return {"error": str(e)}


def main():
    metas, rxns = build_memory("T0")
    model_mets = set(metas.keys())
    # check which metabolites in MET_MAP exist in the model
    mapped_mets = set()
    for bigg, model in MET_MAP.items():
        if model in model_mets:
            mapped_mets.add(bigg)
        else:
            print("UNMAPPED-MODEL-MISSING", bigg, "->", model)

    rows = []
    for bid in BIGG_RXNS:
        d = fetch_reaction(bid)
        if "error" in d:
            print("FETCH_ERR", bid, d["error"])
            continue
        mets = d.get("metabolites", [])
        stoich = {}
        unmapped = []
        for m in mets:
            bm = m.get("bigg_id", "")
            if bm in MET_MAP:
                stoich[MET_MAP[bm]] = stoich.get(MET_MAP[bm], 0.0) + float(m.get("stoichiometry", 0))
            else:
                # skip water/proton/generic, but record unmapped
                unmapped.append(bm)
        # keep only if all non-water reactants map
        real_unmapped = [u for u in unmapped if u not in ("h2o", "h")]
        if real_unmapped:
            rows.append({"bigg_id": bid, "name": d.get("name", ""), "mappable": "no",
                         "unmapped": ",".join(sorted(set(real_unmapped)))})
            continue
        rows.append({"bigg_id": bid, "name": d.get("name", ""), "mappable": "yes",
                     "stoich": json.dumps(stoich), "reversible": str(d.get("reversibility", False))})
        time.sleep(0.05)

    with open("_bigg_reactions.json", "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print("downloaded", len(rows), "reactions; mappable", sum(1 for r in rows if r.get("mappable") == "yes"))
    print("mapped metabolites:", len(mapped_mets))


if __name__ == "__main__":
    main()
