# -*- coding: utf-8 -*-
"""Build H1 = H0 (glucose candidate) + hypothetical OGDH reaction."""
from build_glucose_model import build_memory

BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}


def build_h1(transport="T0"):
    metas, rxns = build_memory(transport)
    rxns["OGDH"] = {
        "sbml_id": "R_OGDH_hypo",
        "name": "2-oxoglutarate dehydrogenase (HYPOTHETICAL_ENGINEERING_REACTION)",
        "reversible": False,
        "lb": 0.0, "ub": 1000.0,
        "stoich": {"akg[c]": -1.0, "coa[c]": -1.0, "nad[c]": -1.0,
                   "succoa[c]": 1.0, "co2[c]": 1.0, "nadh[c]": 1.0},
        "confidence": "", "ec": "1.2.1.105", "pmid": "",
        "subsystem": "HYPOTHETICAL_ENGINEERING_REACTION",
        "gpr": "", "gpr2": "", "protein": "", **BLANK,
    }
    return metas, rxns


if __name__ == "__main__":
    metas, rxns = build_h1("T0")
    print("H1 n_rxn", len(rxns), "OGDH present:", "OGDH" in rxns)
    print("OGDH stoich:", rxns["OGDH"]["stoich"])
