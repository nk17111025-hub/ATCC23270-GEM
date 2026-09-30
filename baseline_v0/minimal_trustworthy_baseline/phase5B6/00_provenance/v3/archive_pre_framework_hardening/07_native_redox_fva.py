# -*- coding: utf-8 -*-
"""FVA for the updated v3 native-redox reactions + update summary."""
from __future__ import annotations

import json

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib
import v3_runtime as V3

OUT = C.OUT
V3DIR = OUT / "00_provenance" / "v3"


def fva(metas, rxns, cond, ov, reactions):
    S, ml, rl, ri, lb, ub = V3.setup_v3(metas, rxns, cond, ov)
    rows = []
    for r in reactions:
        if r not in ri:
            rows.append([r, "", ""])
            continue
        j = ri[r]
        cm = np.zeros(len(rl)); cm[j] = 1.0
        rmin = linprog(cm, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
        cM = np.zeros(len(rl)); cM[j] = -1.0
        rmax = linprog(cM, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
        rows.append([r, f"{rmin.x[j]:.4g}", f"{rmax.x[j]:.4g}"])
    return rows


def main():
    metas, rxns = V3.load_v3("T0")
    # add the two native reactions (already in updated XML, re-add in memory)
    for rid, stoich, ec, gpr in [
        ("NDH2_NATIVE", {"nadh[c]": -1, "h[c]": -1, "q8[c]": -1, "nad[c]": 1, "q8h2[c]": 1},
         "1.6.5.9", "RU820_RS08555"),
        ("NOX_NATIVE", {"nadh[c]": -1, "h[c]": -1, "o2[c]": -0.5, "nad[c]": 1, "h2o[c]": 1},
         "1.6.3.4", "RU820_RS08315"),
    ]:
        rxns[rid] = {"sbml_id": "R_" + rid, "name": rid, "reversible": False, "lb": 0.0,
                     "ub": 1000.0, "stoich": dict(stoich), "confidence": "", "ec": ec,
                     "pmid": "", "subsystem": "Oxidative Phosphorylation", "gpr": gpr,
                     "gpr2": gpr, "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                     "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "",
                     "table1_ub_tsul": ""}

    cond = C.COND
    ref_ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    vref = V3.max_biomass_v3(metas, rxns, cond, ref_ov)
    donor = lib.DONOR_EX[cond]
    ref = {"donor": vref[donor], "o2": vref["Ex_o2[e]"]}
    ov_r0 = {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
             donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
             "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}

    # FVA under R0
    keys = ["NDH2_NATIVE", "NOX_NATIVE", "NADHI", "GAPD1", "GAPD2", "PGK", "PDH",
            "BDGK", "G6PDH2", "CYTBD", "CYTBO3"]
    rows = fva(metas, rxns, cond, ov_r0, keys)
    C.write_tsv(V3DIR / "v3_native_redox_FVA.tsv", ["reaction", "R0_min", "R0_max"], rows)

    # FVA under WT (glucose closed, CO2 fixed) for the two new reactions (must be able to carry flux)
    ov_wt = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    rows_wt = fva(metas, rxns, cond, ov_wt, ["NDH2_NATIVE", "NOX_NATIVE"])
    C.write_tsv(V3DIR / "v3_native_redox_FVA_WT.tsv", ["reaction", "WT_min", "WT_max"], rows_wt)

    # update summary TSV
    C.write_tsv(V3DIR / "v3_native_redox_update.tsv",
                ["item", "value"],
                [["reaction_1", "NDH2_NATIVE"],
                 ["reaction_1_gene", "AFE_1854 -> RU820_RS08555 (WP_012536802.1)"],
                 ["reaction_1_EC", "1.6.5.9"],
                 ["reaction_1_equation", "nadh[c] + h[c] + q8[c] -> nad[c] + q8h2[c]"],
                 ["reaction_1_class", "NATIVE_MODEL_COMPLETION"],
                 ["reaction_1_mass_charge", "BALANCED (residual {} charge 0)"],
                 ["reaction_2", "NOX_NATIVE"],
                 ["reaction_2_gene", "AFE_1803 -> RU820_RS08315 (WP_009568931.1)"],
                 ["reaction_2_EC", "1.6.3.4"],
                 ["reaction_2_equation", "nadh[c] + h[c] + 0.5 o2[c] -> nad[c] + h2o[c]"],
                 ["reaction_2_class", "NATIVE_MODEL_COMPLETION_PROVISIONAL"],
                 ["reaction_2_mass_charge", "BALANCED (residual {} charge 0)"],
                 ["proton_pumping", "NONE (both non-electrogenic)"],
                 ["PMF_generation", "NONE"],
                 ["status", "V3_UPDATED_NATIVE_REDOX"]])

    print(json.dumps({"FVA_R0": rows, "FVA_WT": rows_wt}, indent=2))


if __name__ == "__main__":
    main()
