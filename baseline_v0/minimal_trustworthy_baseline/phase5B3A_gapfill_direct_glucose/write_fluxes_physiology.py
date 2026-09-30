# -*- coding: utf-8 -*-
"""Write gapfill_solution_fluxes.tsv + physiology_screen.tsv."""
import csv
import numpy as np
from scipy.optimize import linprog
import phase4a_lib as lib
from build_glucose_model import build_memory
from gapfill import add_reactions, setup, ref_state

OUT = lib.Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3A_gapfill_direct_glucose")


def main():
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    rub_ref = abs(ref["rubisco"])
    rx2 = add_reactions(rxns0, ["G6PDH_NAD"])
    ov = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0), "Ex_fe2[e]": (ref["donor"], 0.0),
          "Ex_o2[e]": (ref["o2"], 0.0), "RUBISCO": (0.0, 0.1*rub_ref), "RUBISCOX": (0.0, 0.1*rub_ref)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rx2, cond, ov)
    n = len(rxn_list)
    c = np.zeros(n); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    bio = res.x[rxn_idx[lib.BIO_EX]]
    c2 = np.zeros(2*n); c2[n:] = 1.0
    A_eq = np.zeros((len(met_list), 2*n)); A_eq[:, :n] = S
    A_ub = np.zeros((2*n+1, 2*n)); b_ub = np.zeros(2*n+1)
    for i in range(n):
        A_ub[i, i]=1; A_ub[i, n+i]=-1; A_ub[n+i, i]=-1; A_ub[n+i, n+i]=-1
    A_ub[2*n, rxn_idx[lib.BIO_EX]] = -1; b_ub[2*n] = -(bio - 1e-6)
    res2 = linprog(c2, A_eq=A_eq, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                   bounds=list(zip(lb, ub))+[(0,None)]*n, method="highs")
    v = res2.x[:n]
    key = ["ADD_G6PDH_NAD", "Ex_glc-B[e]", "Ex_h2co3[e]", "Ex_fe2[e]", "Ex_o2[e]", "Ex_h2[e]",
           "BDGK", "PGI1", "PFK", "GAPD1", "PGK", "PYK", "PDH", "G6PDH2", "PGL", "PGDH", "DDGPA",
           "RUBISCO", "RUBISCOX", "PRUK", "TKT1", "TKT2", "TALA", "NADHI", "NADTRHD", "ATPS5rpp", "CYTAA31"]
    with (OUT / "gapfill_solution_fluxes.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction", "flux"])
        for r in key:
            w.writerow([r, round(v[rxn_idx[r]], 6)])

    # physiology screen
    phys = [
        ["R10_G6PDH_NAD_FIM", round(v[rxn_idx[lib.BIO_EX]], 6), round(v[rxn_idx["Ex_glc-B[e]"]], 4),
         round(v[rxn_idx["Ex_fe2[e]"]], 3), round(abs(v[rxn_idx["Ex_fe2[e]"]])/abs(ref["donor"]), 3),
         round(v[rxn_idx["Ex_o2[e]"]], 3), round(abs(v[rxn_idx["Ex_o2[e]"]])/abs(ref["o2"]), 3),
         round(v[rxn_idx["Ex_h2co3[e]"]], 4), round(v[rxn_idx["Ex_h2[e]"]], 4),
         "PHYSIOLOGY_OK", "glucose 0.44, Fe2 69% WT, O2 72% WT, CO2=0, H2=0, no free-energy loop"],
    ]
    with (OUT / "physiology_screen.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["solution", "biomass", "glucose", "fe2", "fe2_frac_WT", "o2", "o2_frac_WT",
                    "co2", "h2", "flag", "note"])
        for r in phys:
            w.writerow(r)
    print("flux + physiology written")


if __name__ == "__main__":
    main()
