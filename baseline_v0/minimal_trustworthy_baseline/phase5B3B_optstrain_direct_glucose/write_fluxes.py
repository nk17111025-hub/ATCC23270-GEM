# -*- coding: utf-8 -*-
"""Write solution_fluxes.tsv for the key G6PDH_NAD (R10) solution."""
import csv
import numpy as np
from scipy.optimize import linprog
import phase4a_lib as lib
from build_glucose_model import build_memory
from optstrain import add_reactions, setup, ref_state

OUT = lib.Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3B_optstrain_direct_glucose")


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
    key = ["NON_G6PDH_NAD", "Ex_glc-B[e]", "Ex_h2co3[e]", "Ex_fe2[e]", "Ex_o2[e]", "Ex_h2[e]",
           "BDGK", "PGI1", "PFK", "GAPD1", "PGK", "PYK", "PDH", "G6PDH2", "PGL", "PGDH", "DDGPA",
           "RUBISCO", "RUBISCOX", "PRUK", "TKT1", "TKT2", "TALA", "NADHI", "NADTRHD", "ATPS5rpp", "CYTAA31"]
    with (OUT / "solution_fluxes.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["reaction", "flux"])
        for r in key:
            w.writerow([r, round(v[rxn_idx[r]], 6)])
    print("solution_fluxes written, biomass", round(bio, 6))


if __name__ == "__main__":
    main()
