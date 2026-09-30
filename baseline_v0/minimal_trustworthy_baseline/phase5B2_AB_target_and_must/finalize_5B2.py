# -*- coding: utf-8 -*-
"""pFBA sanity + physiology flags + report data for Phase 5B-2A/2B."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory
from must_lib import setup

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B2_AB_target_and_must")


def min_l1(metas, rxns, cond, override, biomass_min):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    n = len(rxn_list)
    # variables [v (n), t (n)]
    c = np.zeros(2 * n)
    c[n:] = 1.0
    A_eq = np.zeros((len(met_list), 2 * n))
    A_eq[:, :n] = S
    # |v_i| <= t_i
    A_ub = np.zeros((2 * n + 1, 2 * n))
    b_ub = np.zeros(2 * n + 1)
    for i in range(n):
        A_ub[i, i] = 1.0
        A_ub[i, n + i] = -1.0
        A_ub[n + i, i] = -1.0
        A_ub[n + i, n + i] = -1.0
    A_ub[2 * n, rxn_idx[lib.BIO_EX]] = -1.0
    b_ub[2 * n] = -biomass_min
    bnd = list(zip(lb, ub)) + [(0.0, None)] * n
    res = linprog(c, A_eq=A_eq, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub,
                  bounds=bnd, method="highs")
    if not res.success:
        return None
    return res.x[:n]


def main():
    metas, rxns = build_memory("T0")
    priority = ["Ex_glc-B[e]", "BDGK", "PGI1", "PFK", "FBA", "GAPD1", "GAPD2", "PGK", "PYK", "PDH",
                "G6PDH2", "TKT1", "TKT2", "TALA", "PRUK", "RUBISCO", "RUBISCOX", "NADTRHD",
                "CS", "PPC", "MDH", "FUM", "FDH", "MALS", "ACS", "ATPS5rpp"]
    # WT pFBA
    ref = {}
    for cond in ["FIM", "TTM", "TSM"]:
        x = min_l1(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}, 0.052076)
        rxn_idx = {r: j for j, r in enumerate(rxns)}
        ref[cond] = {"donor": x[rxn_idx[lib.DONOR_EX[cond]]], "o2": x[rxn_idx["Ex_o2[e]"]],
                     "co2": x[rxn_idx["Ex_h2co3[e]"]], "bio": x[rxn_idx[lib.BIO_EX]]}

    sanity_rows = []
    for cond in ["FIM", "TTM", "TSM"]:
        donor = lib.DONOR_EX[cond]
        for target, t in [("A", dict(donor_min=0.90, o2_max=1.20, co2_max=0.75)),
                          ("B", dict(donor_min=0.70, o2_max=1.20, co2_max=0.50)),
                          ("C", dict(donor_min=0.50, o2_max=None, co2_max=0.0))]:
            ov = {"Ex_glc-B[e]": (-1.0, 0.0), "Ex_h2co3[e]": (-t["co2_max"] * abs(ref[cond]["co2"]), 0.0),
                  donor: (-1000.0, -t["donor_min"] * abs(ref[cond]["donor"]))}
            if t["o2_max"] is not None:
                ov["Ex_o2[e]"] = (-t["o2_max"] * abs(ref[cond]["o2"]), 1000.0)
            x = min_l1(metas, rxns, cond, ov, 0.052076)
            if x is None:
                continue
            rxn_idx = {r: j for j, r in enumerate(rxns)}
            row = {"condition": cond, "target": target}
            for r in priority:
                row[r] = round(x[rxn_idx[r]], 5) if r in rxn_idx else ""
            row["donor_uptake"] = round(x[rxn_idx[donor]], 5)
            row["o2_uptake"] = round(x[rxn_idx["Ex_o2[e]"]], 5)
            row["co2_uptake"] = round(x[rxn_idx["Ex_h2co3[e]"]], 5)
            row["glucose_uptake"] = round(x[rxn_idx["Ex_glc-B[e]"]], 5)
            sanity_rows.append(row)
    f = ["condition", "target", "glucose_uptake", "donor_uptake", "o2_uptake", "co2_uptake"] + priority
    with (OUT / "pfba_target_sanity.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sanity_rows:
            w.writerow(r)

    # physiology flags
    flags = []
    for r in sanity_rows:
        cond = r["condition"]
        glc = r["glucose_uptake"]
        donor_frac = abs(r["donor_uptake"]) / max(abs(ref[cond]["donor"]), 1e-9)
        o2_frac = abs(r["o2_uptake"]) / max(abs(ref[cond]["o2"]), 1e-9)
        issues = []
        if abs(glc) > 1.0:
            issues.append("high_glucose_flux")
        if donor_frac > 1.05:
            issues.append("donor_above_WT")
        if o2_frac > 1.2:
            issues.append("o2_above_120pct")
        if abs(r["co2_uptake"]) < 1e-5:
            issues.append("co2_near_zero")
        if r.get("NADTRHD", 0) and abs(r["NADTRHD"]) > 50:
            issues.append("high_NADTRHD")
        flag = "PHYSIOLOGY_IMPLAUSIBLE" if ("high_glucose_flux" in issues or "o2_above_120pct" in issues) else ("PHYSIOLOGY_WATCH" if issues else "PHYSIOLOGY_OK")
        flags.append({"condition": cond, "target": r["target"], "glucose_uptake": glc,
                      "donor_frac_WT": round(donor_frac, 3), "o2_frac_WT": round(o2_frac, 3),
                      "co2_uptake": r["co2_uptake"], "issues": ";".join(issues), "flag": flag})
    f = ["condition", "target", "glucose_uptake", "donor_frac_WT", "o2_frac_WT", "co2_uptake", "issues", "flag"]
    with (OUT / "physiology_flags.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in flags:
            w.writerow(r)

    # must pair candidates (minimal: report the single MUST co-occurrence, not a full MILP)
    with (OUT / "must_pair_candidates.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pair", "type", "basis", "note"])
        w.writerow(["Ex_glc-B[e] + GLCtpp", "MUST-ON pair", "glucose uptake requires exchange + transport",
                    "single-reaction MUST co-occurrence; no combinatorial MILP run"])
        w.writerow(["PGI1 + BDGK", "MUST-UP pair", "glucose -> G6P -> F6P requires both",
                    "single-reaction MUST co-occurrence; no combinatorial MILP run"])

    print("finalize 5B2 complete")
    for r in flags:
        print(r["condition"], r["target"], r["flag"], "glc", r["glucose_uptake"], "donor_frac", r["donor_frac_WT"], "co2", r["co2_uptake"])


if __name__ == "__main__":
    main()
