# -*- coding: utf-8 -*-
"""Phase 5B-2C: single-reaction capacity / suppression scan for Target D feasibility."""
import csv
import json
from pathlib import Path

import numpy as np

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B2C_multi_method_route_search")

CANDIDATES = ["NADHI", "NADTRHD", "CYTAA31", "CYTBC1", "CYT1", "CYT2", "CYTA2", "ATPS5rpp",
              "BDGK", "PGI1", "PFK", "FBA", "GAPD1", "GAPD2", "PGK", "PGM1", "ENO", "PYK", "PDH",
              "G6PDH2", "PGL", "PGDH", "DDGPA", "TKT1", "TKT2", "TALA", "RPI", "RPE", "PRUK",
              "RUBISCO", "RUBISCOX", "FBP", "SBPASE", "PPC", "CS", "ICDHyr", "MDH", "FUM", "SUCD",
              "SUCOAS", "ACS", "FDH", "HYD3pp", "Htpp", "GLCtex", "GLCtpp", "Ex_glc-B[e]"]

MULTS = [1.10, 1.25, 1.50]


def setup(metas, rxns, cond, override):
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n); ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]; lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN:
            continue
        rx = rxns[r]
        lv = rx["table1_" + lbc]; uv = rx["table1_" + ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv); ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl; ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol; ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def max_bio(metas, rxns, cond, override):
    from scipy.optimize import linprog
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def ref_state(metas, rxns, cond):
    x = max_bio(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    # need fluxes for donor/o2
    from scipy.optimize import linprog
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    v = res.x
    return {"biomass": v[rxn_idx[lib.BIO_EX]], "donor": v[rxn_idx[lib.DONOR_EX[cond]]],
            "o2": v[rxn_idx["Ex_o2[e]"]], "co2": v[rxn_idx["Ex_h2co3[e]"]]}


def targetD_override(cond, ref, glc_cap, intervention):
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
          donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0)}
    if intervention:
        ov.update(intervention)
    return ov


def apply_intervention(rxns, rxn, kind):
    """Return an override dict for a single-reaction intervention."""
    rx = rxns[rxn]
    lb, ub = rx["lb"], rx["ub"]
    if kind == "open_reverse":
        if lb >= 0:
            return {rxn: (-1000.0, ub)}
        return None  # already reversible
    if kind == "widen_ub_5x":
        return {rxn: (lb, max(ub, 5000.0) if ub > 0 else ub)}
    if kind == "widen_ub_2x":
        return {rxn: (lb, ub * 2.0 if ub > 0 and ub < 1000 else ub)}
    if kind == "close":
        return {rxn: (0.0, 0.0)}
    if kind == "narrow_half":
        return {rxn: (lb * 0.5, ub * 0.5)}
    if kind == "narrow_quarter":
        return {rxn: (lb * 0.25, ub * 0.25)}
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns, cond)
    print("ref:", {k: round(v, 4) for k, v in ref.items()})

    # baseline Target D (no intervention)
    base = max_bio(metas, rxns, cond, targetD_override(cond, ref, 5.0, None))
    print("Target D baseline max biomass:", round(base, 7) if base else "infeasible", "(WT", round(ref["biomass"], 7), ")")

    rows = []
    kinds = ["open_reverse", "widen_ub_5x", "widen_ub_2x", "close", "narrow_half", "narrow_quarter"]
    for rxn in CANDIDATES:
        if rxn not in rxns:
            continue
        for kind in kinds:
            iv = apply_intervention(rxns, rxn, kind)
            if iv is None:
                continue
            bio = max_bio(metas, rxns, cond, targetD_override(cond, ref, 5.0, iv))
            if bio is None:
                continue
            gain = (bio - ref["biomass"]) / ref["biomass"] * 100
            rows.append({"reaction": rxn, "intervention": kind, "max_biomass": round(bio, 7),
                         "biomass_gain_pct": round(gain, 4),
                         "D10": "YES" if bio >= 1.10 * ref["biomass"] - 1e-6 else "no",
                         "D25": "YES" if bio >= 1.25 * ref["biomass"] - 1e-6 else "no",
                         "D50": "YES" if bio >= 1.50 * ref["biomass"] - 1e-6 else "no"})

    # sort by gain
    rows.sort(key=lambda r: -r["biomass_gain_pct"])
    f = ["reaction", "intervention", "max_biomass", "biomass_gain_pct", "D10", "D25", "D50"]
    with (OUT / "single_reaction_capacity_scan.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("=== top 20 single interventions (by gain) ===")
    for r in rows[:20]:
        print(f"  {r['reaction']:12s} {r['intervention']:18s} gain={r['biomass_gain_pct']:8.2f}% D10={r['D10']} D25={r['D25']} D50={r['D50']}")

    # which reactions achieve D10?
    d10 = [r for r in rows if r["D10"] == "YES"]
    print("=== reactions achieving D10 (+10%) ===")
    seen = set()
    for r in d10:
        if r["reaction"] not in seen:
            seen.add(r["reaction"])
            print(f"  {r['reaction']} (best gain {max(x['biomass_gain_pct'] for x in d10 if x['reaction']==r['reaction']):.2f}%)")


if __name__ == "__main__":
    main()
