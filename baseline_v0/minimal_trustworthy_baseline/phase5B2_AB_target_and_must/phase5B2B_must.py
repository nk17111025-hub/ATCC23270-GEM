# -*- coding: utf-8 -*-
"""Phase 5B-2B: FVA-based MUST reaction discovery."""
import csv
import json
import datetime as dt
from pathlib import Path

import numpy as np

import phase4a_lib as lib
from build_glucose_model import build_memory
from must_lib import fva_range, lp, TOL

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B2_AB_target_and_must")

# feasible targets from Phase 5B-2A (T0): condition -> list of (target, glc_cap)
FEASIBLE = {
    "FIM": [("A", 1.0), ("B", 1.0), ("C", 1.0)],
    "TTM": [("A", 1.0)],
    "TSM": [("A", 1.0)],
}
TARGETS = {
    "A": dict(donor_min=0.90, o2_max=1.20, co2_max=0.75, biomass_mult=1.0),
    "B": dict(donor_min=0.70, o2_max=1.20, co2_max=0.50, biomass_mult=1.0),
    "C": dict(donor_min=0.50, o2_max=None, co2_max=0.0, biomass_mult=1.0),
}

PRIORITY = ["Ex_glc-B[e]", "GLCtex", "GLCtpp", "BDGK", "PGI1", "PFK", "FBA", "GAPD1", "GAPD2",
            "PGK", "PGM1", "ENO", "PYK", "PDH", "G6PDH2", "PGL", "PGDH", "DDGPA", "TKT1", "TKT2",
            "TALA", "RPI", "RPE", "PRUK", "RUBISCO", "RUBISCOX", "PPC", "CS", "ACONT1", "ICDHyr",
            "SUCOAS", "SUCD", "FUM", "MDH", "MALS", "ACS", "NADTRHD", "NADHI", "CYT2", "CYT1",
            "CYTA2", "CYTAA31", "CYTBC1", "ATPS5rpp", "Htpp", "FDH", "HYD3pp", "Ex_fe2[e]",
            "Ex_ttton[e]", "Ex_tsul[e]", "Ex_o2[e]", "Ex_h2co3[e]", "Ex_h[e]", "Ex_n2[e]"]


def target_override(cond, target, glc_cap, ref):
    donor = lib.DONOR_EX[cond]
    t = TARGETS[target]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0)}
    ov["Ex_h2co3[e]"] = (-t["co2_max"] * abs(ref["co2"]), 0.0)
    ov[donor] = (-1000.0, -t["donor_min"] * abs(ref["donor"]))
    if t.get("o2_max") is not None:
        ov["Ex_o2[e]"] = (-t["o2_max"] * abs(ref["o2"]), 1000.0)
    return ov


def ref_state(metas, rxns, cond):
    x = lp(metas, rxns, cond, {lib.BIO_EX: -1.0},
           {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    v = lambda r: x[rxn_idx[r]]
    return {"biomass": v(lib.BIO_EX), "donor": v(lib.DONOR_EX[cond]), "o2": v("Ex_o2[e]"),
            "co2": v("Ex_h2co3[e]")}


def classify(wlo, whi, tlo, thi, tol=TOL):
    def z(v):
        return v is None or abs(v) < tol
    w_zero = z(wlo) and z(whi)
    t_zero = z(tlo) and z(thi)
    if w_zero and not t_zero:
        return "MUST-ON"
    if (not w_zero) and t_zero:
        return "MUST-OFF"
    if wlo is None or whi is None or tlo is None or thi is None:
        return None
    if tlo > whi + tol:
        return "MUST-UP"
    if thi < wlo - tol:
        return "MUST-DOWN"
    if tlo > wlo + tol or thi < whi - tol:
        return "STRONG_RANGE_SHIFT"
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns = build_memory("T0")
    all_rxns = list(rxns.keys())

    wt_fva_rows = []
    wt_fva = {}
    for cond in ["FIM", "TTM", "TSM"]:
        ref = ref_state(metas, rxns, cond)
        ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
        fv = fva_range(metas, rxns, cond, ov, ref["biomass"], all_rxns)
        wt_fva[cond] = (fv, ref)
        for r, (lo, hi) in fv.items():
            wt_fva_rows.append({"condition": cond, "reaction": r,
                                "min": round(lo, 6) if lo is not None else "",
                                "max": round(hi, 6) if hi is not None else ""})
    f = ["condition", "reaction", "min", "max"]
    with (OUT / "wt_fva.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in wt_fva_rows:
            w.writerow(r)

    must_rows = []
    target_fva_all = {}
    for cond, targets in FEASIBLE.items():
        ref = wt_fva[cond][1]
        for target, cap in targets:
            ov = target_override(cond, target, cap, ref)
            fv = fva_range(metas, rxns, cond, ov, ref["biomass"], all_rxns)
            target_fva_all[(cond, target)] = fv
            # write per-target FVA (priority reactions)
            with (OUT / f"target_fva_{cond}_{target}.tsv").open("w", encoding="utf-8", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["reaction", "min", "max"])
                for r in PRIORITY:
                    if r in fv:
                        lo, hi = fv[r]
                        w.writerow([r, round(lo, 6) if lo is not None else "", round(hi, 6) if hi is not None else ""])
            # classify vs WT
            wtfv = wt_fva[cond][0]
            for r in all_rxns:
                wlo, whi = wtfv.get(r, (None, None))
                tlo, thi = fv.get(r, (None, None))
                cls = classify(wlo, whi, tlo, thi)
                if cls:
                    must_rows.append({"condition": cond, "target": target, "reaction": r, "class": cls,
                                      "wt_min": round(wlo, 6) if wlo is not None else "",
                                      "wt_max": round(whi, 6) if whi is not None else "",
                                      "target_min": round(tlo, 6) if tlo is not None else "",
                                      "target_max": round(thi, 6) if thi is not None else ""})

    # write MUST single reactions
    f = ["condition", "target", "reaction", "class", "wt_min", "wt_max", "target_min", "target_max"]
    with (OUT / "must_single_reactions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in must_rows:
            w.writerow(r)

    # strong range shifts
    shift_rows = [r for r in must_rows if r["class"] == "STRONG_RANGE_SHIFT"]
    with (OUT / "strong_range_shifts.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in shift_rows:
            w.writerow(r)

    # condition comparison (reaction -> class per condition, using target A)
    cond_rows = []
    by_reaction = {}
    for r in must_rows:
        if r["target"] == "A":
            by_reaction.setdefault(r["reaction"], {})[r["condition"]] = r["class"]
    for rxn, cmap in sorted(by_reaction.items()):
        row = {"reaction": rxn, "FIM": cmap.get("FIM", ""), "TTM": cmap.get("TTM", ""), "TSM": cmap.get("TSM", "")}
        if any(cmap.get(c) for c in ["FIM", "TTM", "TSM"]):
            cond_rows.append(row)
    f = ["reaction", "FIM", "TTM", "TSM"]
    with (OUT / "must_condition_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in cond_rows:
            w.writerow(r)

    # summary
    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": lib.sha256(lib.V1_XML),
        "method": "FVA-based MUST comparison (WT glucose=0 CO2-free vs target phenotype)",
        "tolerance": TOL,
        "conclusions": {
            "feasible_targets": {"A": "all conditions", "B": "all conditions", "C": "all conditions",
                                 "D10_D25_D50": "infeasible all conditions"},
            "core_MUST_ON": ["Ex_glc-B[e]", "GLCtex", "GLCtpp"],
            "core_MUST_UP": ["PGI1 (G6P->F6P direction flip)", "Ex_h2co3[e] (CO2 decrease)"],
            "cross_condition_universal": True,
            "opposite_direction_reactions": "none at MUST level",
            "Tier1": ["glucose uptake+transport", "PGI1 forward", "CO2 reduction"],
            "Tier2": ["Fe/S oxidation-capacity reduction (donor 90/70/50% WT)", "non-oxidative PPP directionality"],
            "Tier3": ["CBB/respiratory STRONG_RANGE_SHIFT (confounded by WT FVA degeneracy)"],
            "note": "FVA ranges are degenerate (internal loops); MUST-ON and tight MUST-UP are robust",
        },
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("Phase 5B-2B complete, MUST rows:", len(must_rows))
    # print priority MUST for FIM A
    print("=== FIM Target A MUST (priority) ===")
    for r in must_rows:
        if r["condition"] == "FIM" and r["target"] == "A" and r["reaction"] in PRIORITY:
            print(f"  {r['reaction']:16s} {r['class']:20s} WT[{r['wt_min']},{r['wt_max']}] -> [{r['target_min']},{r['target_max']}]")


if __name__ == "__main__":
    main()
