# -*- coding: utf-8 -*-
"""Phase 5B-2A: reference states + target phenotype feasibility (minimize CO2 at target)."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib
from build_glucose_model import build_memory
from must_lib import lp

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B2_AB_target_and_must")
CAPS = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
CONDS = ["FIM", "TTM", "TSM"]
FLUX = ["BDGK", "PGI1", "PFK", "FBA", "GAPD1", "GAPD2", "PGK", "PYK", "PDH", "G6PDH2",
        "TKT1", "TKT2", "TALA", "PRUK", "RUBISCO", "RUBISCOX", "NADTRHD", "ATPS5rpp",
        "CS", "PPC", "CYTAA31", "CYTA2", "CYTBC1"]

TARGETS = {
    "A": dict(donor_min=0.90, o2_max=1.20, co2_max=0.75, biomass_mult=1.0),
    "B": dict(donor_min=0.70, o2_max=1.20, co2_max=0.50, biomass_mult=1.0),
    "C": dict(donor_min=0.50, o2_max=None, co2_max=0.0, biomass_mult=1.0),
    "D10": dict(donor_fixed=True, o2_fixed=True, co2_free=True, biomass_mult=1.10),
    "D25": dict(donor_fixed=True, o2_fixed=True, co2_free=True, biomass_mult=1.25),
    "D50": dict(donor_fixed=True, o2_fixed=True, co2_free=True, biomass_mult=1.50),
}


def ref_state(metas, rxns, cond):
    x = lp(metas, rxns, cond, {lib.BIO_EX: -1.0},
           {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    if x is None:
        return None
    rxn_idx = {r: j for j, r in enumerate(rxns)}
    v = lambda r: x[rxn_idx[r]]
    d = {"biomass": v(lib.BIO_EX), "donor": v(lib.DONOR_EX[cond]), "o2": v("Ex_o2[e]"),
         "co2": v("Ex_h2co3[e]"), "glucose": v("Ex_glc-B[e]")}
    for r in FLUX:
        d[r] = v(r) if r in rxn_idx else ""
    return d


def target_lp(metas, rxns, cond, target, glc_cap, ref):
    donor = lib.DONOR_EX[cond]
    t = TARGETS[target]
    ov = {"Ex_glc-B[e]": (-glc_cap, -1e-6)}
    if t.get("co2_free"):
        ov["Ex_h2co3[e]"] = (-1000.0, 0.0)
    else:
        ov["Ex_h2co3[e]"] = (-t["co2_max"] * abs(ref["co2"]), 0.0)
    if t.get("donor_fixed"):
        ov[donor] = (ref["donor"], ref["donor"])
    else:
        ov[donor] = (-1000.0, -t["donor_min"] * abs(ref["donor"]))
    if t.get("o2_fixed"):
        ov["Ex_o2[e]"] = (ref["o2"], ref["o2"])
    elif t.get("o2_max") is not None:
        ov["Ex_o2[e]"] = (-t["o2_max"] * abs(ref["o2"]), 1000.0)
    biomass_min = t["biomass_mult"] * ref["biomass"]
    x = lp(metas, rxns, cond, {"Ex_h2co3[e]": -1.0}, ov, biomass_min=biomass_min)
    return x, biomass_min


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref_rows = []
    ref_data = {}
    for transport in ["T0", "TH"]:
        metas, rxns = build_memory(transport)
        for cond in CONDS:
            r = ref_state(metas, rxns, cond)
            ref_data[(transport, cond)] = r
            row = {"transport": transport, "condition": cond}
            row.update({k: (round(v, 7) if isinstance(v, float) else v) for k, v in r.items()})
            ref_rows.append(row)
    f = ["transport", "condition", "biomass", "donor", "o2", "co2", "glucose"] + FLUX
    with (OUT / "reference_states.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in ref_rows:
            w.writerow(r)

    feas_rows = []
    flux_rows = []
    for transport in ["T0", "TH"]:
        metas, rxns = build_memory(transport)
        for cond in CONDS:
            ref = ref_data[(transport, cond)]
            for target, t in TARGETS.items():
                for cap in CAPS:
                    x, bmin = target_lp(metas, rxns, cond, target, cap, ref)
                    if x is None:
                        feas_rows.append({"transport": transport, "condition": cond, "target": target,
                                          "glucose_cap": cap, "feasible": "no", "biomass": "",
                                          "glucose_uptake": "", "donor_uptake": "", "o2_uptake": "",
                                          "co2_uptake": "", "co2_saving_pct": ""})
                        continue
                    rxn_idx = {r: j for j, r in enumerate(rxns)}
                    v = lambda r: x[rxn_idx[r]]
                    glc = v("Ex_glc-B[e]")
                    co2 = v("Ex_h2co3[e]")
                    saving = (abs(ref["co2"]) - abs(co2)) / abs(ref["co2"]) * 100 if abs(ref["co2"]) > 0 else 0
                    feasible = abs(glc) > 1e-6 and abs(co2) <= t.get("co2_max", 1.0) * abs(ref["co2"]) + 1e-6
                    feas_rows.append({"transport": transport, "condition": cond, "target": target,
                                      "glucose_cap": cap, "feasible": "yes" if feasible else "no",
                                      "biomass": round(v(lib.BIO_EX), 7),
                                      "glucose_uptake": round(glc, 6),
                                      "donor_uptake": round(v(lib.DONOR_EX[cond]), 5),
                                      "o2_uptake": round(v("Ex_o2[e]"), 5),
                                      "co2_uptake": round(co2, 6),
                                      "co2_saving_pct": round(saving, 4)})
                    if feasible:
                        row = {"transport": transport, "condition": cond, "target": target, "glucose_cap": cap}
                        for r in FLUX:
                            row[r] = round(v(r), 6) if r in rxn_idx else ""
                        flux_rows.append(row)
    f = ["transport", "condition", "target", "glucose_cap", "feasible", "biomass", "glucose_uptake",
         "donor_uptake", "o2_uptake", "co2_uptake", "co2_saving_pct"]
    with (OUT / "target_feasibility_matrix.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in feas_rows:
            w.writerow(r)
    f = ["transport", "condition", "target", "glucose_cap"] + FLUX
    with (OUT / "target_flux_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in flux_rows:
            w.writerow(r)

    print("Phase 5B-2A complete")
    feas = {}
    for r in feas_rows:
        if r["feasible"] == "yes":
            feas.setdefault((r["condition"], r["target"]), r["glucose_cap"])
    for k in sorted(feas):
        print(k, "first feasible cap", feas[k])


if __name__ == "__main__":
    main()
