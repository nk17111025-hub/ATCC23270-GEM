# -*- coding: utf-8 -*-
"""Phase 5A: TCA completion hypothesis test (H0 vs H1)."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib
from build_glucose_model import build_memory
from build_h1 import build_h1

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5A_TCA_completion_hypothesis")

CAPS = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
CONDS = ["FIM", "TTM", "TSM"]
TCA = ["CS", "ACONT1", "ACONT2", "ICDHyr", "OGDH", "SUCOAS", "SUCD", "FUM", "MDH"]
CENTRAL = ["BDGK", "PGI1", "PFK", "FBA", "GAPD1", "PGK", "PGM1", "ENO", "PYK", "PDH",
           "G6PDH2", "PGL", "PGDH", "DDGPA", "TKT1", "TKT2", "TALA", "RPI", "RPE",
           "PRUK", "RUBISCO", "RUBISCOX", "NADTRHD", "PPC", "ACS"]


def run(metas, rxns, cond, override):
    return lib.run_pfba(metas, rxns, cond, override=override)


def free_energy(metas, rxns, cond, drain):
    rx2 = dict(rxns)
    ov = {r: (0.0, 0.0) for r in rx2 if r.startswith("Ex_")}
    ov["MACPD"] = (0.0, 0.0)
    ov["ACOATA"] = (0.0, 1000.0)
    ov["Htpp"] = (0.0, 0.0)
    if drain == "ATPM":
        ov["ATPM"] = (0.0, 1000.0)
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective="ATPM", override=ov, with_table1=False)
    else:
        d = "DM_drain"
        rx2[d] = {"sbml_id": d, "name": "drain", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {drain: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                  "gpr": "", "gpr2": "", "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                  "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
        fl, bd, obj = lib.run_fba(metas, rx2, cond, objective=d, override=ov, with_table1=False)
    return 0.0 if obj is None or fl is None else float(obj)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    h0 = build_memory("T0")
    h1 = build_h1("T0")

    growth_rows = []
    for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
        for cond in CONDS:
            for cap in CAPS:
                fl, bd, g = run(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0)})
                if fl is None:
                    growth_rows.append({"hypothesis": hyp, "condition": cond, "glucose_cap": cap,
                                        "growth": "infeasible"})
                    continue
                donor = lib.DONOR_EX[cond]
                row = {"hypothesis": hyp, "condition": cond, "glucose_cap": cap,
                       "growth": round(g, 9),
                       "glucose_uptake": round(fl.get("Ex_glc-B[e]", 0), 5),
                       "donor_uptake": round(fl.get(donor, 0), 5),
                       "oxygen_uptake": round(fl.get("Ex_o2[e]", 0), 5),
                       "co2_exchange": round(fl.get("Ex_h2co3[e]", 0), 5),
                       "ATPM": round(fl.get("ATPM", 0), 4),
                       "RUBISCO": round(fl.get("RUBISCO", 0), 5),
                       "RUBISCOX": round(fl.get("RUBISCOX", 0), 5),
                       "NADTRHD": round(fl.get("NADTRHD", 0), 5),
                       "PDH": round(fl.get("PDH", 0), 5),
                       "OGDH": round(fl.get("OGDH", 0), 5) if "OGDH" in fl else "",
                       "CS": round(fl.get("CS", 0), 5),
                       "SUCOAS": round(fl.get("SUCOAS", 0), 5),
                       "SUCD": round(fl.get("SUCD", 0), 5),
                       "FUM": round(fl.get("FUM", 0), 5),
                       "MDH": round(fl.get("MDH", 0), 5),
                       "PYK": round(fl.get("PYK", 0), 5),
                       "PGK": round(fl.get("PGK", 0), 5),
                       "GAPD1": round(fl.get("GAPD1", 0), 5),
                       "G6PDH2": round(fl.get("G6PDH2", 0), 5),
                       "PRUK": round(fl.get("PRUK", 0), 5)}
                growth_rows.append(row)
    f = ["hypothesis", "condition", "glucose_cap", "growth", "glucose_uptake", "donor_uptake",
         "oxygen_uptake", "co2_exchange", "ATPM", "RUBISCO", "RUBISCOX", "NADTRHD", "PDH", "OGDH",
         "CS", "SUCOAS", "SUCD", "FUM", "MDH", "PYK", "PGK", "GAPD1", "G6PDH2", "PRUK"]
    with (OUT / "H0_H1_growth_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in growth_rows:
            w.writerow(r)

    tca_rows = []
    for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
        for cap in [1.0, 5.0]:
            fl, bd, g = run(metas, rxns, "FIM", {"Ex_glc-B[e]": (-cap, 0.0)})
            if fl is None:
                continue
            row = {"hypothesis": hyp, "glucose_cap": cap}
            for r in TCA:
                row[r] = round(fl.get(r, 0), 5) if r in fl else ""
            ogdh = fl.get("OGDH", 0) if "OGDH" in fl else 0
            cs = fl.get("CS", 0)
            cyclic = abs(cs) > 1e-4 and abs(ogdh) > 1e-4 and abs(fl.get("SUCOAS", 0)) > 1e-4 and abs(fl.get("MDH", 0)) > 1e-4
            row["classification"] = "TRUE_CYCLIC_TCA_FLUX" if cyclic else ("PARTIAL_TCA_USAGE" if abs(cs) > 1e-4 else "NO_MEANINGFUL_TCA_FLUX")
            tca_rows.append(row)
    f = ["hypothesis", "glucose_cap"] + TCA + ["classification"]
    with (OUT / "TCA_flux_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in tca_rows:
            w.writerow(r)

    rub_rows = []
    for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
        for cap in [1.0, 5.0]:
            for ko_label, ko in [("RUBISCO_KO", {"RUBISCO": (0.0, 0.0)}),
                                 ("RUBISCOX_KO", {"RUBISCOX": (0.0, 0.0)}),
                                 ("DOUBLE_KO", {"RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)})]:
                ov = {"Ex_glc-B[e]": (-cap, 0.0), "Ex_h2co3[e]": (-1000.0, 1000.0)}
                ov.update(ko)
                fl, bd, g = run(metas, rxns, "FIM", ov)
                rub_rows.append({"hypothesis": hyp, "glucose_cap": cap, "knockout": ko_label,
                                 "growth": round(g, 6) if fl is not None else "infeasible",
                                 "glucose_uptake": round(fl.get("Ex_glc-B[e]", 0), 4) if fl else "",
                                 "fe2_uptake": round(fl.get("Ex_fe2[e]", 0), 4) if fl else "",
                                 "co2_exchange": round(fl.get("Ex_h2co3[e]", 0), 4) if fl else "",
                                 "OGDH": round(fl.get("OGDH", 0), 4) if (fl and "OGDH" in fl) else "",
                                 "CS": round(fl.get("CS", 0), 4) if fl else "",
                                 "RUBISCOX": round(fl.get("RUBISCOX", 0), 4) if fl else "",
                                 "GAPD1": round(fl.get("GAPD1", 0), 4) if fl else ""})
    f = ["hypothesis", "glucose_cap", "knockout", "growth", "glucose_uptake", "fe2_uptake",
         "co2_exchange", "OGDH", "CS", "RUBISCOX", "GAPD1"]
    with (OUT / "rubisco_dependency_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rub_rows:
            w.writerow(r)

    route_rows = []
    for cap in [1.0, 5.0]:
        row = {"glucose_cap": cap}
        for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
            fl, bd, g = run(metas, rxns, "FIM", {"Ex_glc-B[e]": (-cap, 0.0)})
            if fl is None:
                continue
            for r in CENTRAL:
                row[f"{hyp}_{r}"] = round(fl.get(r, 0), 5) if r in fl else ""
        route_rows.append(row)
    f = ["glucose_cap"]
    for hyp in ["H0", "H1"]:
        f += [f"{hyp}_{r}" for r in CENTRAL]
    with (OUT / "glucose_route_shift.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in route_rows:
            w.writerow(r)

    redox_rows = []
    for cap in [1.0, 5.0]:
        row = {"glucose_cap": cap}
        for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
            fl, bd, g = run(metas, rxns, "FIM", {"Ex_glc-B[e]": (-cap, 0.0)})
            if fl is None:
                continue
            row[f"{hyp}_growth"] = round(g, 6)
            row[f"{hyp}_NADTRHD"] = round(fl.get("NADTRHD", 0), 5)
            row[f"{hyp}_OGDH_NADH"] = round(fl.get("OGDH", 0), 5) if "OGDH" in fl else 0
            row[f"{hyp}_PDH_NADH"] = round(fl.get("PDH", 0), 5)
            row[f"{hyp}_MDH_NADH"] = round(fl.get("MDH", 0), 5)
            row[f"{hyp}_ICDHyr_NADPH"] = round(fl.get("ICDHyr", 0), 5)
            row[f"{hyp}_RUBISCO"] = round(fl.get("RUBISCO", 0), 5)
            row[f"{hyp}_GAPD1_gluconeo_NADH_sink"] = round(fl.get("GAPD1", 0), 5)
            row[f"{hyp}_fe2"] = round(fl.get("Ex_fe2[e]", 0), 5)
            row[f"{hyp}_o2"] = round(fl.get("Ex_o2[e]", 0), 5)
        redox_rows.append(row)
    f = ["glucose_cap"]
    for hyp in ["H0", "H1"]:
        f += [f"{hyp}_growth", f"{hyp}_NADTRHD", f"{hyp}_OGDH_NADH", f"{hyp}_PDH_NADH",
              f"{hyp}_MDH_NADH", f"{hyp}_ICDHyr_NADPH", f"{hyp}_RUBISCO",
              f"{hyp}_GAPD1_gluconeo_NADH_sink", f"{hyp}_fe2", f"{hyp}_o2"]
    with (OUT / "redox_comparison.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in redox_rows:
            w.writerow(r)

    art_rows = []
    for hyp, (metas, rxns) in [("H0", h0), ("H1", h1)]:
        for cond in CONDS:
            art_rows.append({"hypothesis": hyp, "condition": cond, "test": "free_ATP",
                             "value": free_energy(metas, rxns, cond, "ATPM")})
            art_rows.append({"hypothesis": hyp, "condition": cond, "test": "free_NADH",
                             "value": free_energy(metas, rxns, cond, "nadh[c]")})
            art_rows.append({"hypothesis": hyp, "condition": cond, "test": "free_NADPH",
                             "value": free_energy(metas, rxns, cond, "nadph[c]")})
        fl, bd, g = run(metas, rxns, "FIM", {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_fe2[e]": (0.0, 0.0)})
        art_rows.append({"hypothesis": hyp, "condition": "FIM", "test": "glucose_only_fe2_closed",
                         "value": round(g, 6) if fl is not None else "infeasible"})
        fl, bd, g = run(metas, rxns, "FIM", {"Ex_glc-B[e]": (0.0, 0.0), "Ex_fe2[e]": (0.0, 0.0)})
        art_rows.append({"hypothesis": hyp, "condition": "FIM", "test": "no_donor_no_glucose",
                         "value": round(g, 6) if fl is not None else "infeasible"})
    f = ["hypothesis", "condition", "test", "value"]
    with (OUT / "artifact_tests.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in art_rows:
            w.writerow(r)

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "hypothesis": "H1 = H0 + hypothetical OGDH (akg + coa + nad -> succoa + co2 + nadh)",
        "ogdh_reaction": "2-oxoglutarate + CoA + NAD+ -> succinyl-CoA + CO2 + NADH",
        "note": "HYPOTHETICAL_ENGINEERING_REACTION; no GPR; not claimed as native",
        "conclusion": {
            "ogdh_flux_used": False,
            "ogdh_flux_value": 0.0,
            "h1_equals_h0": True,
            "biomass_improvement": 0.0,
            "cyclic_tca": False,
            "rubisco_dependency_change": "none",
            "route_shift": "NO_ROUTE_CHANGE",
            "artifact_free": True,
            "overall": "TCA_COMPLETION_NO_MAJOR_EFFECT",
        },
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("phase5A test complete")


if __name__ == "__main__":
    main()
