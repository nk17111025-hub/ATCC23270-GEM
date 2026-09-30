# -*- coding: utf-8 -*-
"""Phase 4.1G-4A: glucose route diagnosis (no model changes)."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4_1G_glucose_route_diagnosis")

EMP_STEPS = [
    ("glucose -> G6P", "BDGK", "forward", 1.0),
    ("G6P -> F6P", "PGI1", "forward", 0.0),
    ("F6P -> FDP", "PFK", "forward", 1.0),
    ("FDP -> GAP+DHAP", "FBA", "forward", 0.0),
    ("GAP -> 1,3-BPG", "GAPD1", "reverse", -1.0),
    ("1,3-BPG -> 3PG", "PGK", "reverse", -1.0),
    ("3PG -> 2PG", "PGM1", "forward", 0.0),
    ("2PG -> PEP", "ENO", "forward", 0.0),
    ("PEP -> pyruvate", "PYK", "forward", 1.0),
    ("pyruvate -> acetyl-CoA", "PDH", "forward", 1.0),
]


def run_pfba(metas, rxns, override):
    return lib.run_pfba(metas, rxns, "FIM", override=override)


def task1_flux_trace():
    rows = []
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cap in [0.1, 1.0, 5.0]:
            fl, bd, opt = run_pfba(metas, rxns, {"Ex_glc-B[e]": (-cap, 0.0)})
            if fl is None:
                continue
            for r, v in sorted(fl.items(), key=lambda kv: -abs(kv[1])):
                if abs(v) < 1e-3:
                    continue
                rx = rxns[r]
                rows.append({"transport": t, "cap": cap, "reaction_id": r, "name": rx["name"],
                             "equation": lib.equation_str(rx["stoich"]), "flux": round(v, 6),
                             "subsystem": rx["subsystem"], "GPR": rx["gpr"], "confidence": rx["confidence"]})
    f = ["transport", "cap", "reaction_id", "name", "equation", "flux", "subsystem", "GPR", "confidence"]
    with (OUT / "glucose_full_flux_trace.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def task2_lower_emp():
    metas, rxns = build_memory("T0")
    rows = []
    for label, rid, direction, sign in EMP_STEPS:
        rx = rxns[rid]
        base, bd, opt = run_pfba(metas, rxns, {"Ex_glc-B[e]": (-1.0, 0.0)})
        actual = base[rid] if base else None
        row = {"step": label, "reaction_id": rid, "equation": lib.equation_str(rx["stoich"]),
               "lb": rx["lb"], "ub": rx["ub"], "GPR": rx["gpr"], "confidence": rx["confidence"],
               "actual_flux": round(actual, 6) if actual is not None else ""}
        # feasibility of forcing glycolytic direction at increasing magnitudes
        for eps in [1e-6, 1e-4, 1e-2]:
            ov = {"Ex_glc-B[e]": (-1.0, 0.0)}
            if direction == "forward":
                ov[rid] = (eps, max(rx["ub"], eps))
            else:
                ov[rid] = (min(rx["lb"], -eps), -eps)
            fl, _, g = run_pfba(metas, rxns, ov)
            row[f"force_{eps}"] = (round(g, 6) if fl is not None else "infeasible")
        rows.append(row)
    f = ["step", "reaction_id", "equation", "lb", "ub", "GPR", "confidence", "actual_flux",
         "force_1e-06", "force_0.0001", "force_0.01"]
    with (OUT / "lower_EMP_block_diagnosis.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def task3_dependency():
    metas, rxns = build_memory("T0")
    targets = ["RUBISCO", "RUBISCOX", "NADTRHD", "G6PDH2", "PFK", "GAPD1", "GAPD2", "PGK", "PYK", "PDH"]
    rows = []
    for rid in targets:
        # baseline (glucose closed)
        fl0, _, g0 = run_pfba(metas, rxns, {"Ex_glc-B[e]": (0.0, 0.0), rid: (0.0, 0.0)})
        base_g = round(g0, 6) if fl0 is not None else "infeasible"
        # glucose open
        fl1, _, g1 = run_pfba(metas, rxns, {"Ex_glc-B[e]": (-1.0, 0.0), rid: (0.0, 0.0)})
        glc_g = round(g1, 6) if fl1 is not None else "infeasible"
        row = {"reaction_id": rid, "baseline_growth_KO": base_g, "glucose_growth_KO": glc_g}
        if fl1 is not None:
            row.update({"glucose_uptake": round(fl1.get("Ex_glc-B[e]", 0), 4),
                        "fe2_uptake": round(fl1.get("Ex_fe2[e]", 0), 4),
                        "co2_uptake": round(fl1.get("Ex_h2co3[e]", 0), 4),
                        "oxygen_uptake": round(fl1.get("Ex_o2[e]", 0), 4),
                        "RUBISCO": round(fl1.get("RUBISCO", 0), 4),
                        "RUBISCOX": round(fl1.get("RUBISCOX", 0), 4),
                        "NADTRHD": round(fl1.get("NADTRHD", 0), 4),
                        "GAPD1": round(fl1.get("GAPD1", 0), 4),
                        "PGK": round(fl1.get("PGK", 0), 4),
                        "G6PDH2": round(fl1.get("G6PDH2", 0), 4)})
        else:
            row.update({"glucose_uptake": "", "fe2_uptake": "", "co2_uptake": "", "oxygen_uptake": "",
                        "RUBISCO": "", "RUBISCOX": "", "NADTRHD": "", "GAPD1": "", "PGK": "", "G6PDH2": ""})
        # classify
        if base_g == "infeasible":
            cls = "ESSENTIAL_FOR_BASELINE"
        elif glc_g == "infeasible":
            cls = "ESSENTIAL_FOR_GLUCOSE_BENEFIT"
        elif fl1 is not None and g1 < 0.052076 + 1e-6:
            cls = "ESSENTIAL_FOR_GLUCOSE_BENEFIT"
        else:
            cls = "DISPENSABLE"
        row["classification"] = cls
        rows.append(row)
    f = ["reaction_id", "baseline_growth_KO", "glucose_growth_KO", "classification", "glucose_uptake",
         "fe2_uptake", "co2_uptake", "oxygen_uptake", "RUBISCO", "RUBISCOX", "NADTRHD", "GAPD1", "PGK", "G6PDH2"]
    with (OUT / "glucose_dependency_tests.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def task4_routes():
    metas, rxns = build_memory("T0")
    def has(r): return r in rxns
    def rev(r): return rxns[r]["lb"] < 0 if has(r) else False
    routes = [
        ("Entner-Doudoroff (G6P->6PGL->6PG->KDPG->Pyr+GAP)",
         all(has(x) for x in ["G6PDH2", "PGL", "PGDH", "DDGPA"]),
         "G6PDH2" if not rev("G6PDH2") else "COMPLETE"),
        ("oxidative PPP (G6P->6PG->Ru5P+CO2)",
         all(has(x) for x in ["G6PDH2", "PGL", "PGDH"]),
         "PARTIAL (no 6PG dehydrogenase/decarboxylase; PGDH is ED dehydratase)"),
        ("non-oxidative PPP (F6P/GAP -> R5P/X5P/E4P/S7P)",
         all(has(x) for x in ["TKT1", "TKT2", "TALA", "RPI", "RPE"]),
         "COMPLETE"),
        ("glycogen synthesis (G6P->G1P->glycogen)",
         all(has(x) for x in ["PGMT", "PGMT2", "GLCP"]),
         "PARTIAL"),
        ("glycogen degradation (glycogen->G6P)",
         has("GLCP") and rev("GLCP"),
         "PRESENT"),
        ("gluconate pathway", has("GNT") , "ABSENT"),
        ("pyruvate bypass (acetyl-CoA via acetate)",
         all(has(x) for x in ["ACS", "PTA", "ACACT"]),
         "PARTIAL"),
        ("glyoxylate shunt (isocitrate lyase + malate synthase)",
         has("MALS"), "PARTIAL (MALS only, no ICL)"),
    ]
    rows = []
    for name, present, status in routes:
        s = status if present else ("ABSENT" if not present else status)
        if name.startswith("gluconate"):
            s = "ABSENT"
        rows.append({"route": name, "present": "TRUE" if present else "FALSE",
                     "status": s if present else "ABSENT"})
    f = ["route", "present", "status"]
    with (OUT / "alternative_route_inventory.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def task5_rubisco_mechanism():
    metas, rxns = build_memory("T0")
    rows = []
    # dependency on RuBP carboxylase vs oxygenase
    for label, ov in [
        ("RUBISCO_KO_CO2_forced", {"Ex_glc-B[e]": (-1.0, 0.0), "RUBISCO": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}),
        ("RUBISCO_KO_CO2_free", {"Ex_glc-B[e]": (-1.0, 0.0), "RUBISCO": (0.0, 0.0), "Ex_h2co3[e]": (-1000.0, 1000.0)}),
        ("RUBISCOX_KO_CO2_free", {"Ex_glc-B[e]": (-1.0, 0.0), "RUBISCOX": (0.0, 0.0), "Ex_h2co3[e]": (-1000.0, 1000.0)}),
        ("BOTH_KO_CO2_free", {"Ex_glc-B[e]": (-1.0, 0.0), "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0), "Ex_h2co3[e]": (-1000.0, 1000.0)}),
    ]:
        fl, bd, g = run_pfba(metas, rxns, ov)
        rows.append({"case": label,
                     "growth": round(g, 6) if fl is not None else "infeasible",
                     "glucose_uptake": round(fl.get("Ex_glc-B[e]", 0), 4) if fl else "",
                     "co2_uptake": round(fl.get("Ex_h2co3[e]", 0), 4) if fl else "",
                     "RUBISCO": round(fl.get("RUBISCO", 0), 4) if fl else "",
                     "RUBISCOX": round(fl.get("RUBISCOX", 0), 4) if fl else "",
                     "GAPD1": round(fl.get("GAPD1", 0), 4) if fl else "",
                     "PGK": round(fl.get("PGK", 0), 4) if fl else "",
                     "PRUK": round(fl.get("PRUK", 0), 4) if fl else "",
                     "PGLYCP": round(fl.get("PGLYCP", 0), 4) if fl else "",
                     "GLYCH": round(fl.get("GLYCH", 0), 4) if fl else "",
                     "GLYCL": round(fl.get("GLXCL", 0), 4) if fl else "",
                     "GLYCK": round(fl.get("GLYCK", 0), 4) if fl else "",
                     "NADTRHD": round(fl.get("NADTRHD", 0), 4) if fl else ""})
    f = ["case", "growth", "glucose_uptake", "co2_uptake", "RUBISCO", "RUBISCOX", "GAPD1", "PGK",
         "PRUK", "PGLYCP", "GLYCH", "GLYCL", "GLYCK", "NADTRHD"]
    with (OUT / "rubisco_dependency_mechanism.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    task1_flux_trace()
    task2_lower_emp()
    task3_dependency()
    task4_routes()
    task5_rubisco_mechanism()
    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "conclusions": {
            "dominant_glucose_route": "glucose -> G6P -> F6P -> non-oxidative PPP -> Ru5P -> RuBP (PRUK) -> RUBISCO(carboxylase) or RUBISCOX(oxygenase) -> 3PG -> 2PG -> PEP -> pyruvate -> acetyl-CoA",
            "lower_EMP_used": False,
            "first_bottleneck": "F6P is routed to non-oxidative PPP, not PFK/FDP; GAPD1/PGK never carry glycolytic flux",
            "rubisco_dependence": "model requires the RuBP carboxylase/oxygenase system (RUBISCO or RUBISCOX); closing BOTH makes glucose unusable (growth 0)",
            "co2_forced_artifact": "RUBISCO KO + forced CO2 uptake (-2.0) is infeasible because CO2 has no sink",
        },
    }
    (OUT / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("diagnosis complete")


if __name__ == "__main__":
    main()
