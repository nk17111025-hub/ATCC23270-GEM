# -*- coding: utf-8 -*-
"""Parts D-H: baseline regression, glucose sensitivity, route analysis, artifact tests."""
import csv
import json
import datetime as dt
from pathlib import Path

import phase4a_lib as lib
from build_glucose_model import build_memory, sha256

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4_1G_glucose_corrected_model")
FROZEN = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\freeze\minimal_trustworthy_baseline_v1.xml")

CAPS = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0]

FLUX_COLS = ["BDGK", "PFK", "G6PDH2", "PGL", "PGDH", "DDGPA", "RPI", "RPE", "TKT1", "TKT2",
             "TALA", "PRUK", "FBA", "PGI1", "FBP", "GAPD1", "PGK", "PGM1", "ENO", "PYK",
             "PDH", "RUBISCO", "PPC", "CS", "NADTRHD"]
WATCH = ["BDGK", "ACCOAC", "NADTRHD", "MALS", "GLXCL", "GLYCK", "GLXR", "MDH", "FUM"]


def run_case(metas, rxns, cond, override):
    fluxes, bounds, opt = lib.run_pfba(metas, rxns, cond, override=override)
    return fluxes, opt


def base_row(metas, rxns, cond, override):
    fluxes, opt = run_case(metas, rxns, cond, override)
    if fluxes is None:
        return None
    donor = lib.DONOR_EX[cond]
    row = {
        "condition": cond,
        "growth": round(opt, 9),
        "glucose_uptake": round(fluxes.get("Ex_glc-B[e]", 0.0), 6),
        "donor_uptake": round(fluxes.get(donor, 0.0), 6),
        "oxygen_uptake": round(fluxes.get("Ex_o2[e]", 0.0), 6),
        "co2_hco3_uptake": round(fluxes.get("Ex_h2co3[e]", 0.0), 6),
        "ATPM": round(fluxes.get("ATPM", 0.0), 6),
    }
    for r in FLUX_COLS:
        row[r] = round(fluxes.get(r, 0.0), 6) if r in fluxes else ""
    return row


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
    frozen_sha = sha256(FROZEN)
    results = {}

    baseline_rows = []
    sensitivity_rows = []
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cond in ["FIM", "TTM", "TSM"]:
            # baseline (glucose closed)
            br = base_row(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0)})
            if br:
                br["transport"] = t
                br["case"] = "glucose_closed"
                baseline_rows.append(br)
            # HCO3E KO (glucose closed)
            fl, opt = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0), "HCO3E": (0.0, 0.0)})
            hco3e = "infeasible" if fl is None else round(opt, 9)
            baseline_rows.append({"transport": t, "condition": cond, "case": "HCO3E_KO_glucose_closed",
                                  "growth": hco3e, "glucose_uptake": "", "donor_uptake": "", "oxygen_uptake": "",
                                  "co2_hco3_uptake": "", "ATPM": ""})
            # free energy
            baseline_rows.append({"transport": t, "condition": cond, "case": "free_energy_check",
                                  "growth": "", "glucose_uptake": free_energy(metas, rxns, cond, "ATPM"),
                                  "donor_uptake": free_energy(metas, rxns, cond, "nadh[c]"),
                                  "oxygen_uptake": free_energy(metas, rxns, cond, "nadph[c]"),
                                  "co2_hco3_uptake": "", "ATPM": ""})
            # glucose sensitivity
            for cap in CAPS:
                r = base_row(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0)})
                if r:
                    r["transport"] = t
                    r["case"] = f"glucose_cap_{cap}"
                    r["glucose_uptake_cap"] = cap
                    sensitivity_rows.append(r)

    # write baseline regression
    bf = ["transport", "condition", "case", "growth", "glucose_uptake", "donor_uptake",
          "oxygen_uptake", "co2_hco3_uptake", "ATPM"] + FLUX_COLS
    with (OUT / "glucose_baseline_regression.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=bf, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in baseline_rows:
            w.writerow(r)

    sf = ["transport", "condition", "case", "glucose_uptake_cap", "growth", "glucose_uptake",
          "donor_uptake", "oxygen_uptake", "co2_hco3_uptake", "ATPM"] + FLUX_COLS
    with (OUT / "glucose_uptake_sensitivity.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sf, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sensitivity_rows:
            w.writerow(r)

    results["baseline_rows"] = len(baseline_rows)
    results["sensitivity_rows"] = len(sensitivity_rows)

    # route analysis + artifact tests (FIM, caps 0.1/1/5)
    route_rows = []
    artifact_rows = []
    watch_rows = []
    for t in ["T0", "TH"]:
        metas, rxns = build_memory(t)
        for cond in ["FIM", "TTM", "TSM"]:
            for cap in [0.1, 1.0, 5.0]:
                fl, opt_pfba = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0)})
                r = base_row(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0)})
                if not r:
                    continue
                g6p = r.get("BDGK", 0) > 1e-6
                oxppp = r.get("G6PDH2", 0) < -1e-6
                lower_emp = (r.get("GAPD1", 0) < -1e-6) and (r.get("PGK", 0) < -1e-6)
                to_pyr = r.get("PYK", 0) > 1e-6
                to_accoa = r.get("PDH", 0) > 1e-6
                co2_decr = r["co2_hco3_uptake"] > -2.0 + 1e-6
                biomass_inc = r["growth"] > 0.052076387 + 1e-6
                fe_s = r["donor_uptake"] < -1e-6
                rubisco = abs(r.get("RUBISCO", 0)) > 1e-6
                nadtrhd_large = abs(r.get("NADTRHD", 0)) > 5.0
                route_rows.append({
                    "transport": t, "condition": cond, "cap": cap,
                    "Q1_glucose_to_G6P": g6p, "Q2_oxidative_PPP": oxppp, "Q3_lower_EMP_direct": lower_emp,
                    "Q4_to_pyruvate": to_pyr, "Q5_to_acetylCoA": to_accoa,
                    "Q6_inorganic_C_demand_decreases": co2_decr, "Q7_biomass_increases": biomass_inc,
                    "Q8_FeS_active": fe_s, "Q9_RUBISCO_active": rubisco,
                    "Q10_NADTRHD_large": nadtrhd_large,
                    "growth": r["growth"], "glucose_uptake": r["glucose_uptake"],
                    "G6PDH2": r.get("G6PDH2", ""), "GAPD1": r.get("GAPD1", ""), "PGK": r.get("PGK", ""),
                    "RUBISCO": r.get("RUBISCO", ""), "NADTRHD": r.get("NADTRHD", ""),
                })
                # artifact tests
                fl_co2off, opt_co2off = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
                fl_rubko, opt_rubko = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0), "RUBISCO": (0.0, 0.0)})
                fl_hco3ko, opt_hco3ko = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (-cap, 0.0), "HCO3E": (0.0, 0.0)})
                fl_glcoff, opt_glcoff = run_case(metas, rxns, cond, {"Ex_glc-B[e]": (0.0, 0.0)})
                # max internal cycle: max |flux| among non-exchange reactions vs glucose uptake
                internal_max = max((abs(v) for k, v in fl.items() if not k.startswith("Ex_")), default=0.0) if fl else 0.0
                artifact_rows.append({
                    "transport": t, "condition": cond, "cap": cap,
                    "pFBA_growth": r["growth"],
                    "glucose_closed_growth": round(opt_glcoff, 9) if fl_glcoff else "infeasible",
                    "inorganicC_closed_growth": round(opt_co2off, 9) if fl_co2off else "infeasible",
                    "RUBISCO_KO_growth": round(opt_rubko, 9) if fl_rubko else "infeasible",
                    "HCO3E_KO_growth": round(opt_hco3ko, 9) if fl_hco3ko else "infeasible",
                    "free_ATP": free_energy(metas, rxns, cond, "ATPM"),
                    "free_NADH": free_energy(metas, rxns, cond, "nadh[c]"),
                    "free_NADPH": free_energy(metas, rxns, cond, "nadph[c]"),
                    "max_internal_abs_flux": round(internal_max, 4),
                    "glucose_uptake": r["glucose_uptake"],
                    "notes": "",
                })
                if cond == "FIM":
                    wf = {w: round(fl.get(w, 0.0), 6) for w in WATCH} if fl else {}
                    wf.update({"transport": t, "cap": cap})
                    watch_rows.append(wf)

    rf = ["transport", "condition", "cap", "Q1_glucose_to_G6P", "Q2_oxidative_PPP", "Q3_lower_EMP_direct",
          "Q4_to_pyruvate", "Q5_to_acetylCoA", "Q6_inorganic_C_demand_decreases", "Q7_biomass_increases",
          "Q8_FeS_active", "Q9_RUBISCO_active", "Q10_NADTRHD_large", "growth", "glucose_uptake",
          "G6PDH2", "GAPD1", "PGK", "RUBISCO", "NADTRHD"]
    with (OUT / "glucose_flux_route_analysis.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rf, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in route_rows:
            w.writerow(r)

    af = ["transport", "condition", "cap", "pFBA_growth", "glucose_closed_growth", "inorganicC_closed_growth",
          "RUBISCO_KO_growth", "HCO3E_KO_growth", "free_ATP", "free_NADH", "free_NADPH",
          "max_internal_abs_flux", "glucose_uptake", "notes"]
    with (OUT / "glucose_artifact_tests.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=af, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in artifact_rows:
            w.writerow(r)

    wf = ["transport", "cap"] + WATCH
    with (OUT / "glucose_watch_reaction_fluxes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=wf, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in watch_rows:
            w.writerow(r)

    results["route_rows"] = len(route_rows)
    results["artifact_rows"] = len(artifact_rows)
    results["watch_rows"] = len(watch_rows)

    (OUT / "_run_summary.json").write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("baseline", len(baseline_rows), "sensitivity", len(sensitivity_rows),
          "route", len(route_rows), "artifact", len(artifact_rows), "watch", len(watch_rows))


if __name__ == "__main__":
    main()
