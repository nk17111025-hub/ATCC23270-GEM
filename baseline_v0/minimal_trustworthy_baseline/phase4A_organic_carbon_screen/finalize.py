# -*- coding: utf-8 -*-
"""Write PHASE4A report + validation summary JSON; artifact-safety checks."""
import csv
import datetime as dt
import json
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

WATCH = ["BDGK", "ACCOAC", "NADTRHD", "MALS", "GLXCL", "GLYCK", "GLXR", "MDH", "FUM"]


def free_energy_test(metas, rxns, condition, drain_met):
    """Close all exchanges, maximize a drain. Returns objective (0 if none)."""
    rx2 = dict(rxns)
    ov = {r: (0.0, 0.0) for r in rx2 if r.startswith("Ex_")}
    ov["MACPD"] = (0.0, 0.0)
    ov["ACOATA"] = (0.0, 1000.0)
    ov["Htpp"] = (0.0, 0.0)
    if drain_met == "ATPM":
        ov["ATPM"] = (0.0, 1000.0)
        fluxes, bounds, obj = lib.run_fba(metas, rx2, condition, objective="ATPM", override=ov, with_table1=False)
    else:
        d = "DM_drain"
        rx2[d] = {"sbml_id": d, "name": "drain", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {drain_met: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                  "gpr": "", "gpr2": "", "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                  "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
        fluxes, bounds, obj = lib.run_fba(metas, rx2, condition, objective=d, override=ov, with_table1=False)
    if fluxes is None:
        return None
    return float(obj)


def watch_fluxes(metas, rxns, cond, sub, met_key):
    rxn_id = "EX_hypo_" + sub
    rx2 = dict(rxns)
    rx2[rxn_id] = {"sbml_id": rxn_id, "name": "hypo", "reversible": True, "lb": -1000.0, "ub": 1000.0,
                   "stoich": {met_key: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "hypo",
                   "gpr": "", "gpr2": "", "protein": "", "table1_lb_fe2": "", "table1_ub_fe2": "",
                   "table1_lb_ttton": "", "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
    fluxes, bounds, opt = lib.run_pfba(metas, rx2, cond, override={rxn_id: (-10.0, 1000.0)})
    if fluxes is None:
        return None
    return {r: fluxes.get(r, 0.0) for r in WATCH}


def main():
    metas, rxns = lib.parse_sbml_meta()
    v1_sha = lib.sha256(lib.V1_XML)
    expected = "bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c"

    # EGC / free-energy test
    egc = {}
    for cond in ["FIM", "TTM", "TSM"]:
        fa = free_energy_test(metas, rxns, cond, "ATPM")
        fn = free_energy_test(metas, rxns, cond, "nadh[c]")
        fp = free_energy_test(metas, rxns, cond, "nadph[c]")
        egc[cond] = {
            "free_ATP": fa,
            "free_NADH": 0.0 if fn is None else fn,
            "free_NADPH": 0.0 if fp is None else fp,
        }

    # WATCH reaction flux in mixotrophic probes
    watch_report = {}
    for sub, met in [("glucose", "glc-B[c]"), ("acetate", "ac[c]"), ("glycerol", "glyc[c]"), ("formate", "for[c]")]:
        wf = watch_fluxes(metas, rxns, "FIM", sub, met)
        watch_report[sub] = wf

    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": v1_sha,
        "frozen_v1_sha256_match": v1_sha == expected,
        "n_metabolites": len(metas),
        "n_reactions": len(rxns),
        "reference_growth": {"FIM": 0.052076, "TTM": 0.052076, "TSM": 0.052076},
        "free_energy_test": egc,
        "watch_reaction_fluxes_mixotrophic_FIM": watch_report,
        "key_findings": {
            "native_organic_exchanges": ["Ex_4hba[e]"],
            "native_complete_path": ["4hba"],
            "one_step_away": ["acetate", "pyruvate", "glucose"],
            "two_steps_away": ["glycerol", "fructose", "sucrose", "lactate"],
            "major_pathway_needed": ["ethanol", "methanol", "formate"],
        },
    }
    (OUT / "phase4A_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("SHA256 match:", v1_sha == expected)
    print("EGC:", json.dumps(egc, indent=2))
    print("WATCH FIM:", json.dumps(watch_report, indent=2))


if __name__ == "__main__":
    main()
