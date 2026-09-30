# -*- coding: utf-8 -*-
"""Phase 5B-3B: OptStrain-style minimal heterologous reaction search."""
import csv
import json
import itertools
import datetime as dt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import phase4a_lib as lib
from build_glucose_model import build_memory

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase5B3B_optstrain_direct_glucose")

# Non-native candidate database: (id, name, EC, mechanism, reversible, stoich, native_likelihood)
CAND = [
    ("G6PDH_NAD", "glucose-6-phosphate dehydrogenase (NAD)", "1.1.1.363", "ED", False,
     {"g6p-B[c]": -1, "nad[c]": -1, "6pgl[c]": 1, "nadh[c]": 1, "h[c]": 1}, "POSSIBLY_NATIVE"),
    ("PGDH6", "6-phosphogluconate dehydrogenase (decarboxylating)", "1.1.1.44", "PPP", False,
     {"6pgc[c]": -1, "nadp[c]": -1, "ru5p-D[c]": 1, "co2[c]": 1, "nadph[c]": 1}, "POSSIBLY_NATIVE"),
    ("ICL", "isocitrate lyase", "4.1.3.1", "glyoxylate", False,
     {"icit[c]": -1, "glx[c]": 1, "succ[c]": 1}, "POSSIBLY_NATIVE"),
    ("ME_NADP", "malic enzyme (NADP)", "1.1.1.40", "anaplerosis", False,
     {"mal-L[c]": -1, "nadp[c]": -1, "pyr[c]": 1, "co2[c]": 1, "nadph[c]": 1}, "POSSIBLY_NATIVE"),
    ("PC", "pyruvate carboxylase", "6.4.1.1", "anaplerosis", False,
     {"pyr[c]": -1, "co2[c]": -1, "atp[c]": -1, "oaa[c]": 1, "adp[c]": 1, "pi[c]": 1, "h[c]": 1}, "POSSIBLY_NATIVE"),
    ("PEPCK", "PEP carboxykinase (ATP)", "4.1.1.49", "gluconeogenesis", True,
     {"oaa[c]": -1, "atp[c]": -1, "pep[c]": 1, "co2[c]": 1, "adp[c]": 1}, "POSSIBLY_NATIVE"),
    ("PEPS", "PEP synthase", "2.7.9.2", "gluconeogenesis", False,
     {"pyr[c]": -1, "atp[c]": -1, "h2o[c]": -1, "pep[c]": 1, "amp[c]": 1, "pi[c]": 1, "h[c]": 2}, "POSSIBLY_NATIVE"),
    ("PPDK", "pyruvate phosphate dikinase", "2.7.9.1", "gluconeogenesis", False,
     {"pyr[c]": -1, "atp[c]": -1, "pi[c]": -1, "pep[c]": 1, "amp[c]": 1, "ppi[c]": 1}, "POSSIBLY_NATIVE"),
    ("NADHDH", "NADH:ubiquinone oxidoreductase (forward)", "1.6.5.3", "redox", False,
     {"nadh[c]": -1, "q8[c]": -1, "h[c]": -1, "nad[c]": 1, "q8h2[c]": 1}, "POTENTIAL_MODEL_OMISSION"),
    ("NADHox", "NADH oxidase", "1.6.3.4", "redox", False,
     {"nadh[c]": -1, "o2[c]": -1, "h[c]": -1, "nad[c]": 1, "h2o2[c]": 1}, "LIKELY_HETEROLOGOUS"),
    ("PFL", "pyruvate formate lyase", "2.3.1.54", "pyruvate", False,
     {"pyr[c]": -1, "coa[c]": -1, "accoa[c]": 1, "for[c]": 1}, "POSSIBLY_NATIVE"),
    ("GAPDH_NADP_fwd", "glyceraldehyde-3-P dehydrogenase (NADP, oxidative)", "1.2.1.13", "glycolysis", False,
     {"g3p[c]": -1, "nadp[c]": -1, "pi[c]": -1, "13dpg[c]": 1, "nadph[c]": 1, "h[c]": 1}, "MODEL_EQUIVALENT_GAPD2"),
]

CAPS = [0.25, 0.5, 1.0]
BLANK = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
         "table1_lb_tsul": "", "table1_ub_tsul": ""}


def add_reactions(rxns, add_ids):
    rx2 = dict(rxns)
    for rid, name, ec, mech, rev, stoich, nat in CAND:
        if rid not in add_ids:
            continue
        rx2["NON_" + rid] = {"sbml_id": "R_NON_" + rid, "name": name + " (NON_NATIVE_CANDIDATE)",
                             "reversible": rev, "lb": -1000.0 if rev else 0.0, "ub": 1000.0,
                             "stoich": stoich, "confidence": "", "ec": ec, "pmid": "",
                             "subsystem": "NON_NATIVE_CANDIDATE", "gpr": "", "gpr2": "", "protein": "", **BLANK}
    return rx2


def setup(metas, rxns, cond, override):
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n); ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]; lb[j], ub[j] = rx["lb"], rx["ub"]
    lbc, ubc = lib.COND_COL[cond]
    for r in rxn_list:
        if r == lib.BIOMASS_RXN: continue
        rx = rxns[r]
        lv = rx["table1_"+lbc]; uv = rx["table1_"+ubc]
        if lv != "" and uv != "":
            lb[rxn_idx[r]] = float(lv); ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx: lb[rxn_idx[r]] = pl; ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx: lb[rxn_idx[r]] = ol; ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def max_biomass(metas, rxns, cond, override):
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, override)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return res.x[rxn_idx[lib.BIO_EX]]


def ref_state(metas, rxns, cond):
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    S, met_list, rxn_list, rxn_idx, lb, ub = setup(metas, rxns, cond, ov)
    c = np.zeros(len(rxn_list)); c[rxn_idx[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=list(zip(lb, ub)), method="highs")
    v = res.x
    return {"biomass": v[rxn_idx[lib.BIO_EX]], "donor": v[rxn_idx[lib.DONOR_EX[cond]]],
            "o2": v[rxn_idx["Ex_o2[e]"]], "co2": v[rxn_idx["Ex_h2co3[e]"]],
            "rubisco": v[rxn_idx["RUBISCO"]] + v[rxn_idx["RUBISCOX"]]}


def scenario_override(cond, ref, glc_cap, rubisco_limit, ci_mode, added):
    donor = lib.DONOR_EX[cond]
    ov = {"Ex_glc-B[e]": (-glc_cap, 0.0), donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
          "RUBISCO": (0.0, rubisco_limit), "RUBISCOX": (0.0, rubisco_limit)}
    if ci_mode == "low":
        ov["Ex_h2co3[e]"] = (-0.1 * abs(ref["co2"]), 0.0)
    elif ci_mode == "zero":
        ov["Ex_h2co3[e]"] = (0.0, 0.0)
    else:
        ov["Ex_h2co3[e]"] = (-1000.0, 0.0)
    return ov


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metas, rxns0 = build_memory("T0")
    cond = "FIM"
    ref = ref_state(metas, rxns0, cond)
    rub_ref = abs(ref["rubisco"])
    R = {"R25": 0.25 * rub_ref, "R10": 0.10 * rub_ref, "R0": 0.0}

    # reaction database summary
    with (OUT / "reaction_database_summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "name", "EC", "mechanism", "reversible", "native_likelihood", "source"])
        for rid, name, ec, mech, rev, stoich, nat in CAND:
            w.writerow([rid, name, ec, mech, str(rev), nat, "curated universal central-carbon (EC)"])

    # native vs non-native mapping
    with (OUT / "native_vs_non_native_mapping.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["candidate", "native_equivalent", "classification"])
        w.writerow(["GAPDH_NADP_fwd", "GAPD2 (13dpg+NADPH->g3p+NADP, reversible)", "MODEL_EQUIVALENT (excluded as non-native)"])
        w.writerow(["G6PDH_NAD", "G6PDH2 (NADP only, reversible)", "POSSIBLY_NATIVE (NAD isoform missing)"])
        w.writerow(["NADHDH", "NADHI (reverse-ETC only, NAD+->NADH)", "POTENTIAL_MODEL_OMISSION (direction)"])
        w.writerow(["ICL/ME/PC/PEPCK/PEPS/PPDK/PFL/PGDH6", "absent", "TRUE_NETWORK_ADDITION / LIKELY_HETEROLOGOUS"])

    # scenario matrix (no additions)
    matrix = []
    for rid in ["R25", "R10", "R0"]:
        for ci in ["full", "low", "zero"]:
            ov = scenario_override(cond, ref, 1.0, R[rid], ci, [])
            b = max_biomass(metas, rxns0, cond, ov)
            matrix.append({"scenario": rid, "ci": ci, "max_biomass": round(b, 6) if b is not None else "infeasible",
                           "G0": "yes" if b is not None and b >= ref["biomass"] - 1e-6 else "no"})
    with (OUT / "optstrain_scenario_matrix.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["scenario", "ci", "max_biomass", "G0"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in matrix:
            w.writerow(r)

    # OptStrain: minimum non-native set for each scenario
    sols = []
    cand_ids = [c[0] for c in CAND if c[0] != "GAPDH_NADP_fwd"]  # exclude model-equivalent
    for rid in ["R0", "R10", "R25"]:
        for ci in ["zero", "low"]:
            # 1-reaction
            found = []
            for cid in cand_ids:
                rx2 = add_reactions(rxns0, [cid])
                b = max_biomass(metas, rx2, cond, scenario_override(cond, ref, 1.0, R[rid], ci, []))
                if b is not None and b >= ref["biomass"] - 1e-6:
                    found.append((tuple([cid]), b))
            # 2-reaction
            if not found:
                for p in itertools.combinations(cand_ids, 2):
                    rx2 = add_reactions(rxns0, list(p))
                    b = max_biomass(metas, rx2, cond, scenario_override(cond, ref, 1.0, R[rid], ci, []))
                    if b is not None and b >= ref["biomass"] - 1e-6:
                        found.append((p, b))
            # 3-reaction
            if not found:
                for p in itertools.combinations(cand_ids, 3):
                    rx2 = add_reactions(rxns0, list(p))
                    b = max_biomass(metas, rx2, cond, scenario_override(cond, ref, 1.0, R[rid], ci, []))
                    if b is not None and b >= ref["biomass"] - 1e-6:
                        found.append((p, b))
            if not found:
                sols.append({"scenario": rid, "ci": ci, "reactions": "NONE_<=3", "n": None,
                             "max_biomass": "", "mechanism": "", "G10": "no"})
            else:
                # top solutions by biomass
                for p, b in sorted(found, key=lambda x: -x[1])[:8]:
                    mechs = sorted({dict((c[0], c[3]) for c in CAND)[x] for x in p})
                    sols.append({"scenario": rid, "ci": ci, "reactions": "+".join(p), "n": len(p),
                                 "max_biomass": round(b, 6), "mechanism": ",".join(mechs),
                                 "G10": "yes" if b >= 1.10*ref["biomass"] else "no"})

    f = ["scenario", "ci", "reactions", "n", "max_biomass", "mechanism", "G10"]
    with (OUT / "optstrain_minimal_solutions.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=f, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sols:
            w.writerow(r)

    print("=== OptStrain minimal solutions ===")
    for r in sols:
        if r["n"] is not None and r["n"] <= 2:
            print(f"  {r['scenario']} {r['ci']}: {r['reactions']} (n={r['n']}) -> {r['max_biomass']} mech={r['mechanism']}")


if __name__ == "__main__":
    main()
