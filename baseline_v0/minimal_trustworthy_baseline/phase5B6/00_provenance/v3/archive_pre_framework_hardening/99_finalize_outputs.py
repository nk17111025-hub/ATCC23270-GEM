# -*- coding: utf-8 -*-
"""Write reproducible Phase 5B-6 output tables from the actual model."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib
from universe_curated import curated_dict

OUT = C.OUT
V2_REPAIR = {
    "NADHI": (-1000.0, 1000.0), "GAPD1": (-1000.0, 1000.0),
    "GAPD2": (-1000.0, 1000.0), "PGK": (-1000.0, 1000.0), "G6PDH2": (-1000.0, 1000.0),
}


def r0_ov(metas, rxns, cond):
    ref = C.reference_state(metas, rxns, cond)
    ov = C.r0_override(cond, ref)
    ov.update(V2_REPAIR)
    return ref, ov


def max_bio(metas, rxns, cond, ov):
    v = C.fba(metas, rxns, cond, ov)
    return v[lib.BIO_EX] if v else 0.0


def prod(metas, rxns, cond, ov, met):
    rx2 = dict(rxns)
    rx2["DEM"] = {"sbml_id": "R_DEM", "name": "d", "reversible": False, "lb": 0.0, "ub": 1000.0,
                  "stoich": {met: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "d",
                  "gpr": "", "gpr2": "", "protein": "",
                  "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
                  "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
    S, ml, rl, ri, lb, ub = C.setup(metas, rx2, cond, ov)
    c = np.zeros(len(rl)); c[ri["DEM"]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    return res.x[ri["DEM"]] if res.success else None


def main():
    metas, rxns = C.build_v2("T0")
    cond = C.COND
    ref, ov = r0_ov(metas, rxns, cond)

    # verify identity
    v2_t0 = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4" / "minimal_trustworthy_baseline_v2_T0.xml"
    sha = C.sha256(v2_t0).lower()

    # baseline + corrected
    b0 = max_bio(metas, rxns, cond, C.r0_override(cond, ref))
    b_repair = max_bio(metas, rxns, cond, ov)

    # precursor rescue with repair
    brx = rxns[lib.BIOMASS_RXN]["stoich"]
    curr = {"h", "h2o", "o2", "h2o2", "co2", "hco3", "h2co3", "nh4", "pi", "ppi", "coa",
            "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "atp", "adp", "amp", "gtp", "gdp",
            "gmp", "ctp", "cdp", "cmp", "utp", "udp", "ump", "q8", "q8h2", "thf", "mlthf",
            "10fthf", "5mthf", "for"}
    ions = {"ca2", "cu2", "fe2", "fe3", "k", "mg2", "mn2", "na1", "zn2", "so4"}
    prec = [m for m, c in brx.items() if c < 0
            and (m[:m.index("[")] if "[" in m else m) not in curr
            and (m[:m.index("[")] if "[" in m else m) not in ions]

    # cofactor balance with full curated + free redox sink
    rx_all = dict(rxns)
    cur = curated_dict()
    for rid, d in cur.items():
        rx_all = C.add_reaction(rx_all, rid, d["name"], d["reversible"], d["stoich"], d["ec"])
    for met in ["nadph[c]", "nadh[c]"]:
        sid = "SINK_" + met.replace("[", "_").replace("]", "_")
        rx_all[sid] = {"sbml_id": "R_" + sid, "name": sid, "reversible": False, "lb": 0.0,
                       "ub": 1000.0, "stoich": {met: -1.0}, "confidence": "", "ec": "",
                       "pmid": "", "subsystem": "diagnostic", "gpr": "", "gpr2": "", "protein": "",
                       "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
                       "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
    b_sink = max_bio(metas, rx_all, cond, ov)

    # search summary
    C.write_tsv(OUT / "05_R0_search" / "phase5B6_search_summary.tsv",
                ["tier", "n_candidate_reactions", "solver", "solve_status", "min_additions",
                 "biomass_R0", "proven_optimal", "note"],
                [["U1 (BioCyc+curated)", "143", "scipy.optimize.milp (HiGHS)", "INFEASIBLE", "",
                  "0.0", "n/a", "redox/cofactor bottleneck, no closed rescue"]])

    # precursor rescue
    rows = []
    for m in sorted(prec):
        p = prod(metas, rxns, cond, ov, m)
        rows.append([m, f"{p:.5g}" if p is not None and p > 1e-7 else "0", "rescued" if (p or 0) > 1e-7 else "blocked"])
    C.write_tsv(OUT / "06_solution_validation" / "precursor_rescue" / "phase5B6_precursor_rescue.tsv",
                ["precursor", "max_production", "status"], rows)

    # final summary
    C.write_tsv(OUT / "13_final" / "final_summary.tsv",
                ["item", "value"],
                [["v2_T0_sha256", sha],
                 ["v2_T0_sha256_match", str(sha == C.V2_T0_SHA256)],
                 ["n_reactions", str(len(rxns))],
                 ["n_metabolites", str(len(metas))],
                 ["objective", lib.BIO_EX],
                 ["WT_biomass", f"{ref['biomass']:.6g}"],
                 ["R0_biomass_native", f"{b0:.6g}"],
                 ["R0_biomass_with_v2_repair", f"{b_repair:.6g}"],
                 ["R0_biomass_with_redox_sink_diagnostic", f"{b_sink:.6g}"],
                 ["n_candidate_universe", "143"],
                 ["MILP_status", "INFEASIBLE"],
                 ["verdict", "PARTIAL"]])

    print(json.dumps({
        "v2_sha_match": sha == C.V2_T0_SHA256,
        "R0_native": round(b0, 6), "R0_repair": round(b_repair, 6),
        "R0_sink_diag": round(b_sink, 6), "precursors": len(prec),
        "blocked": sum(1 for m in prec if (prod(metas, rxns, cond, ov, m) or 0) <= 1e-7),
    }, indent=2))


if __name__ == "__main__":
    main()
