# -*- coding: utf-8 -*-
"""Framework hardening orchestration with HARD assertions. Never emits PASS unless all pass."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import p5b6_common as C
import fw_common as FW
import v3_model as M
import conditions as CO
import gapfill as G
import qc as Q
import regression2016 as R16

OUT = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B6" / "00_provenance" / "v3" / "framework_hardening"
OUT.mkdir(parents=True, exist_ok=True)
REPRO = C.ROOT / "baseline_v0" / "reproduction_2016"

GROWTH_TOL = 1e-3  # conservative absolute growth delta vs accepted baseline


def main():
    metas, rxns = M.load_v3("T0")

    # 1. identity
    t0_sha = FW.sha256(M.V3_T0).lower()
    th_sha = FW.sha256(M.V3_TH).lower()
    identity_ok = t0_sha == M.V3_T0_SHA and th_sha == M.V3_TH_SHA
    FW.write_tsv(OUT / "authoritative_model_identity.tsv",
                 ["model", "sha256", "expected", "n_met", "n_rxn", "match"],
                 [["v3_T0", t0_sha, M.V3_T0_SHA, len(metas), len(rxns), str(identity_ok)],
                  ["v3_TH", th_sha, M.V3_TH_SHA, "-", "-", str(th_sha == M.V3_TH_SHA)]])

    # 2. round-trip
    m_rt, r_rt, o_rt, sig_rt = M.roundtrip_validate(M.V3_T0)
    m_o, r_o, o_o = M.parse_v3_xml(M.V3_T0)
    sig_o = M.executable_signature(m_o, r_o, o_o)
    roundtrip_ok = sig_rt == sig_o
    FW.write_tsv(OUT / "v3_roundtrip_validation.tsv",
                 ["check", "result"],
                 [["metabolite_ids", "PASS" if len(sig_rt[0]) == len(sig_o[0]) else "FAIL"],
                  ["reaction_ids", "PASS" if len(sig_rt[1]) == len(sig_o[1]) else "FAIL"],
                  ["executable_signature", "PASS" if roundtrip_ok else "FAIL"],
                  ["NDH2_NATIVE_once", "PASS" if list(rxns).count("NDH2_NATIVE") == 1 else "FAIL"],
                  ["NOX_NATIVE_once", "PASS" if list(rxns).count("NOX_NATIVE") == 1 else "FAIL"]])

    # 3. MILP unit tests
    unit = G.run_all_unit_tests()
    unit_pass = sum(1 for _, r, _ in unit if r == "PASS")
    unit_total = len(unit)
    FW.write_tsv(OUT / "milp_unit_tests.tsv", ["test", "result", "detail"], unit)

    # 4. closed zero background + free-loop tests
    czb_status, czb_obj = Q.closed_zero_background_feasible(metas, rxns)
    loop_rows = Q.free_loop_tests(metas, rxns)
    FW.write_tsv(OUT / "free_loop_tests.tsv",
                 ["test", "solver_status", "maximum_flux", "tolerance", "result"], loop_rows)
    # phenotype tests (separate)
    pheno_rows = Q.phenotype_tests(metas, rxns)
    # artifact QC hardened = free-loop + phenotype + closed-zero background
    art_rows = [["CLOSED_ZERO_BACKGROUND_FEASIBILITY", czb_status, f"{czb_obj}", "OPTIMAL",
                 "PASS" if czb_status == "OPTIMAL" else "FAIL"]]
    for r in loop_rows:
        art_rows.append([r[0], r[1], str(r[2]), str(r[3]), r[4]])
    for r in pheno_rows:
        art_rows.append([r[0], r[1], str(r[2]), "phenotype", r[1]])
    FW.write_tsv(OUT / "artifact_qc_hardened.tsv",
                 ["test", "solver_status", "value", "note", "result"], art_rows)

    # 5. basic smoke regression (renamed; NOT the 2016 panel)
    smoke_rows = Q.basic_smoke_regression(metas, rxns)
    smoke_pass = sum(1 for r in smoke_rows if r[4] == "PASS")
    FW.write_tsv(OUT / "basic_v3_smoke_regression.tsv",
                 ["test", "solver_status", "value", "expected", "result"], smoke_rows)

    # 6. 2016 published regression (20 scenarios + 7 validation series)
    scenarios = R16.run_2016(metas, rxns)
    accepted = _read_accepted()
    scen_rows = []
    scen_pass = 0
    for s in scenarios:
        a = accepted[s["idx"]]
        growth_ok = s["growth"] is not None and abs(s["growth"] - float(a["growth_pred"])) < GROWTH_TOL
        donor_ok = s["donor"] is not None and abs(s["donor"] - float(a["donor_pred_magnitude"])) < 1e-3
        o2_ok = (s["o2"] is None and a["o2_pred_magnitude"] == "") or \
                (s["o2"] is not None and a["o2_pred_magnitude"] != "" and
                 abs(s["o2"] - float(a["o2_pred_magnitude"])) < 1e-3)
        status = s["status"]
        result = "PASS" if (status == "optimal" and growth_ok and donor_ok and o2_ok) else "FAIL"
        scen_pass += (1 if result == "PASS" else 0)
        scen_rows.append([s["idx"], s["condition"], str(s["exp_growth"]),
                          a["growth_pred"], f"{s['growth']:.8g}" if s["growth"] is not None else "",
                          f"{abs(s['growth'] - float(a['growth_pred'])):.3g}" if s["growth"] is not None else "",
                          s["donor_type"], str(s["exp_donor"]), a["donor_pred_magnitude"],
                          f"{s['donor']:.6g}" if s["donor"] is not None else "",
                          f"{abs(s['donor'] - float(a['donor_pred_magnitude'])):.3g}" if s["donor"] is not None else "",
                          str(s["exp_o2"] or ""), a["o2_pred_magnitude"],
                          f"{s['o2']:.6g}" if s["o2"] is not None else "",
                          f"{abs(s['o2'] - float(a['o2_pred_magnitude'])):.3g}" if (s["o2"] is not None and a["o2_pred_magnitude"] != "") else "",
                          status, result])
    FW.write_tsv(OUT / "complete_2016_20scenario_regression.tsv",
                 ["scenario_id", "condition", "experimental_growth", "accepted_baseline_growth",
                  "current_v3_growth", "growth_delta", "donor_type", "experimental_donor",
                  "accepted_baseline_donor", "current_v3_donor", "donor_delta", "experimental_o2",
                  "accepted_baseline_o2", "current_v3_o2", "o2_delta", "solver_status", "result"],
                 scen_rows)

    series = R16.validation_series(scenarios)
    accepted_metrics = _read_accepted_metrics()
    ser_rows = []
    ser_pass = 0
    for s in series:
        a = accepted_metrics[s["series"]]
        delta = abs(s["r2_1_minus_sse_sst"] - float(a["r2_coefficient_1_minus_SSE_SST"]))
        result = "PASS" if delta < 1e-3 else "FAIL"
        ser_pass += (1 if result == "PASS" else 0)
        ser_rows.append([s["series"], str(s["n"]), a["r2_coefficient_1_minus_SSE_SST"],
                         f"{s['r2_1_minus_sse_sst']:.8g}", f"{delta:.3g}",
                         a["slope"], f"{s['slope']:.6g}", result])
    FW.write_tsv(OUT / "complete_2016_validation_metrics.tsv",
                 ["series", "n", "accepted_R2_SSE_SST", "current_v3_R2_SSE_SST", "delta",
                  "accepted_slope", "current_v3_slope", "result"], ser_rows)

    mean_r2 = float(np.mean([s["r2_1_minus_sse_sst"] for s in series]))
    accepted_mean_r2 = float(np.mean([float(a["r2_coefficient_1_minus_SSE_SST"]) for a in accepted_metrics.values()]))
    max_growth_delta = max((abs(s["growth"] - float(accepted[s["idx"]]["growth_pred"]))
                            for s in scenarios if s["growth"] is not None), default=0.0)
    summary = {
        "n_scenarios": len(scenarios),
        "n_scenarios_optimal": sum(1 for s in scenarios if s["status"] == "optimal"),
        "n_scenarios_pass": scen_pass,
        "n_validation_series": len(series),
        "mean_R2_SSE_SST": mean_r2,
        "accepted_mean_R2_SSE_SST": accepted_mean_r2,
        "max_growth_delta": max_growth_delta,
        "overall_result": "PASS" if (scen_pass == 20 and ser_pass == 7) else "FAIL",
    }
    FW.write_json(OUT / "complete_2016_regression_summary.json", summary)

    # HARD assertions (before writing PASS)
    loop_solver_statuses = [r[1] for r in loop_rows]
    loop_results = [r[4] for r in loop_rows]
    free_loop_solver_ok = all(s == "OPTIMAL" for s in loop_solver_statuses)
    free_loop_results_ok = all(r == "PASS" for r in loop_results)

    assert len(scenarios) == 20, f"scenarios={len(scenarios)}"
    assert len(series) == 7, f"series={len(series)}"
    assert czb_status == "OPTIMAL", f"closed_zero_background={czb_status}"
    assert free_loop_solver_ok, f"free loop solver statuses: {loop_solver_statuses}"
    assert free_loop_results_ok, f"free loop results: {loop_results}"
    assert not any(r[1] == "INFEASIBLE" and r[4] == "PASS" for r in loop_rows), "INFEASIBLE->PASS detected"

    all_ok = (identity_ok and roundtrip_ok and unit_pass == unit_total and unit_total == 6
              and czb_status == "OPTIMAL" and free_loop_solver_ok and free_loop_results_ok
              and scen_pass == 20 and ser_pass == 7)
    state = {
        "authoritative_model_path": str(M.V3_T0),
        "T0_hash": t0_sha, "TH_hash": th_sha,
        "runtime": "v3_model.py (XML loader)",
        "solver": "scipy.optimize.linprog (HiGHS)",
        "glucose_min_uptake": CO.GLUCOSE_MIN_UPTAKE,
        "basic_smoke_regression": f"{smoke_pass}/4",
        "published_2016_regression": f"{scen_pass}/20",
        "published_2016_validation_series": f"{ser_pass}/7",
        "milp_unit_tests": f"{unit_pass}/{unit_total}",
        "closed_zero_background": czb_status,
        "free_loop_all_optimal": free_loop_solver_ok,
        "free_loop_all_pass": free_loop_results_ok,
        "framework_status": "FRAMEWORK_HARDENING_PASS" if all_ok else "FRAMEWORK_HARDENING_FAIL",
    }
    FW.write_json(OUT / "phase5B6_framework_state.json", state)

    # bug resolution
    bugs = _bug_rows(identity_ok, roundtrip_ok, czb_status, free_loop_solver_ok, free_loop_results_ok,
                     scen_pass, ser_pass, unit_pass)
    FW.write_tsv(OUT / "framework_bug_resolution.tsv",
                 ["bug_id", "description", "old_behavior", "repair", "validation_test", "result", "status"],
                 bugs)

    print(json.dumps(state, indent=2))
    if not all_ok:
        sys.exit(1)


def _read_accepted():
    hdr, rows = FW.read_tsv(REPRO / "table5_20point_reproduction.tsv")
    idx = {h: i for i, h in enumerate(hdr)}
    return {int(r[idx["idx"]]): {h: r[idx[h]] for h in hdr} for r in rows}


def _read_accepted_metrics():
    hdr, rows = FW.read_tsv(REPRO / "seven_validation_metrics.tsv")
    return {r[0]: {h: r[i] for i, h in enumerate(hdr)} for r in rows}


def _bug_rows(identity_ok, roundtrip_ok, czb, loop_opt, loop_pass, scen_pass, ser_pass, unit_pass):
    def r(bug, desc, old, repair, test, result, status):
        return [bug, desc, old, repair, test, result, status]
    rows = [
        r("BLOCKER A", "INFEASIBLE counted as PASS in zero tests", "INFEASIBLE -> PASS",
          "zero_test_result(): only OPTIMAL + |v|<=tol is PASS", "free-loop tests",
          "PASS", "FIXED"),
        r("BLOCKER A2", "free-energy test background infeasible", "maintenance demand blocks zero state",
          "CLOSED_ZERO_BACKGROUND (close donors/carbon, ATPM lb=0)", "closed_zero_background_feasible",
          "PASS" if czb == "OPTIMAL" else "FAIL", "FIXED" if czb == "OPTIMAL" else "NOT_FIXED"),
        r("BLOCKER B", "4-row smoke called complete 2016 regression", "4 rows",
          "reuse Phase 2: 20 scenarios + 7 validation series", "complete_2016_20scenario_regression",
          "PASS" if scen_pass == 20 else "FAIL", "FIXED" if scen_pass == 20 else "NOT_FIXED"),
    ]
    return rows


if __name__ == "__main__":
    main()
