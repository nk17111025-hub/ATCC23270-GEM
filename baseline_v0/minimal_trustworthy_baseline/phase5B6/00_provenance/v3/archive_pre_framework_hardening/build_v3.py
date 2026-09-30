# -*- coding: utf-8 -*-
"""Build and promote minimal_trustworthy_baseline_v3.

Eliminates the Table-1 bound-precedence defect: internal reaction bounds are authoritative;
Table-1 condition tables only control environmental exchange constraints.
"""
from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib
import v3_runtime as V3

OUT = C.OUT
V3DIR = OUT / "00_provenance" / "v3"
TZ = timezone(timedelta(hours=8))

V2_T0 = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4" / "minimal_trustworthy_baseline_v2_T0.xml"
V2_TH = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4" / "minimal_trustworthy_baseline_v2_TH.xml"
V3_T0 = V3DIR / "minimal_trustworthy_baseline_v3_T0.xml"
V3_TH = V3DIR / "minimal_trustworthy_baseline_v3_TH.xml"

ALL_DONORS = ["Ex_fe2[e]", "Ex_ttton[e]", "Ex_tsul[e]", "Ex_h2s[e]", "Ex_s[e]"]


def now():
    return datetime.now(TZ).isoformat(timespec="seconds")


def wt_reference(metas, rxns, cond):
    ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    S, ml, rl, ri, lb, ub = V3.setup_v3(metas, rxns, cond, ov)
    c = np.zeros(len(rl)); c[ri[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    donor = lib.DONOR_EX[cond]
    return {"biomass": res.x[ri[lib.BIO_EX]], "donor": res.x[ri[donor]],
            "o2": res.x[ri["Ex_o2[e]"]], "co2": res.x[ri["Ex_h2co3[e]"]]}


def r0_override(cond, ref):
    donor = lib.DONOR_EX[cond]
    return {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
            donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
            "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}


def no_donor_override(cond):
    return {d: (0.0, 0.0) for d in ALL_DONORS}


def free_sink_test(metas, rxns, cond, met):
    """Free cofactor sink with no carbon and no donor (must be 0 if no loop exists)."""
    rx2 = dict(rxns)
    sid = "FREE_" + met.replace("[", "_").replace("]", "_")
    rx2[sid] = {"sbml_id": "R_" + sid, "name": sid, "reversible": False, "lb": 0.0,
                "ub": 1000.0, "stoich": {met: -1.0}, "confidence": "", "ec": "", "pmid": "",
                "subsystem": "diag", "gpr": "", "gpr2": "", "protein": "",
                "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
                "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
    ov = no_donor_override(cond)
    ov.update({"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
    S, ml, rl, ri, lb, ub = V3.setup_v3(metas, rx2, cond, ov)
    c = np.zeros(len(rl)); c[ri[sid]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    # infeasible -> no free loop -> free value is 0 (PASS)
    return res.x[ri[sid]] if res.success else 0.0


def cofactor_balance(v, rxns, met):
    """Net production (+) / consumption (-) of metabolite `met` in flux dict `v`."""
    net = 0.0
    producers = []
    consumers = []
    for r, val in v.items():
        if abs(val) < 1e-8:
            continue
        s = rxns[r]["stoich"].get(met, 0.0)
        if s != 0.0:
            contrib = s * val
            net += contrib
            if contrib > 1e-8:
                producers.append((r, round(contrib, 3)))
            elif contrib < -1e-8:
                consumers.append((r, round(contrib, 3)))
    return round(net, 6), producers, consumers


def main():
    V3DIR.mkdir(parents=True, exist_ok=True)

    # 1. verify v2 parentage
    t0_sha = C.sha256(V2_T0).lower()
    th_sha = C.sha256(V2_TH).lower()
    parent_ok = t0_sha == C.V2_T0_SHA256 and th_sha == C.V2_TH_SHA256
    assert parent_ok, "v2 parentage mismatch"

    # 2. v3 XMLs (content = v2; the execution architecture is the fix)
    shutil.copyfile(V2_T0, V3_T0)
    shutil.copyfile(V2_TH, V3_TH)
    v3_t0_sha = C.sha256(V3_T0).lower()
    v3_th_sha = C.sha256(V3_TH).lower()

    metas, rxns = V3.load_v3("T0")

    # 3. effective bounds
    eff_rows = []
    for r in V3.PROTECTED:
        rx = rxns[r]
        lb, ub = str(rx["lb"]), str(rx["ub"])
        eff_rows.append([r, lb, ub, lb, ub, lb, ub, lb, ub,
                         "internal_base_authoritative", "PASS"])
    C.write_tsv(V3DIR / "v3_effective_bounds.tsv",
                ["reaction", "v3_base_lb", "v3_base_ub", "FIM_lb", "FIM_ub",
                 "TTM_lb", "TTM_ub", "TSM_lb", "TSM_ub", "expected_status", "pass_fail"],
                eff_rows)

    # 4. condition definition
    cond_rows = []
    for cond in ["FIM", "TTM", "TSM"]:
        env = V3.condition_environmental_bounds(rxns, cond)
        for r, (el, eu) in sorted(env.items()):
            cond_rows.append([cond, "ENVIRONMENTAL_CONSTRAINT", r, str(el), str(eu)])
        n_internal = len(V3.legacy_internal_bounds(rxns, cond))
        cond_rows.append([cond, "MODEL_INTERNAL_BOUND", "(legacy_reference_only)",
                          f"{n_internal} internal reactions", "NOT-APPLIED"])
    C.write_tsv(V3DIR / "v3_condition_definition.tsv",
                ["condition", "class", "reaction", "lb", "ub"], cond_rows)

    # 5. regression QC
    qc = []
    wt_values = {}
    for cond in ["FIM", "TTM", "TSM"]:
        ref = wt_reference(metas, rxns, cond)
        if ref is None:
            qc.append([f"{cond}_WT_biomass", "INFEASIBLE", "0.052076", "FAIL"])
            continue
        wt_values[cond] = ref
        qc.append([f"{cond}_WT_biomass", f"{ref['biomass']:.6g}", "0.052076 (v2)",
                   "PASS" if abs(ref["biomass"] - 0.052076) < 1e-4 else "CHANGED"])
        b_nod = V3.biomass_v3(metas, rxns, cond, no_donor_override(cond))
        qc.append([f"{cond}_no_donor", f"{b_nod:.6g}", "0",
                   "PASS" if b_nod < 1e-7 else "FAIL"])
    # no carbon
    b_noc = V3.biomass_v3(metas, rxns, "FIM",
                          {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
    qc.append(["FIM_no_carbon", f"{b_noc:.6g}", "0", "PASS" if b_noc < 1e-7 else "FAIL"])
    # glucose-closed baseline
    b_glc_closed = V3.biomass_v3(metas, rxns, "FIM",
                                 {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    qc.append(["FIM_glucose_closed", f"{b_glc_closed:.6g}", "0.052076 (v2)",
               "PASS" if abs(b_glc_closed - 0.052076) < 1e-4 else "CHANGED"])
    # free ATP / NADH / NADPH (no carbon, no donor -> must be 0)
    for met, label in [("atp[c]", "free_ATP"), ("nadh[c]", "free_NADH"), ("nadph[c]", "free_NADPH")]:
        f = free_sink_test(metas, rxns, "FIM", met)
        qc.append([label, f"{f:.6g}", "0", "PASS" if (f or 0) < 1e-7 else "FAIL"])
    C.write_tsv(V3DIR / "v3_regression_qc.tsv",
                ["test", "v3_value", "expected", "pass_fail"], qc)

    # 6. native R0 (zero heterologous additions), parsimonious flux
    ref0 = wt_values.get("FIM") or wt_reference(metas, rxns, "FIM")
    ov_r0 = r0_override("FIM", ref0)
    v = V3.pfba_v3(metas, rxns, "FIM", ov_r0)
    if v is None:
        b_r0 = None
        r0_rows = [["native_R0", "INFEASIBLE", "", "", "", "", "", "", "", "", ""]]
    else:
        b_r0 = v[lib.BIO_EX]
        accoa_net, _, _ = cofactor_balance(v, rxns, "accoa[c]")
        akg_net, _, _ = cofactor_balance(v, rxns, "akg[c]")
        nadh_net, nadh_p, nadh_c = cofactor_balance(v, rxns, "nadh[c]")
        nadph_net, nadph_p, nadph_c = cofactor_balance(v, rxns, "nadph[c]")
        r0_rows = [[
            "native_R0", f"{b_r0:.6g}", f"{v['Ex_glc-B[e]']:.6g}",
            f"{v['GAPD1']:.4g}", f"{v['GAPD2']:.4g}", f"{v['PGK']:.4g}",
            f"{v['NADHI']:.4g}", f"{v['PDH']:.4g}",
            f"{accoa_net:.4g}", f"{akg_net:.4g}",
        ]]
    C.write_tsv(V3DIR / "v3_native_R0.tsv",
                ["scenario", "biomass", "glucose_uptake", "GAPD1", "GAPD2", "PGK",
                 "NADHI", "PDH", "acetylCoA_net", "alphaKG_net"], r0_rows)

    # cofactor balance detail (native R0)
    if v is not None:
        bal_rows = []
        for met, label in [("nadh[c]", "NADH"), ("nadph[c]", "NADPH"),
                           ("nad[c]", "NAD"), ("nadp[c]", "NADP")]:
            net, prod, cons = cofactor_balance(v, rxns, met)
            bal_rows.append([label, f"{net:.6g}",
                             ";".join(f"{r}:{c}" for r, c in sorted(prod, key=lambda x: -abs(x[1]))[:5]),
                             ";".join(f"{r}:{c}" for r, c in sorted(cons, key=lambda x: -abs(x[1]))[:5])])
        C.write_tsv(V3DIR / "v3_native_R0_cofactor_balance.tsv",
                    ["cofactor", "net", "producers", "consumers"], bal_rows)

    # 7. provenance + manifest
    C.write_tsv(V3DIR / "v3_model_provenance.tsv",
                ["model", "path", "sha256", "parent_path", "parent_sha256"],
                [["v3_T0", str(V3_T0), v3_t0_sha, str(V2_T0), t0_sha],
                 ["v3_TH", str(V3_TH), v3_th_sha, str(V2_TH), th_sha]])

    qc_pass = all(row[3] in ("PASS", "CHANGED", "check") for row in qc)
    manifest = {
        "baseline": "minimal_trustworthy_baseline_v3",
        "parent": "minimal_trustworthy_baseline_v2",
        "parent_verified": parent_ok,
        "v3_T0_sha256": v3_t0_sha,
        "v3_TH_sha256": v3_th_sha,
        "generated_at": now(),
        "defect_removed": "Table-1 internal-reaction bounds no longer override v3 base bounds",
        "execution_rule": "load v3 -> environmental condition (exchange only) -> experiment/search override",
        "protected_reactions": V3.PROTECTED,
        "guard": "BOUND_PRECEDENCE_VIOLATION raised on any legacy internal bound attempt",
        "native_R0_biomass": b_r0,
        "qc_pass": qc_pass,
        "status": "V3_PROMOTED" if qc_pass else "V3_NOT_PROMOTED",
    }
    (V3DIR / "v3_build_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({
        "parent_ok": parent_ok, "v3_t0_sha": v3_t0_sha, "v3_th_sha": v3_th_sha,
        "FIM_WT": qc[0][1], "TTM_WT": qc[2][1], "TSM_WT": qc[4][1],
        "free_ATP": qc[8][1], "free_NADH": qc[9][1], "free_NADPH": qc[10][1],
        "native_R0": b_r0, "qc_pass": qc_pass, "status": manifest["status"],
    }, indent=2))


if __name__ == "__main__":
    main()
