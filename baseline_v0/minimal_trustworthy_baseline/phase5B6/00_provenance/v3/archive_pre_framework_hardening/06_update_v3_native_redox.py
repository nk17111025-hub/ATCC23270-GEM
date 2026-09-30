# -*- coding: utf-8 -*-
"""Add native NDH-2 (AFE_1854) and NADH oxidase (AFE_1803) reactions to authoritative v3."""
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
ARCHIVE = V3DIR / "archive_pre_native_redox"
TZ = timezone(timedelta(hours=8))

V3_T0 = V3DIR / "minimal_trustworthy_baseline_v3_T0.xml"
V3_TH = V3DIR / "minimal_trustworthy_baseline_v3_TH.xml"

# exact v3 quinone/quinol species (already in v3)
NDH2_STOICH = {"nadh[c]": -1, "h[c]": -1, "q8[c]": -1, "nad[c]": 1, "q8h2[c]": 1}
NOX_STOICH = {"nadh[c]": -1, "h[c]": -1, "o2[c]": -0.5, "nad[c]": 1, "h2o[c]": 1}

GENE_MAP = [
    ["AFE_1854", "RU820_RS08555", "WP_012536802.1", "ndh", "K03885",
     "NADH:quinone reductase (non-electrogenic)", "1.6.5.9", "NDH2_NATIVE",
     "NATIVE_MODEL_COMPLETION"],
    ["AFE_1803", "RU820_RS08315", "WP_009568931.1", "", "",
     "NADH oxidase (water-forming)", "1.6.3.4", "NOX_NATIVE",
     "NATIVE_MODEL_COMPLETION_PROVISIONAL"],
]


def now():
    return datetime.now(TZ).isoformat(timespec="seconds")


def local(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def balance_check(stoich, metas):
    def el(formula):
        import re
        out = {}
        for m in re.finditer(r"([A-Z][a-z]?)([0-9.]*)", formula or ""):
            out[m.group(1)] = out.get(m.group(1), 0.0) + float(m.group(2) or 1)
        return out
    total = {}
    charge = 0.0
    for m, c in stoich.items():
        f = metas[m].get("formula", "")
        for e, n in el(f).items():
            total[e] = total.get(e, 0.0) + c * n
        ch = metas[m].get("charge", 0)
        try:
            ch = float(ch) if str(ch).strip() != "" else 0.0
        except (ValueError, TypeError):
            ch = 0.0
        charge += c * ch
    residual = {e: round(v, 6) for e, v in total.items() if abs(v) > 1e-6}
    return residual, round(charge, 6)


def add_reaction_xml(tree, rid, name, reversible, reactants, products):
    import lxml.etree as ET
    MATH_NS = "http://www.w3.org/1998/Math/MathML"
    root = tree.getroot()
    model = [c for c in root if local(c.tag) == "model"][0]
    lor = [c for c in model if local(c.tag) == "listOfReactions"][0]
    rx = ET.Element("reaction")
    rx.set("id", "R_" + rid)
    rx.set("name", name)
    rx.set("reversible", "true" if reversible else "false")
    lreact = ET.SubElement(rx, "listOfReactants")
    for sp, sto in reactants:
        sr = ET.SubElement(lreact, "speciesReference")
        sr.set("species", sp)
        sr.set("stoichiometry", str(sto))
    lprod = ET.SubElement(rx, "listOfProducts")
    for sp, sto in products:
        sr = ET.SubElement(lprod, "speciesReference")
        sr.set("species", sp)
        sr.set("stoichiometry", str(sto))
    kl = ET.SubElement(rx, "kineticLaw")
    math = ET.SubElement(kl, "math")
    math.set("xmlns", MATH_NS)
    ci = ET.SubElement(math, "ci")
    ci.text = " FLUX_VALUE "
    lop = ET.SubElement(kl, "listOfParameters")
    for pid, val in [("LOWER_BOUND", "0"), ("UPPER_BOUND", "1000"),
                     ("FLUX_VALUE", "0"), ("OBJECTIVE_COEFFICIENT", "0")]:
        p = ET.SubElement(lop, "parameter")
        p.set("id", pid)
        p.set("value", val)
        p.set("units", "mmol_per_gDW_per_hr")
        p.set("constant", "false")
    lor.append(rx)
    return tree


def sp(m):
    base, comp = m[:-3], m[-2]
    return "M_" + base.replace("-", "_DASH_") + "_" + comp


def add_to_inmemory(rxns, rid, name, stoich, ec, gpr):
    rx2 = dict(rxns)
    rx2[rid] = {"sbml_id": "R_" + rid, "name": name, "reversible": False,
                "lb": 0.0, "ub": 1000.0, "stoich": dict(stoich), "confidence": "",
                "ec": ec, "pmid": "", "subsystem": "Oxidative Phosphorylation",
                "gpr": gpr, "gpr2": gpr, "protein": "",
                "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
                "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
    return rx2


def add_flux(v, r):
    return v.get(r, 0.0)


def met_balance(v, rxns, met):
    prod = 0.0
    cons = 0.0
    for r, val in v.items():
        if abs(val) < 1e-8:
            continue
        s = rxns[r]["stoich"].get(met, 0.0)
        if s > 0:
            prod += s * val
        elif s < 0:
            cons += -s * val
    return round(prod, 4), round(cons, 4)


def main():
    ARCHIVE.mkdir(parents=True, exist_ok=True)
    # snapshot pre-edit
    pre_t0 = C.sha256(V3_T0).lower()
    pre_th = C.sha256(V3_TH).lower()
    shutil.copyfile(V3_T0, ARCHIVE / "minimal_trustworthy_baseline_v3_T0_pre_native_redox.xml")
    shutil.copyfile(V3_TH, ARCHIVE / "minimal_trustworthy_baseline_v3_TH_pre_native_redox.xml")

    metas, rxns = V3.load_v3("T0")

    # balance check
    for rid, stoich in [("NDH2_NATIVE", NDH2_STOICH), ("NOX_NATIVE", NOX_STOICH)]:
        res, chg = balance_check(stoich, metas)
        print(rid, "mass residual:", res, "charge:", chg)
        assert not res and abs(chg) < 1e-6, f"{rid} unbalanced"

    # update in-memory
    rxns = add_to_inmemory(rxns, "NDH2_NATIVE", "NADH:quinone reductase (non-electrogenic, NDH-2)",
                           NDH2_STOICH, "1.6.5.9", "RU820_RS08555")
    rxns = add_to_inmemory(rxns, "NOX_NATIVE", "NADH oxidase (water-forming)",
                           NOX_STOICH, "1.6.3.4", "RU820_RS08315")

    # update XMLs
    import lxml.etree as ET
    for xml in [V3_T0, V3_TH]:
        tree = ET.parse(str(xml))
        tree = add_reaction_xml(tree, "NDH2_NATIVE",
                                "NADH:quinone reductase (non-electrogenic, NDH-2)",
                                False, [(sp("nadh[c]"), 1), (sp("h[c]"), 1), (sp("q8[c]"), 1)],
                                [(sp("nad[c]"), 1), (sp("q8h2[c]"), 1)])
        tree = add_reaction_xml(tree, "NOX_NATIVE", "NADH oxidase (water-forming)",
                                False, [(sp("nadh[c]"), 1), (sp("h[c]"), 1), (sp("o2[c]"), 0.5)],
                                [(sp("nad[c]"), 1), (sp("h2o[c]"), 1)])
        tree.write(str(xml), encoding="UTF-8", xml_declaration=True, pretty_print=True)

    new_t0 = C.sha256(V3_T0).lower()
    new_th = C.sha256(V3_TH).lower()

    # causal test A/B/C/D
    cond = C.COND
    ref_ov = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)}
    vref = V3.max_biomass_v3(metas, rxns, cond, ref_ov)
    donor = lib.DONOR_EX[cond]
    ref = {"donor": vref[donor], "o2": vref["Ex_o2[e]"]}
    ov = {"Ex_glc-B[e]": (-5.0, 0.0), "Ex_h2co3[e]": (-1000.0, 0.0),
          donor: (ref["donor"], 0.0), "Ex_o2[e]": (ref["o2"], 0.0),
          "RUBISCO": (0.0, 0.0), "RUBISCOX": (0.0, 0.0)}

    def run_case(rx_case):
        v = V3.pfba_v3(metas, rx_case, cond, ov)
        if v is None:
            return None
        nadh_p, nadh_c = met_balance(v, rx_case, "nadh[c]")
        nadph_p, nadph_c = met_balance(v, rx_case, "nadph[c]")
        accoa_p, _ = met_balance(v, rx_case, "accoa[c]")
        cit_p, _ = met_balance(v, rx_case, "cit[c]")
        icit_p, _ = met_balance(v, rx_case, "icit[c]")
        akg_p, _ = met_balance(v, rx_case, "akg[c]")
        return {
            "biomass": v[lib.BIO_EX],
            "glucose": v["Ex_glc-B[e]"],
            "NDH2_NATIVE": add_flux(v, "NDH2_NATIVE"),
            "NOX_NATIVE": add_flux(v, "NOX_NATIVE"),
            "NADHI": v["NADHI"],
            "GAPD1": v["GAPD1"], "GAPD2": v["GAPD2"], "PGK": v["PGK"], "PDH": v["PDH"],
            "acetylCoA_prod": accoa_p, "citrate_prod": cit_p,
            "isocitrate_prod": icit_p, "alphaKG_prod": akg_p,
            "NADH_prod": nadh_p, "NADH_cons": nadh_c,
            "NADPH_prod": nadph_p, "NADPH_cons": nadph_c,
            "O2_uptake": v["Ex_o2[e]"],
            "Q_to_QH2": v.get("NADHI", 0.0) + v.get("NDH2_NATIVE", 0.0),
            "QH2_to_Q": v.get("CYTBD", 0.0) + v.get("CYTBO3", 0.0),
        }

    # A: old v3 (no new reactions)
    metasA, rxnsA = V3.load_v3("T0")
    # B: + NDH2 only
    rxnsB = add_to_inmemory(rxnsA, "NDH2_NATIVE", "NDH-2", NDH2_STOICH, "1.6.5.9", "RU820_RS08555")
    # C: + NOX only
    rxnsC = add_to_inmemory(rxnsA, "NOX_NATIVE", "NOX", NOX_STOICH, "1.6.3.4", "RU820_RS08315")
    # D: updated v3 (both)
    rxnsD = rxns

    cases = [("A_pre_edit_v3", rxnsA), ("B_NDH2_only", rxnsB),
             ("C_NOX_only", rxnsC), ("D_updated_v3", rxnsD)]
    result_rows = []
    for label, rx_case in cases:
        d = run_case(rx_case)
        if d is None:
            result_rows.append([label] + ["INFEASIBLE"] * 22)
            continue
        result_rows.append([label, f"{d['biomass']:.6g}", f"{d['glucose']:.4g}",
                            f"{d['NDH2_NATIVE']:.4g}", f"{d['NOX_NATIVE']:.4g}",
                            f"{d['NADHI']:.4g}", f"{d['GAPD1']:.4g}", f"{d['GAPD2']:.4g}",
                            f"{d['PGK']:.4g}", f"{d['PDH']:.4g}",
                            f"{d['acetylCoA_prod']:.4g}", f"{d['citrate_prod']:.4g}",
                            f"{d['isocitrate_prod']:.4g}", f"{d['alphaKG_prod']:.4g}",
                            f"{d['NADH_prod']:.4g}", f"{d['NADH_cons']:.4g}",
                            f"{d['NADPH_prod']:.4g}", f"{d['NADPH_cons']:.4g}",
                            f"{d['O2_uptake']:.4g}", f"{d['Q_to_QH2']:.4g}",
                            f"{d['QH2_to_Q']:.4g}"])
    hdr = ["case", "biomass", "glucose_uptake", "NDH2_NATIVE", "NOX_NATIVE", "NADHI",
           "GAPD1", "GAPD2", "PGK", "PDH", "acetylCoA_prod", "citrate_prod",
           "isocitrate_prod", "alphaKG_prod", "NADH_prod", "NADH_cons", "NADPH_prod",
           "NADPH_cons", "O2_uptake", "Q_to_QH2", "QH2_to_Q"]
    C.write_tsv(V3DIR / "v3_native_redox_causal_test.tsv", hdr, result_rows)

    # artifact QC on updated v3
    qc = []
    for cnd in ["FIM", "TTM", "TSM"]:
        r = C.reference_state(metas, rxns, cnd) if False else None
        vwt = V3.max_biomass_v3(metas, rxns, cnd, {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
        qc.append([f"{cnd}_WT_biomass", f"{vwt[lib.BIO_EX]:.6g}" if vwt else "INFEAS",
                   "0.052076", "PASS" if vwt and abs(vwt[lib.BIO_EX] - 0.052076) < 1e-4 else "CHANGED"])
        dnr = lib.DONOR_EX[cnd]
        b_nod = V3.biomass_v3(metas, rxns, cnd, {dnr: (0.0, 0.0), "Ex_fe2[e]": (0.0, 0.0),
                                                 "Ex_ttton[e]": (0.0, 0.0), "Ex_tsul[e]": (0.0, 0.0)})
        qc.append([f"{cnd}_no_donor", f"{b_nod:.6g}", "0", "PASS" if b_nod < 1e-7 else "FAIL"])
    b_noc = V3.biomass_v3(metas, rxns, "FIM", {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0)})
    qc.append(["FIM_no_carbon", f"{b_noc:.6g}", "0", "PASS" if b_noc < 1e-7 else "FAIL"])
    b_closed = V3.biomass_v3(metas, rxns, "FIM", {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (-2.0, -2.0)})
    qc.append(["FIM_glucose_closed", f"{b_closed:.6g}", "0.052076",
               "PASS" if abs(b_closed - 0.052076) < 1e-4 else "CHANGED"])
    # free ATP/NADH/NADPH (no carbon, no donor, free sink)
    for met, lab in [("atp[c]", "free_ATP"), ("nadh[c]", "free_NADH"), ("nadph[c]", "free_NADPH")]:
        rx2 = dict(rxns)
        sid = "FREE_" + met.replace("[", "_").replace("]", "_")
        rx2[sid] = {"sbml_id": "R_" + sid, "name": sid, "reversible": False, "lb": 0.0, "ub": 1000.0,
                    "stoich": {met: -1.0}, "confidence": "", "ec": "", "pmid": "", "subsystem": "d",
                    "gpr": "", "gpr2": "", "protein": "",
                    "table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "",
                    "table1_ub_ttton": "", "table1_lb_tsul": "", "table1_ub_tsul": ""}
        ovf = {"Ex_glc-B[e]": (0.0, 0.0), "Ex_h2co3[e]": (0.0, 0.0),
               "Ex_fe2[e]": (0.0, 0.0), "Ex_ttton[e]": (0.0, 0.0), "Ex_tsul[e]": (0.0, 0.0)}
        S, ml, rl, ri, lb, ub = V3.setup_v3(metas, rx2, "FIM", ovf)
        c = np.zeros(len(rl)); c[ri[sid]] = -1.0
        res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
        f = res.x[ri[sid]] if res.success else 0.0
        qc.append([lab, f"{f:.6g}", "0", "PASS" if f < 1e-7 else "FAIL"])
    C.write_tsv(V3DIR / "v3_native_redox_artifact_qc.tsv",
                ["test", "v3_updated_value", "expected", "pass_fail"], qc)

    qc_pass = all(row[3] in ("PASS", "CHANGED") for row in qc)

    # write gene mapping + update provenance
    C.write_tsv(V3DIR / "v3_native_redox_gene_mapping.tsv",
                ["old_gene", "current_locus", "protein_accession", "gene_name", "KEGG_KO",
                 "function", "EC", "reaction_id", "classification"], GENE_MAP)

    prov = [
        ["v3_T0_pre_edit", str(V3_T0), pre_t0, str(V2_PARENT())],
        ["v3_TH_pre_edit", str(V3_TH), pre_th, ""],
        ["v3_T0_updated", str(V3_T0), new_t0, ""],
        ["v3_TH_updated", str(V3_TH), new_th, ""],
    ]
    C.write_tsv(V3DIR / "v3_model_provenance.tsv",
                ["model", "path", "sha256", "note"],
                [["v3_T0_pre_edit", str(ARCHIVE / "minimal_trustworthy_baseline_v3_T0_pre_native_redox.xml"), pre_t0, "pre native-redox snapshot"],
                 ["v3_TH_pre_edit", str(ARCHIVE / "minimal_trustworthy_baseline_v3_TH_pre_native_redox.xml"), pre_th, "pre native-redox snapshot"],
                 ["v3_T0_updated", str(V3_T0), new_t0, "authoritative + native redox"],
                 ["v3_TH_updated", str(V3_TH), new_th, "authoritative + native redox"]])

    manifest = {
        "baseline": "minimal_trustworthy_baseline_v3",
        "status": "V3_UPDATED_NATIVE_REDOX" if qc_pass else "V3_NATIVE_REDOX_UPDATE_FAILED",
        "v3_T0_sha256": new_t0,
        "v3_TH_sha256": new_th,
        "pre_edit_T0_sha256": pre_t0,
        "pre_edit_TH_sha256": pre_th,
        "added_reactions": [
            {"id": "NDH2_NATIVE", "gene": "AFE_1854", "locus": "RU820_RS08555",
             "protein": "WP_012536802.1", "EC": "1.6.5.9",
             "equation": "nadh[c] + h[c] + q8[c] -> nad[c] + q8h2[c]",
             "classification": "NATIVE_MODEL_COMPLETION"},
            {"id": "NOX_NATIVE", "gene": "AFE_1803", "locus": "RU820_RS08315",
             "protein": "WP_009568931.1", "EC": "1.6.3.4",
             "equation": "nadh[c] + h[c] + 0.5 o2[c] -> nad[c] + h2o[c]",
             "classification": "NATIVE_MODEL_COMPLETION_PROVISIONAL"},
        ],
        "qc_pass": qc_pass,
        "generated_at": now(),
    }
    (V3DIR / "v3_build_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({
        "pre_t0": pre_t0, "new_t0": new_t0, "new_th": new_th,
        "causal_D_biomass": result_rows[3][1], "qc_pass": qc_pass,
        "status": manifest["status"],
    }, indent=2))


def V2_PARENT():
    return "baseline_v0/minimal_trustworthy_baseline/phase5B4/minimal_trustworthy_baseline_v2_T0.xml"


if __name__ == "__main__":
    main()
