# -*- coding: utf-8 -*-
"""Part B/C: build evidence-corrected glucose models (T0, TH) + write XML + pre-edit table."""
import csv
import hashlib
import json
from pathlib import Path

import lxml.etree as ET

import phase4a_lib as lib

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
FROZEN = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase4_1G_glucose_corrected_model"

NS = "{http://www.sbml.org/sbml/level2}"
MATH_NS = "http://www.w3.org/1998/Math/MathML"
XHTML_NS = "http://www.w3.org/1999/xhtml"


def local(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


# Part B pre-edit table (reaction_id, equation, old_lb, old_ub, old_GPR, gene, evidence, proposed_change, reason)
CORRECTIONS = [
    ["BDGK", "atp[c] + glc-B[c] -> adp[c] + g6p-B[c] + h[c]", -1000.0, 1000.0, "",
     "AFE_2841/RU820_RS13140 (candidate only)", "ROK-family annotation only; no assay",
     "set lb 0 (forward-only)", "prevent reverse ATP-generating glucose formation; glucokinase activity unresolved"],
    ["G6PDH2", "6pgl[c] + h[c] + nadph[c] -> g6p-B[c] + nadp[c]", 0.0, 1000.0, "AFE_2025",
     "AFE_2025/RU820_RS09360 (zwf, EC 1.1.1.49)", "G6PDH annotation (EC + domains); no flux assay",
     "set lb -1000 (open oxidative direction)", "reaction identity is G6PDH (oxidative PPP); reductive-only bound is a reconstruction constraint"],
    ["GAPD1", "13dpg[c] + h[c] + nadh[c] -> g3p[c] + nad[c] + pi[c]", 0.0, 1000.0, "AFE_3251",
     "AFE_3251/RU820_RS15125 (gap, EC 1.2.1.12)", "type-I GAPDH annotation; reversible enzyme",
     "set lb -1000 (open glycolytic direction)", "type-I GAPDH is reversible; one-way gluconeogenic bound is a reconstruction constraint"],
    ["GAPD2", "13dpg[c] + h[c] + nadph[c] -> g3p[c] + nadp[c] + pi[c]", 0.0, 1000.0, "AFE_3251",
     "AFE_3251/RU820_RS15125 (gap, EC 1.2.1.13)", "same gene, NADP-dependent GAPDH",
     "set lb -1000 (open glycolytic direction)", "same reasoning as GAPD1"],
    ["PGK", "3pg[c] + atp[c] -> 13dpg[c] + adp[c]", 0.0, 1000.0, "AFE_3250",
     "AFE_3250/RU820_RS15120 (pgk, EC 2.7.2.3)", "Swiss-Prot reviewed; reversible reaction; GO glycolysis+gluconeogenesis",
     "set lb -1000 (open glycolytic direction)", "reviewed bidirectional PGK; one-way gluconeogenic bound is a reconstruction constraint"],
    ["PFK", "atp[c] + f6p-B[c] -> adp[c] + fdp[c] + h[c]", 0.0, 1000.0, "AFE_1807",
     "AFE_1807/RU820_RS08325 (pfkB)", "Wang 2012: functional PfkB (heteroexpression + mutant)",
     "no change (GPR already correct)", "experimentally supported PfkB; already forward glycolysis direction"],
]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_xml(transport: str):
    """transport = 'T0' (facilitated) or 'TH' (H+ symport). Returns modified XML tree."""
    tree = ET.parse(str(FROZEN))
    root = tree.getroot()
    model = [ch for ch in root if local(ch.tag) == "model"][0]

    # 1. bound corrections
    bound_map = {
        "R_BDGK": 0.0,
        "R_G6PDH2": -1000.0,
        "R_GAPD1": -1000.0,
        "R_GAPD2": -1000.0,
        "R_PGK": -1000.0,
    }
    for r in model.iter():
        if local(r.tag) == "reaction" and r.get("id") in bound_map:
            for p in r.iter():
                if local(p.tag) == "parameter" and p.get("id") == "LOWER_BOUND":
                    p.set("value", str(bound_map[r.get("id")]))

    # 2. add species glc-B[e], glc-B[p]
    los = [ch for ch in model if local(ch.tag) == "listOfSpecies"][0]
    template = None
    for s in los:
        if s.get("id") == "M_glc_DASH_B_c":
            template = s
            break
    for spid, comp, nm in [("M_glc_DASH_B_e", "e", "beta-D-Glucose"), ("M_glc_DASH_B_p", "p", "beta-D-Glucose")]:
        ns = ET.fromstring(ET.tostring(template))
        ns.set("id", spid)
        ns.set("compartment", comp)
        ns.set("name", nm)
        los.append(ns)

    # 3. add reactions
    lor = [ch for ch in model if local(ch.tag) == "listOfReactions"][0]

    def new_reaction(rid, name, reversible, reactants, products, lb, ub):
        rx = ET.Element("reaction")
        rx.set("id", rid)
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
        for pid, val in [("LOWER_BOUND", str(lb)), ("UPPER_BOUND", str(ub)),
                         ("FLUX_VALUE", "0"), ("OBJECTIVE_COEFFICIENT", "0")]:
            p = ET.SubElement(lop, "parameter")
            p.set("id", pid)
            p.set("value", val)
            p.set("units", "mmol_per_gDW_per_hr")
            p.set("constant", "false")
        return rx

    # exchange
    lor.append(new_reaction("R_Ex_glc_DASH_B_LSQBKT_e_RSQBKT_", "beta-D-Glucose exchange (GPR unresolved)",
                            True, [("M_glc_DASH_B_e", 1)], [], -1000.0, 1000.0))
    # outer porin
    lor.append(new_reaction("R_GLCtex", "glucose outer-membrane porin (GPR unresolved)",
                            True, [("M_glc_DASH_B_e", 1)], [("M_glc_DASH_B_p", 1)], -1000.0, 1000.0))
    # inner transporter
    if transport == "T0":
        lor.append(new_reaction("R_GLCtpp", "glucose inner-membrane facilitated transport (GPR unresolved)",
                                True, [("M_glc_DASH_B_p", 1)], [("M_glc_DASH_B_c", 1)], -1000.0, 1000.0))
    else:
        lor.append(new_reaction("R_GLCtpp", "glucose inner-membrane H+ symport (GPR unresolved)",
                                True, [("M_glc_DASH_B_p", 1), ("M_h_p", 1)],
                                [("M_glc_DASH_B_c", 1), ("M_h_c", 1)], -1000.0, 1000.0))
    return tree


def build_memory(transport: str):
    """In-memory corrected model dicts for FBA."""
    metas, rxns = lib.parse_sbml_meta()
    # bound corrections
    rxns["BDGK"]["lb"] = 0.0
    for rid in ["G6PDH2", "GAPD1", "GAPD2", "PGK"]:
        rxns[rid]["lb"] = -1000.0
    # add metabolites
    for key, comp in [("glc-B[e]", "e"), ("glc-B[p]", "p")]:
        metas[key] = {"name": "beta-D-Glucose", "formula": "C6H12O6", "charge": 0,
                      "location": comp, "kegg": "", "pubchem": "", "base": "glc-B", "compartment": comp}
    blank = {"table1_lb_fe2": "", "table1_ub_fe2": "", "table1_lb_ttton": "", "table1_ub_ttton": "",
             "table1_lb_tsul": "", "table1_ub_tsul": ""}
    rxns["Ex_glc-B[e]"] = {"sbml_id": "R_Ex_glc_DASH_B_LSQBKT_e_RSQBKT_", "name": "beta-D-Glucose exchange",
                           "reversible": True, "lb": -1000.0, "ub": 1000.0, "stoich": {"glc-B[e]": -1.0},
                           "confidence": "", "ec": "", "pmid": "", "subsystem": "hypothetical glucose uptake",
                           "gpr": "", "gpr2": "", "protein": "", **blank}
    rxns["GLCtex"] = {"sbml_id": "R_GLCtex", "name": "glucose outer porin", "reversible": True,
                      "lb": -1000.0, "ub": 1000.0, "stoich": {"glc-B[e]": -1.0, "glc-B[p]": 1.0},
                      "confidence": "", "ec": "", "pmid": "", "subsystem": "hypothetical glucose uptake",
                      "gpr": "", "gpr2": "", "protein": "", **blank}
    if transport == "T0":
        rxns["GLCtpp"] = {"sbml_id": "R_GLCtpp", "name": "glucose facilitated inner transport",
                          "reversible": True, "lb": -1000.0, "ub": 1000.0,
                          "stoich": {"glc-B[p]": -1.0, "glc-B[c]": 1.0},
                          "confidence": "", "ec": "", "pmid": "", "subsystem": "hypothetical glucose uptake",
                          "gpr": "", "gpr2": "", "protein": "", **blank}
    else:
        rxns["GLCtpp"] = {"sbml_id": "R_GLCtpp", "name": "glucose H+ symport inner transport",
                          "reversible": True, "lb": -1000.0, "ub": 1000.0,
                          "stoich": {"glc-B[p]": -1.0, "h[p]": -1.0, "glc-B[c]": 1.0, "h[c]": 1.0},
                          "confidence": "", "ec": "", "pmid": "", "subsystem": "hypothetical glucose uptake",
                          "gpr": "", "gpr2": "", "protein": "", **blank}
    return metas, rxns


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frozen_sha = sha256(FROZEN)

    # write corrections TSV
    cf = ["reaction_id", "equation", "old_lb", "old_ub", "old_GPR", "relevant_gene",
          "evidence", "proposed_change", "reason"]
    with (OUT / "glucose_model_corrections.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cf)
        for row in CORRECTIONS:
            w.writerow(row)
        # transport additions
        w.writerow(["Ex_glc-B[e] / GLCtex / GLCtpp", "glc-B[e] <-> glc-B[p] <-> glc-B[c]", "absent", "absent", "",
                    "unresolved transporter", "ATCC 23270 assimilates exogenous glucose (Wang 2012); transporter gene unresolved",
                    "ADD provisional native glucose uptake (two energy variants T0/TH); no transporter GPR", "glucose uptake is a model omission; mechanism/energetics unresolved"])

    # write XML files
    for transport in ["T0", "TH"]:
        tree = build_xml(transport)
        outp = OUT / f"minimal_trustworthy_baseline_v1_glucose_corrected_{transport}.xml"
        tree.write(str(outp), encoding="UTF-8", xml_declaration=True, pretty_print=True)
        print(transport, "written", outp.name, "sha", sha256(outp))

    # validate in-memory models parse and report
    for transport in ["T0", "TH"]:
        metas, rxns = build_memory(transport)
        print(transport, "n_met", len(metas), "n_rxn", len(rxns), "glc-B[e] present:", "glc-B[e]" in metas,
              "Ex present:", "Ex_glc-B[e]" in rxns, "GLCtpp present:", "GLCtpp" in rxns)

    summary = {"frozen_v1_sha256": frozen_sha,
               "frozen_sha_match": frozen_sha == "bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c",
               "corrections": CORRECTIONS}
    (OUT / "_build_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("frozen sha", frozen_sha, "match", summary["frozen_sha_match"])


if __name__ == "__main__":
    main()
