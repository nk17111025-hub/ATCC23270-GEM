# -*- coding: utf-8 -*-
"""Authoritative v3 XML loader + round-trip validation + idempotent mutation."""
from __future__ import annotations

import sys
from pathlib import Path

import lxml.etree as ET
import numpy as np

import p5b6_common as C
import fw_common as FW

V3DIR = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B6" / "00_provenance" / "v3"
V3_T0 = V3DIR / "minimal_trustworthy_baseline_v3_T0.xml"
V3_TH = V3DIR / "minimal_trustworthy_baseline_v3_TH.xml"

V3_T0_SHA = "5ba88e6dd31d375ea0acbd94442e190b0de9f67d2e61169f8b73380d57becb77"
V3_TH_SHA = "4c5c5474c43d215bea6d4b24f999f1f21eab3a75dea10a997e4ffeaf47fbe1b2"

EXPECTED_RXNS = ["NDH2_NATIVE", "NOX_NATIVE"]
BIOMASS_EX = "Ex_bio[e]"


def local(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def decode_species(sid):
    comp = sid[-1]
    base = sid[2:-2].replace("_DASH_", "-")
    return f"{base}[{comp}]"


def decode_reaction(rid):
    return rid[2:].replace("_LSQBKT_", "[").replace("_RSQBKT_", "]").replace("_DASH_", "-")


def encode_species(met):
    base = met[: met.index("[")]
    comp = met[met.index("[") + 1 : -1]
    return "M_" + base.replace("-", "_DASH_") + "_" + comp


def encode_reaction(r):
    return "R_" + r.replace("[", "_LSQBKT_").replace("]", "_RSQBKT_").replace("-", "_DASH_")


def parse_v3_xml(path: Path):
    tree = ET.parse(str(path))
    root = tree.getroot()
    model = [c for c in root if local(c.tag) == "model"][0]

    metabolites = {}
    for s in [c for c in model if local(c.tag) == "listOfSpecies"][0]:
        if local(s.tag) != "species":
            continue
        sid = s.get("id")
        name = decode_species(sid)
        comp = s.get("compartment")
        formula = None
        charge = None
        for p in s.iter():
            if local(p.tag) == "p":
                txt = (p.text or "").strip()
                if txt.startswith("FORMULA:"):
                    formula = txt.split(":", 1)[1].strip()
                elif txt.startswith("CHARGE:"):
                    try:
                        charge = float(txt.split(":", 1)[1].strip())
                    except ValueError:
                        charge = None
        metabolites[name] = {"name": s.get("name") or "", "compartment": comp,
                             "formula": formula, "charge": charge}

    reactions = {}
    objective = None
    for r in [c for c in model if local(c.tag) == "listOfReactions"][0]:
        if local(r.tag) != "reaction":
            continue
        rid = decode_reaction(r.get("id"))
        reversible = r.get("reversible") == "true"
        lb = ub = 0.0
        objc = 0.0
        for p in r.iter():
            if local(p.tag) == "parameter":
                if p.get("id") == "LOWER_BOUND":
                    lb = float(p.get("value"))
                elif p.get("id") == "UPPER_BOUND":
                    ub = float(p.get("value"))
                elif p.get("id") == "OBJECTIVE_COEFFICIENT":
                    objc = float(p.get("value"))
        stoich = {}
        for sr in r.iter():
            if local(sr.tag) != "speciesReference":
                continue
            sp = decode_species(sr.get("species"))
            coef = float(sr.get("stoichiometry"))
            parent = local(sr.getparent().tag)
            if parent == "listOfReactants":
                stoich[sp] = stoich.get(sp, 0.0) - coef
            else:
                stoich[sp] = stoich.get(sp, 0.0) + coef
        reactions[rid] = {"name": r.get("name") or "", "reversible": reversible,
                          "lb": lb, "ub": ub, "stoich": {m: c for m, c in stoich.items() if abs(c) > 1e-12},
                          "objective_coefficient": objc}
        if objc != 0.0:
            objective = rid
    return metabolites, reactions, objective


def load_v3(transport="T0"):
    path = V3_T0 if transport == "T0" else V3_TH
    expect_sha = V3_T0_SHA if transport == "T0" else V3_TH_SHA
    actual = FW.sha256(path)
    assert actual == expect_sha, f"v3 {transport} hash mismatch: {actual} != {expect_sha}"
    metas, rxns, obj = parse_v3_xml(path)
    # identity assertions
    assert len(set(rxns)) == len(rxns), "duplicate reaction IDs"
    assert len(set(metas)) == len(metas), "duplicate species IDs"
    for r in EXPECTED_RXNS:
        assert r in rxns, f"{r} missing from v3 {transport}"
    assert obj is not None, "no objective reaction"
    return metas, rxns


def executable_signature(metas, rxns, objective):
    """Return a canonical tuple for round-trip comparison."""
    met_sig = [(m, metas[m]["compartment"]) for m in sorted(metas)]
    rxn_sig = []
    for r in sorted(rxns):
        rx = rxns[r]
        s = tuple(sorted(rx["stoich"].items()))
        rxn_sig.append((r, rx["reversible"], rx["lb"], rx["ub"], s))
    return (met_sig, rxn_sig, objective)


def roundtrip_validate(path: Path):
    """write already exists; reload and compare against a freshly parsed copy."""
    metas1, rxns1, obj1 = parse_v3_xml(path)
    sig1 = executable_signature(metas1, rxns1, obj1)
    return metas1, rxns1, obj1, sig1
