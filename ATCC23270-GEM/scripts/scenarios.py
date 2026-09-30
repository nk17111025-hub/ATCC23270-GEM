# -*- coding: utf-8 -*-
"""Deterministic FIM/TTM/TSM scenario and perturbation application for ATCC23270-GEM.

Recovers the validated historical scenario semantics:

    load SBML (base bounds authoritative)
      -> apply condition exchange bounds (environment only)
      -> apply perturbation override (WT / no_donor / no_carbon / glucose_closed)

Internal reaction bounds are NEVER overridden by a condition table.
The objective is the biomass exchange (Ex_bio[e]).

Pure-Python FBA (scipy.linprog) so the new repository has no COBRApy/solver dependency
beyond scipy + lxml.
"""
from __future__ import annotations

import json
from pathlib import Path

import lxml.etree as ET
import numpy as np
from scipy.optimize import linprog


SBML_NS = "{http://www.sbml.org/sbml/level2}"
CONFIG = Path(__file__).resolve().parent.parent / "config" / "scenarios.json"


def _local(tag: str) -> str:
    return tag.split("}", 1)[1] if "}" in tag else tag


def encode_reaction(name: str) -> str:
    """'Ex_bio[e]' -> 'Ex_bio_LSQBKT_e_RSQBKT_' (SBML id fragment)."""
    return name.replace("[", "_LSQBKT_").replace("]", "_RSQBKT_").replace("-", "_DASH_")


def decode_reaction(rid: str) -> str:
    """'Ex_bio_LSQBKT_e_RSQBKT_' -> 'Ex_bio[e]'."""
    return rid[2:].replace("_LSQBKT_", "[").replace("_RSQBKT_", "]").replace("_DASH_", "-")


def decode_species(sid: str) -> str:
    comp = sid[-1]
    base = sid[2:-2].replace("_DASH_", "-")
    return f"{base}[{comp}]"


def load_model(path: str | Path) -> dict:
    """Parse an SBML Level 2 model into {reactions, species, objective}."""
    tree = ET.parse(str(path))
    root = tree.getroot()
    model = next(c for c in root if _local(c.tag) == "model")

    species = {}
    for sp in model.iter():
        if _local(sp.tag) != "species":
            continue
        species[sp.get("id")] = sp.get("compartment")

    reactions = {}
    objective = None
    for r in model.iter():
        if _local(r.tag) != "reaction":
            continue
        sid = r.get("id")
        lb = ub = 0.0
        objc = 0.0
        for p in r.iter():
            if _local(p.tag) == "parameter":
                if p.get("id") == "LOWER_BOUND":
                    lb = float(p.get("value"))
                elif p.get("id") == "UPPER_BOUND":
                    ub = float(p.get("value"))
                elif p.get("id") == "OBJECTIVE_COEFFICIENT":
                    objc = float(p.get("value"))
        stoich = {}
        for sr in r.iter():
            if _local(sr.tag) != "speciesReference":
                continue
            met = decode_species(sr.get("species"))
            coef = float(sr.get("stoichiometry"))
            parent = _local(sr.getparent().tag)
            if parent == "listOfReactants":
                stoich[met] = stoich.get(met, 0.0) - coef
            else:
                stoich[met] = stoich.get(met, 0.0) + coef
        reactions[sid] = {
            "name": r.get("name") or "",
            "reversible": r.get("reversible") == "true",
            "lb": lb,
            "ub": ub,
            "stoich": {m: c for m, c in stoich.items() if abs(c) > 1e-12},
        }
        if objc != 0.0:
            objective = sid
    return {"reactions": reactions, "species": species, "objective": objective}


def _build_matrices(model: dict):
    rxns = model["reactions"]
    rxn_list = list(rxns)
    met_list = sorted({m for r in rxn_list for m in rxns[r]["stoich"]})
    met_idx = {m: i for i, m in enumerate(met_list)}
    S = np.zeros((len(met_list), len(rxn_list)))
    for j, r in enumerate(rxn_list):
        for m, c in rxns[r]["stoich"].items():
            if m in met_idx:
                S[met_idx[m], j] = c
    return S, met_list, rxn_list, met_idx


def apply_scenario(model: dict, condition: str, perturbation: str | None = None,
                   config: dict | None = None) -> dict:
    """Return (S, met_list, rxn_list, rxn_idx, lb, ub, objective_id) for one scenario."""
    if config is None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
    rxns = model["reactions"]
    S, met_list, rxn_list, met_idx = _build_matrices(model)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        lb[j], ub[j] = rxns[r]["lb"], rxns[r]["ub"]

    # 1. condition exchange bounds (environment only)
    cond_bounds = config["conditions"][condition]["exchange_bounds"]
    for ex_name, (el, eu) in cond_bounds.items():
        sbml_id = "R_" + encode_reaction(ex_name)
        if sbml_id in rxn_idx:
            lb[rxn_idx[sbml_id]] = float(el)
            ub[rxn_idx[sbml_id]] = float(eu)

    # 2. perturbation override
    if perturbation == "no_donor":
        donor = config["conditions"][condition]["donor_exchange"]
        sbml_id = "R_" + encode_reaction(donor)
        if sbml_id in rxn_idx:
            lb[rxn_idx[sbml_id]] = 0.0
            ub[rxn_idx[sbml_id]] = 0.0
    elif perturbation:
        ov = config["perturbations"][perturbation]
        for ex_name, (ol, ou) in ov.items():
            sbml_id = "R_" + encode_reaction(ex_name)
            if sbml_id in rxn_idx:
                lb[rxn_idx[sbml_id]] = float(ol)
                ub[rxn_idx[sbml_id]] = float(ou)

    return S, met_list, rxn_list, rxn_idx, lb, ub


def optimize(model: dict, condition: str, perturbation: str | None = None,
             config: dict | None = None) -> tuple:
    """Maximize biomass. Return (status, objective_value, flux_dict_or_None)."""
    S, met_list, rxn_list, rxn_idx, lb, ub = apply_scenario(
        model, condition, perturbation, config)
    objective = model["objective"]
    c = np.zeros(len(rxn_list))
    if objective in rxn_idx:
        c[rxn_idx[objective]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)),
                  bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        status = "infeasible" if res.status == 2 else "unbounded" if res.status == 3 else f"error_{res.status}"
        return status, None, None
    fluxes = {r: res.x[rxn_idx[r]] for r in rxn_list}
    return "optimal", float(-res.fun), fluxes


def scenario_key(condition: str, perturbation: str) -> str:
    return f"{condition}_{perturbation}"
