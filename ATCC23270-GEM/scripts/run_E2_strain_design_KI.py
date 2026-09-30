#!/usr/bin/env python
"""OptStrain-like, KI-only search for E2 using the installed StrainDesign API.

The authoritative host and E1 pool are read-only inputs. Candidate reactions
exist in the search model with source-supported bounds; StrainDesign's ki_cost
decides which additions are enabled in each design.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import cobra

ROOT = Path(__file__).resolve().parents[1]
HOST = ROOT.parent / "before final" / "updatedv3.1_verified" / "updatedv3.1.xml"
EXPECTED_SHA = "026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3"
POOL = ROOT / "results/engineering/E1_donor_reaction_pool/candidate_reactions.tsv"
DONOR_DIR = ROOT / "external/reference/donor_models"
OUT = ROOT / "results/engineering/E2_strain_design_KI"
STRAINDESIGN = ROOT / "external/straindesign"
DONORS = {
    "iML1515": DONOR_DIR / "iML1515.json",
    "iJN1463": DONOR_DIR / "iJN1463.xml",
    "iCN1361": DONOR_DIR / "iCN1361.xml",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def decode_id(s: str) -> str:
    return s.replace("_LSQBKT_", "[").replace("_RSQBKT_", "]").replace("_DASH_", "-")


def resolve_rxn(model, requested: str):
    decoded = decode_id(requested)
    for rxn in model.reactions:
        if rxn.id == requested or decode_id(rxn.id) == decoded:
            return rxn
    raise KeyError(f"Reaction not found: {requested}")


def load_e1_helpers():
    path = ROOT / "scripts/build_E1_donor_reaction_pool.py"
    spec = importlib.util.spec_from_file_location("e1_helpers", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_models(e1):
    os.chdir(ROOT)
    host, _ = cobra.io.validate_sbml_model("../before final/updatedv3.1_verified/updatedv3.1.xml")
    if host is None:
        raise RuntimeError("Failed to read authoritative host")
    donors = {}
    for name, path in DONORS.items():
        rel = path.relative_to(ROOT).as_posix()
        donors[name] = cobra.io.load_json_model(rel) if path.suffix == ".json" else cobra.io.read_sbml_model(rel)
    alias_models = [("host", host), *donors.items()]
    alias_index = e1.build_alias_index(alias_models)
    alias_canonical, _ = e1.build_alias_canonical(alias_models, alias_index)
    return host, donors, alias_index, alias_canonical


def get_source_map(rows):
    source_map = defaultdict(list)
    for row in rows:
        cid = row["candidate_id"]
        for item in row["source_reactions"].split(";"):
            donor, rid = item.split(":", 1)
            orientation = "forward"
            for entry in row.get("orientation_to_canonical", "").split(";"):
                if entry.startswith(item + "="):
                    orientation = entry.rsplit("=", 1)[1]
                    break
            source_map[cid].append((donor, rid, orientation))
    return source_map


def add_candidate_reactions(host, donors, alias_index, alias_canonical, rows):
    """Add candidate source chemistry, with irreversible opposite directions split."""
    h = load_e1_helpers()
    host_met_by_key = defaultdict(list)
    for met in host.metabolites:
        key, _ = h.metabolite_key(met, "host", alias_index, alias_canonical)
        comp = h.compartment_class(met)
        host_met_by_key[(key, comp)].append(met)
    comp_ids = {h.compartment_class(m): m.compartment for m in host.metabolites}
    src_map = get_source_map(rows)
    candidate_cost = {}
    provenance = {}
    variant_meta = {}
    for row in rows:
        cid = row["candidate_id"]
        srcs = []
        for donor, rid, orient in src_map[cid]:
            rxn = donors[donor].reactions.get_by_id(rid)
            # Transform each donor interval to canonical candidate orientation.
            if orient == "reverse":
                lb, ub = -rxn.upper_bound, -rxn.lower_bound
            else:
                lb, ub = rxn.lower_bound, rxn.upper_bound
            srcs.append((donor, rxn, orient, lb, ub))
        has_fwd = any(ub > 0 for _, _, _, _, ub in srcs)
        has_rev = any(lb < 0 for _, _, _, lb, _ in srcs)
        supports_both = any(lb < 0 and ub > 0 for _, _, _, lb, ub in srcs)
        split = has_fwd and has_rev and not supports_both
        variants = [("fwd", [x for x in srcs if x[4] > 0], 0.0, max(x[4] for x in srcs if x[4] > 0))] if split else []
        if split:
            variants.append(("rev", [x for x in srcs if x[3] < 0], min(x[3] for x in srcs if x[3] < 0), 0.0))
        else:
            variants = [("", srcs, min([x[3] for x in srcs] + [0.0]), max([x[4] for x in srcs] + [0.0]))]
        for suffix, contributors, lb, ub in variants:
            rid = cid + ("__" + suffix if suffix else "")
            exemplar = contributors[0][1]
            new = cobra.Reaction(rid, name=f"E1 candidate {cid}" + (f" ({suffix})" if suffix else ""))
            new.lower_bound, new.upper_bound = lb, ub
            # A row can have several donor sources for the same chemistry. Use
            # one normalized source stoichiometry; contributors only support
            # the available direction/bounds and provenance.
            stoich = defaultdict(float)
            donor, source, orient, _, _ = contributors[0]
            sign = -1.0 if orient == "reverse" else 1.0
            for met, coeff in source.metabolites.items():
                key, _ = h.metabolite_key(met, donor, alias_index, alias_canonical)
                cls = h.compartment_class(met)
                matches = host_met_by_key.get((key, cls), [])
                if matches:
                    target = matches[0]
                else:
                    comp = comp_ids.get(cls, met.compartment or "c")
                    safe = re.sub(r"[^A-Za-z0-9_]+", "_", key)
                    mid = "E1M_" + hashlib.sha1((safe + "@" + cls).encode()).hexdigest()[:14] + "_" + comp
                    if mid not in host.metabolites:
                        host.add_metabolites([cobra.Metabolite(mid, name=met.name or key, compartment=comp)])
                    target = host.metabolites.get_by_id(mid)
                    host_met_by_key[(key, cls)].append(target)
                stoich[target] += sign * coeff
            new.add_metabolites({m: v for m, v in stoich.items() if abs(v) > 1e-12})
            new.annotation["e1_candidate_id"] = cid
            host.add_reactions([new])
            candidate_cost[rid] = 1.0
            donor_ids = sorted({f"{d}:{rx.id}" for d, rx, *_ in contributors})
            provenance[rid] = {"candidate_id": cid, "source_reactions": donor_ids, "variant": suffix or "combined"}
            variant_meta[rid] = {"candidate_id": cid, "lower_bound": lb, "upper_bound": ub}
    return candidate_cost, provenance, variant_meta


def apply_condition(model, config, relaxed=True):
    fim = config["conditions"]["FIM"]["exchange_bounds"]
    for raw, bounds in fim.items():
        resolve_rxn(model, raw).bounds = tuple(bounds)
    for raw, bounds in config["perturbations"]["WT"].items():
        resolve_rxn(model, raw).bounds = tuple(bounds)
    if relaxed:
        resolve_rxn(model, "Ex_h2co3[e]").bounds = (-2.0, 0.0)


def fba_summary(model, ids):
    sol = model.optimize()
    if sol.status != "optimal":
        return {"status": sol.status}
    return {"status": sol.status, "biomass": float(sol.fluxes[ids["biomass"].id]),
            "glucose_uptake": float(sol.fluxes[ids["glucose"].id]),
            "rubisco_flux": float(sol.fluxes[ids["rubisco"].id]),
            "co2_flux": float(sol.fluxes[ids["co2"].id]),
            "fe2_flux": float(sol.fluxes[ids["fe2"].id]), "o2_flux": float(sol.fluxes[ids["o2"].id])}


def scenario_constraints(ids, base, factor=None):
    cons = [f"{ids['biomass'].id} >= {0.8 * base['biomass']:.12g}",
            f"{ids['glucose'].id} <= {-0.25 * abs(base['glucose_uptake']):.12g}",
            f"{ids['fe2'].id} <= -0.01", f"{ids['rubisco'].id} >= 0"]
    if factor is not None:
        cons.append(f"{ids['rubisco'].id} <= {factor * base['rubisco_flux']:.12g}")
    return cons


def set_cbb_scenario_bounds(model, ids, base, factor):
    rubisco = model.reactions.get_by_id(ids["rubisco"].id)
    rubisco.bounds = (0.0, rubisco.upper_bound if factor is None else factor * base["rubisco_flux"])


def sd_solve(model, candidates, constraints, max_cost, time_limit, max_solutions=3):
    sys.path.insert(0, str(STRAINDESIGN))
    from straindesign import SDModule, compute_strain_designs
    # The default GLPK precheck becomes numerically unreliable on this merged
    # network; compute_strain_designs rechecks feasibility using selected SCIP.
    module = SDModule(model, "protect", constraints=constraints, skip_checks=True)
    start = time.time()
    sols = compute_strain_designs(model, sd_modules=[module], solver="scip", ki_cost=candidates,
                                  ko_cost={}, max_cost=max_cost, max_solutions=max_solutions,
                                  solution_approach="best", compress=True, time_limit=time_limit,
                                  skip_preprocessing_fvas=True, seed=1)
    return sols, time.time() - start


def materialize_candidate_bounds(model, candidate_ids, active):
    for rid in candidate_ids:
        r = model.reactions.get_by_id(rid)
        # In the no-KI validation, candidates are physically unavailable.
        if not active:
            r.bounds = (0.0, 0.0)


def summarize_solutions(sols, provenance):
    out = []
    for i, design in enumerate(sols.get_reaction_sd(), 1):
        selected = sorted(k for k, v in design.items() if k in provenance and float(v) > 0)
        out.append({"solution_id": f"S{i:02d}", "reaction_ids": selected,
                    "candidate_ids": sorted({provenance[r]["candidate_id"] for r in selected}),
                    "source_reactions": sorted({s for r in selected for s in provenance[r]["source_reactions"]})})
    return out


def smoke_sd_ki_semantics():
    """Tiny positive control: a protect phenotype must select one KI."""
    sys.path.insert(0, str(STRAINDESIGN))
    from straindesign import SDModule, compute_strain_designs
    model = cobra.Model("E2_KI_smoke")
    a = cobra.Metabolite("a_c", compartment="c")
    b = cobra.Metabolite("b_c", compartment="c")
    source = cobra.Reaction("SRC_A"); source.add_metabolites({a: 1}); source.bounds = (0, 10)
    ki = cobra.Reaction("E1R_SMOKE_KI"); ki.add_metabolites({a: -1, b: 1}); ki.bounds = (0, 10)
    demand = cobra.Reaction("DM_B"); demand.add_metabolites({b: -1}); demand.bounds = (0, 10)
    model.add_reactions([source, ki, demand]); model.objective = demand
    with model:
        ki.bounds = (0, 0)
        model.add_cons_vars(model.problem.Constraint(demand.flux_expression, lb=1))
        no_ki_status = model.optimize().status
    module = SDModule(model, "protect", constraints=["DM_B >= 1"], skip_checks=True)
    sols = compute_strain_designs(model, sd_modules=[module], solver="scip", ki_cost={ki.id: 1}, ko_cost={},
                                  max_cost=1, max_solutions=1, solution_approach="best", compress=False,
                                  time_limit=30, seed=1)
    designs = sols.get_reaction_sd()
    selected = [rxn for design in designs for rxn, val in design.items() if rxn == ki.id and float(val) > 0]
    return {"no_KI_status": no_ki_status, "solver_status": str(sols.status), "selected_KI": sorted(set(selected)),
            "passed": no_ki_status == "infeasible" and sols.status == "optimal" and ki.id in selected}


def phenotype_fba(model, ids, constraints, candidate_bounds=None):
    """Maximize biomass under the same linear phenotype constraints."""
    with model:
        for rid, bounds in (candidate_bounds or {}).items():
            model.reactions.get_by_id(rid).bounds = bounds
        for text in constraints:
            m = re.fullmatch(r"(.+?)\s*(<=|>=|=)\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)", text)
            if not m:
                raise ValueError(f"Unsupported phenotype constraint: {text}")
            rid, op, value = m.groups()
            rxn = model.reactions.get_by_id(rid)
            val = float(value)
            lb, ub = (None, val) if op == "<=" else (val, None) if op == ">=" else (val, val)
            model.add_cons_vars(model.problem.Constraint(rxn.flux_expression, lb=lb, ub=ub, name="E2_" + str(len(model.solver.constraints))))
        sol = model.optimize()
        if sol.status != "optimal":
            return {"status": sol.status}
        return {"status": sol.status, "biomass": float(sol.fluxes[ids["biomass"].id]),
                "glucose_uptake": float(sol.fluxes[ids["glucose"].id]),
                "rubisco_flux": float(sol.fluxes[ids["rubisco"].id]),
                "co2_flux": float(sol.fluxes[ids["co2"].id]), "fe2_flux": float(sol.fluxes[ids["fe2"].id]),
                "o2_flux": float(sol.fluxes[ids["o2"].id])}


def glucose_off_max_biomass(model, ids, candidate_bounds=None):
    """Maximum biomass without glucose under the current scenario's CBB/Fe bounds."""
    glucose_off_constraints = [f"{ids['fe2'].id} <= -0.01", f"{ids['rubisco'].id} >= 0"]
    with model:
        model.reactions.get_by_id(ids["glucose"].id).bounds = (0.0, 0.0)
        for rid, bounds in (candidate_bounds or {}).items():
            model.reactions.get_by_id(rid).bounds = bounds
        for text in glucose_off_constraints:
            m = re.fullmatch(r"(.+?)\s*(<=|>=|=)\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)", text)
            rid, op, value = m.groups()
            val = float(value)
            lb, ub = (None, val) if op == "<=" else (val, None)
            rxn = model.reactions.get_by_id(rid)
            model.add_cons_vars(model.problem.Constraint(rxn.flux_expression, lb=lb, ub=ub,
                                                         name="E2_off_" + str(len(model.solver.constraints))))
        sol = model.optimize()
        if sol.status == "infeasible":
            return {"status": "infeasible", "max_biomass": 0.0}
        if sol.status != "optimal":
            return {"status": sol.status, "max_biomass": None}
        return {"status": "optimal", "max_biomass": float(sol.fluxes[ids["biomass"].id]),
                "glucose_uptake": float(sol.fluxes[ids["glucose"].id]),
                "rubisco_flux": float(sol.fluxes[ids["rubisco"].id]),
                "co2_flux": float(sol.fluxes[ids["co2"].id]), "fe2_flux": float(sol.fluxes[ids["fe2"].id]),
                "o2_flux": float(sol.fluxes[ids["o2"].id])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="run baseline checks and one bounded KI formulation test")
    ap.add_argument("--run", action="store_true", help="run scenario costs 1 through 6")
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument("--no-resume", action="store_false", dest="resume")
    ap.add_argument("--time-limit", type=int, default=900)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "checkpoints").mkdir(exist_ok=True)
    (OUT / "logs").mkdir(exist_ok=True)
    observed = sha256(HOST)
    if observed != EXPECTED_SHA:
        raise SystemExit(f"Host checksum mismatch: {observed}")
    e1 = load_e1_helpers()
    rows = read_tsv(POOL)
    if len(rows) != 2424:
        raise SystemExit(f"Expected 2424 E1 candidates, observed {len(rows)}")
    host, donors, alias_index, alias_canonical = load_models(e1)
    config = json.loads((ROOT / "config/scenarios.json").read_text(encoding="utf-8"))
    ids = {k: resolve_rxn(host, v) for k, v in {"biomass": "Ex_bio[e]", "glucose": "Ex_glc-B[e]",
          "rubisco": "RUBISCO", "co2": "Ex_h2co3[e]", "fe2": "Ex_fe2[e]", "o2": "Ex_o2[e]"}.items()}
    host.objective = ids["biomass"]
    # Record original FIM + WT phenotype; search baseline has finite CO2 availability.
    original = host.copy(); apply_condition(original, config, relaxed=False)
    original_baseline = fba_summary(original, ids)
    search_host = host.copy(); apply_condition(search_host, config, relaxed=True)
    initial_search_baseline = fba_summary(search_host, ids)
    # Bound the search to at most the baseline glucose supply and oxygen uptake.
    # Negative exchange flux denotes uptake in this model.
    search_glucose = resolve_rxn(search_host, "Ex_glc-B[e]")
    search_glucose.bounds = (max(search_glucose.lower_bound, initial_search_baseline["glucose_uptake"]), 0.0)
    search_o2 = resolve_rxn(search_host, "Ex_o2[e]")
    search_o2.upper_bound = min(search_o2.upper_bound, 0.0)
    base = fba_summary(search_host, ids)
    if base.get("status") != "optimal" or base["biomass"] <= 0 or base["glucose_uptake"] >= 0:
        raise SystemExit(f"Unexpected relaxed glucose baseline: {base}")
    candidates, provenance, variant_meta = add_candidate_reactions(search_host, donors, alias_index, alias_canonical, rows)
    config_path = ROOT / "config/scenarios.json"
    fingerprint = hashlib.sha256((observed + sha256(POOL) + sha256(config_path) + sha256(Path(__file__)) +
                                  "|" + "|".join(sha256(p) for p in DONORS.values())).encode()).hexdigest()
    # Candidate reactions are present as StrainDesign KI targets. The smoke no-KI
    # model closes all of them and establishes that the 50% phenotype needs KIs.
    report_base = {"host_path": str(HOST), "host_sha256": observed, "e1_candidate_rows": len(rows),
                   "donor_files": {k: str(v) for k, v in DONORS.items()}, "search_candidate_reaction_variants": len(candidates),
                   "baseline_original_FIM_WT": original_baseline, "baseline_relaxed_CO2": base,
                   "initial_relaxed_baseline_before_uptake_caps": initial_search_baseline,
                   "solver_seed": 1,
                   "search_bound_policy": {"glucose": [search_glucose.lower_bound, search_glucose.upper_bound],
                                           "oxygen": [search_o2.lower_bound, search_o2.upper_bound],
                                           "basis": "glucose uptake capped at pre-cap relaxed baseline; oxygen exchange cannot produce O2"},
                   "reference_rubisco_flux": base["rubisco_flux"], "reference_biomass": base["biomass"]}
    atomic_json(OUT / "checkpoints/input_manifest.json", report_base)
    if args.smoke:
        cons = scenario_constraints(ids, base, 0.5)
        no_ki_bounds = {rid: (0.0, 0.0) for rid in candidates}
        no_ki_model = search_host.copy(); set_cbb_scenario_bounds(no_ki_model, ids, base, 0.5)
        no_ki_result = phenotype_fba(no_ki_model, ids, cons, no_ki_bounds)
        infeasible = no_ki_result.get("status") != "optimal"
        full_pool_model = search_host.copy(); set_cbb_scenario_bounds(full_pool_model, ids, base, 0.5)
        full_pool_result = phenotype_fba(full_pool_model, ids, cons)
        full_pool_off = glucose_off_max_biomass(full_pool_model, ids)
        viability_threshold = 0.8 * base["biomass"]
        glucose_use_proof = full_pool_off.get("max_biomass") is not None and full_pool_off["max_biomass"] < viability_threshold
        sys.path.insert(0, str(STRAINDESIGN))
        from straindesign import fba as sd_fba
        sd_full_pool_result = sd_fba(full_pool_model, constraints=cons, solver="scip")
        semantics = smoke_sd_ki_semantics()
        # Actual cost-1 search on the full E1 optional pool verifies that KI
        # switches are formulated by StrainDesign and not free in the phenotype.
        try:
            if full_pool_result.get("status") != "optimal":
                result = {"no_KI_phenotype_status": no_ki_result.get("status"), "no_KI_phenotype": no_ki_result,
                          "no_KI_infeasible": infeasible, "full_pool_phenotype": full_pool_result,
                          "full_pool_glucose_off": full_pool_off, "viability_threshold": viability_threshold,
                          "glucose_use_proof": glucose_use_proof,
                          "straindesign_full_pool_status": str(sd_full_pool_result.status),
                          "KI_semantics_positive_control": semantics,
                          "test_SD_status": "not_run_full_pool_proven_infeasible", "solutions_found": 0,
                          "ki_cost_entries": len(candidates), "ko_cost": {}, "input_fingerprint": fingerprint}
            else:
                smoke_model = search_host.copy(); set_cbb_scenario_bounds(smoke_model, ids, base, 0.5)
                sols, elapsed = sd_solve(smoke_model, candidates, cons, 1, min(args.time_limit, 90), 1)
                status = getattr(sols, "status", "unknown")
                designs = summarize_solutions(sols, provenance)
                result = {"no_KI_phenotype_status": no_ki_result.get("status"), "no_KI_phenotype": no_ki_result,
                          "no_KI_infeasible": infeasible, "full_pool_phenotype": full_pool_result,
                          "full_pool_glucose_off": full_pool_off, "viability_threshold": viability_threshold,
                          "glucose_use_proof": glucose_use_proof,
                          "straindesign_full_pool_status": str(sd_full_pool_result.status),
                          "KI_semantics_positive_control": semantics,
                          "test_SD_status": str(status), "test_runtime_seconds": elapsed,
                          "ki_cost_entries": len(candidates), "ko_cost": {}, "solutions_found": len(designs), "solutions": designs,
                          "api": "SDModule('protect') + compute_strain_designs(solver='scip')",
                          "input_fingerprint": fingerprint}
        except Exception as exc:
            result = {"no_KI_phenotype_status": no_ki_result.get("status"), "no_KI_phenotype": no_ki_result,
                      "no_KI_infeasible": infeasible, "full_pool_phenotype": full_pool_result,
                      "full_pool_glucose_off": full_pool_off, "viability_threshold": viability_threshold,
                      "glucose_use_proof": glucose_use_proof,
                      "straindesign_full_pool_status": str(sd_full_pool_result.status),
                      "KI_semantics_positive_control": semantics,
                      "test_SD_status": "error", "error": repr(exc), "input_fingerprint": fingerprint}
        result["smoke_pass"] = bool(infeasible and full_pool_result.get("status") == "optimal" and
                                    glucose_use_proof and semantics.get("passed"))
        atomic_json(OUT / "checkpoints/smoke.json", result)
        print(json.dumps({**report_base, "smoke": result}, indent=2))
        return
    if not args.run:
        print(json.dumps(report_base, indent=2)); return

    scenario_defs = [("reference", None), ("rubisco_75", .75), ("rubisco_50", .50), ("rubisco_25", .25)]
    result_rows, ki_rows = [], []
    for scenario, factor in scenario_defs:
        constraints = scenario_constraints(ids, base, factor)
        if scenario == "reference":
            reference_model = search_host.copy()
            set_cbb_scenario_bounds(reference_model, ids, base, None)
            phenotype = phenotype_fba(reference_model, ids, constraints, {rid: (0.0, 0.0) for rid in candidates})
            off = glucose_off_max_biomass(reference_model, ids, {rid: (0.0, 0.0) for rid in candidates})
            full_pool_off = glucose_off_max_biomass(reference_model, ids)
            viability_threshold = 0.8 * base["biomass"]
            data = {"scenario": scenario, "constraint_factor": None, "max_cost": 0, "status": phenotype["status"],
                    "runtime_seconds": 0.0, "solutions_found": 1, "phenotype": phenotype,
                    "glucose_off_max_biomass": off, "full_pool_glucose_off": full_pool_off,
                    "glucose_viability_threshold": viability_threshold,
                    "glucose_use_proof_pass": full_pool_off.get("max_biomass") is not None and full_pool_off["max_biomass"] < viability_threshold,
                    "phenotype_constraints": constraints, "input_fingerprint": fingerprint, "solutions": []}
            data["host_sha256"] = observed
            atomic_json(OUT / "checkpoints/reference_cost0.json", data)
            result_rows.append({"scenario": scenario, "constraint_factor": None, "rubisco_cap": None, "max_cost": 0, "status": data["status"],
                                "runtime_seconds": 0.0, "solutions_found": 1, **{k: phenotype.get(k) for k in ("biomass", "glucose_uptake", "rubisco_flux", "co2_flux", "fe2_flux", "o2_flux")},
                                "full_pool_glucose_off_max_biomass": full_pool_off.get("max_biomass"),
                                "glucose_viability_threshold": viability_threshold,
                                "glucose_use_proof_pass": data["glucose_use_proof_pass"],
                                "solution_glucose_off_max_biomass": off.get("max_biomass"),
                                "host_sha256": observed,
                                "biomass_glucose_on_off_delta": phenotype.get("biomass", 0) - (off.get("max_biomass") or 0)})
            ki_rows.append({"scenario": scenario, "minimum_KI_count": 0, "solution_id": "S00", "E1_candidate_ids": "",
                            "donor_source_reactions": "", "solver_status": data["status"], "search_reaction_ids": "",
                            "max_cost_tested": 0, "host_sha256": observed, **{k: phenotype.get(k) for k in
                            ("biomass", "glucose_uptake", "rubisco_flux", "co2_flux", "fe2_flux", "o2_flux")},
                            "glucose_off_max_biomass": off.get("max_biomass"),
                            "biomass_glucose_on_off_delta": phenotype.get("biomass", 0) - (off.get("max_biomass") or 0)})
            write_tsv(OUT / "scenario_results.tsv", result_rows); write_tsv(OUT / "minimal_KI_sets.tsv", ki_rows)
            continue
        preflight_model = search_host.copy(); set_cbb_scenario_bounds(preflight_model, ids, base, factor)
        full_pool_phenotype = phenotype_fba(preflight_model, ids, constraints)
        full_pool_off = glucose_off_max_biomass(preflight_model, ids)
        viability_threshold = 0.8 * base["biomass"]
        off_biomass = full_pool_off.get("max_biomass")
        off_proof_pass = full_pool_off.get("status") in {"optimal", "infeasible"} and off_biomass is not None and off_biomass < viability_threshold
        if full_pool_phenotype.get("status") != "optimal" or not off_proof_pass:
            preflight_status = ("infeasible_full_pool_lp" if full_pool_phenotype.get("status") != "optimal"
                                else "preflight_failed_glucose_off_not_below_threshold")
            for cost in range(1, 7):
                data = {"scenario": scenario, "constraint_factor": factor, "max_cost": cost,
                        "status": preflight_status, "runtime_seconds": 0.0,
                        "solutions_found": 0, "solutions": [], "phenotype_constraints": constraints,
                        "full_pool_phenotype": full_pool_phenotype, "full_pool_glucose_off": full_pool_off,
                        "glucose_viability_threshold": viability_threshold, "glucose_use_proof_pass": off_proof_pass,
                        "preflight_error": ("Full-pool glucose-off maximum is not below the biomass threshold; monotonic glucose-use proof is unavailable. Requires a coupled protect/suppress formulation." if not off_proof_pass else None),
                        "phenotype": None, "input_fingerprint": fingerprint, "host_sha256": observed, "error": None}
                atomic_json(OUT / "checkpoints" / f"{scenario}_cost{cost}.json", data)
                preflight_time = datetime.now(timezone.utc).isoformat()
                append_jsonl(OUT / "logs/run_log.jsonl", {"scenario": scenario, "max_cost": cost,
                             "started_at_utc": preflight_time, "ended_at_utc": preflight_time,
                             "solver_status": preflight_status, "runtime_seconds": 0.0,
                             "exception": None, "checkpoint": str(OUT / "checkpoints" / f"{scenario}_cost{cost}.json")})
                result_rows.append({"scenario": scenario, "constraint_factor": factor, "max_cost": cost,
                                    "status": data["status"], "runtime_seconds": 0.0, "solutions_found": 0,
                                    "biomass": None, "glucose_uptake": None, "rubisco_flux": None,
                                    "co2_flux": None, "fe2_flux": None, "o2_flux": None,
                                    "rubisco_cap": factor * base["rubisco_flux"],
                                    "full_pool_glucose_off_max_biomass": off_biomass,
                                    "glucose_viability_threshold": viability_threshold,
                                    "glucose_use_proof_pass": off_proof_pass, "host_sha256": observed})
            write_tsv(OUT / "scenario_results.tsv", result_rows); write_tsv(OUT / "minimal_KI_sets.tsv", ki_rows)
            continue
        found_min = False
        for cost in range(1, 7):
            cp = OUT / "checkpoints" / f"{scenario}_cost{cost}.json"
            data = None
            if args.resume and cp.exists():
                cached = json.loads(cp.read_text(encoding="utf-8"))
                retryable = {"time_limit", "time_limit_w_sol", "time_limit_w_sols", "time_limit_no_sol", "error", "validation_error"}
                if cached.get("input_fingerprint") == fingerprint and cached.get("status") not in retryable:
                    data = cached
            if data is None:
                started_at = datetime.now(timezone.utc).isoformat()
                active_model = search_host.copy()
                set_cbb_scenario_bounds(active_model, ids, base, factor)
                t0 = time.time()
                try:
                    sols, runtime = sd_solve(active_model, candidates, constraints, cost, args.time_limit, 3)
                    status = str(getattr(sols, "status", "unknown"))
                    designs = summarize_solutions(sols, provenance)
                    phenotypes = []
                    off_phenotypes = []
                    if designs:
                        for design in designs:
                            selected = set(design["reaction_ids"])
                            candidate_bounds = {rid: (active_model.reactions.get_by_id(rid).lower_bound,
                                                      active_model.reactions.get_by_id(rid).upper_bound) if rid in selected else (0.0, 0.0)
                                                for rid in candidates}
                            phenotypes.append(phenotype_fba(active_model, ids, constraints, candidate_bounds))
                            off_phenotypes.append(glucose_off_max_biomass(active_model, ids, candidate_bounds))
                    phenotype = phenotypes[0] if phenotypes else None
                    for design, pheno, off_pheno in zip(designs, phenotypes, off_phenotypes):
                        design["phenotype"] = pheno
                        design["glucose_off_max_biomass"] = off_pheno
                        design["biomass_glucose_on_off_delta"] = (pheno.get("biomass", 0) - off_pheno.get("max_biomass", 0)
                                                                  if pheno.get("status") == "optimal" and off_pheno.get("max_biomass") is not None else None)
                        design["phenotype_valid"] = (pheno.get("status") == "optimal" and
                                                      off_pheno.get("status") in {"optimal", "infeasible"} and
                                                      off_pheno.get("max_biomass") is not None and
                                                      off_pheno["max_biomass"] < 0.8 * base["biomass"])
                    valid_designs = [d for d in designs if d.get("phenotype_valid")]
                    if status == "optimal" and valid_designs:
                        for design in valid_designs:
                            design["selected_KI_count"] = len(design["reaction_ids"])
                        min_count = min(d["selected_KI_count"] for d in valid_designs)
                        if min_count > cost:
                            status = "validation_error"
                            valid_designs = []
                        else:
                            valid_designs = [d for d in valid_designs if d["selected_KI_count"] == min_count]
                            for design in valid_designs:
                                design["minimum_KI_count"] = min_count
                    if status == "optimal" and designs and not valid_designs:
                        status = "validation_error"
                    designs = valid_designs if status in {"optimal", "validation_error"} else designs
                    phenotype = designs[0].get("phenotype") if designs else phenotype
                    data = {"scenario": scenario, "constraint_factor": factor, "max_cost": cost, "status": status,
                            "runtime_seconds": runtime, "solutions_found": len(designs), "solutions": designs,
                            "phenotype_constraints": constraints, "full_pool_phenotype": full_pool_phenotype,
                            "full_pool_glucose_off": full_pool_off, "glucose_viability_threshold": viability_threshold,
                            "glucose_use_proof_pass": off_proof_pass,
                            "phenotype": phenotype, "input_fingerprint": fingerprint, "host_sha256": observed, "error": None}
                except Exception as exc:
                    data = {"scenario": scenario, "constraint_factor": factor, "max_cost": cost,
                            "status": "error", "runtime_seconds": time.time()-t0, "solutions_found": 0,
                            "solutions": [], "phenotype_constraints": constraints, "phenotype": None,
                            "input_fingerprint": fingerprint, "error": repr(exc)}
                atomic_json(cp, data)
                append_jsonl(OUT / "logs/run_log.jsonl", {"scenario": scenario, "max_cost": cost,
                             "started_at_utc": started_at, "ended_at_utc": datetime.now(timezone.utc).isoformat(),
                             "solver_status": data.get("status"), "runtime_seconds": data.get("runtime_seconds"),
                             "exception": data.get("error"), "checkpoint": str(cp)})
            rr = {k: data.get(k) for k in ("scenario", "constraint_factor", "max_cost", "status", "runtime_seconds", "solutions_found")}
            for k in ("biomass", "glucose_uptake", "rubisco_flux", "co2_flux", "fe2_flux", "o2_flux"):
                rr[k] = (data.get("phenotype") or {}).get(k)
            rr["host_sha256"] = observed
            rr["rubisco_cap"] = None if factor is None else factor * base["rubisco_flux"]
            rr["full_pool_glucose_off_max_biomass"] = (data.get("full_pool_glucose_off") or {}).get("max_biomass")
            rr["glucose_viability_threshold"] = data.get("glucose_viability_threshold", 0.8 * base["biomass"])
            rr["glucose_use_proof_pass"] = data.get("glucose_use_proof_pass")
            first_design = (data.get("solutions") or [{}])[0]
            rr["solution_glucose_off_max_biomass"] = (first_design.get("glucose_off_max_biomass") or {}).get("max_biomass")
            rr["biomass_glucose_on_off_delta"] = first_design.get("biomass_glucose_on_off_delta")
            result_rows.append(rr)
            write_tsv(OUT / "scenario_results.tsv", result_rows)
            write_tsv(OUT / "minimal_KI_sets.tsv", ki_rows)
            if data["solutions_found"] and data["status"] == "optimal" and not found_min:
                found_min = True
                for design in data["solutions"]:
                    ki_rows.append({"scenario": scenario, "minimum_KI_count": len(design["reaction_ids"]),
                                    "solution_id": design["solution_id"], "E1_candidate_ids": ";".join(design["candidate_ids"]),
                                    "donor_source_reactions": ";".join(design["source_reactions"]), "solver_status": data["status"],
                                    "search_reaction_ids": ";".join(design["reaction_ids"]), "max_cost_tested": cost,
                                    "host_sha256": observed, **{k: (design.get("phenotype") or {}).get(k) for k in
                                    ("biomass", "glucose_uptake", "rubisco_flux", "co2_flux", "fe2_flux", "o2_flux")},
                                    "glucose_off_max_biomass": (design.get("glucose_off_max_biomass") or {}).get("max_biomass"),
                                    "biomass_glucose_on_off_delta": design.get("biomass_glucose_on_off_delta")})
                break
            if data["status"] not in {"optimal", "infeasible"}:
                break
        # Rewrite checkpoint-backed tabular summaries after every scenario.
        write_tsv(OUT / "scenario_results.tsv", result_rows)
        write_tsv(OUT / "minimal_KI_sets.tsv", ki_rows)
    write_report(report_base, result_rows, ki_rows, provenance)
    print(json.dumps({"scenarios": result_rows, "minimum_sets": ki_rows}, indent=2))


def write_tsv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def append_jsonl(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")


def write_report(manifest, result_rows, ki_rows, provenance):
    lines = ["# E2 KI-only StrainDesign search", "", f"- Authoritative host: `{manifest['host_path']}`",
             f"- Host SHA256: `{manifest['host_sha256']}`", f"- E1 candidates loaded: {manifest['e1_candidate_rows']}",
             f"- Search candidate direction variants: {manifest['search_candidate_reaction_variants']}",
             f"- Original FIM+WT baseline: `{json.dumps(manifest['baseline_original_FIM_WT'])}`",
             f"- Relaxed CO2 search baseline: `{json.dumps(manifest['baseline_relaxed_CO2'])}`", "",
             f"StrainDesign was called through `SDModule('protect')` and `compute_strain_designs` with SCIP, explicit `ko_cost={{}}`, unit `ki_cost` for each optional source-supported candidate reaction, and fixed seed `{manifest.get('solver_seed', 1)}`. Opposite one-way donor directions are separate cost-1 variants; a reversible source supports one reversible KI. Candidate sources and bounds are reconstructed from the three local donor GEMs. The search caps glucose uptake at the relaxed host baseline and disallows O2 secretion. A small positive Fe2 uptake guard preserves the Fe oxidation condition without setting a high Fe2 fraction. The model-default GLPK precheck was bypassed because it reported infeasible on the large network; StrainDesign rechecks feasibility with selected SCIP.", "",
             "## Search results", "", "| Scenario | Max cost | Status | Runtime (s) | Solutions | Rubisco cap | Biomass | Glucose | Rubisco | CO2 | Fe2 | O2 | Full-pool glucose-off max | Glucose threshold | Proof | KI-set glucose-off max | On-off delta |",
             "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|"]
    for r in result_rows:
        lines.append("| " + " | ".join(str(r.get(k, "")) for k in
                      ("scenario", "max_cost", "status", "runtime_seconds", "solutions_found", "rubisco_cap", "biomass",
                       "glucose_uptake", "rubisco_flux", "co2_flux", "fe2_flux", "o2_flux", "full_pool_glucose_off_max_biomass",
                       "glucose_viability_threshold", "glucose_use_proof_pass", "solution_glucose_off_max_biomass",
                       "biomass_glucose_on_off_delta")) + " |")
    for scenario in ("rubisco_75", "rubisco_50", "rubisco_25"):
        states = [r for r in result_rows if r.get("scenario") == scenario]
        selected = [r for r in ki_rows if r.get("scenario") == scenario]
        if selected:
            lines.append(f"- `{scenario}`: feasible; minimum confirmed at {selected[0]['minimum_KI_count']} KIs by an optimal StrainDesign solve.")
            first = states[0] if states else {}
            lines.append(f"  Full-pool glucose-off maximum biomass={first.get('full_pool_glucose_off_max_biomass')}, viability threshold={first.get('glucose_viability_threshold')}; all selected solutions are checked with glucose closed.")
        elif states and all(r.get("status") in {"infeasible", "infeasible_full_pool_lp"} for r in states):
            lines.append(f"- `{scenario}`: infeasible through 6 KIs (full-pool LP proof where status is `infeasible_full_pool_lp`).")
        elif states and all(r.get("status") == "preflight_failed_glucose_off_not_below_threshold" for r in states):
            lines.append(f"- `{scenario}`: preflight failed; full-pool glucose-off biomass did not fall below the viability threshold, so glucose dependence was not proven. A coupled protect/suppress formulation would be needed.")
        elif states:
            lines.append(f"- `{scenario}`: unresolved; last completed status `{states[-1].get('status')}`. A timeout is not proof of infeasibility or minimum size.")
    lines += ["", "## Minimum KI sets", ""]
    if ki_rows:
        for row in ki_rows:
            lines.append(f"- **{row['scenario']} ({row['minimum_KI_count']} KI):** {row['E1_candidate_ids']} — {row['donor_source_reactions']} (status `{row['solver_status']}`).")
            lines.append(f"  Phenotype: biomass={row.get('biomass')}, glucose={row.get('glucose_uptake')}, Rubisco={row.get('rubisco_flux')}, CO2={row.get('co2_flux')}, Fe2={row.get('fe2_flux')}, O2={row.get('o2_flux')}; glucose-off max biomass={row.get('glucose_off_max_biomass')}, biomass on-off delta={row.get('biomass_glucose_on_off_delta')}.")
    else:
        lines.append("No minimum set was confirmed in the completed runs.")
    lines += ["", "Scenario checkpoints under `checkpoints/` are authoritative for interrupted or unresolved runs. `time_limit` and `error` statuses are retained and must not be interpreted as proven infeasibility.", ""]
    (OUT / "E2_report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
