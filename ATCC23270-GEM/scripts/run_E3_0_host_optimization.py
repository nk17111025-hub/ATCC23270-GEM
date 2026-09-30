#!/usr/bin/env python
"""E3.0 whole-host KO/regulatory optimization, independently per E2R3 branch."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import re
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib.metadata import version as package_version
from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import add_pfba

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/engineering/E3_0_host_optimization"
HOST = ROOT.parent / "before final/updatedv3.1_verified/updatedv3.1.xml"
POOL = ROOT / "results/engineering/E1_donor_reaction_pool/candidate_reactions.tsv"
E2_PATH = ROOT / "scripts/run_E2_strain_design_KI.py"
E2R_PATH = ROOT / "scripts/run_E2R_resource_matched_KI.py"
E2R3_PATH = ROOT / "scripts/run_E2R3_two_stage_full_pool.py"
SUPPORT_JSON = ROOT / "results/engineering/E2R3_two_stage_full_pool/checkpoints/support_sets.json"
E2R3_MANIFEST = ROOT / "results/engineering/E2R3_two_stage_full_pool/checkpoints/input_manifest.json"
SD_PATCH_FILES = (ROOT / "external/straindesign/straindesign/solver_interface.py",
                  ROOT / "external/straindesign/straindesign/strainDesignMILP.py")
EXPECTED_HOST = "026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3"
EXPECTED_SUPPORT = {"S01": {"E1R01016", "E1R01493", "E1R01642"},
                    "S03": {"E1R00732", "E1R00741"}}
IDS = {"biomass": "Ex_bio[e]", "glucose": "Ex_glc-B[e]", "rubisco": "RUBISCO",
       "co2": "Ex_h2co3[e]", "fe2": "Ex_fe2[e]", "o2": "Ex_o2[e]"}
BRANCHES = ("S01", "S03")
MODES = {
    "CBB_min": ("rubisco", "minimize"),
    "Fe2_efficiency": ("fe2", "maximize"),
    "O2_efficiency": ("o2", "maximize"),
    "CO2_efficiency": ("co2", "maximize"),
}
BUDGETS = (1, 2, 4, 6)
TOL = 1e-8
MAX_TIME = 300
MAX_SOLUTIONS = 5


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def save_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def load_ids(model, e2):
    return {k: e2.resolve_rxn(model, v) for k, v in IDS.items()}


def load_supports():
    payload = json.loads(SUPPORT_JSON.read_text(encoding="utf-8"))
    matches = {branch: None for branch in BRANCHES}
    # Support checkpoint stores strategy rows. The human-readable TSV is the
    # source for exact E1 IDs; verify these against deterministic checkpoint data.
    support_tsv = ROOT / "results/engineering/E2R3_two_stage_full_pool/stage1_support_sets.tsv"
    with support_tsv.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sid = row.get("support_set_id")
            if sid in matches:
                matches[sid] = row
    out = {}
    for branch, expected in EXPECTED_SUPPORT.items():
        row = matches.get(branch)
        if not row:
            raise RuntimeError(f"E2R3 support set missing: {branch}")
        ids = {x for x in row.get("E1_candidate_IDs", row.get("E1_candidate_ids", "")).split(";") if x}
        if ids != expected:
            raise RuntimeError(f"E2R3 {branch} candidates mismatch: {sorted(ids)}")
        search_ids = {x for x in row.get("search_reaction_IDs", row.get("search_reaction_ids", "")).split(";") if x}
        if not search_ids:
            # Older E2R3 TSV header spelling fallback; do not infer by candidate ID.
            search_ids = ids.copy()
        out[branch] = {"candidate_ids": sorted(ids), "search_ids": sorted(search_ids),
                       "source_reactions": row.get("donor_source_reactions", row.get("donor/source_reaction_ids", "")),
                       "checkpoint_fingerprint": payload.get("fingerprint")}
    return out


def apply_environment(model, ids, ref):
    ids["fe2"].bounds = (-float(ref["Fe_WT"]), 0.0)
    ids["o2"].bounds = (-float(ref["O2_WT"]), 0.0)
    ids["glucose"].bounds = (-0.5, 0.0)
    ids["co2"].bounds = (-2.0, 0.0)
    ids["rubisco"].lower_bound = 0.0
    ids["rubisco"].upper_bound = min(ids["rubisco"].upper_bound, float(ref["Stage2"]["rubisco_max"]))
    model.add_cons_vars(model.problem.Constraint(ids["biomass"].flux_expression,
        lb=float(ref["Stage2"]["biomass_min"]), name="E3_biomass_min"))
    model.add_cons_vars(model.problem.Constraint(ids["glucose"].flux_expression,
        ub=-0.05, name="E3_glucose_min"))
    model.add_cons_vars(model.problem.Constraint(ids["fe2"].flux_expression,
        ub=-0.01, name="E3_trace_fe2"))


def is_pseudo(rxn, biomass_id, fixed_ids):
    if rxn.id in fixed_ids or rxn.id == biomass_id or rxn.boundary:
        return True
    name = f"{rxn.id} {rxn.name}".lower()
    return bool(re.search(r"(^|[_\s])(biomass|demand|sink|exchange|pseudo)([_\s]|$)", name))


def read_support_model(branch, e2, e2r, config, rows, donors, aliases_i, aliases_c, ref):
    host, _ = cobra.io.validate_sbml_model("../before final/updatedv3.1_verified/updatedv3.1.xml")
    if host is None:
        raise RuntimeError("Could not load authoritative host")
    e2r.native_condition(host, config, "FIM", e2)
    ids = load_ids(host, e2)
    native_host_ids = {r.id for r in host.reactions}
    host.notes["e2r3_fe_wt"] = float(ref["Fe_WT"])
    host.notes["e2r3_o2_wt"] = float(ref["O2_WT"])
    candidates, provenance, variants = e2.add_candidate_reactions(host, donors, aliases_i, aliases_c, rows)
    support = load_supports()[branch]
    if not set(support["search_ids"]).issubset(candidates):
        raise RuntimeError(f"E2R3 {branch} search reaction IDs could not be reconstructed")
    for rid in candidates:
        reaction = host.reactions.get_by_id(rid)
        if rid not in support["search_ids"]:
            reaction.bounds = (0.0, 0.0)
    apply_environment(host, ids, ref)
    return host, ids, support, candidates, provenance, variants, native_host_ids


def branch_reference(model, ids):
    model.objective = ids["biomass"]
    optimum = model.optimize()
    if optimum.status != "optimal":
        raise RuntimeError(f"branch biomass reference failed: {optimum.status}")
    with model:
        model.objective = ids["biomass"]
        pf = cobra.flux_analysis.pfba(model, fraction_of_optimum=1.0)
    return optimum, pf


def native_reaction_ids(model, ids, fixed_ids):
    return sorted(r.id for r in model.reactions if not is_pseudo(r, ids["biomass"].id, fixed_ids))


def generate_reg_candidates(model, ids, native_ids, ref_fluxes, bmin, rcap):
    # Flux variability is evaluated only inside the branch's viable Stage-2 phenotype.
    with model:
        model.add_cons_vars(model.problem.Constraint(ids["biomass"].flux_expression, lb=bmin,
                          name="E3_reg_fva_biomass"))
        model.add_cons_vars(model.problem.Constraint(ids["glucose"].flux_expression, ub=-0.05,
                          name="E3_reg_fva_glucose"))
        ranges = flux_variability_analysis(model, reaction_list=native_ids, fraction_of_optimum=0.0,
                                           loopless=None)
    reg = {}
    details = []
    for rid in native_ids:
        v = float(ref_fluxes[rid]); lo = float(ranges.loc[rid, "minimum"]); hi = float(ranges.loc[rid, "maximum"])
        directions = []
        if v > TOL and hi > v + TOL:
            target = 1.5 * v if 1.5 * v <= hi - TOL else v + 0.8 * (hi - v)
            if target > v + TOL:
                directions.append((f"{rid} >= {target:.12g}", "UPREG", target))
        elif v < -TOL and lo < v - TOL:
            target = 1.5 * v if 1.5 * v >= lo + TOL else v - 0.8 * (v - lo)
            if target < v - TOL:
                directions.append((f"{rid} <= {target:.12g}", "UPREG", target))
        elif abs(v) <= TOL:
            span = hi - lo
            if hi > TOL and span > TOL:
                target = max(TOL * 10, 0.1 * max(0.0, hi))
                if target < hi - TOL:
                    directions.append((f"{rid} >= {target:.12g}", "ACTIVATION_UPREG", target))
            if lo < -TOL and span > TOL:
                target = min(-TOL * 10, 0.1 * min(0.0, lo))
                if target > lo + TOL:
                    directions.append((f"{rid} <= {target:.12g}", "ACTIVATION_UPREG", target))
        for constraint, kind, target in directions:
            reg[constraint] = 1.0
            details.append({"native_reaction_id": rid, "intervention_type": kind,
                "reference_flux": v, "fva_minimum": lo, "fva_maximum": hi,
                "regulatory_constraint": constraint, "intervention_cost": 1.0})
    return reg, ranges, details


def candidate_rows(branch, model, ids, native_ids, reg_details, fva, ref_fluxes):
    rows = []
    for rid in native_ids:
        rows.append({"branch": branch, "native_reaction_id": rid, "intervention_type": "KO",
            "reference_flux": ref_fluxes[rid], "fva_minimum": float(fva.loc[rid,"minimum"]),
            "fva_maximum": float(fva.loc[rid,"maximum"]), "regulatory_constraint": "",
            "intervention_cost": 1.0})
    rows.extend({"branch": branch, **x} for x in reg_details)
    return rows


def sd_api():
    # E2 helper resolves the project's vendored StrainDesign package and imports it.
    e2 = load_module(E2_PATH, "e2_helpers_e3")
    sys.path.insert(0, str(e2.STRAINDESIGN))
    from straindesign import SDModule, compute_strain_designs
    return SDModule, compute_strain_designs


def solve_mode(model, ids, mode, ko_cost, reg_cost, max_cost, time_limit,
               compress=True, skip_preprocessing_fvas=False):
    SDModule, compute_strain_designs = sd_api()
    rxn_key, sense = MODES[mode]
    constraint_strs = [f"{ids['biomass'].id} >= {model.notes['e3_bmin']:.12g}",
        f"{ids['glucose'].id} <= -0.05", f"{ids['fe2'].id} <= -0.01",
        f"{ids['rubisco'].id} >= 0", f"{ids['rubisco'].id} <= {model.notes['e3_r10']:.12g}"]
    module = SDModule(model, "robustknock", inner_objective={ids["biomass"].id: 1.0},
        inner_opt_sense="maximize", outer_objective={ids[rxn_key].id: 1.0},
        outer_opt_sense=sense, constraints=constraint_strs)
    return compute_strain_designs(model, sd_modules=[module], solution_approach="best",
        max_cost=max_cost, max_solutions=MAX_SOLUTIONS, ko_cost=dict(ko_cost), reg_cost=dict(reg_cost),
        solver="scip", seed=1, time_limit=time_limit, compress=compress,
        skip_preprocessing_fvas=skip_preprocessing_fvas)


def reconstruct(model, ids, support, candidates, sd, mode):
    # Apply interventions once, fix the biomass-optimal face, and calculate
    # all four worst-case coordinates so cross-mode Pareto comparisons use
    # consistent phenotype metrics.
    trial = model.copy()
    tids = {key: trial.reactions.get_by_id(reaction.id) for key, reaction in ids.items()}
    for rid, value in sd.items():
        if value == -1:
            trial.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        elif value is True and isinstance(rid, str) and (" >= " in rid or " <= " in rid or " = " in rid):
            # Only selected regulatory intervention strings are encoded True.
            left, op, right = re.split(r"\s*(>=|<=|=)\s*", rid, maxsplit=1)
            expression = trial.reactions.get_by_id(left).flux_expression
            val = float(right)
            c = (trial.problem.Constraint(expression, lb=val, name=f"E3reg_{len(trial.solver.constraints)}")
                 if op == ">=" else trial.problem.Constraint(expression, ub=val, name=f"E3reg_{len(trial.solver.constraints)}")
                 if op == "<=" else trial.problem.Constraint(expression, lb=val, ub=val, name=f"E3reg_{len(trial.solver.constraints)}"))
            trial.add_cons_vars(c)
    trial.objective = tids["biomass"]
    inner = trial.optimize()
    if inner.status != "optimal":
        return trial, inner, None
    optimum = float(inner.fluxes[tids["biomass"].id])
    trial.add_cons_vars(trial.problem.Constraint(tids["biomass"].flux_expression,
        lb=max(float(model.notes["e3_bmin"]), optimum - 1e-8), name="E3_inner_optimum_floor"))
    trial.objective = tids["biomass"]
    representative = cobra.flux_analysis.pfba(trial, fraction_of_optimum=1.0)
    robust = {}
    failures = []
    for key in ("rubisco", "fe2", "o2", "co2"):
        trial.objective = tids[key]
        trial.objective_direction = "max" if key == "rubisco" else "min"
        solution = trial.optimize()
        if solution.status == "optimal":
            robust[key] = float(solution.fluxes[tids[key].id])
        else:
            failures.append(f"{key}:{solution.status}")
    return trial, representative, optimum, robust, failures


def phenotype(sol, ids):
    if sol.status != "optimal":
        return {"status": sol.status}
    f = sol.fluxes
    return {"status": sol.status, "biomass": float(f[ids["biomass"].id]),
        "glucose": float(f[ids["glucose"].id]), "rubisco": float(f[ids["rubisco"].id]),
        "fe2": float(f[ids["fe2"].id]), "o2": float(f[ids["o2"].id]),
        "co2_h2co3": float(f[ids["co2"].id])}


def pareto(rows):
    # All values converted to desirability: biomass higher, all resource/use counts lower.
    kept = []
    for i, row in enumerate(rows):
        if row.get("solver_status") not in ("optimal", "time_limit_w_sols", "time_limit", "baseline_reference") or row.get("biomass") is None:
            continue
        if row.get("reconstruction_status") not in ("optimal", "baseline_reference"):
            continue
        a = [int(row["total_intervention_cost"]), float(row["rubisco"]), abs(float(row["fe2"])),
             abs(float(row["o2"])), abs(float(row["co2_h2co3"])), -float(row["biomass"])]
        dominated = False
        for j, other in enumerate(rows):
            if i == j or other.get("solver_status") not in ("optimal", "time_limit_w_sols", "time_limit", "baseline_reference") or other.get("biomass") is None:
                continue
            if other.get("reconstruction_status") not in ("optimal", "baseline_reference"):
                continue
            b = [int(other["total_intervention_cost"]), float(other["rubisco"]), abs(float(other["fe2"])),
                 abs(float(other["o2"])), abs(float(other["co2_h2co3"])), -float(other["biomass"])]
            if all(y <= x + 1e-12 for x,y in zip(a,b)) and any(y < x - 1e-12 for x,y in zip(a,b)):
                dominated = True; break
        if not dominated:
            kept.append(row)
    return kept


def branch_job(branch, e2, e2r, config, rows, donors, aliases_i, aliases_c, ref, time_limit, prepare_only=False):
    bdir = OUT / branch
    (bdir / "checkpoints").mkdir(parents=True, exist_ok=True)
    (bdir / "logs").mkdir(exist_ok=True)
    model, ids, support, candidates, provenance, variants, native_host_ids = read_support_model(branch,e2,e2r,config,rows,donors,aliases_i,aliases_c,ref)
    bmin=float(ref["Stage2"]["biomass_min"]); rcap=float(ref["Stage2"]["rubisco_max"])
    model.notes["e3_bmin"]=bmin; model.notes["e3_r10"]=rcap
    optimum, pf = branch_reference(model,ids)
    expected_bio={"S01":0.0795062,"S03":0.0800899}[branch]
    if abs(float(optimum.fluxes[ids["biomass"].id])-expected_bio)>2e-5:
        raise RuntimeError(f"{branch} branch reference discrepancy: {optimum.fluxes[ids['biomass'].id]}")
    fixed=set(support["search_ids"])
    native_ids=[rid for rid in native_reaction_ids(model,ids,fixed) if rid in native_host_ids]
    ref_fluxes={rid:float(pf.fluxes[rid]) for rid in native_ids}
    reg_cost, fva, reg_details=generate_reg_candidates(model,ids,native_ids,ref_fluxes,bmin,rcap)
    ko_cost={rid:1.0 for rid in native_ids}
    candidate_table=candidate_rows(branch,model,ids,native_ids,reg_details,fva,ref_fluxes)
    write_tsv(bdir/"intervention_candidates.tsv",candidate_table)
    base={"branch":branch,"biomass_optimum":float(optimum.fluxes[ids["biomass"].id]),
        "reference_pfba_biomass":float(pf.fluxes[ids["biomass"].id]),
        "glucose":float(pf.fluxes[ids["glucose"].id]),"rubisco":float(pf.fluxes[ids["rubisco"].id]),
        "fe2":float(pf.fluxes[ids["fe2"].id]),"o2":float(pf.fluxes[ids["o2"].id]),
        "co2_h2co3":float(pf.fluxes[ids["co2"].id]),"KO_candidate_count":len(ko_cost),
        "UPREG_candidate_count":sum(x["intervention_type"]=="UPREG" for x in reg_details),
        "ACTIVATION_UPREG_candidate_count":sum(x["intervention_type"]=="ACTIVATION_UPREG" for x in reg_details),
        "fixed_E1_candidate_ids":";".join(support["candidate_ids"]),
        "fixed_search_reaction_ids":";".join(support["search_ids"])}
    _, ref_pf, ref_bio, ref_robust, ref_failures = reconstruct(model,ids,support,candidates,{},"CBB_min")
    base.update({f"worst_case_{key}":value for key,value in ref_robust.items()})
    base["worst_case_validation_failures"]=";".join(ref_failures)
    write_tsv(bdir/"branch_reference.tsv",[base])
    fpdata={"branch":branch,"host":sha(HOST),"pool":sha(POOL),"E2R3_fingerprint":ref["fingerprint"],
        "support":support,"script":sha(Path(__file__)),"bmin":bmin,"r10":rcap,"ko_ids":native_ids,
        "reg_candidates":sorted(reg_cost),"modes":MODES,"budgets":BUDGETS,"seed":1,"time_limit":time_limit,
        "straindesign_patch_sha256":{p.name:sha(p) for p in SD_PATCH_FILES},
        "solver_preprocessing_fallback":{"primary":{"compress":True,"skip_preprocessing_fvas":False},
                                         "on_error":{"compress":False,"skip_preprocessing_fvas":True}}}
    fp=hashlib.sha256(json.dumps(fpdata,sort_keys=True).encode()).hexdigest()
    save_json(bdir/"checkpoints/input_manifest.json",{"fingerprint":fp,**fpdata})
    if prepare_only:
        return {"branch":branch,"fingerprint":fp,"KO_candidates":len(ko_cost),"reg_candidates":len(reg_cost)}
    all_rows=[{"branch":branch,"objective_mode":"BASELINE_REFERENCE","max_intervention_cost":0,
        "solver_status":"baseline_reference","solution_id":f"{branch}_BASELINE_C0",
        "total_intervention_cost":0,"KO_count":0,"OE_count":0,"KO_reaction_ids":"",
        "regulatory_constraints":"","biomass_inner_optimum":ref_bio,
        "biomass":float(ref_pf.fluxes[ids["biomass"].id]),
        "glucose":float(ref_pf.fluxes[ids["glucose"].id]),
        "rubisco":ref_robust.get("rubisco"),"fe2":ref_robust.get("fe2"),
        "o2":ref_robust.get("o2"),"co2_h2co3":ref_robust.get("co2"),
        "robust_validation_failures":";".join(ref_failures),
        "reconstruction_status":"optimal" if not ref_failures else "robust_validation_failed",
        "input_fingerprint":fp}]
    sd_path=bdir/"checkpoints/optimization_results.json"
    cached={}
    if sd_path.exists():
        try:
            old=json.loads(sd_path.read_text(encoding="utf-8"))
            if old.get("fingerprint")==fp:
                cached=old.get("results",{})
                for entry in cached.values():
                    if entry.get("status") not in ("error", "unknown", ""):
                        all_rows.extend(entry.get("rows",[]))
        except Exception: pass
    for mode,(rxn_key,sense) in MODES.items():
        for budget in BUDGETS:
            key=f"{mode}|{budget}"
            if key in cached and cached[key].get("status") not in ("error", "unknown", ""):
                continue
            started=time.time(); status="error"; records=[]; error=""; attempts=[]
            try:
                result=None; setting=None
                for setting in ({"compress":True,"skip_preprocessing_fvas":False},
                                {"compress":False,"skip_preprocessing_fvas":True}):
                    attempt_started=time.time()
                    try:
                        result=solve_mode(model,ids,mode,ko_cost,reg_cost,budget,time_limit,**setting)
                        attempt_status=str(getattr(result,"status","unknown"))
                        attempts.append({**setting,"solver":"scip","status":attempt_status,
                                         "runtime_seconds":time.time()-attempt_started})
                        if attempt_status not in ("error","unknown"):
                            break
                    except Exception as attempt_exc:
                        attempts.append({**setting,"solver":"scip","status":"error",
                                         "error":f"{type(attempt_exc).__name__}: {attempt_exc}",
                                         "traceback":traceback.format_exc(),
                                         "runtime_seconds":time.time()-attempt_started})
                        result=None
                if result is None:
                    raise RuntimeError("Both StrainDesign SCIP preprocessing modes failed")
                status=str(getattr(result,"status","unknown"))
                designs=getattr(result,"reaction_sd",[])
                costs=getattr(result,"sd_cost",[])
                if status=="time_limit_w_sols" and not designs:
                    status="time_limit"
                for idx,sd in enumerate(designs[:MAX_SOLUTIONS]):
                    ko_ids=sorted(r for r,v in sd.items() if v == -1 and r in ko_cost)
                    reg_ids=sorted(k for k,v in sd.items() if v is True and k in reg_cost)
                    trial,sol,bio_opt,robust,robust_failures=reconstruct(model,ids,support,candidates,sd,mode)
                    ph=phenotype(sol,ids)
                    if ph.get("status")=="optimal":
                        for metric_key,value in robust.items(): ph["co2_h2co3" if metric_key=="co2" else metric_key]=value
                        if robust_failures: ph["status"]="robust_validation_failed"
                    cost=float(costs[idx]) if idx<len(costs) else len(ko_ids)+len(reg_ids)
                    records.append({"branch":branch,"objective_mode":mode,"max_intervention_cost":budget,
                        "solver_status":status,"solution_id":f"{branch}_{mode}_C{budget}_S{idx+1}",
                        "total_intervention_cost":int(round(cost)),"KO_count":len(ko_ids),"OE_count":len(reg_ids),
                        "KO_reaction_ids":";".join(ko_ids),"regulatory_constraints":";".join(reg_ids),
                        "biomass_inner_optimum":bio_opt,
                        "robust_validation_failures":";".join(robust_failures),
                        "solver_attempts":json.dumps(attempts,sort_keys=True),
                        **{k:v for k,v in ph.items() if k!="status"},"reconstruction_status":ph["status"],
                        "runtime_seconds":time.time()-started,"input_fingerprint":fp})
                if not records:
                    records=[{"branch":branch,"objective_mode":mode,"max_intervention_cost":budget,
                        "solver_status":status,"solution_id":"","total_intervention_cost":"",
                        "KO_count":"","OE_count":"","KO_reaction_ids":"","regulatory_constraints":"",
                        "solver_attempts":json.dumps(attempts,sort_keys=True),
                            "runtime_seconds":time.time()-started,"input_fingerprint":fp}]
            except Exception as exc:
                error=f"{type(exc).__name__}: {exc}"; status="error"
                records=[{"branch":branch,"objective_mode":mode,"max_intervention_cost":budget,
                    "solver_status":status,"error":error,"traceback":traceback.format_exc(),
                    "solver_attempts":json.dumps(attempts,sort_keys=True),"runtime_seconds":time.time()-started,
                    "input_fingerprint":fp}]
            cached[key]={"status":status,"error":error,"rows":records}
            save_json(sd_path,{"fingerprint":fp,"results":cached})
            all_rows.extend(records)
            write_tsv(bdir/"optimization_results.tsv",all_rows)
    # Budget 8 is a single conditional extension: require a measured objective gain at budget 6.
    for mode in MODES:
        r6=[r for r in all_rows if r.get("objective_mode")==mode and r.get("max_intervention_cost")==6 and r.get("reconstruction_status")=="optimal"]
        r4=[r for r in all_rows if r.get("objective_mode")==mode and r.get("max_intervention_cost")==4 and r.get("reconstruction_status")=="optimal"]
        key=f"{mode}|8"
        if key in cached and cached[key].get("status") not in ("error", "unknown", ""):
            continue
        if r6 and r4 and mode+"|6" in cached:
            metric="rubisco" if mode=="CBB_min" else {"Fe2_efficiency":"fe2","O2_efficiency":"o2","CO2_efficiency":"co2_h2co3"}[mode]
            score=lambda r: float(r[metric]) if mode=="CBB_min" else abs(float(r[metric]))
            improvement=min(score(r) for r in r4)-min(score(r) for r in r6)
            # For all objectives a lower absolute metric is the requested improvement.
            if improvement > max(1e-7,1e-3*max(1.0,abs(min(score(r) for r in r4)))) and key not in cached:
                try:
                    started=time.time(); result=None; attempts=[]
                    for setting in ({"compress":True,"skip_preprocessing_fvas":False},
                                    {"compress":False,"skip_preprocessing_fvas":True}):
                        attempt_started=time.time()
                        try:
                            result=solve_mode(model,ids,mode,ko_cost,reg_cost,8,time_limit,**setting)
                            attempt_status=str(getattr(result,"status","unknown"))
                            attempts.append({**setting,"solver":"scip","status":attempt_status,
                                             "runtime_seconds":time.time()-attempt_started})
                            if attempt_status not in ("error","unknown"):
                                break
                        except Exception as attempt_exc:
                            attempts.append({**setting,"solver":"scip","status":"error",
                                             "error":f"{type(attempt_exc).__name__}: {attempt_exc}",
                                             "traceback":traceback.format_exc(),
                                             "runtime_seconds":time.time()-attempt_started})
                            result=None
                    if result is None:
                        raise RuntimeError("Both StrainDesign SCIP preprocessing modes failed")
                    status=str(getattr(result,"status","unknown")); records=[]
                    designs=getattr(result,"reaction_sd",[])
                    if status=="time_limit_w_sols" and not designs:
                        status="time_limit"
                    for idx,sd in enumerate(designs[:MAX_SOLUTIONS]):
                        ko_ids=sorted(r for r,v in sd.items() if v == -1 and r in ko_cost)
                        reg_ids=sorted(k for k,v in sd.items() if v is True and k in reg_cost)
                        trial,sol,bio_opt,robust,robust_failures=reconstruct(model,ids,support,candidates,sd,mode); ph=phenotype(sol,ids)
                        if ph.get("status")=="optimal":
                            for metric_key,value in robust.items(): ph["co2_h2co3" if metric_key=="co2" else metric_key]=value
                            if robust_failures: ph["status"]="robust_validation_failed"
                        records.append({"branch":branch,"objective_mode":mode,"max_intervention_cost":8,
                            "solver_status":status,"solution_id":f"{branch}_{mode}_C8_S{idx+1}",
                            "total_intervention_cost":len(ko_ids)+len(reg_ids),"KO_count":len(ko_ids),"OE_count":len(reg_ids),
                            "KO_reaction_ids":";".join(ko_ids),"regulatory_constraints":";".join(reg_ids),
                            "biomass_inner_optimum":bio_opt,
                            "robust_validation_failures":";".join(robust_failures),
                            "solver_attempts":json.dumps(attempts,sort_keys=True),
                            **{k:v for k,v in ph.items() if k!="status"},"reconstruction_status":ph["status"],
                            "runtime_seconds":time.time()-started,"input_fingerprint":fp})
                    if not records:
                        records=[{"branch":branch,"objective_mode":mode,"max_intervention_cost":8,
                                  "solver_status":status,"solver_attempts":json.dumps(attempts,sort_keys=True),
                                  "runtime_seconds":time.time()-started,"input_fingerprint":fp}]
                    cached[key]={"status":status,"rows":records}; save_json(sd_path,{"fingerprint":fp,"results":cached})
                    all_rows.extend(records)
                except Exception as exc:
                    records=[{"branch":branch,"objective_mode":mode,"max_intervention_cost":8,
                        "solver_status":"error","error":f"{type(exc).__name__}: {exc}",
                        "traceback":traceback.format_exc(),"solver_attempts":json.dumps(attempts,sort_keys=True),
                        "input_fingerprint":fp}]
                    cached[key]={"status":"error","rows":records}; save_json(sd_path,{"fingerprint":fp,"results":cached}); all_rows.extend(records)
    write_tsv(bdir/"optimization_results.tsv",all_rows)
    front=pareto(all_rows)
    # Glucose-off only for final nondominated reconstructed designs.
    final=[]
    for row in front:
        design_key=row.get("solution_id")
        sdrow=next((x for x in all_rows if x.get("solution_id")==design_key),None)
        if not sdrow: continue
        ko_ids=[x for x in sdrow.get("KO_reaction_ids","").split(";") if x]
        reg_ids=[x for x in sdrow.get("regulatory_constraints","").split(";") if x]
        trial=model.copy()
        tids={key:trial.reactions.get_by_id(reaction.id) for key,reaction in ids.items()}
        for rid in ko_ids: trial.reactions.get_by_id(rid).bounds=(0,0)
        for rule in reg_ids:
            left,op,right=re.split(r"\s*(>=|<=|=)\s*",rule,maxsplit=1); expr=trial.reactions.get_by_id(left).flux_expression; val=float(right)
            cons=trial.problem.Constraint(expr,lb=val,name=f"E3final_{len(trial.solver.constraints)}") if op==">=" else trial.problem.Constraint(expr,ub=val,name=f"E3final_{len(trial.solver.constraints)}")
            trial.add_cons_vars(cons)
        # Biomass-on state is the already reconstructed result stored in row.
        with trial:
            tids["glucose"].bounds=(0,0)
            trial.remove_cons_vars([c for c in trial.solver.constraints if c.name in ("E3_glucose_min","E3_biomass_min")])
            trial.objective=tids["biomass"]; off=trial.optimize()
        row["glucose_off_status"]=off.status
        row["glucose_off_biomass"]=float(off.fluxes[tids["biomass"].id]) if off.status=="optimal" else None
        row["glucose_contribution_delta"]=(float(row.get("biomass",0))-float(row["glucose_off_biomass"])
            if row["glucose_off_status"]=="optimal" and row.get("biomass") is not None else None)
        final.append(row)
    write_tsv(bdir/"pareto_designs.tsv",final)
    return {"branch":branch,"fingerprint":fp,"KO_candidates":len(ko_cost),"reg_candidates":len(reg_cost),"pareto_designs":len(final)}


def aggregate_outputs(common):
    OUT.mkdir(parents=True,exist_ok=True); (OUT/"comparison").mkdir(exist_ok=True)
    save_json(OUT/"E3_0_manifest.json",common)
    combined=[]
    for name in BRANCHES:
        path=OUT/name/"pareto_designs.tsv"
        if path.exists():
            with path.open(encoding="utf-8-sig",newline="") as f:
                combined.extend(csv.DictReader(f,delimiter="\t"))
    write_tsv(OUT/"comparison/S01_vs_S03.tsv",combined)
    write_tsv(OUT/"comparison/nondominated_combined.tsv",pareto(combined))
    def summary(branch):
        path=OUT/branch/"pareto_designs.tsv"
        if not path.exists():
            return f"### {branch}\n\nNo completed Pareto output is available.\n"
        with path.open(encoding="utf-8-sig",newline="") as stream:
            designs=list(csv.DictReader(stream,delimiter="\t"))
        used=[]
        for row in designs:
            used.extend(x for x in row.get("KO_reaction_ids","").split(";") if x)
            used.extend(x for x in row.get("regulatory_constraints","").split(";") if x)
        repeat=sorted((x,used.count(x)) for x in set(used) if used.count(x)>1)
        recurring="; ".join(f"{x} ({n})" for x,n in repeat) or "none observed"
        return f"### {branch}\n\nPareto designs retained: {len(designs)}. Repeated interventions: {recurring}.\n"
    report=["# E3.0 whole-host optimization", "",
        f"Host SHA256: `{common['host_sha256']}`  ",
        f"E1 pool SHA256: `{common['E1_pool_sha256']}`  ",
        f"E2R3 input fingerprint: `{common['E2R3_input_fingerprint']}`  ",
        f"E3.0 script SHA256: `{common['E3_script_sha256']}`  ",
        f"Solver: {common['solver']}; seed {common['seed']}; time limit {common['time_limit_seconds']} s per solve.", "",
        "Unit-cost reaction KOs and branch-specific regulatory flux constraints were searched together using RobustKnock. "
        "Budgets: 1, 2, 4, 6; budget 8 is considered only after a measurable gain at budget 6. "
        "The objective modes minimize Rubisco or reduce Fe2, O2, and inorganic-carbon uptake. Regulatory thresholds are model-level flux hypotheses.", "",
        summary("S01"), "", summary("S03"), "",
        f"Combined nondominated designs retained: {len(pareto(combined))}. See the branch and comparison TSV files for phenotypes and interventions.", ""]
    (OUT/"E3_0_report.md").write_text("\n".join(report),encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch",choices=["S01","S03","all"],default="all")
    parser.add_argument("--prepare-only",action="store_true")
    parser.add_argument("--time-limit",type=int,default=MAX_TIME)
    parser.add_argument("--aggregate-only",action="store_true",help="write common manifest/comparison after independent branch jobs finish")
    args=parser.parse_args(); os.chdir(ROOT)
    if sha(HOST)!=EXPECTED_HOST: raise SystemExit(f"Authoritative host SHA mismatch: {sha(HOST)}")
    e2=load_module(E2_PATH,"e2_helpers_e3_main"); e2r=load_module(E2R_PATH,"e2r_helpers_e3_main")
    e1=e2.load_e1_helpers(); host,donors,aliases_i,aliases_c=e2.load_models(e1)
    cfgpath=ROOT/"config/scenarios.json"; config=json.loads(cfgpath.read_text(encoding="utf-8"))
    rows=e2.read_tsv(POOL)
    if len(rows)!=2424: raise SystemExit(f"E1 row count expected 2424, got {len(rows)}")
    if sha(E2R3_PATH)!=json.loads(E2R3_MANIFEST.read_text(encoding="utf-8"))["inputs"]["script"]:
        # E2R3 manifest records the exact input-stage script hash; refuse stale support context.
        raise SystemExit("E2R3 script fingerprint differs from the completed E2R3 manifest")
    ref=json.loads(E2R3_MANIFEST.read_text(encoding="utf-8"))
    supports=load_supports()
    common={"created_utc":datetime.now(timezone.utc).isoformat(),"host_path":str(HOST),"host_sha256":sha(HOST),
        "E1_pool_sha256":sha(POOL),"E2R3_input_fingerprint":ref["fingerprint"],"fixed_sets":supports,
        "donor_model_sha256":{name:sha(path) for name,path in e2.DONORS.items()},
        "helper_script_sha256":{"E2":sha(E2_PATH),"E2R":sha(E2R_PATH),"E2R3":sha(E2R3_PATH)},
        "package_versions":{},"solver":"StrainDesign RobustKnock / SCIP","seed":1,"time_limit_seconds":args.time_limit,
        "straindesign_patch_sha256":{p.name:sha(p) for p in SD_PATCH_FILES},
        "max_solutions":MAX_SOLUTIONS,"budgets":list(BUDGETS),"modes":MODES,
        "environment":ref["Stage2"],"E3_script_sha256":sha(Path(__file__))}
    for pkg in ("cobra","straindesign","optlang","pyscipopt","scipy"):
        try: common["package_versions"][pkg]=package_version(pkg)
        except Exception: common["package_versions"][pkg]="not_installed"
    OUT.mkdir(parents=True,exist_ok=True)
    if args.aggregate_only:
        aggregate_outputs(common)
        return
    if args.branch=="all":
        aggregate_outputs(common)
    branches=BRANCHES if args.branch=="all" else (args.branch,)
    results=[]
    for branch in branches:
        try:
            results.append(branch_job(branch,e2,e2r,config,rows,donors,aliases_i,aliases_c,ref,args.time_limit,args.prepare_only))
        except Exception as exc:
            result={"branch":branch,"status":"error","error":f"{type(exc).__name__}: {exc}"}
            results.append(result); save_json(OUT/branch/"checkpoints/branch_error.json",result)
    if args.branch=="all":
        save_json(OUT/"checkpoints_run_summary.json",{"results":results})
        if not args.prepare_only: aggregate_outputs(common)
    print(json.dumps(results,indent=2))


if __name__=="__main__": main()
