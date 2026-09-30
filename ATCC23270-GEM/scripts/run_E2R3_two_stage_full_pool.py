#!/usr/bin/env python
"""E2R3: characterize the Stage-1 full-pool rescue and validate compact sets.

All reactions are added to a derived in-memory host copy via the established
E2 candidate constructor. The authoritative host and prior E1/E2/E2R outputs
are read-only inputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version as package_version
from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.parsimonious import add_pfba

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/engineering/E2R3_two_stage_full_pool"
E2_PATH = ROOT / "scripts/run_E2_strain_design_KI.py"
E2R_PATH = ROOT / "scripts/run_E2R_resource_matched_KI.py"
HOST = ROOT.parent / "before final" / "updatedv3.1_verified" / "updatedv3.1.xml"
POOL = ROOT / "results/engineering/E1_donor_reaction_pool/candidate_reactions.tsv"
EXPECTED_HOST = "026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3"
EXPECTED_ROWS = 2424
IDS = {"biomass": "Ex_bio[e]", "glucose": "Ex_glc-B[e]", "rubisco": "RUBISCO",
       "co2": "Ex_h2co3[e]", "fe2": "Ex_fe2[e]", "o2": "Ex_o2[e]"}
TOL = 1e-8
BIOMASS_MULTIPLIER = 1.5
GLUCOSE_CAP = 0.5
GLUCOSE_MIN = 0.05
CO2_BOUNDS = (-2.0, 0.0)
STAGE1_FRACTION = 0.30
STAGE2_FRACTION = 0.10
MAX_SUPPORTS = 5


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def fingerprint_tree(path: Path) -> str:
    h = hashlib.sha256()
    for item in sorted(path.rglob("*.py")):
        h.update(item.relative_to(path).as_posix().encode())
        h.update(item.read_bytes())
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


def add_constraints(model, ids, bmin: float | None, rubisco_cap: float | None, glucose_on: bool = True):
    """Apply the common E2R3 environment to this derived model copy."""
    ids["fe2"].lower_bound = -model.notes["e2r3_fe_wt"]
    ids["fe2"].upper_bound = 0
    ids["o2"].lower_bound = -model.notes["e2r3_o2_wt"]
    ids["o2"].upper_bound = 0
    ids["glucose"].bounds = (-GLUCOSE_CAP, 0) if glucose_on else (0, 0)
    ids["co2"].bounds = CO2_BOUNDS
    ids["rubisco"].lower_bound = 0
    if rubisco_cap is not None:
        ids["rubisco"].upper_bound = min(ids["rubisco"].upper_bound, rubisco_cap)
    if bmin is not None:
        bconstraint = model.problem.Constraint(ids["biomass"].flux_expression, lb=bmin,
                                                name=f"E2R3_growth_{len(model.solver.constraints)}")
        model.add_cons_vars(bconstraint)
    if glucose_on:
        model.add_cons_vars(model.problem.Constraint(ids["glucose"].flux_expression, ub=-GLUCOSE_MIN,
                    name=f"E2R3_glucose_{len(model.solver.constraints)}"))
    model.add_cons_vars(model.problem.Constraint(ids["fe2"].flux_expression, ub=-0.01,
                    name=f"E2R3_trace_fe_{len(model.solver.constraints)}"))


def flux_summary(solution, ids):
    if solution.status != "optimal":
        return {"status": solution.status}
    f = solution.fluxes
    return {"status": "optimal", "biomass": float(f[ids["biomass"].id]),
            "glucose": float(f[ids["glucose"].id]), "rubisco": float(f[ids["rubisco"].id]),
            "co2": float(f[ids["co2"].id]), "fe2": float(f[ids["fe2"].id]),
            "o2": float(f[ids["o2"].id])}


def phenotype(model, ids, candidate_ids, selected, bmin, rcap, glucose_on=True, parsimonious=False):
    """Run independently reconstructed FBA/pFBA with only selected KIs open."""
    selected = set(selected)
    with model:
        for rid in candidate_ids:
            if rid not in selected:
                model.reactions.get_by_id(rid).bounds = (0, 0)
        # For glucose-off, remove the with-glucose biomass target and uptake
        # floor; retain all other resource and Rubisco constraints.
        add_constraints(model, ids, bmin if glucose_on else None, rcap, glucose_on)
        model.objective = ids["biomass"]
        if not parsimonious:
            sol = model.optimize()
        else:
            # The explicit phenotype constraints retain the 1.5x WT floor;
            # minimize total flux without imposing extra growth maximization.
            model.objective = model.problem.Objective(0, direction="max")
            add_pfba(model, fraction_of_optimum=0)
            sol = model.optimize()
        return flux_summary(sol, ids)


def reaction_metadata(candidates, provenance, donors, variants):
    out = {}
    for rid in candidates:
        rxn = None
        sources = provenance[rid]["source_reactions"]
        source_rxn = None
        for source in sources:
            donor, source_id = source.split(":", 1)
            try:
                source_rxn = donors[donor].reactions.get_by_id(source_id)
                break
            except KeyError:
                pass
        out[rid] = {"candidate_id": provenance[rid]["candidate_id"],
                    "source_reactions": sources,
                    "direction": provenance[rid].get("variant", "combined"),
                    "lower_bound":float(variants[rid]["lower_bound"]),
                    "upper_bound":float(variants[rid]["upper_bound"]),
                    "equation": source_rxn.reaction if source_rxn else ""}
    return out


def candidate_fva(model, ids, candidates, bmin, rcap, fingerprint, chunk_size=120, resume=True):
    """FVA in restartable chunks over all states satisfying Stage-1 phenotype."""
    checkpoint=OUT/"checkpoints/stage1_fva_chunks.json"
    completed={}
    if resume and checkpoint.exists():
        try:
            old=json.loads(checkpoint.read_text(encoding="utf-8"))
            if old.get("fingerprint")==fingerprint:
                completed=old.get("chunks",{})
        except Exception:
            completed={}
    candidate_ids=list(candidates)
    chunks=[candidate_ids[i:i+chunk_size] for i in range(0,len(candidate_ids),chunk_size)]
    for index, chunk in enumerate(chunks):
        key=str(index)
        if key in completed:
            continue
        with model:
            add_constraints(model, ids, bmin, rcap, True)
            model.objective=ids["biomass"]
            frame=flux_variability_analysis(model,reaction_list=chunk,
                fraction_of_optimum=0.0,processes=1,loopless=None)
        completed[key]={rid:{"minimum":float(frame.at[rid,"minimum"]),
                             "maximum":float(frame.at[rid,"maximum"])} for rid in chunk}
        save_json(checkpoint,{"fingerprint":fingerprint,"chunk_size":chunk_size,
                              "chunk_count":len(chunks),"completed_chunks":len(completed),
                              "chunks":completed})
        with (OUT/"logs/run_log.jsonl").open("a",encoding="utf-8") as log:
            log.write(json.dumps({"event":"fva_chunk_complete","chunk":index,
                "chunk_count":len(chunks),"reaction_count":len(chunk),
                "fingerprint":fingerprint,"utc":datetime.now(timezone.utc).isoformat()})+"\n")
    merged={}
    for chunk_result in completed.values(): merged.update(chunk_result)
    missing=set(candidates)-set(merged)
    if missing: raise RuntimeError(f"FVA checkpoint is missing {len(missing)} reaction IDs")
    import pandas as pd
    return pd.DataFrame.from_dict(merged,orient="index")


def classify_candidates(model, ids, candidates, metadata, pfba_flux, fva):
    rows = []
    for rid in candidates:
        lo, hi = float(fva.at[rid, "minimum"]), float(fva.at[rid, "maximum"])
        pf = float(pfba_flux.get(rid, 0.0))
        mandatory = lo > TOL or hi < -TOL
        active = abs(pf) > TOL
        if abs(lo) <= TOL and abs(hi) <= TOL:
            category = "inactive_stage1"
        elif mandatory:
            category = "mandatory_under_stage1"
        elif active:
            category = "active_in_pfba"
        else:
            category = "optional_stage1"
        rxn = model.reactions.get_by_id(rid)
        rows.append({"E1_candidate_id": metadata[rid]["candidate_id"], "search_reaction_id": rid,
          "donor_source_reaction": ";".join(metadata[rid]["source_reactions"]),
          "donor_direction_variant": metadata[rid]["direction"],
          "donor_supported_lower_bound":metadata[rid]["lower_bound"],
          "donor_supported_upper_bound":metadata[rid]["upper_bound"], "pfba_flux": pf,
          "FVA_minimum": lo, "FVA_maximum": hi,
          "active_in_pfba":active,"mandatory_under_stage1":mandatory,
          "classification": category, "direction": "none" if abs(lo)<=TOL and abs(hi)<=TOL else "forward" if lo >= -TOL else "reverse" if hi <= TOL else "both",
          "reaction_equation": rxn.reaction})
    return rows


def prune_seed(model, ids, candidates, seed, bmin, rcap, order, strategy):
    """Deterministic single-deletion pruning to an irreducible set."""
    kept = set(seed)
    for rid in order(kept):
        if rid not in kept:
            continue
        trial = kept - {rid}
        result = phenotype(model, ids, candidates, trial, bmin, rcap)
        if result["status"] == "optimal":
            kept = trial
    # Single-deletion feasibility is monotone under further reaction removal:
    # a reaction that cannot be removed now cannot become removable later.
    # The sequential sweep therefore ends irreducible without a duplicate pass.
    ph = phenotype(model, ids, candidates, kept, bmin, rcap)
    if ph["status"] != "optimal":
        return None
    return {"reaction_ids": sorted(kept), "strategy": strategy, "phenotype": ph}


def forced_pfba_seed(model, ids, candidates, mandatory, optional_id, fva, bmin, rcap):
    """Return a deterministic alternative active set requiring one optional flux.

    Forcing an FVA-capable reaction before pFBA can expose routes requiring
    multiple reactions that were inactive in the original parsimonious state.
    """
    lo=float(fva.at[optional_id,"minimum"]); hi=float(fva.at[optional_id,"maximum"])
    if hi > 2*TOL:
        direction=1; extreme=hi
    elif lo < -2*TOL:
        direction=-1; extreme=abs(lo)
    else:
        return None
    epsilon=min(1e-3,extreme/2.0)
    if epsilon <= TOL:
        return None
    with model:
        add_constraints(model,ids,bmin,rcap,True)
        rxn=model.reactions.get_by_id(optional_id)
        if direction>0:
            model.add_cons_vars(model.problem.Constraint(rxn.flux_expression,lb=epsilon,
                 name="E2R3_optional_force_"+hashlib.sha1(optional_id.encode()).hexdigest()[:12]))
        else:
            model.add_cons_vars(model.problem.Constraint(rxn.flux_expression,ub=-epsilon,
                 name="E2R3_optional_force_"+hashlib.sha1(optional_id.encode()).hexdigest()[:12]))
        model.objective=model.problem.Objective(0,direction="max")
        add_pfba(model,fraction_of_optimum=0)
        sol=model.optimize()
        if sol.status!="optimal": return None
        flux=sol.fluxes
        seed={rid for rid in candidates if abs(float(flux.get(rid,0)))>TOL}
    seed.update(mandatory)
    return {"reaction_ids":sorted(seed),"forced_reaction":optional_id,
            "forced_direction":"forward" if direction>0 else "reverse","forced_flux_minimum":epsilon}


def build_support_sets(model, ids, candidates, metadata, class_rows, pfba_flux, bmin, r30, fva, fingerprint, resume=True):
    mandatory = {x["search_reaction_id"] for x in class_rows if x["classification"] == "mandatory_under_stage1"}
    active = {rid for rid in candidates if abs(float(pfba_flux.get(rid, 0))) > TOL}
    base = mandatory | active
    strategy_list = [
      ("low_pfba_flux_first", lambda s: sorted(s, key=lambda x:(abs(float(pfba_flux.get(x,0))),x))),
      ("high_pfba_flux_first", lambda s: sorted(s, key=lambda x:(-abs(float(pfba_flux.get(x,0))),x))),
      ("candidate_id_ascending", lambda s: sorted(s)),
      ("candidate_id_descending", lambda s: sorted(s, reverse=True)),
    ]
    # Optional reactions with non-zero FVA span can enter alternative seeds.
    # Try them in descending span, then stable reaction-ID order, stopping when
    # the target number of distinct validated sets has been reached.
    optional = [rid for rid in candidates if rid not in base and
                (abs(float(fva.at[rid,"minimum"])) > 2*TOL or abs(float(fva.at[rid,"maximum"])) > 2*TOL)]
    optional.sort(key=lambda rid:(-max(abs(float(fva.at[rid,"minimum"])),abs(float(fva.at[rid,"maximum"]))),rid))
    optional=optional[:32]
    proposals = [(base, label, order) for label, order in strategy_list]
    proposals.extend((None, f"forced_pfba_optional_{rid}",
                      lambda s, optional_id=rid: sorted(s, key=lambda x:(x==optional_id,abs(float(pfba_flux.get(x,0))),x)),rid)
                     for rid in optional)
    checkpoint=OUT/"checkpoints/support_sets.json"
    supports=[]; strategy_log=[]
    if resume and checkpoint.exists():
        try:
            old=json.loads(checkpoint.read_text(encoding="utf-8"))
            if old.get("fingerprint")==fingerprint:
                supports=old.get("supports",[]); strategy_log=old.get("strategies",[])
        except Exception:
            pass
    seen={tuple(x["reaction_ids"]) for x in supports}
    finished={x["strategy"] for x in strategy_log}
    normalized=[]
    for item in proposals:
        normalized.append((*item,None) if len(item)==3 else item)
    for seed, label, order, forced_id in normalized:
        if len(supports) >= MAX_SUPPORTS:
            break
        if label in finished:
            continue
        if forced_id is not None:
            alternative=forced_pfba_seed(model,ids,candidates,mandatory,forced_id,fva,bmin,r30)
            if alternative is None:
                strategy_log.append({"strategy":label,"status":"forced_pfba_infeasible_or_no_signal"})
                save_json(checkpoint,{"fingerprint":fingerprint,"terminal":False,"supports":supports,"strategies":strategy_log})
                continue
            seed=alternative["reaction_ids"]
        pruned=prune_seed(model,ids,candidates,seed,bmin,r30,order,label)
        if pruned is None:
            strategy_log.append({"strategy":label,"status":"seed_infeasible"})
            save_json(checkpoint,{"fingerprint":fingerprint,"terminal":False,"supports":supports,"strategies":strategy_log})
            continue
        key=tuple(pruned["reaction_ids"])
        if key in seen:
            strategy_log.append({"strategy":label,"status":"duplicate","support_set_id":next(x["support_set_id"] for x in supports if tuple(x["reaction_ids"])==key)})
            save_json(checkpoint,{"fingerprint":fingerprint,"terminal":False,"supports":supports,"strategies":strategy_log})
            continue
        validation=phenotype(model,ids,candidates,pruned["reaction_ids"],bmin,r30)
        if validation["status"]!="optimal":
            strategy_log.append({"strategy":label,"status":"independent_validation_failed"})
            save_json(checkpoint,{"fingerprint":fingerprint,"terminal":False,"supports":supports,"strategies":strategy_log})
            continue
        # prune_seed has already tested each single deletion. This independent
        # reconstruction below is the required final Stage-1 validation.
        off=phenotype(model,ids,candidates,pruned["reaction_ids"],bmin,r30,glucose_on=False)
        idx=len(supports)+1; sid=f"S{idx:02d}"
        mandatory_selected=sorted(set(pruned["reaction_ids"]) & mandatory)
        support={**pruned,"support_set_id":sid,"validation":validation,"glucose_off":off,
                 "delta_biomass":validation["biomass"]-off["biomass"] if off["status"]=="optimal" else None,
                 "mandatory_reaction_ids":mandatory_selected}
        supports.append(support); seen.add(key)
        strategy_log.append({"strategy":label,"status":"validated","support_set_id":sid,"reaction_count":len(key)})
        save_json(checkpoint,{"fingerprint":fingerprint,"terminal":False,"supports":supports,"strategies":strategy_log})
    save_json(checkpoint,{"fingerprint":fingerprint,"terminal":True,"supports":supports,"strategies":strategy_log,
                          "mandatory_reaction_ids":sorted(mandatory)})
    return supports, strategy_log, mandatory


def map_support(support, metadata):
    ids=support["reaction_ids"]
    return {"E1_candidate_ids": sorted({metadata[r]["candidate_id"] for r in ids}),
            "donor_source_reactions": sorted({s for r in ids for s in metadata[r]["source_reactions"]})}


def support_rows(supports, metadata):
    rows=[]
    for x in supports:
        ph=x["validation"]; off=x["glucose_off"]; trace=map_support(x,metadata)
        rows.append({"support_set_id":x["support_set_id"],"reaction_count":len(x["reaction_ids"]),
          "search_variant_count":len(x["reaction_ids"]),"E1_candidate_count":len(trace["E1_candidate_ids"]),
          "E1_candidate_ids":";".join(trace["E1_candidate_ids"]),"search_reaction_ids":";".join(x["reaction_ids"]),
          "donor_source_reactions":";".join(trace["donor_source_reactions"]),"generation_pruning_strategy":x["strategy"],
          "Stage1_biomass":ph["biomass"],"Stage1_glucose_uptake":ph["glucose"],"Stage1_Rubisco_flux":ph["rubisco"],
          "Stage1_CO2_flux":ph["co2"],"Stage1_Fe2_flux":ph["fe2"],"Stage1_O2_flux":ph["o2"],
          "Stage1_glucose_off_biomass":off.get("biomass"),"Stage1_delta_biomass":x["delta_biomass"],
          "Stage1_glucose_contribution_pass":off.get("status")=="optimal" and x["delta_biomass"] is not None and x["delta_biomass"]>TOL,
          "Stage1_validation_status":ph["status"],"mandatory_search_reaction_ids":";".join(x["mandatory_reaction_ids"])})
    return rows


def run_stage2(model, ids, candidates, supports, bmin, r10):
    rows=[]; results={}
    for x in supports:
        ph=phenotype(model,ids,candidates,x["reaction_ids"],bmin,r10)
        if ph["status"] not in ("optimal","infeasible"):
            raise RuntimeError(f"Stage 2 unresolved for {x['support_set_id']}: {ph['status']}")
        off=phenotype(model,ids,candidates,x["reaction_ids"],bmin,r10,glucose_on=False) if ph["status"]=="optimal" else {"status":"not_run"}
        if ph["status"]=="optimal" and off["status"] not in ("optimal","infeasible"):
            raise RuntimeError(f"Stage 2 glucose-off unresolved for {x['support_set_id']}: {off['status']}")
        passed=ph["status"]=="optimal"
        cls="PASS_BOTH" if passed else "PASS_30_ONLY"
        results[x["support_set_id"]]={"phenotype":ph,"glucose_off":off,"classification":cls}
        rows.append({"support_set_id":x["support_set_id"],"reaction_count":len(x["reaction_ids"]),
          "Stage1_status":"optimal","Stage2_status":ph["status"],"Stage2_biomass":ph.get("biomass"),
          "Stage2_glucose_uptake":ph.get("glucose"),"Stage2_Rubisco_flux":ph.get("rubisco"),
          "Stage2_CO2_flux":ph.get("co2"),"Stage2_Fe2_flux":ph.get("fe2"),"Stage2_O2_flux":ph.get("o2"),
          "Stage2_glucose_off_biomass":off.get("biomass"),
          "Stage2_delta_biomass":ph.get("biomass",0)-off.get("biomass",0) if passed and off.get("status")=="optimal" else None,
          "final_classification":cls})
        save_json(OUT/"checkpoints/stage2_validation_partial.json",{"completed_support_sets":len(rows),"results":results})
    return rows,results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only",action="store_true",help="reproduce baselines, full-pool FBA, pFBA, and FVA")
    parser.add_argument("--run",action="store_true",help="run full Stage-1 characterization and two-stage support analysis")
    parser.add_argument("--resume",action="store_true",help="reuse terminal checkpoints matching the current input fingerprint")
    args=parser.parse_args()
    os.chdir(ROOT)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"checkpoints").mkdir(exist_ok=True); (OUT/"logs").mkdir(exist_ok=True)
    e2=load_module(E2_PATH,"e2_helpers")
    e2r=load_module(E2R_PATH,"e2r_helpers")
    host_hash=sha(HOST)
    if host_hash!=EXPECTED_HOST: raise SystemExit(f"Authoritative host SHA256 mismatch: {host_hash}")
    e1_rows=e2.read_tsv(POOL)
    if len(e1_rows)!=EXPECTED_ROWS: raise SystemExit(f"Expected {EXPECTED_ROWS} E1 rows; found {len(e1_rows)}")
    e1_hash=sha(POOL)
    helper=e2.load_e1_helpers(); host,donors,alias_i,alias_c=e2.load_models(helper)
    cfgpath=ROOT/"config/scenarios.json"; config=json.loads(cfgpath.read_text(encoding="utf-8"))
    # WT glucose-closed FIM, with established inorganic-carbon supply fixed.
    wt=host.copy(); e2r.native_condition(wt,config,"FIM",e2); wid=e2r.ids_for(wt,e2)
    wid["glucose"].bounds=(0,0); wid["co2"].bounds=(-2,-2); wt.objective=wid["biomass"]
    best=wt.optimize()
    if best.status!="optimal": raise RuntimeError(f"Glucose-closed WT FBA: {best.status}")
    bwt=float(best.fluxes[wid["biomass"].id])
    with wt:
        wt.objective=wid["biomass"]
        # pFBA at the exact WT maximum is the established resource reference.
        wt_psol=cobra.flux_analysis.pfba(wt,fraction_of_optimum=1.0)
    wt_flux=e2r.fluxes(wt_psol.fluxes,wid); fe_wt=abs(wt_flux["fe2"]); o2_wt=abs(wt_flux["o2"])
    expected={"B_WT":0.0520763871,"Fe_WT":163.781529,"O2_WT":38.481543,"R_REF":0.791344958}
    for key,val in (("B_WT",bwt),("Fe_WT",fe_wt),("O2_WT",o2_wt)):
        target=expected[key]; tol=max(1e-6,abs(target)*1e-6)
        if abs(val-target)>tol: raise SystemExit(f"E2R baseline discrepancy: {key}={val:.12g}, expected≈{target}")
    search=host.copy(); e2r.native_condition(search,config,"FIM",e2)
    ids=e2r.ids_for(search,e2); search.notes["e2r3_fe_wt"]=fe_wt; search.notes["e2r3_o2_wt"]=o2_wt
    ids["fe2"].bounds=(-fe_wt,0); ids["o2"].bounds=(-o2_wt,0); ids["co2"].bounds=CO2_BOUNDS
    ids["glucose"].bounds=(-GLUCOSE_CAP,0)
    mixed=e2.fba_summary(search,ids)
    if mixed.get("status")!="optimal": raise RuntimeError(f"Resource-matched mixed reference: {mixed}")
    rref=float(mixed["rubisco_flux"])
    if abs(rref-expected["R_REF"])>max(1e-6,expected["R_REF"]*1e-6):
        raise SystemExit(f"E2R baseline discrepancy: R_REF={rref:.12g}, expected≈{expected['R_REF']}")
    bmin=BIOMASS_MULTIPLIER*bwt; r30=STAGE1_FRACTION*rref; r10=STAGE2_FRACTION*rref
    candidates,provenance,variants=e2.add_candidate_reactions(search,donors,alias_i,alias_c,e1_rows)
    metadata=reaction_metadata(candidates,provenance,donors,variants)
    donor_hashes={name:sha(path) for name,path in e2.DONORS.items()}
    sd_path=e2.STRAINDESIGN/"straindesign"
    package_versions={}
    for pkg in ("cobra","straindesign","optlang","swiglpk","pyscipopt","scipy"):
        try: package_versions[pkg]=package_version(pkg)
        except Exception: package_versions[pkg]="not_installed"
    solver_interface=getattr(search.solver.interface,"__name__",str(search.solver.interface))
    solver_info={"interface":solver_interface,
                 "solver_name":str(getattr(search.solver,"name","unknown")),
                 "backend_version":package_versions.get("swiglpk") if "glpk" in solver_interface.lower() else package_versions.get("scipy"),
                 "configured_tolerances":str(getattr(search.solver.configuration,"tolerances",{}))}
    fpdata={"host":host_hash,"E1_pool":e1_hash,"E2_script":sha(E2_PATH),"E2R_script":sha(E2R_PATH),
      "script":sha(Path(__file__)),"config":sha(cfgpath),"donors":donor_hashes,
      "thresholds":{"bmin":bmin,"r30":r30,"r10":r10,"Fe_WT":fe_wt,"O2_WT":o2_wt,
                    "glucose_cap":GLUCOSE_CAP,"glucose_min":GLUCOSE_MIN,"CO2":CO2_BOUNDS}}
    fingerprint=hashlib.sha256(json.dumps(fpdata,sort_keys=True).encode()).hexdigest()
    base={"host_sha256":host_hash,"E1_pool_sha256":e1_hash,"E1_candidate_rows":len(e1_rows),
      "donor_model_sha256":json.dumps(donor_hashes,sort_keys=True),"B_WT":bwt,
      "WT_pFBA_Fe2_flux":wt_flux["fe2"],"WT_pFBA_Fe2_uptake_magnitude":fe_wt,
      "WT_pFBA_O2_flux":wt_flux["o2"],"WT_pFBA_O2_uptake_magnitude":o2_wt,
      "WT_pFBA_Rubisco_flux":wt_flux["rubisco"],"WT_pFBA_CO2_flux":wt_flux["co2"],
      "mixed_no_KI_biomass":mixed["biomass"],"mixed_no_KI_glucose_flux":mixed["glucose_uptake"],
      "R_REF":rref,"mixed_no_KI_Rubisco_flux":mixed["rubisco_flux"],"mixed_no_KI_CO2_flux":mixed["co2_flux"],
      "mixed_no_KI_Fe2_flux":mixed["fe2_flux"],"mixed_no_KI_O2_flux":mixed["o2_flux"],
      "Stage1_biomass_threshold":bmin,"Stage1_Rubisco_cap":r30,"Stage2_Rubisco_cap":r10,
      "glucose_bounds":f"[-{GLUCOSE_CAP},0]","glucose_minimum_use":GLUCOSE_MIN,
      "Fe2_bounds":f"[-{fe_wt},0]","O2_bounds":f"[-{o2_wt},0]","CO2_bounds":str(CO2_BOUNDS),
      "input_fingerprint":fingerprint}
    write_tsv(OUT/"baseline_summary.tsv",[base])
    manifest={"created_utc":datetime.now(timezone.utc).isoformat(),"fingerprint":fingerprint,
      "inputs":fpdata,"host_path":str(HOST),"candidate_variants":len(candidates),"B_WT":bwt,
      "Fe_WT":fe_wt,"O2_WT":o2_wt,"mixed_reference":mixed,"R_REF":rref,
      "Stage1":{"biomass_min":bmin,"rubisco_max":r30,"glucose_min_uptake":GLUCOSE_MIN},
      "Stage2":{"biomass_min":bmin,"rubisco_max":r10,"glucose_min_uptake":GLUCOSE_MIN},
      "software_versions":package_versions,"solver":solver_info,
      "tolerance":TOL,"pfba":"WT: maximize biomass then COBRApy pFBA at exact optimum. Stage1: retain explicit phenotype constraints including biomass >=1.5xB_WT; objective=0 and COBRApy add_pfba(fraction_of_optimum=0) minimize sum of forward/reverse flux without forcing maximum biomass.",
      "fva":"COBRApy flux_variability_analysis; fraction_of_optimum=0; exact explicit Stage1 phenotype constraints; loopless=None"}
    save_json(OUT/"checkpoints/input_manifest.json",manifest)
    (OUT/"logs/run_log.jsonl").touch(exist_ok=True)
    if not args.run and not args.prepare_only:
        print(json.dumps(manifest,indent=2)); return
    # Positive control validates that the reused StrainDesign KI semantics still work.
    control=e2.smoke_sd_ki_semantics()
    if not control.get("passed"): raise RuntimeError(f"KI semantics positive control failed: {control}")
    save_json(OUT/"checkpoints/KI_semantics_positive_control.json",control)
    b30=bmin
    checkpoint=OUT/"checkpoints/stage1_full_pool.json"
    old_full=None
    if args.resume and checkpoint.exists():
        try:
            old_full=json.loads(checkpoint.read_text(encoding="utf-8"))
            if old_full.get("fingerprint")!=fingerprint or not old_full.get("terminal"):
                old_full=None
        except Exception:
            old_full=None
    if old_full and old_full.get("status")=="STAGE1_FULL_POOL_INFEASIBLE":
        write_report(manifest,"STAGE1_FULL_POOL_INFEASIBLE",old_full.get("phenotype",{}),[],[],[],[],[],None)
        return
    t0=time.time()
    if old_full and old_full.get("status")=="feasible":
        full=old_full["phenotype_fba"]
    else:
        full=phenotype(search,ids,candidates,candidates,b30,r30,True,False)
        if full["status"]!="optimal":
            terminal=full["status"]=="infeasible"
            status="STAGE1_FULL_POOL_INFEASIBLE" if terminal else "error"
            save_json(checkpoint,{"fingerprint":fingerprint,"terminal":terminal,"status":status,
                 "phenotype":full,"runtime_seconds":time.time()-t0})
            if terminal:
                write_report(manifest,status,full,[],[],[],[],[],None)
            else:
                raise RuntimeError(f"Stage1 full-pool FBA unresolved: {full['status']}")
            return
    # Threshold-conditioned pFBA with full Stage-1 constraints enforced.
    parsimonious=phenotype(search,ids,candidates,candidates,b30,r30,True,True)
    if parsimonious["status"]!="optimal": raise RuntimeError(f"Stage1 pFBA failed: {parsimonious}")
    with search:
        add_constraints(search,ids,b30,r30,True)
        search.objective=ids["biomass"]
        sol=search.optimize()
        if sol.status!="optimal": raise RuntimeError(f"Stage1 full-pool solution changed: {sol.status}")
        flux=sol.fluxes.copy()
        # Replace with threshold-conditioned parsimonious fluxes by a second solve.
        search.objective=search.problem.Objective(0,direction="max")
        add_pfba(search,fraction_of_optimum=0)
        psol=search.optimize()
        if psol.status!="optimal": raise RuntimeError(f"Stage1 pFBA refit failed: {psol.status}")
        pflux=psol.fluxes.copy()
    fva=candidate_fva(search,ids,candidates,b30,r30,fingerprint,resume=args.resume)
    class_rows=classify_candidates(search,ids,candidates,metadata,pflux,fva)
    full_flux_rows=[]
    for rxn in search.reactions:
        v=float(flux.get(rxn.id,0)); pv=float(pflux.get(rxn.id,0))
        if rxn.id in metadata or abs(v)>TOL or abs(pv)>TOL:
            meta=metadata.get(rxn.id,{})
            full_flux_rows.append({"reaction_id":rxn.id,"reaction_type":"E1_candidate" if meta else "host",
              "E1_candidate_id":meta.get("candidate_id",""),"donor_source_reaction":";".join(meta.get("source_reactions",[])),
              "reaction_equation":rxn.reaction,"FBA_flux":v,"pFBA_flux":pv,
              "FVA_minimum":float(fva.at[rxn.id,"minimum"]) if rxn.id in fva.index else None,
              "FVA_maximum":float(fva.at[rxn.id,"maximum"]) if rxn.id in fva.index else None})
    write_tsv(OUT/"stage1_full_pool_flux.tsv",full_flux_rows)
    write_tsv(OUT/"stage1_candidate_classification.tsv",class_rows)
    save_json(checkpoint,{"fingerprint":fingerprint,"terminal":True,"status":"feasible",
      "phenotype_fba":full,"phenotype_pfba":parsimonious,"stage1_flux":{k:float(flux[k]) for k in flux.index if abs(float(flux[k]))>TOL},
      "stage1_pfba_flux":{k:float(pflux[k]) for k in pflux.index if abs(float(pflux[k]))>TOL},
      "runtime_seconds":time.time()-t0})
    mandatory=[r["search_reaction_id"] for r in class_rows if r["classification"]=="mandatory_under_stage1"]
    active=[rid for rid in candidates if abs(float(pflux.get(rid,0)))>TOL]
    fva_state={rid:{"minimum":float(fva.at[rid,"minimum"]),"maximum":float(fva.at[rid,"maximum"])} for rid in candidates}
    supports,strategy_log,mandatory=build_support_sets(search,ids,candidates,metadata,class_rows,pflux,b30,r30,fva,fingerprint,args.resume)
    write_tsv(OUT/"stage1_support_sets.tsv",support_rows(supports,metadata))
    stage2_rows,stage2=run_stage2(search,ids,candidates,supports,b30,r10)
    write_tsv(OUT/"stage2_validation.tsv",stage2_rows)
    final=[]
    for x in supports:
        second=stage2[x["support_set_id"]]; stage1=x["validation"]; off1=x["glucose_off"]
        trace=map_support(x,metadata); cls=second["classification"]
        final.append({"priority_group":"A_PASS_BOTH" if cls=="PASS_BOTH" else "B_PASS_30_ONLY",
          "support_set_id":x["support_set_id"],"reaction_count":len(x["reaction_ids"]),
          "search_variant_count":len(x["reaction_ids"]),"E1_candidate_count":len(trace["E1_candidate_ids"]),
          "E1_candidate_ids":";".join(trace["E1_candidate_ids"]),"search_reaction_ids":";".join(x["reaction_ids"]),
          "donor_source_reactions":";".join(trace["donor_source_reactions"]),
          "Stage1_biomass":stage1["biomass"],"Stage1_glucose":stage1["glucose"],"Stage1_Rubisco":stage1["rubisco"],
          "Stage1_CO2":stage1["co2"],"Stage1_Fe2":stage1["fe2"],"Stage1_O2":stage1["o2"],
          "Stage2_status":second["phenotype"]["status"],"Stage2_biomass":second["phenotype"].get("biomass"),
          "Stage2_glucose":second["phenotype"].get("glucose"),"Stage2_Rubisco":second["phenotype"].get("rubisco"),
          "Stage2_CO2":second["phenotype"].get("co2"),"Stage2_Fe2":second["phenotype"].get("fe2"),"Stage2_O2":second["phenotype"].get("o2"),
          "glucose_contribution_stage1":off1.get("status")=="optimal" and x["delta_biomass"]>TOL,
          "Stage1_glucose_off_biomass":off1.get("biomass"),"Stage1_delta_biomass":x["delta_biomass"],
          "Stage2_glucose_off_biomass":second["glucose_off"].get("biomass"),
          "Stage2_delta_biomass":second["phenotype"].get("biomass",0)-second["glucose_off"].get("biomass",0) if second["phenotype"]["status"]=="optimal" and second["glucose_off"].get("status")=="optimal" else None,
          "mandatory_under_stage1_FVA":";".join(sorted(set(x["reaction_ids"]) & set(mandatory))),
          "optional_or_pfba_selected":";".join(sorted(set(x["reaction_ids"])-set(mandatory))),
          "final_classification":cls})
    def priority_key(row):
        stage1_margin=float(row["Stage1_biomass"])-bmin
        if row["priority_group"].startswith("A"):
            margin=min(stage1_margin,float(row["Stage2_biomass"])-bmin)
            glucose_ok=bool(row["glucose_contribution_stage1"] and row["Stage2_delta_biomass"] is not None and row["Stage2_delta_biomass"]>TOL)
            delta=min(float(row["Stage1_delta_biomass"] or 0),float(row["Stage2_delta_biomass"] or 0))
        else:
            margin=stage1_margin
            glucose_ok=bool(row["glucose_contribution_stage1"])
            delta=float(row["Stage1_delta_biomass"] or 0)
        return (0 if row["priority_group"].startswith("A") else 1,row["reaction_count"],-margin,
                0 if glucose_ok else 1,-delta)
    final.sort(key=priority_key)
    write_tsv(OUT/"final_prioritized_sets.tsv",final)
    stage2_cp=OUT/"checkpoints/stage2_validation.json"
    save_json(stage2_cp,{"fingerprint":fingerprint,"terminal":True,"results":stage2})
    with (OUT/"logs/run_log.jsonl").open("a",encoding="utf-8") as log:
        log.write(json.dumps({"utc":datetime.now(timezone.utc).isoformat(),"fingerprint":fingerprint,
          "stage1_full_pool":full,"stage1_pfba":parsimonious,"support_count":len(supports),
          "stage2":stage2_rows},allow_nan=False,sort_keys=True)+"\n")
    write_report(manifest,"feasible",full,parsimonious,class_rows,supports,stage2_rows,final,strategy_log)
    print(json.dumps({"Stage1_status":"feasible","mandatory_count":len(mandatory),"pfba_active_count":len(active),
      "candidate_classifications":len(class_rows),"support_sets":len(supports),
      "PASS_BOTH":sum(x["priority_group"].startswith("A") for x in final),
      "PASS_30_ONLY":sum(x["priority_group"].startswith("B") for x in final)},indent=2))


def write_report(manifest,stage1_status,full,pfba_result,class_rows,supports,stage2,final,strategy_log):
    if stage1_status=="STAGE1_FULL_POOL_INFEASIBLE":
        body=["# E2R3 两阶段全池筛选","",f"Stage 1 全候选池不可行（{full['status']}），按流程停止；未运行 pFBA、FVA 或 Stage 2。"]
    else:
        cats={}
        for row in class_rows: cats[row["classification"]]=cats.get(row["classification"],0)+1
        body=["# E2R3 两阶段全池筛选","",f"- 权威宿主 SHA256：{manifest['inputs']['host']}",
          f"- E1 pool SHA256：{manifest['inputs']['E1_pool']}；候选行：{EXPECTED_ROWS}；构建后方向变体：{manifest['candidate_variants']}。",
          f"- WT glucose-closed FIM：B_WT={manifest['B_WT']:.12g}；pFBA Fe2={manifest['Fe_WT']:.12g}、O2={manifest['O2_WT']:.12g}。",
          f"- 资源匹配混合零 KI：biomass={manifest['mixed_reference']['biomass']:.12g}；glucose={manifest['mixed_reference']['glucose_uptake']:.12g}；R_REF={manifest['R_REF']:.12g}。",
          f"- 两阶段共同 biomass 下限={manifest['Stage1']['biomass_min']:.12g}；Rubisco cap：Stage 1={manifest['Stage1']['rubisco_max']:.12g}（30%），Stage 2={manifest['Stage2']['rubisco_max']:.12g}（10%）。",
          f"- Fe2 bounds=[-{manifest['Fe_WT']:.12g},0] 且 uptake≥0.01；O2 bounds=[-{manifest['O2_WT']:.12g},0]；glucose bounds=[-0.5,0] 且 uptake≥0.05；CO2/H2CO3=[-2,0]。",
          f"- Stage 1 全池 FBA：{full['status']}，biomass={full.get('biomass')}；pFBA：{pfba_result.get('status')}，biomass={pfba_result.get('biomass')}。",
          f"- FVA 分类（阈值容差 {TOL:g}）："+", ".join(f"{k}={v}" for k,v in sorted(cats.items())),
          "- Stage-1 pFBA 保留 biomass 下限和全部其余表型约束，以零主目标最小化总通量；它是一个简约通量代表。FVA mandatory 表示区间不含零。对 FVA 可活动范围最大的 32 个 optional variants，逐个强制方向一致的非零 flux 后再做阈值条件 pFBA，以发现需要多个共同反应的备选路线。",
          "- 每个支持集经 Stage 1 独立 FBA 重建，并逐反应删除验证不可再单删；该性质不代表全局 KI 数最小。", "",
          "- reaction_count 与优先排序按 donor-supported search reaction variants 计数；E1_candidate_count 按去重后的 E1 candidate_id 计数，两个数值均保存在支持集表。",
          "## Stage-1 validated irreducible support sets",""]
        for x in supports:
            ph=x["validation"]; off=x["glucose_off"]
            delta="NA" if x["delta_biomass"] is None else f"{x['delta_biomass']:.8g}"
            body.append(f"- **{x['support_set_id']}**（{len(x['reaction_ids'])} search variants；{x['strategy']}）：{'; '.join(x['reaction_ids'])}。B={ph['biomass']:.8g}，glucose={ph['glucose']:.8g}，Rubisco={ph['rubisco']:.8g}，Fe2={ph['fe2']:.8g}，O2={ph['o2']:.8g}；glucose-off B={off.get('biomass')}，ΔB={delta}。")
        body += ["","## Stage 2: same-set 10% challenge",""]
        for x in stage2:
            body.append(f"- {x['support_set_id']}：{x['final_classification']}；status={x['Stage2_status']}；B={x['Stage2_biomass']}，glucose={x['Stage2_glucose_uptake']}，Rubisco={x['Stage2_Rubisco_flux']}，glucose-off B={x['Stage2_glucose_off_biomass']}，ΔB={x['Stage2_delta_biomass']}。Stage 2 使用同一反应集合，未增加 KI。")
        body += ["","## 优先级",""]
        body += [f"- {x['priority_group']} / {x['support_set_id']}（{x['reaction_count']} reactions）：Stage1 B={x['Stage1_biomass']:.8g} / Rubisco={x['Stage1_Rubisco']:.8g}；Stage2 B={x['Stage2_biomass']} / Rubisco={x['Stage2_Rubisco']}；Stage1 glucose contribution={x['glucose_contribution_stage1']}。" for x in final]
        if not final: body.append("- 未产生通过 Stage 1 独立验证的支持集。")
        body += ["","## Reproducibility",f"- Input fingerprint：`{manifest['fingerprint']}`",
          f"- Donor model SHA256：`{manifest['inputs']['donors']}`",
          f"- Package versions：`{json.dumps(manifest['software_versions'],sort_keys=True)}`",
          f"- pFBA/FVA 设置：{manifest['pfba']}；{manifest['fva']}。",
          "- FVA 使用精确 Stage-1 表型约束（biomass 下限、Rubisco 上限、glucose 最低摄取、Fe2 trace、Fe2/O2 caps、CO2 bounds），fraction_of_optimum=0。",
          "- 结果、输入指纹检查点及单独 E2R3 日志均保存在本目录；未触碰 E1/E2/E2R 结果。"]
    (OUT/"E2R3_report.md").write_text("\n".join(body)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
