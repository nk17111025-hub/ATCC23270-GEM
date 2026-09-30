#!/usr/bin/env python
"""Resource-matched E2R; reuses E2 candidate construction and StrainDesign."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, os, re, time, sys
from importlib.metadata import version as package_version
from datetime import datetime, timezone
from pathlib import Path
import cobra
from cobra.flux_analysis import pfba

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/engineering/E2R_resource_matched_KI"
E2_PATH=ROOT/"scripts/run_E2_strain_design_KI.py"
HOST=ROOT.parent/"before final"/"updatedv3.1_verified"/"updatedv3.1.xml"
POOL=ROOT/"results/engineering/E1_donor_reaction_pool/candidate_reactions.tsv"
EXPECTED="026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3"
IDS={"biomass":"Ex_bio[e]","glucose":"Ex_glc-B[e]","rubisco":"RUBISCO",
     "co2":"Ex_h2co3[e]","fe2":"Ex_fe2[e]","o2":"Ex_o2[e]"}
BS=(1.,1.25,1.5); RS=(1.,.9,.75,.5,.25,.1,0.)
RETRY={"time_limit","error","validation_error"}

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def tree_fingerprint(directory):
    h=hashlib.sha256()
    for p in sorted(Path(directory).rglob("*.py")):
        h.update(p.relative_to(directory).as_posix().encode()); h.update(p.read_bytes())
    return h.hexdigest()

def write_tsv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    cols=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def save_json(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8"); tmp.replace(path)

def load_config(e2):
    os.chdir(ROOT); helper=e2.load_e1_helpers()
    host,donors,ai,ac=e2.load_models(helper)
    cfgpath=ROOT/"config/scenarios.json"; cfg=json.loads(cfgpath.read_text(encoding="utf-8"))
    return host,donors,ai,ac,cfg,cfgpath

def ids_for(model,e2):
    return {k:e2.resolve_rxn(model,v) for k,v in IDS.items()}

def native_condition(model,cfg,name,e2):
    for rid,b in cfg["conditions"][name]["exchange_bounds"].items(): e2.resolve_rxn(model,rid).bounds=tuple(b)
    # Established WT perturbation is glucose opening; use for FIM only.
    if name=="FIM":
        for rid,b in cfg.get("perturbations",{}).get("WT",{}).items(): e2.resolve_rxn(model,rid).bounds=tuple(b)

def fluxes(f,ids):
    return {"biomass":float(f[ids["biomass"].id]),"glucose":float(f[ids["glucose"].id]),
            "rubisco":float(f[ids["rubisco"].id]),"co2":float(f[ids["co2"].id]),
            "fe2":float(f[ids["fe2"].id]),"o2":float(f[ids["o2"].id])}

def constraints(ids,bmin,rcap):
    return [f"{ids['biomass'].id} >= {bmin:.12g}",f"{ids['glucose'].id} <= -0.05",
            f"{ids['fe2'].id} <= -0.01",f"{ids['rubisco'].id} >= 0",f"{ids['rubisco'].id} <= {rcap:.12g}"]

def bounds_closed(cands): return {x:(0.,0.) for x in cands}
def bounds_selected(m,cands,selected):
    s=set(selected)
    return {x:((m.reactions.get_by_id(x).lower_bound,m.reactions.get_by_id(x).upper_bound) if x in s else (0.,0.)) for x in cands}

def solve_fba(e2,m,ids,cons,candidate_bounds=None,off=False):
    # E2 helper enforces the scenario constraints using reversible model context.
    # Off phenotype comparison measures maximum biomass with glucose shut, retaining
    # Rubisco/resource caps, and omits only the biomass target and glucose floor.
    use=[x for x in cons if not (off and (x.startswith(ids["biomass"].id+" >=") or x.startswith(ids["glucose"].id+" <=")))]
    with m:
        if off: ids["glucose"].bounds=(0,0)
        if off:
            for rid,b in (candidate_bounds or {}).items(): m.reactions.get_by_id(rid).bounds=b
            for text in use:
                mt=re.fullmatch(r"(.+?)\s*(<=|>=|=)\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",text)
                rid,op,val=mt.groups(); val=float(val); rx=m.reactions.get_by_id(rid)
                lb,ub=(None,val) if op=="<=" else (val,None) if op==">=" else (val,val)
                m.add_cons_vars(m.problem.Constraint(rx.flux_expression,lb=lb,ub=ub,name="E2R_off_"+str(len(m.solver.constraints))))
            sol=m.optimize()
            if sol.status=="infeasible": return {"status":"infeasible","biomass":0.}
            return {"status":sol.status,**(fluxes(sol.fluxes,ids) if sol.status=="optimal" else {})}
        result=e2.phenotype_fba(m,ids,use,candidate_bounds)
        if result.get("status")!="optimal": return result
        return {"status":"optimal","biomass":result["biomass"],
                "glucose":result["glucose_uptake"],"rubisco":result["rubisco_flux"],
                "co2":result["co2_flux"],"fe2":result["fe2_flux"],"o2":result["o2_flux"]}

def sid(b,r):
    rv=format(float(r*100),".6g") if r else "0"
    if "." in rv: rv=rv.rstrip("0").rstrip(".")
    rb=rv.replace(".","p")
    return "B"+format(float(b),".6g").replace(".","p")+"_R"+rb

def run_one(defn,search,ids,cands,prov,rref,bwt,fp,limit,e2,resume=True):
    name=defn["scenario_id"]; b=defn["biomass_multiplier"]; r=defn["rubisco_fraction"]
    bmin=b*bwt; rcap=max(0,r*rref); cons=constraints(ids,bmin,rcap)
    cp=OUT/"checkpoints"/(name+".json")
    if resume and cp.exists():
        try:
            old=json.loads(cp.read_text(encoding="utf-8"))
            if old.get("fingerprint")==fp and old.get("status") not in RETRY: return old
        except Exception: pass
    t0=time.time(); closed=bounds_closed(cands)
    no=solve_fba(e2,search,ids,cons,closed)
    if no.get("status") not in ("optimal","infeasible"):
        data={"scenario_id":name,"stage":defn["stage"],"biomass_multiplier":b,
          "rubisco_fraction":r,"biomass_threshold":bmin,"rubisco_cap":rcap,
          "zero_KI_status":no.get("status"),"full_pool_status":"not_run",
          "status":"error","solver_status":"error","minimum_KI_count":None,"solutions":[],
          "runtime_seconds":time.time()-t0,"fingerprint":fp,"host_sha256":EXPECTED,
          "error":"No-KI phenotype LP returned neither optimal nor proven infeasible."}
        save_json(cp,data); return data
    if no.get("status")=="optimal":
        off=solve_fba(e2,search,ids,cons,closed,True)
        sols=[{"solution_id":"S00","reaction_ids":[],"candidate_ids":[],"source_reactions":[],
               "selected_KI_count":0,"phenotype":no,"glucose_off":off,
               "delta_biomass":no["biomass"]-off.get("biomass",0)}]
        data={"status":"zero_KI_feasible","zero_KI_status":"optimal","full_pool_status":"not_needed",
              "solver_status":"zero_KI_feasible","minimum_KI_count":0,"solutions":sols}
    else:
        full=solve_fba(e2,search,ids,cons)
        if full.get("status")!="optimal":
            final="infeasible_full_pool" if full.get("status")=="infeasible" else "error"
            data={"status":final,"zero_KI_status":no.get("status"),
                  "full_pool_status":full.get("status"),"solver_status":"not_run","minimum_KI_count":None,"solutions":[]}
            data.update({k:defn[k] for k in ("scenario_id","stage","biomass_multiplier","rubisco_fraction")})
            data.update(biomass_threshold=bmin,rubisco_cap=rcap,runtime_seconds=time.time()-t0,fingerprint=fp)
            save_json(cp,data); return data
        try:
            sd,dt=e2.sd_solve(search.copy(),cands,cons,6,limit,4)
            raw=str(getattr(sd,"status","unknown")).lower()
            found=e2.summarize_solutions(sd,prov); sols=[]
            for d in found:
                cb=bounds_selected(search,cands,d["reaction_ids"])
                ph=solve_fba(e2,search,ids,cons,cb)
                off=solve_fba(e2,search,ids,cons,cb,True)
                d.update(selected_KI_count=len(d["reaction_ids"]),phenotype=ph,glucose_off=off,
                         delta_biomass=(ph.get("biomass")-off.get("biomass")) if ph.get("status")=="optimal" and off.get("biomass") is not None else None)
                sols.append(d)
            if "time_limit" in raw: status="time_limit"
            elif raw=="infeasible": status="infeasible_within_6_KI"
            elif raw!="optimal": status="error"
            elif sols:
                mn=min(d["selected_KI_count"] for d in sols)
                sols=[d for d in sols if d["selected_KI_count"]==mn]
                for d in sols:
                    d["design_validated"]=(d["phenotype"].get("status")=="optimal" and
                      d.get("delta_biomass") is not None and d["delta_biomass"]>1e-7)
                rejected=[d for d in sols if not d["design_validated"]]
                sols=[d for d in sols if d["design_validated"]]
                status="optimal" if sols else "validation_error"
            else: status="validation_error"
            data={"status":status,"solver_status":status,"zero_KI_status":no.get("status"),
                  "full_pool_status":"feasible","minimum_KI_count":min((d["selected_KI_count"] for d in sols),default=None) if status=="optimal" else None,
                  "solutions":sols,"sd_runtime_seconds":dt}
            if "rejected" in locals(): data["rejected_solutions"]=rejected
        except Exception as ex:
            data={"status":"error","solver_status":"error","zero_KI_status":no.get("status"),
                  "full_pool_status":"feasible","minimum_KI_count":None,"solutions":[],"error":repr(ex)}
    data.update({k:defn[k] for k in ("scenario_id","stage","biomass_multiplier","rubisco_fraction")})
    data.update(biomass_threshold=bmin,rubisco_cap=rcap,runtime_seconds=time.time()-t0,fingerprint=fp,host_sha256=EXPECTED)
    save_json(cp,data)
    with (OUT/"logs/run_log.jsonl").open("a",encoding="utf-8") as f: f.write(json.dumps(data,sort_keys=True)+"\n")
    return data

def row_of(d):
    sol=(d.get("solutions") or [{}])[0]; ph=sol.get("phenotype") or {}; off=sol.get("glucose_off") or {}
    return {"scenario_id":d["scenario_id"],"stage":d["stage"],"biomass_multiplier":d["biomass_multiplier"],
      "biomass_threshold":d.get("biomass_threshold"),"rubisco_fraction":d["rubisco_fraction"],"rubisco_absolute_cap":d.get("rubisco_cap"),
      "zero_KI_status":d.get("zero_KI_status"),"full_pool_preflight_status":d.get("full_pool_status"),
      "solver_status":d.get("solver_status"),"runtime_seconds":d.get("runtime_seconds"),
      "host_sha256":EXPECTED,
      "minimum_KI_count":d.get("minimum_KI_count"),"optimal_alternatives_retained":len(d.get("solutions",[])),
      "biomass":ph.get("biomass"),"glucose_uptake":ph.get("glucose"),"glucose_off_biomass":off.get("biomass"),
      "delta_biomass":sol.get("delta_biomass"),"rubisco_flux":ph.get("rubisco"),"external_CO2_flux":ph.get("co2"),
      "Fe2_flux":ph.get("fe2"),"O2_flux":ph.get("o2"),
      "glucose_contribution_validated":sol.get("design_validated",sol.get("delta_biomass") is not None and sol["delta_biomass"]>1e-7)}

def materialize(defs,data):
    scenarios=[]; sets=[]
    for x in defs:
        d=data.get(x["scenario_id"])
        if not d: continue
        scenarios.append(row_of(d))
        if d.get("solver_status") not in ("optimal","zero_KI_feasible"): continue
        for s in d.get("solutions",[]):
            ph=s.get("phenotype") or {}; off=s.get("glucose_off") or {}
            sets.append({"scenario_id":d["scenario_id"],"biomass_multiplier":d["biomass_multiplier"],
              "rubisco_fraction":d["rubisco_fraction"],"minimum_KI_count":d.get("minimum_KI_count"),
              "host_sha256":EXPECTED,
              "solution_id":s["solution_id"],"E1_candidate_ids":";".join(s["candidate_ids"]),
              "search_reaction_ids":";".join(s["reaction_ids"]),"donor_source_reaction_ids":";".join(s["source_reactions"]),
              "rubisco_absolute_cap":d.get("rubisco_cap"),
              "biomass":ph.get("biomass"),"glucose_uptake":ph.get("glucose"),"glucose_off_biomass":off.get("biomass"),
              "delta_biomass":s.get("delta_biomass"),"rubisco_flux":ph.get("rubisco"),"CO2_flux":ph.get("co2"),
              "Fe2_flux":ph.get("fe2"),"O2_flux":ph.get("o2"),"solver_status":d.get("solver_status"),
              "glucose_contribution_validated":s.get("delta_biomass") is not None and s["delta_biomass"]>1e-7,
              "design_validated":s.get("design_validated",s.get("delta_biomass") is not None and s["delta_biomass"]>1e-7),
              "design_FBA_valid":ph.get("status")=="optimal"})
    eligible=[x for x in sets if x["design_FBA_valid"] and x["design_validated"] and
              x["solver_status"] in ("optimal","zero_KI_feasible")]
    frontier=[]
    biomass_tol=1e-6
    for x in eligible:
        b=float(x["biomass"]); r=float(x["rubisco_absolute_cap"]); k=int(x["minimum_KI_count"] or 0)
        if not any(float(y["biomass"])>=b-biomass_tol and float(y["rubisco_absolute_cap"])<=r+1e-9 and int(y["minimum_KI_count"] or 0)<=k and
          (float(y["biomass"])>b+biomass_tol or float(y["rubisco_absolute_cap"])<r-1e-9 or int(y["minimum_KI_count"] or 0)<k) for y in eligible if y is not x): frontier.append(x)
    # Collapse repeated scenario labels that encode the same engineering point.
    # Prefer the most stringent biomass target when several labels produce the
    # same KI set, count, biomass (1e-6 tolerance), and Rubisco cap.
    unique=[]
    for x in sorted(frontier,key=lambda q:(-float(q["biomass_multiplier"]),q["scenario_id"])):
        same=next((y for y in unique if y["E1_candidate_ids"]==x["E1_candidate_ids"] and
          int(y["minimum_KI_count"] or 0)==int(x["minimum_KI_count"] or 0) and
          abs(float(y["biomass"])-float(x["biomass"]))<=1e-6 and
          abs(float(y["rubisco_absolute_cap"])-float(x["rubisco_absolute_cap"]))<=1e-9),None)
        if same is None: unique.append(x)
    frontier=unique
    write_tsv(OUT/"scenario_results.tsv",scenarios); write_tsv(OUT/"minimal_KI_sets.tsv",sets); write_tsv(OUT/"pareto_frontier.tsv",frontier)
    return scenarios,sets,frontier

def refinements(rows):
    look={(float(x["biomass_multiplier"]),float(x["rubisco_fraction"])):x for x in rows}; cand=[]
    good={"optimal","zero_KI_feasible"}
    proven=good|{"infeasible_full_pool","infeasible_within_6_KI"}
    for b in BS:
      for hi,lo in zip(RS[:-1],RS[1:]):
        a,z=look.get((b,hi)),look.get((b,lo))
        if a and z and a["solver_status"] in proven and z["solver_status"] in proven and ((a["solver_status"] in good)!=(z["solver_status"] in good) or
          (a["minimum_KI_count"] is not None and z["minimum_KI_count"] is not None and a["minimum_KI_count"]!=z["minimum_KI_count"])):
            feasible=(a["solver_status"] in good)!=(z["solver_status"] in good)
            jump=a["minimum_KI_count"] is not None and z["minimum_KI_count"] is not None and abs(a["minimum_KI_count"]-z["minimum_KI_count"])>=2
            cand.append((0 if feasible else 1 if jump else 2,b,(hi+lo)/2,
                         "feasibility_boundary" if feasible else "KI_jump" if jump else "KI_or_pareto_transition",
                         a["scenario_id"],z["scenario_id"]))
    for hi,lo in zip(BS[:-1],BS[1:]):
      for r in RS:
        a,z=look.get((hi,r)),look.get((lo,r))
        if a and z and a["solver_status"] in proven and z["solver_status"] in proven and ((a["solver_status"] in good)!=(z["solver_status"] in good) or
          (a["minimum_KI_count"] is not None and z["minimum_KI_count"] is not None and abs(a["minimum_KI_count"]-z["minimum_KI_count"])>=1)):
            feasible=(a["solver_status"] in good)!=(z["solver_status"] in good)
            jump=a["minimum_KI_count"] is not None and z["minimum_KI_count"] is not None and abs(a["minimum_KI_count"]-z["minimum_KI_count"])>=2
            cand.append((0 if feasible else 1 if jump else 2,(hi+lo)/2,r,
                         "feasibility_boundary" if feasible else "KI_jump" if jump else "KI_or_pareto_transition",
                         a["scenario_id"],z["scenario_id"]))
    out=[]; seen=set()
    for priority,b,r,reason,left,right in sorted(cand):
      if (b,r) not in seen:
        seen.add((b,r)); out.append({"scenario_id":sid(b,r),"stage":"refined","biomass_multiplier":b,"rubisco_fraction":r,
                                     "reason":reason,"adjacent_scenario_ids":left+";"+right,"priority":priority})
      if len(out)==12: break
    return out

def ttm_tsm(search,cfg,cands,front,e2,o2_wt,b_wt):
    if not front:
        (OUT/"ttm_tsm_validation.tsv").write_text("scenario_id\tcondition\tvalidation_mode\tstatus\tbiomass\n",encoding="utf-8")
        return []
    valid=[x for x in front if x.get("design_validated")]
    if not valid:
        (OUT/"ttm_tsm_validation.tsv").write_text("scenario_id\tcondition\tvalidation_mode\tstatus\tbiomass\n",encoding="utf-8")
        return []
    high=max(valid,key=lambda x:float(x["biomass"]))
    knee_candidates=[x for x in valid if float(x["biomass"])>=1.5*b_wt]
    knee=min(knee_candidates,key=lambda x:float(x["rubisco_absolute_cap"])) if knee_candidates else None
    strongest=min(valid,key=lambda x:float(x["rubisco_absolute_cap"]))
    chosen={x["scenario_id"]:x for x in (high,knee,strongest) if x is not None}
    out=[]
    for x in chosen.values():
      for cond,donor in (("TTM","Ex_ttton[e]"),("TSM","Ex_tsul[e]")):
        for mode in ("native_Rubisco_unrestricted","strict_FIM_cap","donor_closed_control"):
          m=search.copy()
          for rid,b in cfg["conditions"][cond]["exchange_bounds"].items(): e2.resolve_rxn(m,rid).bounds=tuple(b)
          ids=ids_for(m,e2); ids["fe2"].bounds=(-.01,0); ids["o2"].bounds=(-o2_wt,0)
          ids["glucose"].bounds=(-.5,0); ids["co2"].bounds=(-2,-2)
          dr=e2.resolve_rxn(m,donor); dr.bounds=(min(dr.lower_bound,-1000),0)
          selected=set(filter(None,x["search_reaction_ids"].split(";")))
          for rid in cands:
            rx=m.reactions.get_by_id(rid)
            if rid not in selected: rx.bounds=(0,0)
          if mode=="strict_FIM_cap":
            ids["rubisco"].upper_bound=min(ids["rubisco"].upper_bound,float(x["rubisco_absolute_cap"]))
          if mode=="donor_closed_control":
            dr.bounds=(0,0)
          else:
            m.add_cons_vars(m.problem.Constraint(dr.flux_expression,ub=-.01,name="E2R_sulfur_donor_uptake"))
          m.objective=ids["biomass"]; sol=m.optimize()
          out.append({"scenario_id":x["scenario_id"],"minimum_KI_count":x["minimum_KI_count"],
            "condition":cond,"validation_mode":mode,"rubisco_cap":float(x["rubisco_absolute_cap"]) if mode=="strict_FIM_cap" else "unrestricted",
            "status":sol.status,"biomass":float(sol.fluxes[ids["biomass"].id]) if sol.status=="optimal" else None,
            "Fe2_bounds":"[-0.01,0] trace allowance","O2_lower_bound":-o2_wt,
            "sulfur_donor":donor,"sulfur_donor_flux":float(sol.fluxes[dr.id]) if sol.status=="optimal" else None,
            "glucose_flux":float(sol.fluxes[ids["glucose"].id]) if sol.status=="optimal" else None,
            "CO2_flux":float(sol.fluxes[ids["co2"].id]) if sol.status=="optimal" else None,
            "rubisco_flux":float(sol.fluxes[ids["rubisco"].id]) if sol.status=="optimal" else None})
    write_tsv(OUT/"ttm_tsm_validation.tsv",out); return out

def make_report(manifest,rows,sets,front,refine,tvals):
    coarse=[x for x in rows if x["stage"]=="coarse"]
    refined=[x for x in rows if x["stage"]=="refined"]
    count=lambda records,status:sum(x.get("solver_status")==status for x in records)
    nonzero=sum(1 for x in rows if x.get("minimum_KI_count") not in (None,0))
    boundary=[]
    good={"optimal","zero_KI_feasible"}
    cmap={(float(x["biomass_multiplier"]),float(x["rubisco_fraction"])):x for x in coarse}
    for b in BS:
      for hi,lo in zip(RS[:-1],RS[1:]):
        a,z=cmap.get((b,hi)),cmap.get((b,lo))
        if a and z and a["solver_status"] in good and z["solver_status"] not in good and z["solver_status"] in {"infeasible_full_pool","infeasible_within_6_KI"}:
          boundary.append(f"biomass ×{b:g}: Rubisco {hi*100:g}% → {lo*100:g}%")
    for hi,lo in zip(BS[:-1],BS[1:]):
      for r in RS:
        a,z=cmap.get((hi,r)),cmap.get((lo,r))
        if a and z and a["solver_status"] in good and z["solver_status"] not in good and z["solver_status"] in {"infeasible_full_pool","infeasible_within_6_KI"}:
          boundary.append(f"Rubisco {r*100:g}%: biomass ×{hi:g} → ×{lo:g}")
    threshold=1.5*float(manifest["B_WT"])
    knee_pool=[x for x in front if float(x["biomass"])>=threshold]
    knee=min(knee_pool,key=lambda x:float(x["rubisco_absolute_cap"])) if knee_pool else None
    nextcap=None
    if knee:
      stronger=[x for x in front if float(x["rubisco_absolute_cap"])<float(knee["rubisco_absolute_cap"])-1e-10]
      if stronger: nextcap=max(stronger,key=lambda x:float(x["rubisco_absolute_cap"]))
    all_kis={int(x["minimum_KI_count"] or 0) for x in sets}
    total_infeasible=count(rows,"infeasible_within_6_KI")
    full_pool_ok=sum(x.get("solver_status")=="infeasible_within_6_KI" and x.get("full_pool_preflight_status")=="feasible" for x in rows)
    lines=["# E2R 资源匹配 KI 搜索","",f"- 权威宿主 SHA256：{EXPECTED}",
      f"- WT glucose-closed FIM pFBA：B_WT={manifest['B_WT']:.9g}；Fe2={manifest['WT']['fe2']:.9g}；O2={manifest['WT']['o2']:.9g}；Rubisco={manifest['WT']['rubisco']:.9g}；CO2={manifest['WT']['co2']:.9g}。",
      f"- 混合无 KI：biomass={manifest['mixed']['biomass']:.9g}；glucose={manifest['mixed']['glucose_uptake']:.9g}；R_REF={manifest['R_REF']:.9g}；CO2={manifest['mixed']['co2_flux']:.9g}；Fe2={manifest['mixed']['fe2_flux']:.9g}；O2={manifest['mixed']['o2_flux']:.9g}。",
      f"- 资源上限：Fe2≥{-manifest['Fe_WT']:.9g}、O2≥{-manifest['O2_WT']:.9g}；CO2/H2CO3为[-2,0]（以当前 FIM 的已定义无机碳供应上限为范围）；glucose 为[-0.5,0]且至少摄取0.05。历史配置 perturbation 中名为 WT 的状态开启了 glucose；本报告的 B_WT 明确使用 glucose-closed FIM。",
      "- FIM 情景还要求至少 0.01 Fe2 uptake，以保留铁氧化条件。",
      "- StrainDesign KI-only：显式 ko_cost={}，单位 ki_cost，max_cost=6，best，SCIP，seed=1；已通过正控。","",
      "## 结果概览",
      f"- 粗网格：{len(coarse)} 个；无 KI 可行 {count(coarse,'zero_KI_feasible')} 个；六 KI 内搜索确认无解 {count(coarse,'infeasible_within_6_KI')} 个；全池 LP 不可行 {count(coarse,'infeasible_full_pool')} 个；超时 {count(coarse,'time_limit')} 个。",
      f"- 细化：{len(refined)} 个；无 KI 可行 {count(refined,'zero_KI_feasible')} 个；六 KI 内无解 {count(refined,'infeasible_within_6_KI')} 个；超时 {count(refined,'time_limit')} 个。",
      f"- 当前共确认 {nonzero} 个非零 KI 最小解。{total_infeasible} 个场景的 StrainDesign 搜索在六 KI 内确认无解，其中 {full_pool_ok} 个全池 LP 可行；这些场景的设计若存在，需超过六个 KI。",
      f"- 主要边界：{'; '.join(boundary) if boundary else '粗网格未发现已证实的相邻可行性边界。'}",
      "- 六 KI 内搜索无解仅表示当前上限内无解；E2R 不据此推断更高 KI 数仍不可行。","",
      "## 场景摘要","","| 场景 | 阶段 | B倍数 | Rubisco比例 | 状态 | KI数 | B | glucose | glucose-off B | ΔB | Rubisco | CO2 | Fe2 | O2 |",
      "|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for x in rows: lines.append("| "+" | ".join(str(x.get(k,"")) for k in ("scenario_id","stage","biomass_multiplier","rubisco_fraction","solver_status","minimum_KI_count","biomass","glucose_uptake","glucose_off_biomass","delta_biomass","rubisco_flux","external_CO2_flux","Fe2_flux","O2_flux"))+" |")
    lines += ["",f"## 自适应细化","",f"根据粗网格相邻可行性或 KI 数变化增加 {len(refine)} 个场景（最多 12）。"]
    lines += [f"- {x['scenario_id']}：biomass×{x['biomass_multiplier']}，Rubisco {x['rubisco_fraction']*100:g}%" for x in refine]
    lines += ["","## 最小 KI 组合"]
    lines += [f"- {x['scenario_id']} / {x['minimum_KI_count']} KI / {x['solution_id']}：候选 {x['E1_candidate_ids']}，供体 {x['donor_source_reaction_ids']}；B={x['biomass']}，glucose={x['glucose_uptake']}，glucose-off B={x['glucose_off_biomass']}，ΔB={x['delta_biomass']}，贡献验证={x['glucose_contribution_validated']}。" for x in sets]
    lines += ["","## Pareto 前沿与候选拐点",f"去重后非支配工程点 {len(front)} 个。"]
    lines += [f"- {x['scenario_id']}：{x['minimum_KI_count']} KI，B={x['biomass']}，Rubisco cap={x['rubisco_absolute_cap']}，实际 Rubisco={x['rubisco_flux']}，glucose-off B={x['glucose_off_biomass']}，ΔB={x['delta_biomass']}，CO2={x['CO2_flux']}。" for x in front]
    if front:
      simplest=min(front,key=lambda x:(int(x["minimum_KI_count"]),float(x["rubisco_absolute_cap"])))
      strongest=min(front,key=lambda x:float(x["rubisco_absolute_cap"]))
      lines.append(f"- 最低复杂度前沿点：{simplest['scenario_id']}（{simplest['minimum_KI_count']} KI）；最强 CBB 限制前沿点：{strongest['scenario_id']}（Rubisco cap={strongest['rubisco_absolute_cap']}）。")
    if knee:
      lines.append(f"- 高增长拐点：{knee['scenario_id']}，Rubisco cap={knee['rubisco_absolute_cap']}，biomass={knee['biomass']}（目标阈值 {threshold:.6g}）。")
      if nextcap:
        loss=float(knee["biomass"])-float(nextcap["biomass"])
        lines.append(f"- 再加强 Rubisco 限制至 {nextcap['rubisco_absolute_cap']} 后，biomass 降至 {nextcap['biomass']}，减少 {loss:.6g}（{loss/max(float(knee['biomass']),1e-12)*100:.2f}%）。")
    if all_kis and all(k==0 for k in all_kis):
      lines.append("- KI 数量拐点：本次有效解均为 0 KI，暂未出现可比较的 KI 复杂度拐点。")
    lines += ["","## TTM/TSM 固定 KI 验证",
      "抽查高增长点、满足 1.5×B_WT 的 Rubisco 限制拐点和前沿中最强 Rubisco 限制点（相同点合并）。native_Rubisco_unrestricted 使用 TTM/TSM 原条件，不施加 FIM Rubisco cap；strict_FIM_cap 额外施加 FIM 场景 cap；donor_closed_control 关闭对应硫供体。三种模式都使用 Fe2 微量营养范围、FIM O2 资源上限、glucose 上限和固定 CO2。"]
    lines += [f"- {x['scenario_id']} / {x['condition']} / {x['validation_mode']}：{x['status']}，biomass={x['biomass']}，Rubisco={x.get('rubisco_flux')}，CO2={x.get('CO2_flux')}，{x['sulfur_donor']} flux={x['sulfur_donor_flux']}。" for x in tvals]
    if not tvals: lines.append("- 当前没有通过 glucose-off 验证的 FIM 设计，跳过 TTM/TSM 固定 KI 验证。")
    unresolved=[x for x in rows if x["solver_status"] in ("time_limit","error","validation_error")]
    lines += ["",f"未解决场景 {len(unresolved)} 个；超时均保留为未解决，未记为不可行。",""]
    (OUT/"E2R_report.md").write_text("\n".join(lines),encoding="utf-8")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--baseline",action="store_true",help="references and KI semantics check only")
    ap.add_argument("--run",action="store_true"); ap.add_argument("--time-limit",type=int,default=900)
    ap.add_argument("--no-resume",action="store_true"); args=ap.parse_args()
    global E2
    E2=module(E2_PATH,"e2_helpers")
    for p in (OUT,OUT/"checkpoints",OUT/"logs"): p.mkdir(parents=True,exist_ok=True)
    observed=sha(HOST)
    if observed!=EXPECTED: raise SystemExit(f"Host SHA256 mismatch: {observed}")
    rows=E2.read_tsv(POOL)
    if len(rows)!=2424: raise SystemExit(f"Expected 2424 E1 rows; found {len(rows)}")
    host,donors,ai,ac,cfg,cfgpath=load_config(E2)
    wt=host.copy(); native_condition(wt,cfg,"FIM",E2); wi=ids_for(wt,E2)
    wi["glucose"].bounds=(0,0); wi["co2"].bounds=(-2,-2)
    wt.objective=wi["biomass"]
    wtres=wt.optimize()
    if wtres.status!="optimal": raise RuntimeError(f"WT FIM biomass FBA status {wtres.status}")
    bwt=float(wtres.fluxes[wi["biomass"].id])
    wt.objective=wi["biomass"]
    psol=pfba(wt,fraction_of_optimum=1.0)
    wtvals=fluxes(psol.fluxes,wi); fe_wt=abs(wtvals["fe2"]); o2_wt=abs(wtvals["o2"])
    search=host.copy(); native_condition(search,cfg,"FIM",E2); ids=ids_for(search,E2)
    ids["fe2"].bounds=(-fe_wt,0); ids["o2"].bounds=(-o2_wt,0)
    ids["co2"].bounds=(-2,0); ids["glucose"].bounds=(-.5,0)
    mixed=E2.fba_summary(search,ids); rref=mixed["rubisco_flux"]
    if abs(rref)<1e-8: raise RuntimeError("R_REF≈0; percentage cap grid cannot be applied")
    cands,prov,variants=E2.add_candidate_reactions(search,donors,ai,ac,rows)
    sem=E2.smoke_sd_ki_semantics()
    if not sem.get("passed"): raise RuntimeError(f"KI semantics positive control failed: {sem}")
    base=[{"host_sha256":observed,"B_WT":bwt,"WT_pFBA_Fe2_flux":wtvals["fe2"],"WT_pFBA_O2_flux":wtvals["o2"],
      "WT_pFBA_Rubisco_flux":wtvals["rubisco"],"WT_pFBA_CO2_flux":wtvals["co2"],
      "mixed_no_KI_biomass":mixed["biomass"],"mixed_no_KI_glucose_flux":mixed["glucose_uptake"],
      "R_REF":rref,"mixed_no_KI_CO2_flux":mixed["co2_flux"],"mixed_no_KI_Fe2_flux":mixed["fe2_flux"],
      "mixed_no_KI_O2_flux":mixed["o2_flux"],"glucose_cap":.5,"glucose_minimum_use":.05,
      "Fe2_lower_bound":-fe_wt,"O2_lower_bound":-o2_wt,"CO2_bounds":"[-2,0]","minimum_Fe2_uptake_guard":0.01}]
    write_tsv(OUT/"baseline_summary.tsv",base)
    save_json(OUT/"checkpoints/KI_semantics_positive_control.json",sem)
    donorsrc={k:sha(v) for k,v in E2.DONORS.items()}
    sd_tree=tree_fingerprint(E2.STRAINDESIGN/"straindesign")
    try: sd_version=package_version("straindesign")
    except Exception: sd_version="unknown"
    definition={"B":BS,"R":RS,"Fe_WT":fe_wt,"O2_WT":o2_wt,"R_REF":rref,"glucose":[-.5,0],"co2":[-2,0]}
    fp=hashlib.sha256(json.dumps({"host":observed,"pool":sha(POOL),"config":sha(cfgpath),
      "script":sha(Path(__file__)),"E2_helper":sha(E2_PATH),"StrainDesign_source":sd_tree,
      "StrainDesign_version":sd_version,"donors":donorsrc,"scenario_definition":definition},sort_keys=True).encode()).hexdigest()
    manifest={"host":str(HOST),"host_sha256":observed,"E1_rows":len(rows),"candidate_variants":len(cands),
      "WT":wtvals,"B_WT":bwt,"Fe_WT":fe_wt,"O2_WT":o2_wt,"mixed":{k:mixed[k] for k in ("biomass","glucose_uptake","rubisco_flux","co2_flux","fe2_flux","o2_flux")},
      "R_REF":rref,"fingerprint":fp,"KI_positive_control":sem,"E2_helper_sha256":sha(E2_PATH),
      "StrainDesign_source_sha256":sd_tree,"StrainDesign_version":sd_version,
      "donor_files":{k:str(v) for k,v in E2.DONORS.items()},
      "solver":"StrainDesign SCIP, KI-only, ko_cost empty, unit costs, best max_cost=6, seed 1"}
    save_json(OUT/"checkpoints/input_manifest.json",manifest)
    if args.baseline and not args.run:
        print(json.dumps(manifest,indent=2)); return
    if not args.run: print(json.dumps(manifest,indent=2)); return
    defs=[{"scenario_id":sid(b,r),"stage":"coarse","biomass_multiplier":b,"rubisco_fraction":r} for b in BS for r in RS]
    data={}
    for d in defs:
        data[d["scenario_id"]]=run_one(d,search,ids,cands,prov,rref,bwt,fp,args.time_limit,E2,not args.no_resume)
        write_tsv(OUT/"scenario_results.tsv",[row_of(data[x["scenario_id"]]) for x in defs if x["scenario_id"] in data])
    coarse=[row_of(data[x["scenario_id"]]) for x in defs]
    refine=refinements(coarse); write_tsv(OUT/"refinement_scenarios.tsv",refine)
    alldefs=defs+refine
    for d in refine:
        if d["scenario_id"] not in data or args.no_resume:
            data[d["scenario_id"]]=run_one(d,search,ids,cands,prov,rref,bwt,fp,args.time_limit,E2,not args.no_resume)
    scenarios,sets,front=materialize(alldefs,data)
    tvals=ttm_tsm(search,cfg,cands,front,E2,o2_wt,bwt)
    make_report(manifest,scenarios,sets,front,refine,tvals)
    print(json.dumps({"scenarios":len(scenarios),"refinements":len(refine),"KI_sets":len(sets),"pareto":len(front),"TTM_TSM":tvals},indent=2))

if __name__=="__main__": main()
