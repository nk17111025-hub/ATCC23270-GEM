#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv, datetime as dt, hashlib, json, platform, subprocess, sys
from pathlib import Path
import xlrd
ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
PAPER = ROOT / "02_论文" / "02_基因组与代谢模型" / "2016主要模型GEM.pdf"
OUT = ROOT / "baseline_v0"
VENV_PY = ROOT / "Archive_历史版本与废弃文件" / "A0" / "运行环境" / ".venv" / "Scripts" / "python.exe"
def now_iso():
    return dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds")
def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()
for name in ["original", "environment", "scenario_definition", "qc", "reproduction_2016", "fba", "pfba", "fva", "sensitivity", "logs"]:
    (OUT / name).mkdir(parents=True, exist_ok=True)
SOURCES = [
    ("2016 iMC507 SBML model (frozen network)", ORIG / "mmc3.xml"),
    ("2016 iMC507 reaction table (Supplementary Table 1, mmc1.xls)", ORIG / "mmc1.xls"),
    ("2016 supplementary Word file (Supplementary Tables/Figures, mmc2.docx)", ORIG / "mmc2.docx"),
    ("2016 Campodonico et al. paper (Metabolic Engineering Communications 3:84-96)", PAPER),
]
manifest = []
for role, src in SOURCES:
    dst = OUT / "original" / src.name
    data = src.read_bytes()
    dst.write_bytes(data)
    sh = sha256(src)
    fh = sha256(dst)
    manifest.append(dict(role=role, original_source_path=str(src.relative_to(ROOT)), frozen_copy_path=str(dst.relative_to(ROOT)), size_bytes=src.stat().st_size, sha256=sh, frozen_sha256=fh, match=(sh == fh)))
fields = ["role", "original_source_path", "frozen_copy_path", "size_bytes", "sha256", "frozen_sha256", "match"]
with (OUT / "input_manifest.tsv").open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in manifest:
        w.writerow(r)
meta = dict(generated_at=now_iso(), python_executable=str(VENV_PY), python_version=platform.python_version(), python_version_full=sys.version, platform=platform.platform(), machine=platform.machine(), processor=platform.processor())
import cobra, optlang, libsbml, openpyxl
meta["cobra_version"] = cobra.__version__
meta["optlang_version"] = optlang.__version__
meta["python_libsbml_version"] = libsbml.getLibSBMLVersionString()
meta["xlrd_version"] = xlrd.__version__
meta["openpyxl_version"] = openpyxl.__version__
try:
    import swiglpk
    meta["swiglpk_version"] = getattr(swiglpk, "__version__", "n/a")
    meta["glpk_version"] = swiglpk.glp_version()
except Exception as e:
    meta["swiglpk_error"] = repr(e)
try:
    from cobra.util.solver import solvers
    meta["cobra_solvers_available"] = sorted(solvers.keys())
except Exception as e:
    meta["cobra_solvers_error"] = repr(e)
try:
    cfg = cobra.Configuration()
    meta["cobra_config"] = dict(solver=str(cfg.solver), tolerance=str(cfg.tolerance), lower_bound=str(cfg.lower_bound), upper_bound=str(cfg.upper_bound), processes=str(cfg.processes))
except Exception as e:
    meta["cobra_config_error"] = repr(e)
meta["solver_tolerances"] = dict(status="UNKNOWN", note="optlang GLPK interface does not expose GLPK primal/dual feasibility tolerance or MILP gap tolerance through cobra.Configuration(). Left UNKNOWN until a direct GLPK probe in Phase 1.")
(OUT / "environment" / "environment_snapshot.txt").write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
proc = subprocess.run([str(VENV_PY), "-m", "pip", "freeze"], capture_output=True, text=True, encoding="utf-8", errors="replace")
(OUT / "environment" / "requirements_locked.txt").write_text(proc.stdout or "", encoding="utf-8")
book = xlrd.open_workbook(str(ORIG / "mmc1.xls"), on_demand=True)
s = book.sheet_by_name("Table 1")
table1 = []
for r in range(2, s.nrows):
    rid = s.cell_value(r, 0)
    if rid == "" or rid is None:
        continue
    def cv(c):
        v = s.cell_value(r, c)
        if v == "":
            return ""
        if isinstance(v, float) and v.is_integer():
            return int(v)
        return v
    table1.append(dict(reaction_id=str(rid), reaction_name=str(cv(1)), equation=str(cv(2)), confidence_level=cv(3), ec_number=str(cv(4)), pmid=cv(5), subsystem=str(cv(6)), gene_reaction_association=str(cv(7)), gene_protein_reaction_association=str(cv(8)), protein_reaction_association=str(cv(9)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15), source_row=r + 1))
book.release_resources()
gpr_fields = ["reaction_id", "reaction_name", "equation", "confidence_level", "ec_number", "pmid", "subsystem", "gene_reaction_association", "gene_protein_reaction_association", "protein_reaction_association", "lb_fe2", "ub_fe2", "lb_ttton", "ub_ttton", "lb_tsul", "ub_tsul", "source_row"]
with (OUT / "scenario_definition" / "gpr_2016_reference.tsv").open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=gpr_fields, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in table1:
        w.writerow(r)
exchange = [r for r in table1 if r["reaction_id"].startswith("Ex_")]
biomass = [r for r in table1 if not r["reaction_id"].startswith("Ex_") and ("biomass" in r["reaction_id"].lower() or r["subsystem"] == "Exchange")]
transport = [r for r in table1 if r not in exchange and r not in biomass and "Transport" in r["subsystem"]]
metabolic = [r for r in table1 if r not in exchange and r not in biomass and r not in transport]
classification = dict(total=len(table1), exchange=len(exchange), biomass_objective=len(biomass), transport=len(transport), metabolic=len(metabolic), non_exchange=len(table1) - len(exchange), exchange_ids=sorted(r["reaction_id"] for r in exchange))
(OUT / "qc" / "reaction_classification.json").write_text(json.dumps(classification, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
book = xlrd.open_workbook(str(ORIG / "mmc1.xls"), on_demand=True)
s5 = book.sheet_by_name("Table 5")
points = []
def isnum(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)
for r in range(3, s5.nrows):
    def cv5(c):
        v = s5.cell_value(r, c)
        return v if v != "" else None
    fe2_mu, fe2_up, fe2_o2, fe2_co2 = cv5(0), cv5(1), cv5(2), cv5(3)
    s4_mu, s4_up, s4_co2 = cv5(5), cv5(6), cv5(7)
    s2_mu, s2_up, s2_co2 = cv5(9), cv5(10), cv5(11)
    if isnum(fe2_mu):
        points.append(dict(condition="FIM", electron_donor="fe2", mu_h_minus1=fe2_mu, electron_donor_uptake_mmol_gDW_h=fe2_up, o2_uptake_mmol_gDW_h=fe2_o2, co2_consumption_mmol_gDW_h=fe2_co2))
    if isnum(s4_mu):
        points.append(dict(condition="TTM", electron_donor="s4o6", mu_h_minus1=s4_mu, electron_donor_uptake_mmol_gDW_h=s4_up, o2_uptake_mmol_gDW_h=None, co2_consumption_mmol_gDW_h=s4_co2))
    if isnum(s2_mu):
        points.append(dict(condition="TSM", electron_donor="s2o3", mu_h_minus1=s2_mu, electron_donor_uptake_mmol_gDW_h=s2_up, o2_uptake_mmol_gDW_h=None, co2_consumption_mmol_gDW_h=s2_co2))
book.release_resources()
exchange_rows = [r for r in table1 if r["reaction_id"].startswith("Ex_")]
def bounds(r, cond):
    return dict(lower_bound=r["lb_" + cond], upper_bound=r["ub_" + cond])
minimal_media = [dict(reaction_id=r["reaction_id"], equation=r["equation"], FIM=bounds(r, "fe2"), TTM=bounds(r, "ttton"), TSM=bounds(r, "tsul")) for r in exchange_rows]
parameters = []
def add(pid, meaning, *, value=None, units="", lower_bound=None, upper_bound=None, conditions=None, reported=False, reported_text="", model_default=False, model_default_text="", evidence_status="UNKNOWN", source_file="", source_section="", notes=""):
    parameters.append(dict(parameter_or_reaction_id=pid, biological_meaning=meaning, value=value, units=units, lower_bound=lower_bound, upper_bound=upper_bound, conditions=conditions, reported_in_paper=reported, reported_value=reported_text, model_default=model_default, model_default_value=model_default_text, evidence_status=evidence_status, source_file=source_file, source_section=source_section, notes=notes))
add("Ex_h2co3[e]", "inorganic carbon uptake (CO2/H2CO3) - unique carbon source", value=-2, units="mmol/gDW/h", lower_bound=-2, upper_bound=-2, conditions=["FIM", "TTM", "TSM"], reported=True, reported_text="h2co3 uptake rate set at 2 mmol/gDW/h (Methods 2.4; Fig captions)", model_default=True, model_default_text="mmc1.xls Table 1 Ex_h2co3[e] lb=ub=-2 for fe2/ttton/tsul", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 cols 11-16; Methods 2.4", notes="REPORTED and MODEL_DEFAULT agree at -2/-2.")
add("Ex_o2[e]", "oxygen uptake (aerobic respiration)", lower_bound=-1000, upper_bound=1000, conditions=["FIM", "TTM", "TSM"], reported=True, reported_text="aerobic conditions; external O2 exchange set to -1000/1000 (Methods 2.4)", model_default=True, model_default_text="mmc1.xls Table 1 Ex_o2[e] -1000/1000", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 cols 11-16; Methods 2.4")
add("Ex_fe2[e]", "ferrous iron electron donor (FIM)", lower_bound=-1000, upper_bound=1000, conditions=["FIM"], reported=True, reported_text="fe2 as electron donor under FIM (Methods 2.4)", model_default=True, model_default_text="mmc1.xls Table 1 Ex_fe2[e] -1000/1000 (FIM)", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 cols 11-12; Methods 2.4")
add("Ex_ttton[e]", "tetrathionate electron donor (TTM)", lower_bound=-1000, upper_bound=0, conditions=["TTM"], reported=True, reported_text="tetrathionate as electron donor under TTM (Methods 2.4)", model_default=True, model_default_text="mmc1.xls Table 1 Ex_ttton[e] -1000/0 (TTM)", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 cols 13-14; Methods 2.4")
add("Ex_tsul[e]", "thiosulfate electron donor (TSM)", lower_bound=-1000, upper_bound=0, conditions=["TSM"], reported=True, reported_text="thiosulfate as electron donor under TSM (Methods 2.4)", model_default=True, model_default_text="mmc1.xls Table 1 Ex_tsul[e] -1000/0 (TSM)", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 cols 15-16; Methods 2.4")
add("CYT2", "cytochrome c/rusticyanin complex (Fe2 oxidation)", lower_bound=0, upper_bound=0, conditions=["TTM", "TSM"], reported=True, reported_text="CYT2 reaction bounds set to 0 for TTM/TSM (Methods 2.4)", model_default=True, model_default_text="mmc1.xls Table 1 CYT2 0/0 for ttton and tsul", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 row CYT2; Methods 2.4")
add("Ex_n2[e]", "dinitrogen uptake (nitrogen source)", lower_bound=-1000, upper_bound=1000, conditions=["FIM", "TTM", "TSM"], reported=False, model_default=True, model_default_text="mmc1.xls Table 1 Ex_n2[e] -1000/1000 (all three)", evidence_status="MODEL_DEFAULT", source_file="mmc1.xls", source_section="Table 1 row Ex_n2[e]", notes="Nitrogen source not explicitly stated in paper Methods; model uses N2 (open) and NH4 (closed), consistent with NITF nitrogen fixation.")
add("Ex_nh4[e]", "ammonia uptake", lower_bound=0, upper_bound=0, conditions=["FIM", "TTM", "TSM"], reported=False, model_default=True, model_default_text="mmc1.xls Table 1 Ex_nh4[e] 0/0 (all three)", evidence_status="MODEL_DEFAULT", source_file="mmc1.xls", source_section="Table 1 row Ex_nh4[e]", notes="NH4 closed in all three minimal media; nitrogen comes via N2 fixation.")
add("Ex_so4[e]", "sulfate exchange (export of oxidized S)", lower_bound=0, upper_bound=0, conditions=["FIM"], reported=False, model_default=True, model_default_text="mmc1.xls Table 1 Ex_so4[e] 0/0 (FIM), 0/1000 (TTM, TSM)", evidence_status="MODEL_DEFAULT", source_file="mmc1.xls", source_section="Table 1 row Ex_so4[e]", notes="Sulfate export allowed only under TTM/TSM (S oxidation produces sulfate).")
add("ATPM (NGAM)", "non-growth associated maintenance (ATP hydrolysis)", value=3.475, units="mmol/gDW/h", lower_bound=3.475, upper_bound=3.475, conditions=["FIM", "TTM", "TSM"], reported=True, reported_text="NGAM set as 2.5% of GAM (Feist et al., 2006); GAM=139 -> 3.475", model_default=True, model_default_text="mmc1.xls Table 1 ATPM lb=ub=3.475", evidence_status="REPORTED", source_file="mmc1.xls; 2016 paper", source_section="Table 1 row ATPM; Methods 2.5", notes="NGAM value 3.475 is INFERRED from GAM*2.5%, but fixed in model as ATPM=3.475.")
add("Afe_biomass_mc507_WT_139p0M (GAM)", "biomass objective function with growth-associated maintenance", value=139, units="mmol ATP/gDW", conditions=["FIM", "TTM", "TSM"], reported=True, reported_text="biomass reaction labelled with 139 GAM estimate; GAM fitted by genetic algorithm", model_default=True, model_default_text="mmc3.xml / mmc1.xls biomass reaction contains 139.1 atp and 139 h2o", evidence_status="REPORTED", source_file="mmc3.xml; mmc1.xls; 2016 paper", source_section="Table 1 row Afe_biomass_mc507_WT_139p0M; Methods 2.5; Supp Note 3")
add("Objective", "maximize biomass exchange", value=1.0, units="coefficient", conditions=["FIM", "TTM", "TSM"], reported=True, reported_text="FBA run to maximize flux through BOF (Methods; Fig captions)", model_default=True, model_default_text="mmc3.xml legacy OBJECTIVE_COEFFICIENT=1.0 on R_Ex_bio_LSQBKT_e_RSQBKT_", evidence_status="REPORTED", source_file="mmc3.xml; 2016 paper", source_section="Methods 2.3; mmc3.xml objective parameter")
scenario = dict(schema_version="1.0", title="2016 iMC507 scenario constraints (traceable run-condition definition)", generated_at=now_iso(), note="Defines 2016 model run conditions, not a wet-lab medium recipe. Values not explicitly reported are not guessed.", source_files=dict(sbml=str((ORIG / "mmc3.xml").relative_to(ROOT)), reaction_table=str((ORIG / "mmc1.xls").relative_to(ROOT)), supplementary_word=str((ORIG / "mmc2.docx").relative_to(ROOT)), paper=str(PAPER.relative_to(ROOT))), model_scope=dict(total_reactions_sbml=615, exchange_reactions=28, non_exchange_reactions=587, metabolites=573, genes=507, compartments=["c", "e", "p"], objective_reaction="Ex_bio[e] (R_Ex_bio_LSQBKT_e_RSQBKT_)"), minimal_media_from_table1=minimal_media, experimental_points_table5=points, parameters=parameters)
(OUT / "scenario_definition" / "2016_scenario_constraints.json").write_text(json.dumps(scenario, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
log = ["generated_at=" + now_iso(), "python=" + sys.version.split()[0], "total_reactions=" + str(classification["total"]), "exchange=" + str(classification["exchange"]), "biomass_objective=" + str(classification["biomass_objective"]), "transport=" + str(classification["transport"]), "metabolic=" + str(classification["metabolic"]), "non_exchange=" + str(classification["non_exchange"]), "experimental_points_table5=" + str(len(points))]
(OUT / "logs" / "phase0_freeze.log").write_text("\n".join(log) + "\n", encoding="utf-8")
print("DONE")
print("manifest_match_all=", all(r["match"] for r in manifest))
print(classification)
print("points=", len(points))
