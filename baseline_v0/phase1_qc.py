#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 structural & numerical QC of frozen 2016 iMC507 baseline."""
from __future__ import annotations
import csv
import datetime as dt
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
import libsbml
import xlrd
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import find_blocked_reactions

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
ORIG = ROOT / "01_原始数据" / "02_2016_iMC507原始模型"
FROZEN_SBML = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ORIG / "mmc1.xls"
OUT = ROOT / "baseline_v0" / "qc"
LOGS = ROOT / "baseline_v0" / "logs"
TOL = 1e-9

def now_iso():
    return dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds")

def enc(s):
    return s.replace("[", "_LSQBKT_").replace("]", "_RSQBKT_").replace("-", "_DASH_")

def read_table1():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 1")
    rows = []
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
        rows.append(dict(reaction_id=str(rid), reaction_name=str(cv(1)), equation=str(cv(2)), confidence=cv(3), ec=str(cv(4)), pmid=cv(5), subsystem=str(cv(6)), gra=str(cv(7)), gpra=str(cv(8)), pra=str(cv(9)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15), source_row=r + 1))
    book.release_resources()
    return rows

def read_table2():
    book = xlrd.open_workbook(str(MMC1), on_demand=True)
    s = book.sheet_by_name("Table 2")
    met = {}
    for r in range(1, s.nrows):
        abbr = s.cell_value(r, 0)
        if abbr == "" or abbr is None:
            continue
        formula = s.cell_value(r, 2)
        charge = s.cell_value(r, 3)
        met[str(abbr)] = dict(formula=(str(formula) if formula != "" else ""), charge=(charge if charge != "" else None))
    book.release_resources()
    return met

def load_model():
    raw = FROZEN_SBML.read_bytes()
    doc = libsbml.readSBMLFromString(raw.decode("utf-8", "replace"))
    model = _sbml_to_model(doc)
    return doc, model

ARROW = re.compile(r"(<=>|<->|=>|->|→|↔)")

def parse_equation(text):
    t = text.strip()
    m = ARROW.search(t)
    if not m:
        return None, {}, {}
    arrow = m.group(1)
    left, right = t[: m.start()], t[m.end():]
    reversible = arrow in ("<=>", "<->", "↔")
    return reversible, parse_side(left), parse_side(right)

TERM = re.compile(r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+(.+?)\s*$")

def parse_side(side):
    out = defaultdict(float)
    for part in side.split("+"):
        part = part.strip()
        if not part:
            continue
        m = TERM.match(part)
        if m:
            coef = float(m.group(1))
            met = m.group(2).strip()
        else:
            coef = 1.0
            met = part.strip()
        out[met] += coef
    return dict(out)

ELEM = re.compile(r"([A-Z][a-z]?)(\d*)")

def parse_formula(formula):
    counts = defaultdict(int)
    for m in ELEM.finditer(formula):
        el = m.group(1)
        n = int(m.group(2)) if m.group(2) else 1
        counts[el] += n
    return dict(counts)

def element_balance(left, right, met_map):
    els = defaultdict(float)
    assessable = True
    for side, sign in ((left, -1.0), (right, +1.0)):
        for met, coef in side.items():
            info = met_map.get(met)
            if info is None or info["formula"] == "":
                assessable = False
                continue
            for el, n in parse_formula(info["formula"]).items():
                els[el] += sign * coef * n
    if not assessable:
        return None, assessable
    imbalance = {el: round(v, 9) for el, v in els.items() if abs(v) > 1e-6}
    return imbalance, assessable

def charge_balance(left, right, met_map):
    total = 0.0
    assessable = True
    for side, sign in ((left, -1.0), (right, +1.0)):
        for met, coef in side.items():
            info = met_map.get(met)
            if info is None or info["charge"] is None:
                assessable = False
                continue
            total += sign * coef * float(info["charge"])
    if not assessable:
        return None, assessable
    return round(total, 9), assessable

def reaction_type(r, subsystem=None):
    rid = r.id
    if rid.startswith("Ex_"):
        return "exchange"
    if rid == "Afe_biomass_mc507_WT_139p0M":
        return "biomass"
    if subsystem is not None and "Transport" in subsystem:
        return "transport"
    if len(r.products) == 0:
        return "demand"
    if len(r.reactants) == 0:
        return "sink"
    return "internal"

def build_t1_map(rows):
    return {r["reaction_id"]: r for r in rows}

def apply_condition_bounds(model, rows, cond, log_entries):
    cobra_ids = {r.id for r in model.reactions}
    applied = 0
    unmapped = []
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in cobra_ids:
            unmapped.append(r["reaction_id"])
            continue
        if not (cid.startswith("Ex_") or cid == "CYT2"):
            continue
        lb = r["lb_" + cond]
        ub = r["ub_" + cond]
        rx = model.reactions.get_by_id(cid)
        rx.lower_bound = float(lb)
        rx.upper_bound = float(ub)
        applied += 1
    log_entries.append("condition={} applied_bounds={} unmapped={}".format(cond, applied, len(unmapped)))
    return applied, unmapped

def set_objective_biomass(model):
    bio = model.reactions.get_by_id("Ex_bio_LSQBKT_e_RSQBKT_")
    model.objective = bio

def max_sv_residual(model, fluxes):
    import cobra
    import numpy as np
    S = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    order = [r.id for r in model.reactions]
    v = [float(fluxes.get(rid, 0.0)) for rid in order]
    resid = S @ v
    return float(np.max(np.abs(resid)))

def solve_record(model, label):
    sol = model.optimize()
    fluxes = sol.fluxes
    resid = max_sv_residual(model, fluxes)
    return dict(label=label, status=sol.status, objective=sol.objective_value, max_sv_residual=resid)

def classify_blocked(rid, t1map):
    info = t1map.get(rid)
    if info is None:
        return "unknown"
    if info["reaction_id"].startswith("Ex_"):
        return "expected condition-specific block"
    sub = info["subsystem"].lower()
    if "transport" in sub or "exchange" in sub:
        return "expected condition-specific block"
    if any(k in rid.lower() for k in ("cyt2", "suldo", "sccr", "tsqoc", "h2s")):
        return "potentially biologically meaningful"
    return "potentially suspicious"

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    log = []
    issues = []
    doc, model = load_model()
    rows = read_table1()
    met_map = read_table2()
    t1map = build_t1_map(rows)

    log.append("libsbml_errors={} level={} version={} model_id={!r}".format(doc.getNumErrors(), doc.getLevel(), doc.getVersion(), doc.getModel().getId()))
    for i in range(doc.getNumErrors()):
        e = doc.getError(i)
        log.append("libsbml_error_{}: {} {} {}".format(i, e.getSeverityAsString(), e.getLine(), e.getMessage()))

    n_rxn = len(model.reactions)
    n_met = len(model.metabolites)
    n_gene = len(model.genes)
    comps = list(model.compartments)
    log.append("cobra reactions={} metabolites={} genes={} compartments={}".format(n_rxn, n_met, n_gene, comps))

    struct_rows = []
    type_counts = Counter()
    exchange_ids = []
    fixed_ids = []
    empty_stoich = []
    lb_gt_ub = []
    for r in model.reactions:
        t1 = t1map.get(r.id)
        typ = reaction_type(r, subsystem=(t1["subsystem"] if t1 else ""))
        type_counts[typ] += 1
        if typ == "exchange":
            exchange_ids.append(r.id)
        if r.lower_bound == r.upper_bound:
            fixed_ids.append(r.id)
        if len(r.metabolites) == 0:
            empty_stoich.append(r.id)
        if r.lower_bound > r.upper_bound:
            lb_gt_ub.append(r.id)
        struct_rows.append(dict(reaction_id=r.id, table1_id=(t1["reaction_id"] if t1 else ""), reaction_name=r.name, type=typ, subsystem=(t1["subsystem"] if t1 else ""), reversible=r.reversibility, lower_bound=r.lower_bound, upper_bound=r.upper_bound, fixed=(r.lower_bound == r.upper_bound), empty_stoichiometry=(len(r.metabolites) == 0), lb_gt_ub=(r.lower_bound > r.upper_bound), is_duplicate=False))

    rids = [r.getId() for r in doc.getModel().getListOfReactions()]
    sids = [sp.getId() for sp in doc.getModel().getListOfSpecies()]
    dup_rxn = [k for k, v in Counter(rids).items() if v > 1]
    dup_met = [k for k, v in Counter(sids).items() if v > 1]
    log.append("duplicate_reaction_ids={} duplicate_metabolite_ids={}".format(dup_rxn, dup_met))

    expected = dict(total=615, metabolites=573, compartments=3, exchange=28, non_exchange=587)
    actual = dict(total=n_rxn, metabolites=n_met, compartments=len(comps), exchange=len(exchange_ids), non_exchange=n_rxn - len(exchange_ids))
    for k in expected:
        if expected[k] != actual[k]:
            issues.append(dict(check="structural_count_" + k, expected=expected[k], actual=actual[k], detail="count mismatch vs Phase 0"))

    log.append("objective_expression={}".format(model.objective.expression))
    log.append("objective_direction={}".format(model.objective.direction))

    with (OUT / "structural_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["reaction_id", "table1_id", "reaction_name", "type", "subsystem", "reversible", "lower_bound", "upper_bound", "fixed", "empty_stoichiometry", "lb_gt_ub", "is_duplicate"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in struct_rows:
            w.writerow(r)

    bal_rows = []
    bal_counts = Counter()
    for r in rows:
        reversible, left, right = parse_equation(r["equation"])
        if reversible is None:
            overall = "NOT_ASSESSABLE"
            mass_status = "NOT_ASSESSABLE"
            charge_status = "NOT_ASSESSABLE"
            elem_imb = ""
            charge_imb = ""
        elif r["reaction_id"].startswith("Ex_"):
            mass_status = "NOT_ASSESSABLE"
            charge_status = "NOT_ASSESSABLE"
            overall = "NOT_ASSESSABLE"
            elem_imb = ""
            charge_imb = ""
        else:
            elem_imb, mass_assessable = element_balance(left, right, met_map)
            charge_imb, charge_assessable = charge_balance(left, right, met_map)
            mass_status = ("BALANCED" if (mass_assessable and not elem_imb) else ("MASS_IMBALANCED" if mass_assessable else "NOT_ASSESSABLE"))
            charge_status = ("BALANCED" if (charge_assessable and charge_imb == 0) else ("CHARGE_IMBALANCED" if charge_assessable else "NOT_ASSESSABLE"))
            if mass_status == "NOT_ASSESSABLE" or charge_status == "NOT_ASSESSABLE":
                overall = "NOT_ASSESSABLE"
            elif mass_status == "MASS_IMBALANCED" and charge_status == "CHARGE_IMBALANCED":
                overall = "MASS_AND_CHARGE_IMBALANCED"
            elif mass_status == "MASS_IMBALANCED":
                overall = "MASS_IMBALANCED"
            elif charge_status == "CHARGE_IMBALANCED":
                overall = "CHARGE_IMBALANCED"
            else:
                overall = "BALANCED"
        bal_counts[overall] += 1
        typ = "exchange" if r["reaction_id"].startswith("Ex_") else ("biomass" if "biomass" in r["reaction_id"].lower() else ("transport" if "Transport" in r["subsystem"] else "internal"))
        bal_rows.append(dict(reaction_id=r["reaction_id"], reaction_name=r["reaction_name"], type=typ, subsystem=r["subsystem"], equation=r["equation"], mass_status=mass_status, charge_status=charge_status, overall_status=overall, element_imbalance=(json.dumps(elem_imb, ensure_ascii=False) if isinstance(elem_imb, dict) else ""), charge_imbalance=("" if charge_imb is None else str(charge_imb)), note=("boundary exchange" if r["reaction_id"].startswith("Ex_") else "")))

    with (OUT / "reaction_balance_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["reaction_id", "reaction_name", "type", "subsystem", "equation", "mass_status", "charge_status", "overall_status", "element_imbalance", "charge_imbalance", "note"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in bal_rows:
            w.writerow(r)

    gpr_rows = []
    gene_counter = Counter()
    non_empty = 0
    empty = 0
    malformed = 0
    for r in rows:
        gpr = r["gra"]
        genes = re.findall(r"AFE_\d+", gpr)
        for g in genes:
            gene_counter[g] += 1
        mal = False
        if gpr.strip():
            non_empty += 1
            if gpr.count("(") != gpr.count(")"):
                mal = True
            if re.search(r"\b(and|or)\s+(and|or)\b", gpr, re.I):
                mal = True
        else:
            empty += 1
        if mal:
            malformed += 1
        gpr_rows.append(dict(reaction_id=r["reaction_id"], gpr=gpr, gene_count=len(genes), unique_genes=",".join(sorted(set(genes))), malformed_boolean=mal, note=("empty" if not gpr.strip() else "")))

    with (OUT / "gpr_reference_qc.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["reaction_id", "gpr", "gene_count", "unique_genes", "malformed_boolean", "note"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in gpr_rows:
            w.writerow(r)

    dead_rows = []
    for met in model.metabolites:
        can_produce = False
        can_consume = False
        producers = []
        consumers = []
        for rxn in met.reactions:
            coeff = rxn.get_coefficient(met)
            if coeff > 0 and rxn.upper_bound > 0:
                can_produce = True
                producers.append(rxn.id)
            if coeff < 0 and rxn.upper_bound > 0:
                can_consume = True
                consumers.append(rxn.id)
            if coeff > 0 and rxn.lower_bound < 0:
                can_consume = True
                consumers.append(rxn.id)
            if coeff < 0 and rxn.lower_bound < 0:
                can_produce = True
                producers.append(rxn.id)
        role = ""
        if not can_produce:
            role = "no_producer"
        elif not can_consume:
            role = "no_consumer"
        if role:
            dead_rows.append(dict(metabolite_id=met.id, compartment=met.compartment, role=role, producers=",".join(sorted(set(producers))), consumers=",".join(sorted(set(consumers))), is_boundary=any(r.id.startswith("Ex_") for r in met.reactions)))

    with (OUT / "dead_end_metabolites.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["metabolite_id", "compartment", "role", "producers", "consumers", "is_boundary"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in dead_rows:
            w.writerow(r)

    conds = ["fe2", "ttton", "tsul"]
    blocked = {}
    for cond in conds:
        cp = model.copy()
        apply_condition_bounds(cp, rows, cond, log)
        set_objective_biomass(cp)
        blocked_ids = find_blocked_reactions(cp, zero_cutoff=1e-7)
        blocked[cond] = blocked_ids
        log.append("blocked_{}={}".format(cond, len(blocked_ids)))
        rows_out = []
        for rid in sorted(blocked_ids):
            rx = cp.reactions.get_by_id(rid)
            info = t1map.get(rid)
            rows_out.append(dict(reaction_id=rid, table1_id=(info["reaction_id"] if info else ""), name=rx.name, type=reaction_type(rx, subsystem=(info["subsystem"] if info else "")), subsystem=(info["subsystem"] if info else ""), classification=classify_blocked(rid, t1map)))
        cond_label = "FIM" if cond == "fe2" else ("TTM" if cond == "ttton" else "TSM")
        with (OUT / "blocked_reactions_{}.tsv".format(cond_label)).open("w", encoding="utf-8", newline="") as f:
            fields = ["reaction_id", "table1_id", "name", "type", "subsystem", "classification"]
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows_out:
                w.writerow(r)

    num_rows = []
    base = model.copy()
    set_objective_biomass(base)
    rec = solve_record(base, "baseline_untouched")
    num_rows.append(dict(condition="baseline_untouched", **rec))
    for cond in conds:
        cp = model.copy()
        apply_condition_bounds(cp, rows, cond, log)
        set_objective_biomass(cp)
        cond_label = "FIM" if cond == "fe2" else ("TTM" if cond == "ttton" else "TSM")
        rec = solve_record(cp, cond_label)
        num_rows.append(dict(condition=cond_label, **rec))

    exact_rows = []
    for label, cp0 in [("baseline_untouched", base), ("FIM", model.copy())]:
        if label == "FIM":
            apply_condition_bounds(cp0, rows, "fe2", log)
        set_objective_biomass(cp0)
        try:
            cp0.solver = "glpk_exact"
            sol = cp0.optimize()
            exact_rows.append(dict(condition=label, status=sol.status, objective=sol.objective_value))
        except Exception as exc:
            exact_rows.append(dict(condition=label, status="ERROR", objective=str(exc)))
        cp0.solver = "glpk"

    with (OUT / "numerical_feasibility.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["condition", "label", "status", "objective", "max_sv_residual"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in num_rows:
            w.writerow(r)

    cycle_rows = []
    c_cycle = model.copy()
    for r in c_cycle.reactions:
        if r.id.startswith("Ex_"):
            r.lower_bound = 0
            r.upper_bound = 0
    from cobra.flux_analysis import flux_variability_analysis
    fva = flux_variability_analysis(c_cycle, fraction_of_optimum=1.0)
    internal_nonzero = []
    for rid, row in fva.iterrows():
        if abs(row["minimum"]) > 1e-7 or abs(row["maximum"]) > 1e-7:
            internal_nonzero.append(rid)
    cycle_rows.append(dict(test_id="internal_cycles_all_exchanges_closed", description="FVA with all Ex_ bounds=0; internal reactions able to carry nonzero flux", status=("POSITIVE" if internal_nonzero else "NEGATIVE"), result=";".join(sorted(internal_nonzero)), notes=""))

    c_atp = model.copy()
    for r in c_atp.reactions:
        if r.id.startswith("Ex_"):
            r.lower_bound = 0
            r.upper_bound = 0
    atpm = c_atp.reactions.get_by_id("ATPM")
    atpm.lower_bound = 0
    atpm.upper_bound = 1000
    c_atp.objective = atpm
    sol_atp = c_atp.optimize()
    cycle_rows.append(dict(test_id="free_atp_generation", description="All Ex_=0, ATPM relaxed to [0,1000], maximize ATPM", status=("POSITIVE" if sol_atp.objective_value > 1e-7 else "NEGATIVE"), result=str(sol_atp.objective_value), notes="ATPM bounds changed to [0,1000] in diagnostic copy only"))

    for name, metid in [("free_nadh_generation", "nadh_c"), ("free_nadph_generation", "nadph_c")]:
        c = model.copy()
        for r in c.reactions:
            if r.id.startswith("Ex_"):
                r.lower_bound = 0
                r.upper_bound = 0
        drain = cobra.Reaction("TMP_drain_" + name)
        c.add_reactions([drain])
        drain.add_metabolites({metid: -1.0})
        drain.lower_bound = 0
        drain.upper_bound = 1000
        c.objective = drain
        sol = c.optimize()
        cycle_rows.append(dict(test_id=name, description="All Ex_=0, temporary demand {} -> [] added, maximize".format(metid), status=("POSITIVE" if sol.objective_value > 1e-7 else "NEGATIVE"), result=str(sol.objective_value), notes="temporary demand reaction added to diagnostic copy only"))

    with (OUT / "cycle_diagnostics.tsv").open("w", encoding="utf-8", newline="") as f:
        fields = ["test_id", "description", "status", "result", "notes"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in cycle_rows:
            w.writerow(r)

    probe_lines = []
    probe_lines.append("solver_probe " + now_iso())
    probe_lines.append("cobra_version=" + cobra.__version__)
    cfg = cobra.Configuration()
    probe_lines.append("cobra_tolerance=" + str(cfg.tolerance))
    probe_lines.append("cobra_solver=" + str(cfg.solver))
    import swiglpk
    probe_lines.append("glpk_version=" + swiglpk.glp_version())
    from cobra.util.solver import solvers as solver_list
    probe_lines.append("solvers_available=" + str(sorted(solver_list.keys())))
    for name, ctor in [("glp_smcp", swiglpk.glp_smcp), ("glp_iocp", swiglpk.glp_iocp), ("glp_iptcp", swiglpk.glp_iptcp)]:
        try:
            p = ctor()
            d = {}
            for attr in ("tol_bnd", "tol_dj", "tol_piv", "mip_gap", "tol_obj"):
                try:
                    d[attr] = float(getattr(p, attr))
                except Exception:
                    pass
            probe_lines.append(name + "=" + str(d))
        except Exception as exc:
            probe_lines.append(name + "=ERROR " + repr(exc))
    probe_lines.append("glpk_exact_available=" + str("glpk_exact" in solver_list))
    with (OUT / "solver_probe.txt").open("w", encoding="utf-8") as f:
        f.write("\n".join(probe_lines) + "\n")

    summary = dict(generated_at=now_iso(), structural=dict(reactions=n_rxn, metabolites=n_met, genes=n_gene, compartments=comps, type_counts=dict(type_counts), exchange_ids=sorted(exchange_ids), fixed_ids=fixed_ids, empty_stoichiometry=empty_stoich, lb_gt_ub=lb_gt_ub, duplicate_reactions=dup_rxn, duplicate_metabolites=dup_met), balance_counts=dict(bal_counts), gpr=dict(non_empty=non_empty, empty=empty, malformed=malformed, unique_genes=len(gene_counter)), dead_end_metabolites=len(dead_rows), blocked={cond: len(blocked[cond]) for cond in conds}, numerical=num_rows, exact=exact_rows, cycles=cycle_rows, issues=issues)
    (OUT / "phase1_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    log.append("type_counts=" + str(dict(type_counts)))
    log.append("balance_counts=" + str(dict(bal_counts)))
    log.append("gpr non_empty={} empty={} malformed={} unique_genes={}".format(non_empty, empty, malformed, len(gene_counter)))
    log.append("dead_end_metabolites=" + str(len(dead_rows)))
    log.append("numerical=" + str(num_rows))
    log.append("exact=" + str(exact_rows))
    log.append("cycles_status=" + str([(c["test_id"], c["status"]) for c in cycle_rows]))
    (LOGS / "phase1_qc.log").write_text("\n".join(log) + "\n", encoding="utf-8")

    print("PHASE1 QC done")
    print("type_counts", dict(type_counts))
    print("balance_counts", dict(bal_counts))
    print("gpr", non_empty, empty, malformed, len(gene_counter))
    print("dead_end", len(dead_rows))
    print("blocked", {cond: len(blocked[cond]) for cond in conds})
    print("numerical", num_rows)
    print("exact", exact_rows)
    print("cycles", [(c["test_id"], c["status"]) for c in cycle_rows])

if __name__ == "__main__":
    main()
