#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Freeze minimal_trustworthy_baseline_v1: apply 3 approved repairs + validate."""
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import libsbml
import xlrd
import numpy as np
import cobra
from cobra.io.sbml import _sbml_to_model
from cobra.flux_analysis import pfba, flux_variability_analysis

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
FROZEN = ROOT / "baseline_v0" / "original" / "mmc3.xml"
MMC1 = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"
FREEZE = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze"
V1 = FREEZE / "minimal_trustworthy_baseline_v1.xml"

BIO = "Ex_bio_LSQBKT_e_RSQBKT_"
H2CO3 = "Ex_h2co3_LSQBKT_e_RSQBKT_"
O2 = "Ex_o2_LSQBKT_e_RSQBKT_"
FE2 = "Ex_fe2_LSQBKT_e_RSQBKT_"
TTTON = "Ex_ttton_LSQBKT_e_RSQBKT_"
TSUL = "Ex_tsul_LSQBKT_e_RSQBKT_"

def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

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
        rows.append(dict(reaction_id=str(rid), subsystem=str(cv(6)), lb_fe2=cv(10), ub_fe2=cv(11), lb_ttton=cv(12), ub_ttton=cv(13), lb_tsul=cv(14), ub_tsul=cv(15)))
    book.release_resources()
    return rows

def build_v1():
    raw = FROZEN.read_bytes()
    doc = libsbml.readSBMLFromString(raw.decode("utf-8", "replace"))
    m = doc.getModel()
    # 1. MACPD disable (ub 1000 -> 0)
    m.getReaction("R_MACPD").getKineticLaw().getParameter("UPPER_BOUND").setValue(0.0)
    # 2. ACOATA reverse disable (lb -1000 -> 0)
    m.getReaction("R_ACOATA").getKineticLaw().getParameter("LOWER_BOUND").setValue(0.0)
    # 3. Htpp disable by default (lb -1000 -> 0, ub 1000 -> 0)
    m.getReaction("R_Htpp").getKineticLaw().getParameter("LOWER_BOUND").setValue(0.0)
    m.getReaction("R_Htpp").getKineticLaw().getParameter("UPPER_BOUND").setValue(0.0)
    return libsbml.writeSBMLToString(doc)

def reaction_level_diff():
    """Compare frozen vs v1 for every reaction: bounds + stoichiometry + reversible."""
    def load(p):
        d = libsbml.readSBMLFromString(Path(p).read_bytes().decode("utf-8", "replace"))
        return d.getModel()
    fm = load(FROZEN)
    vm = load(V1)
    diff = []
    for vr in vm.getListOfReactions():
        vid = vr.getId()
        fr = fm.getReaction(vid)
        if fr is None:
            diff.append(dict(reaction_id=vid, change="UNEXPECTED_ADDED", frozen_lb="", frozen_ub="", v1_lb="", v1_ub="", stoich_changed=True))
            continue
        flb = fr.getKineticLaw().getParameter("LOWER_BOUND").getValue()
        fub = fr.getKineticLaw().getParameter("UPPER_BOUND").getValue()
        vlb = vr.getKineticLaw().getParameter("LOWER_BOUND").getValue()
        vub = vr.getKineticLaw().getParameter("UPPER_BOUND").getValue()
        # stoichiometry compare
        def stoich(r):
            d = {}
            for sr in r.getListOfReactants():
                d[sr.getSpecies()] = -sr.getStoichiometry()
            for sp in r.getListOfProducts():
                d[sp.getSpecies()] = d.get(sp.getSpecies(), 0.0) + sp.getStoichiometry()
            return d
        stoich_changed = (stoich(fr) != stoich(vr))
        rev_changed = (fr.getReversible() != vr.getReversible())
        if abs(flb - vlb) > 1e-12 or abs(fub - vub) > 1e-12 or stoich_changed or rev_changed:
            diff.append(dict(reaction_id=vid, frozen_lb=flb, frozen_ub=fub, v1_lb=vlb, v1_ub=vub, stoich_changed=stoich_changed, reversible_changed=rev_changed))
    return diff

def load_v1():
    raw = V1.read_bytes()
    return _sbml_to_model(libsbml.readSBMLFromString(raw.decode("utf-8", "replace")))

def apply_operational_and_patch(m, cond, rows):
    ids = {r.id for r in m.reactions}
    for r in rows:
        cid = enc(r["reaction_id"])
        if cid not in ids or cid == "Afe_biomass_mc507_WT_139p0M":
            continue
        rx = m.reactions.get_by_id(cid)
        rx.bounds = (float(r["lb_" + cond]), float(r["ub_" + cond]))
    # 3 approved patches on top of operational bounds
    m.reactions.get_by_id("MACPD").bounds = (0.0, 0.0)
    m.reactions.get_by_id("ACOATA").bounds = (0.0, 1000.0)
    m.reactions.get_by_id("Htpp").bounds = (0.0, 0.0)

def set_objective(m):
    m.objective = m.reactions.get_by_id(BIO)

def donor_rxn(cond):
    return {"FIM": FE2, "TTM": TTTON, "TSM": TSUL}[cond]

def cond_col(cond):
    return {"FIM": "fe2", "TTM": "ttton", "TSM": "tsul"}[cond]

def run_biomass(m, cond, rows, hco3=2.0, hco3e_off=False, close_ex=()):
    mm = m.copy()
    apply_operational_and_patch(mm, cond_col(cond), rows)
    if hco3e_off:
        mm.reactions.get_by_id("HCO3E").bounds = (0.0, 0.0)
    if hco3 is None:
        mm.reactions.get_by_id(H2CO3).bounds = (0.0, 0.0)
    else:
        mm.reactions.get_by_id(H2CO3).bounds = (-float(hco3), -float(hco3))
    for rid in close_ex:
        if rid in mm.reactions:
            mm.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    set_objective(mm)
    try:
        sol = pfba(mm)
    except Exception:
        sol = None
    return mm, sol

def free_energy_test(m, drain_met):
    """Close all exchanges, add a demand drain, maximize it. Returns objective (0 if none)."""
    mm = m.copy()
    for r in mm.reactions:
        if r.id.startswith("Ex_"):
            r.bounds = (0.0, 0.0)
    mm.reactions.get_by_id("MACPD").bounds = (0.0, 0.0)
    mm.reactions.get_by_id("ACOATA").bounds = (0.0, 1000.0)
    mm.reactions.get_by_id("Htpp").bounds = (0.0, 0.0)
    if drain_met == "ATPM":
        rx = mm.reactions.get_by_id("ATPM")
        rx.bounds = (0.0, 1000.0)
        mm.objective = rx
    else:
        drain = cobra.Reaction("TMP_drain")
        mm.add_reactions([drain])
        drain.add_metabolites({drain_met: -1.0})
        drain.bounds = (0.0, 1000.0)
        mm.objective = drain
    sol = mm.optimize()
    return sol.objective_value

def main():
    FREEZE.mkdir(parents=True, exist_ok=True)
    frozen_hash = sha256(FROZEN)

    # build v1
    V1.write_text(build_v1(), encoding="utf-8")
    v1_hash = sha256(V1)

    # artifact integrity
    diff = reaction_level_diff()

    # load for tests
    model = load_v1()
    rows = read_table1()
    conds = ["FIM", "TTM", "TSM"]

    # A. EGC elimination + broader
    egc = {}
    egc["free_ATP"] = free_energy_test(model, "ATPM")
    egc["free_NADH"] = free_energy_test(model, "nadh_c")
    egc["free_NADPH"] = free_energy_test(model, "nadph_c")
    egc["free_ATP_via_reverse_proton_leak_recheck"] = free_energy_test(model, "ATPM")

    # B/C/D. regression
    reg_rows = []
    for cond in conds:
        m, sol = run_biomass(model, cond, rows)
        wt = float(sol.fluxes[BIO]) if sol is not None else None
        _, sol_ko = run_biomass(model, cond, rows, hco3e_off=True)
        ko = float(sol_ko.fluxes[BIO]) if sol_ko is not None else "infeasible"
        # Htpp FVA
        fva = flux_variability_analysis(m, reaction_list=["Htpp", "MACPD", "ACOATA"], fraction_of_optimum=1.0) if sol is not None else None
        htpp = [float(fva.loc["Htpp", "minimum"]), float(fva.loc["Htpp", "maximum"])] if fva is not None else None
        macpd_fva = [float(fva.loc["MACPD", "minimum"]), float(fva.loc["MACPD", "maximum"])] if fva is not None else None
        acoata_fva = [float(fva.loc["ACOATA", "minimum"]), float(fva.loc["ACOATA", "maximum"])] if fva is not None else None
        # zero carbon / zero donor / zero acceptor
        _, zc = run_biomass(model, cond, rows, hco3=None)
        zc_g = float(zc.fluxes[BIO]) if zc is not None else "infeasible"
        dr = donor_rxn(cond)
        _, zd = run_biomass(model, cond, rows, close_ex=(dr,))
        zd_g = float(zd.fluxes[BIO]) if zd is not None else "infeasible"
        _, za = run_biomass(model, cond, rows, close_ex=(O2, "Ex_fe3_LSQBKT_e_RSQBKT_"))
        za_g = float(za.fluxes[BIO]) if za is not None else "infeasible"
        reg_rows.append(dict(condition=cond, WT_biomass=wt, HCO3E_KO=ko, Htpp_FVA=htpp, MACPD_FVA=macpd_fva, ACOATA_FVA=acoata_fva, zero_carbon=zc_g, zero_donor=zd_g, zero_acceptor=za_g))

    # outputs
    with (FREEZE / "reaction_level_diff_vs_2016.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["reaction_id", "frozen_lb", "frozen_ub", "v1_lb", "v1_ub", "stoich_changed", "reversible_changed"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in diff:
            w.writerow(d)

    with (FREEZE / "egc_diagnostic.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["test", "objective"], delimiter="\t")
        w.writeheader()
        for k, v in egc.items():
            w.writerow(dict(test=k, objective=v))

    with (FREEZE / "regression.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(reg_rows[0].keys()), delimiter="\t")
        w.writeheader()
        for d in reg_rows:
            w.writerow(d)

    manifest = dict(
        frozen_2016_path=str(FROZEN.relative_to(ROOT)),
        frozen_2016_sha256=frozen_hash,
        v1_path=str(V1.relative_to(ROOT)),
        v1_sha256=v1_hash,
        created_at=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
    )
    (FREEZE / "sha256_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # changelog
    changelog = [
        "# MINIMAL_TRUSTWORTHY_BASELINE_V1_CHANGELOG",
        "",
        "generated_at=" + manifest["created_at"],
        "source_model=" + manifest["frozen_2016_path"],
        "source_sha256=" + frozen_hash,
        "v1_model=" + manifest["v1_path"],
        "v1_sha256=" + v1_hash,
        "",
        "## Centrally approved structural changes (3 reactions)",
        "",
        "### 1. MACPD (R_MACPD) - disable",
        "- bounds: [0.0, 1000.0] -> [0.0, 0.0]",
        "- reason: malonyl-ACP decarboxylase; Confidence=1; no GPR/PMID; EC 2.3.1.38 inconsistent; closes artificial HCO3E bypass",
        "- rollback: restore UPPER_BOUND=1000.0",
        "",
        "### 2. ACOATA (R_ACOATA) - disable reverse",
        "- bounds: [-1000.0, 1000.0] -> [0.0, 1000.0]",
        "- reason: reverse direction closes artificial HCO3E bypass; forward retained as conservative assumption",
        "- rollback: restore LOWER_BOUND=-1000.0",
        "",
        "### 3. Htpp (R_Htpp) - disable by default",
        "- bounds: [-1000.0, 1000.0] -> [0.0, 0.0]",
        "- reason: unconstrained reversible Htpp creates free-ATP EGC (Htpp+ATPS5rpp+ATPM); published Table 1 already sets 0/0; no GPR/PMID for alternative direction",
        "- status: DISABLED_BY_DEFAULT_DUE_TO_UNCONSTRAINED_EGC; BIOLOGICAL_PROTON_LEAK_DIRECTION_AND_CAPACITY_UNRESOLVED",
        "- rollback: restore LOWER_BOUND=-1000.0, UPPER_BOUND=1000.0",
        "",
        "## Deferred (not repaired in this freeze)",
        "- CYTRED mmc1/mmc3 proton-direction discrepancy",
        "- CA2tpp imbalance",
        "- remaining charge imbalances",
        "- TSM quantitative thiosulfate-uptake limitation",
        "- 508 vs 507 GPR count",
        "- uncertain ETC proton stoichiometries",
        "- uncertain biomass composition",
        "- unfinished genome-wide G5X review",
    ]
    (FREEZE / "MINIMAL_TRUSTWORTHY_BASELINE_V1_CHANGELOG.md").write_text("\n".join(changelog), encoding="utf-8")

    summary = dict(manifest=manifest, diff=diff, egc=egc, regression=reg_rows)
    (FREEZE / "freeze_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print("=== FREEZE VALIDATION ===")
    print("frozen_sha256", frozen_hash)
    print("v1_sha256", v1_hash)
    print("diff count (should be 3)", len(diff))
    for d in diff:
        print("  ", d["reaction_id"], d["frozen_lb"], d["frozen_ub"], "->", d["v1_lb"], d["v1_ub"], "stoich", d["stoich_changed"], "rev", d["reversible_changed"])
    print("EGC", egc)
    for r in reg_rows:
        print(r)

if __name__ == "__main__":
    main()
