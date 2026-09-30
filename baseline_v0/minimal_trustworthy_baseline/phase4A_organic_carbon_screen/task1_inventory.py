# -*- coding: utf-8 -*-
import csv
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

ORG_METS = set()
for s in lib.SUBSTRATE_SPECS:
    if s["c_met"]:
        ORG_METS.add(s["c_met"])
    for a in s["alt_c_met"]:
        ORG_METS.add(a)
    if s["ex_met"]:
        ORG_METS.add(s["ex_met"])
    if s["p_met"]:
        ORG_METS.add(s["p_met"])
for k, v in lib.CENTRAL_TARGETS.items():
    ORG_METS.update(v)
ORG_METS.update([
    "fdp[c]", "3pg[c]", "2pg[c]", "pep[c]", "actACP[c]", "acACP[c]", "malcoa[c]",
    "succoa[c]", "6pgc[c]", "6pgl[c]", "2ddg6p[c]", "glyclt[c]", "2pglyc[c]",
    "gcald[c]", "2h3oppan[c]", "g1p[c]", "g1p-B[c]", "mal-B[c]", "malttr[c]",
    "maltttr[c]", "glycogen[c]", "glycogendeb[c]", "man6p[c]", "gam6p[c]", "gam1p[c]",
])


def role_of(rxn, rx):
    rid = rxn
    name = rx["name"].lower()
    sub = rx["subsystem"].lower()
    stoich = rx["stoich"]
    if rid.startswith("Ex_"):
        return "exchange"
    if "transport" in sub or ("tex" in rid.lower() or "tpp" in rid.lower()):
        return "transport"
    if "kinase" in name or "phosphor" in name or rx["ec"].startswith("2.7."):
        return "phosphorylation"
    if any(x in name for x in ["dehydrogenase", "oxidase", "reductase", "dioxygenase", "epimerase", "isomerase"]):
        if any(k in stoich for k in ["nadh[c]", "nadph[c]", "nad[c]", "nadp[c]", "fadh2[c]", "fad[c]", "o2[c]", "h2o2[c]"]):
            return "oxidation_reduction"
        return "isomerization"
    if "coa" in stoich and "coa[c]" in stoich:
        return "activation_to_CoA"
    if "co2[c]" in stoich or "hco3[c]" in stoich or "h2co3[c]" in stoich:
        return "carboxylation_decarboxylation"
    return "conversion_central"


def main():
    metas, rxns = lib.parse_sbml_meta()
    precursors = set(lib.biomass_precursors(rxns))

    # reachability maps (as-is, reversible-aware)
    edges = lib.graph_edges(metas, rxns, all_reversible=False)
    exchange_mets = set()
    for r, rx in rxns.items():
        if r.startswith("Ex_"):
            for m, c in rx["stoich"].items():
                exchange_mets.add(m)
    # from exchange: BFS forward from exchange metabolites
    from_ex = set(exchange_mets)
    stack = list(exchange_mets)
    while stack:
        cur = stack.pop()
        for (nxt, _r, _rev) in edges.get(cur, []):
            if nxt not in from_ex:
                from_ex.add(nxt)
                stack.append(nxt)
    # to biomass: reverse BFS from precursors
    rev_edges = {m: [] for m in metas}
    for m in metas:
        for (nxt, r, rev) in edges.get(m, []):
            rev_edges.setdefault(nxt, []).append((m, r, rev))
    to_bio = set(precursors)
    stack = list(precursors)
    while stack:
        cur = stack.pop()
        for (prv, _r, _rev) in rev_edges.get(cur, []):
            if prv not in to_bio:
                to_bio.add(prv)
                stack.append(prv)

    rows = []
    seen = set()
    for r, rx in sorted(rxns.items()):
        org = sorted(set(rx["stoich"]) & ORG_METS)
        is_ex = r.startswith("Ex_")
        is_transport = ("transport" in rx["subsystem"].lower()) or ("tex" in r.lower() or "tpp" in r.lower())
        if not org and not is_ex and not is_transport:
            continue
        if r in seen:
            continue
        seen.add(r)
        substrate_or_met = ";".join(org) if org else (";".join(sorted(rx["stoich"])) )
        comps = set()
        for m in org:
            comps.add(metas.get(m, {}).get("compartment", ""))
        role = role_of(r, rx)
        conn_ex = "TRUE" if any(m in from_ex for m in org) else "FALSE"
        conn_bio = "TRUE" if any(m in to_bio for m in org) else "FALSE"
        rows.append({
            "substrate_or_metabolite": substrate_or_met,
            "compartment": ";".join(sorted(comps)),
            "reaction_id": r,
            "reaction_equation": lib.equation_str(rx["stoich"]),
            "reaction_role": role,
            "lower_bound": rx["lb"],
            "upper_bound": rx["ub"],
            "confidence": rx["confidence"],
            "GPR": rx["gpr"],
            "EC": rx["ec"],
            "PMID": rx["pmid"],
            "subsystem": rx["subsystem"],
            "currently_connected_to_exchange": conn_ex,
            "currently_connected_to_biomass_network": conn_bio,
            "notes": "",
        })

    fields = ["substrate_or_metabolite", "compartment", "reaction_id", "reaction_equation",
              "reaction_role", "lower_bound", "upper_bound", "confidence", "GPR", "EC",
              "PMID", "subsystem", "currently_connected_to_exchange",
              "currently_connected_to_biomass_network", "notes"]
    out = OUT / "organic_carbon_entry_inventory.tsv"
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote", out, "rows", len(rows))
    print("exchange_mets:", sorted(exchange_mets))
    print("organic exchange reactions:", [r for r in rxns if r.startswith("Ex_") and any(m in ORG_METS for m in rxns[r]["stoich"])])


if __name__ == "__main__":
    main()
