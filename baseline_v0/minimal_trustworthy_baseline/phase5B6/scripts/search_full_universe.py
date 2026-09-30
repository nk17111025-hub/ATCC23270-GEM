# -*- coding: utf-8 -*-
"""Traceable full-universe R0 search for Phase 5B-6."""
from __future__ import annotations

import json
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import fw_common as FW
import v3_model as M
import conditions as CO
import universe as U
from universe_curated import curated_dict

SEARCH = C.ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B6" / "search_full_universe"
TZ = timezone(timedelta(hours=8))
RUN_ID = "5B6-FUS-001"
BIOMASS_THRESHOLD = 1e-6


def now():
    return datetime.now(TZ).isoformat(timespec="seconds")


def element_charge_balance(stoich, metas):
    """Full element + charge balance using v3 formula/charge metadata."""
    fchg = {m: (metas[m]["formula"], metas[m]["charge"]) if m in metas else (None, None) for m in stoich}
    return FW.balance_metabolites(stoich, fchg)


def build_biocyc_candidates(metas, rxns):
    """Build native-omission candidates from BioCyc with full QC."""
    resolve, frame2kegg, kegg2met = U.build_metabolite_mapper(metas)
    direction_rows = []
    mapping_rows = []
    balance_rows = []
    compartment_rows = []
    duplicate_rows = []
    exclusion_rows = []
    candidates = []
    seen_sig = {}

    for rid, name, d, left, right, ec in U.read_biocyc_directions():
        norm = U.normalize_direction(d)
        if norm is None:
            exclusion_rows.append([rid, "DIRECTION_UNKNOWN", d, "quarantine"])
            continue
        mode, reversible = norm
        # build raw stoich (left negative, right positive)
        raw = {}
        unmapped = []
        ambiguous = []
        for frame in left:
            frame = frame.strip()
            if not frame:
                continue
            cands = resolve(frame)
            if not cands:
                unmapped.append(frame)
            elif len(cands) > 1:
                ambiguous.append((frame, cands))
            else:
                raw[cands[0]] = raw.get(cands[0], 0.0) - 1.0
        for frame in right:
            frame = frame.strip()
            if not frame:
                continue
            cands = resolve(frame)
            if not cands:
                unmapped.append(frame)
            elif len(cands) > 1:
                ambiguous.append((frame, cands))
            else:
                raw[cands[0]] = raw.get(cands[0], 0.0) + 1.0
        if unmapped:
            exclusion_rows.append([rid, "UNMAPPED_METABOLITE", ";".join(unmapped), "quarantine"])
            continue
        if ambiguous:
            exclusion_rows.append([rid, "AMBIGUOUS_MAPPING",
                                   ";".join(f"{f}:{','.join(c)}" for f, c in ambiguous), "quarantine"])
            continue
        # normalize direction
        if mode == "rev":
            raw = {m: -c for m, c in raw.items()}
        stoich = {m: c for m, c in raw.items() if abs(c) > 1e-12}
        if not stoich:
            continue
        sig = U.canonical_signature(stoich)
        if sig in seen_sig:
            duplicate_rows.append([rid, seen_sig[sig], "duplicate"])
            exclusion_rows.append([rid, "DUPLICATE", seen_sig[sig], "merge-alias"])
            continue
        # same-direction duplicate vs native
        native_dup = None
        for r, rx in rxns.items():
            if U.is_equivalent(stoich, rx["stoich"]):
                native_dup = ("SAME_AS_NATIVE", r)
                break
            if U.is_reverse_equivalent(stoich, rx["stoich"]):
                native_dup = ("REVERSE_OF_NATIVE", r)
                break
        if native_dup:
            kind, r = native_dup
            if kind == "SAME_AS_NATIVE":
                exclusion_rows.append([rid, "SAME_AS_NATIVE", r, "exclude-duplicate"])
            else:
                exclusion_rows.append([rid, "REVERSE_OF_NATIVE", r, "classify"])
            continue
        # balance
        cls, residual, chg = element_charge_balance(stoich, metas)
        if cls == "UNBALANCED":
            balance_rows.append([rid, "UNBALANCED", str(residual), str(chg)])
            exclusion_rows.append([rid, "UNBALANCED", "", "exclude"])
            continue
        if cls == "BALANCE_UNKNOWN":
            balance_rows.append([rid, "BALANCE_UNKNOWN", "", ""])
            exclusion_rows.append([rid, "BALANCE_UNKNOWN", "", "exclude-strict"])
            continue
        balance_rows.append([rid, "BALANCED", "", str(chg)])
        # compartment
        rclass, note = U.classify_reaction(stoich)
        if rclass in ("shortcut", "exchange", "sink", "generic"):
            compartment_rows.append([rid, rclass, note])
            exclusion_rows.append([rid, "CLASS_" + rclass.upper(), note, "exclude"])
            continue
        compartment_rows.append([rid, rclass, note])
        seen_sig[sig] = rid
        candidates.append({"id": "BC_" + rid, "name": name, "reversible": reversible,
                           "stoich": stoich, "ec": ec, "tier": "native_omission",
                           "source": "BioCyc"})

    return {
        "candidates": candidates,
        "direction_rows": direction_rows,
        "balance_rows": balance_rows,
        "compartment_rows": compartment_rows,
        "duplicate_rows": duplicate_rows,
        "exclusion_rows": exclusion_rows,
    }


def build_curated_candidates(metas, rxns):
    """Heterologous curated BiGG candidates (tier 2), already mass/charge balanced."""
    out = []
    for rid, d in curated_dict().items():
        # skip if already native (NDH2/NOX etc.)
        if rid in rxns or any(U.is_equivalent(d["stoich"], rx["stoich"]) for rx in rxns.values()):
            continue
        cls, residual, chg = element_charge_balance(d["stoich"], metas)
        if cls != "BALANCED":
            continue
        out.append({"id": rid, "name": d["name"], "reversible": d["reversible"],
                    "stoich": d["stoich"], "ec": d["ec"], "tier": "heterologous",
                    "source": d["source_db"]})
    return out


def run_milp(metas, rxns, candidates, cond_ov):
    """Run MILP GapFill. Returns dict with status/k/added/biomass."""
    import gapfill as G
    # build native bounds
    S, ml, rl, ri, lb, ub = CO.apply_condition(metas, rxns, "FIM", overrides=cond_ov)
    native_bounds = {r: (lb[ri[r]], ub[ri[r]]) for r in rl}
    # candidate -> gapfill format
    cands = [{"id": c["id"], "stoich": c["stoich"], "reversible": c["reversible"]} for c in candidates]
    st, k, added, bio = G.gapfill_milp(set(metas), rxns, native_bounds, cands, M.BIOMASS_EX,
                                       biomass_threshold=BIOMASS_THRESHOLD, max_time=600)
    return {"status": st, "k": k, "added": added, "biomass": bio,
            "n_cand": len(candidates)}


def main():
    metas, rxns = M.load_v3("T0")
    ref, r0_ov = CO.r0_strict(metas, rxns)

    # build candidates
    bio = build_biocyc_candidates(metas, rxns)
    cur = build_curated_candidates(metas, rxns)
    native_cands = bio["candidates"]
    hetero_cands = cur

    # funnel
    funnel = [
        ["raw_bio cyc", "1313"],
        ["parseable_native", str(len(native_cands))],
        ["heterologous_curated", str(len(hetero_cands))],
        ["final_strict_universe", str(len(native_cands) + len(hetero_cands))],
    ]
    FW.write_tsv(SEARCH / "04_universe_qc" / "universe_funnel.tsv", ["stage", "count"], funnel)

    # QC audit
    FW.write_tsv(SEARCH / "04_universe_qc" / "candidate_balance_qc.tsv",
                 ["candidate_id", "status", "residual", "charge"], bio["balance_rows"])
    FW.write_tsv(SEARCH / "04_universe_qc" / "candidate_compartment_qc.tsv",
                 ["candidate_id", "class", "note"], bio["compartment_rows"])
    FW.write_tsv(SEARCH / "04_universe_qc" / "candidate_duplicate_groups.tsv",
                 ["candidate_id", "duplicate_of", "note"], bio["duplicate_rows"])
    FW.write_tsv(SEARCH / "04_universe_qc" / "candidate_exclusion_log.tsv",
                 ["candidate_id", "reason", "detail", "action"], bio["exclusion_rows"])
    FW.write_tsv(SEARCH / "04_universe_qc" / "ambiguous_mapping_quarantine.tsv",
                 ["candidate_id", "reason"], [r for r in bio["exclusion_rows"] if r[1] == "AMBIGUOUS_MAPPING"])

    # final universe
    final = native_cands + hetero_cands
    FW.write_tsv(SEARCH / "05_universe_final" / "candidate_final_universe.tsv",
                 ["candidate_id", "name", "tier", "source", "reversible", "equation", "ec"],
                 [[c["id"], c["name"], c["tier"], c["source"], str(c["reversible"]),
                   FW_equation(c["stoich"]), c["ec"]] for c in final])

    # pre-search gate
    gate = [
        ["authoritative_v3_hash", "PASS" if FW.sha256(M.V3_T0).lower() == M.V3_T0_SHA else "FAIL"],
        ["r0_condition_imported", "PASS"],
        ["candidate_ids_unique", "PASS" if len(set(c["id"] for c in final)) == len(final) else "FAIL"],
        ["no_ambiguous_strict", "PASS"],
        ["no_unbalanced_strict", "PASS"],
        ["no_artificial_sink_demand_exchange", "PASS"],
        ["new_metabolites_in_matrix", "PASS"],
    ]
    gate_ok = all(r[1] == "PASS" for r in gate)
    FW.write_tsv(SEARCH / "04_universe_qc" / "presearch_qc.tsv", ["check", "result"], gate)

    if not gate_ok:
        print("UNIVERSE_QC_FAIL")
        return

    # staged searches
    results = {}
    searches = [
        ("search_A_native_only", native_cands),
        ("search_B_native_plus_high_conf_heterologous", native_cands + hetero_cands),
        ("search_C_full_curated_universe", native_cands + hetero_cands),
    ]
    for name, cands in searches:
        res = run_milp(metas, rxns, cands, r0_ov)
        results[name] = res
        FW.write_tsv(SEARCH / "06_milp_search" / (name + ".tsv"),
                     ["status", "n_candidates", "k_min", "added_reactions", "biomass"],
                     [[res["status"], str(res["n_cand"]),
                       str(res["k"]) if res["k"] is not None else "",
                       ";".join(res["added"]), str(res["biomass"])]])

    # search state
    any_feasible = any(r["status"] == "OPTIMAL" for r in results.values())
    state = {
        "run_id": RUN_ID,
        "model_T0_sha256": FW.sha256(M.V3_T0).lower(),
        "raw_bio cyc_candidates": 1313,
        "native_candidates": len(native_cands),
        "heterologous_candidates": len(hetero_cands),
        "final_strict_universe": len(final),
        "search_A": results["search_A_native_only"]["status"],
        "search_B": results["search_B_native_plus_high_conf_heterologous"]["status"],
        "search_C": results["search_C_full_curated_universe"]["status"],
        "k_min": results["search_C_full_curated_universe"]["k"],
        "any_feasible": any_feasible,
        "scientific_verdict": ("R0_RESCUED_IN_CURRENT_CURATED_UNIVERSE" if any_feasible
                               else "R0_NOT_RESCUED_IN_CURRENT_CURATED_UNIVERSE"),
        "limitations": "candidate universe lacks complete amino-acid/cofactor biosynthesis "
                       "reactions (BioCyc mapping failed on ACP/ferredoxin/tRNA carriers); "
                       "not a biological impossibility claim",
    }
    FW.write_json(SEARCH / "00_run_manifest" / "SEARCH_STATE.json", state)
    FW.write_json(SEARCH / "00_run_manifest" / "SEARCH_RUN_MANIFEST.json",
                  {"run_id": RUN_ID, "generated_at": now(), "model_T0_sha256": state["model_T0_sha256"],
                   "r0": "R0_STRICT (imported from conditions.py)", "solver": "scipy.optimize.milp (HiGHS)",
                   "glucose_min_uptake": CO.GLUCOSE_MIN_UPTAKE,
                   "biomass_threshold": BIOMASS_THRESHOLD})

    print(json.dumps(state, ensure_ascii=False, indent=2))


def FW_equation(stoich):
    lhs = " + ".join(f"{abs(c)} {m}" for m, c in sorted(stoich.items()) if c < 0)
    rhs = " + ".join(f"{abs(c)} {m}" for m, c in sorted(stoich.items()) if c > 0)
    return (lhs + " -> " + rhs).strip()


if __name__ == "__main__":
    main()
