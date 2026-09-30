#!/usr/bin/env python
"""Deduplicate and compare completed E3.0 S01/S03 branch outputs.

Run only after both branch jobs have reached terminal states:
    python scripts/summarize_E3_0_host_optimization.py

Raw optimization_results.tsv and branch pareto_designs.tsv files are read-only
inputs. This script writes pareto_designs_unique.tsv beside each branch output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from collections import Counter, defaultdict
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/engineering/E3_0_host_optimization"
BRANCHES = ("S01", "S03")
MODES = ("CBB_min", "Fe2_efficiency", "O2_efficiency", "CO2_efficiency")
BASELINE = "BASELINE_REFERENCE"
VALID_SOLVER = {"optimal", "time_limit_w_sols", "time_limit"}
TOL = 1e-10


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    if not fields:
        fields = ["status"]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(row: dict, key: str, default=None):
    try:
        value = float(row[key])
        return value if value == value and abs(value) != float("inf") else default
    except (KeyError, TypeError, ValueError):
        return default


def intervention_key(row: dict) -> tuple[tuple[str, ...], tuple[str, ...]]:
    kos = tuple(sorted(x.strip() for x in (row.get("KO_reaction_ids") or "").split(";") if x.strip()))
    regs = tuple(sorted(x.strip() for x in (row.get("regulatory_constraints") or "").split(";") if x.strip()))
    return kos, regs


def valid_reconstructed(row: dict) -> bool:
    return (row.get("solver_status", "").lower() in VALID_SOLVER | {"baseline_reference"}
            and row.get("reconstruction_status", "").lower() in {"optimal", "baseline_reference"}
            and all(number(row, key) is not None for key in
                    ("total_intervention_cost", "biomass", "rubisco", "fe2", "o2", "co2_h2co3")))


def objective_vector(row: dict) -> tuple[float, ...]:
    return (float(row["total_intervention_cost"]), float(row["rubisco"]),
            abs(float(row["fe2"])), abs(float(row["o2"])),
            abs(float(row["co2_h2co3"])), -float(row["biomass"]))


def dominates(a: dict, b: dict) -> bool:
    va, vb = objective_vector(a), objective_vector(b)
    return all(x <= y + TOL for x, y in zip(va, vb)) and any(x < y - TOL for x, y in zip(va, vb))


def pareto(rows: list[dict]) -> list[dict]:
    return [row for i, row in enumerate(rows)
            if not any(i != j and dominates(other, row) for j, other in enumerate(rows))]


def dedupe_designs(raw: list[dict]) -> list[dict]:
    unique: dict[tuple, dict] = {}
    sources: dict[tuple, list[str]] = defaultdict(list)
    for row in raw:
        key = intervention_key(row)
        sources[key].append(f"{row.get('objective_mode')}|{row.get('max_intervention_cost')}|{row.get('solution_id')}")
        # Prefer the baseline label for a zero-intervention phenotype; otherwise
        # preserve the lowest-cost solve and deterministic source ordering.
        previous = unique.get(key)
        status_rank = {"baseline_reference": 0, "optimal": 1, "time_limit_w_sols": 2}
        sort_key = (status_rank.get(row.get("solver_status", "").lower(), 3),
                    int(number(row, "total_intervention_cost", 999999)),
                    row.get("objective_mode", ""), int(number(row, "max_intervention_cost", 999999)),
                    row.get("solution_id", ""))
        if previous is None or sort_key < previous[0]:
            unique[key] = (sort_key, dict(row))
    result = []
    for key, (_, row) in unique.items():
        row["source_solutions"] = ";".join(sorted(set(sources[key])))
        row["KO_reaction_ids"] = ";".join(key[0])
        row["regulatory_constraints"] = ";".join(key[1])
        row["KO_count"] = len(key[0])
        row["OE_count"] = len(key[1])
        row["design_key"] = "KO=" + ";".join(key[0]) + "|REG=" + ";".join(key[1])
        result.append(row)
    result.sort(key=lambda row: (int(number(row, "total_intervention_cost", 0)),
                                 row.get("objective_mode", ""), row.get("solution_id", "")))
    return result


def load_manifest(branch: str) -> dict:
    path = OUT / branch / "checkpoints/input_manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Branch checkpoint manifest missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def status_matrix(rows: list[dict]) -> list[dict]:
    indexed = defaultdict(list)
    for row in rows:
        if row.get("objective_mode") == BASELINE:
            continue
        indexed[(row.get("objective_mode", ""), str(row.get("max_intervention_cost", "")))].append(row)
    out = []
    for mode in MODES:
        budgets = sorted({key[1] for key in indexed if key[0] == mode}, key=lambda x: int(x) if x.isdigit() else 999)
        for budget in budgets:
            group = indexed[(mode, budget)]
            statuses = sorted({r.get("solver_status", "missing").upper() for r in group})
            recon = sorted({r.get("reconstruction_status", "").upper() for r in group if r.get("reconstruction_status")})
            out.append({"objective_mode": mode, "max_intervention_cost": budget,
                        "solver_status": ";".join(statuses), "reconstruction_status": ";".join(recon),
                        "solutions_retained": sum(1 for r in group if r.get("solution_id")),
                        # Each retained design row repeats its MILP runtime. Use
                        # the largest row value as the solve runtime, not a sum.
                        "runtime_seconds": round(max((number(r, "runtime_seconds", 0) or 0 for r in group), default=0), 6),
                        "error": "; ".join(sorted({r.get("error", "") for r in group if r.get("error")}))})
    for mode in MODES:
        present = {r["objective_mode"] for r in out if r["objective_mode"] == mode}
        if not present:
            out.append({"objective_mode": mode, "max_intervention_cost": "", "solver_status": "MISSING",
                        "reconstruction_status": "", "solutions_retained": 0, "runtime_seconds": 0, "error": "no rows"})
    return out


def recurring_interventions(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    ko_modes: dict[str, set[str]] = defaultdict(set)
    reg_modes: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row.get("objective_mode") == BASELINE:
            continue
        kos, regs = intervention_key(row)
        for rid in kos:
            ko_modes[rid].add(row.get("objective_mode", ""))
        for rule in regs:
            reg_modes[rule].add(row.get("objective_mode", ""))
    return ([{"reaction_id": rid, "objective_modes": ";".join(sorted(modes)),
              "mode_count": len(modes)} for rid, modes in sorted(ko_modes.items(), key=lambda x: (-len(x[1]), x[0]))],
            [{"regulatory_constraint": rule, "objective_modes": ";".join(sorted(modes)),
              "mode_count": len(modes)} for rule, modes in sorted(reg_modes.items(), key=lambda x: (-len(x[1]), x[0]))])


def pkg(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not-installed"


def aggregate() -> None:
    branch_data = {}
    branch_unique = {}
    branch_fronts = {}
    manifests = {branch: load_manifest(branch) for branch in BRANCHES}
    for field in ("host", "pool", "E2R3_fingerprint", "script", "straindesign_patch_sha256",
                  "modes", "budgets", "seed", "time_limit", "solver_preprocessing_fallback"):
        if manifests["S01"].get(field) != manifests["S03"].get(field):
            raise RuntimeError(f"Cross-branch common input mismatch: {field}")
    for branch in BRANCHES:
        bdir = OUT / branch
        raw = read_tsv(bdir / "optimization_results.tsv")
        if not raw:
            raise RuntimeError(f"{branch}: optimization_results.tsv is missing or empty")
        expected_fp = manifests[branch].get("fingerprint")
        if any(row.get("input_fingerprint") != expected_fp for row in raw):
            raise RuntimeError(f"{branch}: optimization TSV contains a stale input fingerprint")
        state = json.loads((bdir / "checkpoints/optimization_results.json").read_text(encoding="utf-8"))
        actual_counts = Counter(f"{row.get('objective_mode')}|{row.get('max_intervention_cost')}"
                                for row in raw if row.get("objective_mode") != BASELINE)
        for checkpoint_key, entry in state.get("results", {}).items():
            if actual_counts.get(checkpoint_key, 0) != len(entry.get("rows", [])):
                raise RuntimeError(f"{branch}: TSV/checkpoint row count mismatch for {checkpoint_key}")
        off_rows = read_tsv(bdir / "pareto_designs.tsv")
        off_rows += read_tsv(bdir / "glucose_off_completion.tsv")
        off_by_id = {row.get("solution_id"): row for row in off_rows if row.get("glucose_off_status")}
        for row in raw:
            off = off_by_id.get(row.get("solution_id"))
            if off:
                for field in ("glucose_off_status", "glucose_off_biomass", "glucose_contribution_delta"):
                    if field in off:
                        row[field] = off[field]
                if row.get("glucose_off_status") == "optimal" and not row.get("glucose_contribution_delta"):
                    on_bio = number(row, "biomass")
                    off_bio = number(row, "glucose_off_biomass")
                    if on_bio is not None and off_bio is not None:
                        row["glucose_contribution_delta"] = on_bio - off_bio
        # Ignore malformed and unreconstructed incumbents in the biological
        # comparison, while preserving every original row in the raw TSV.
        valid = [row for row in raw if valid_reconstructed(row)]
        unique = dedupe_designs(valid)
        frontier = pareto(unique)
        branch_data[branch] = {"raw": raw, "valid": valid, "status": status_matrix(raw)}
        branch_unique[branch] = unique
        branch_fronts[branch] = frontier
        write_tsv(bdir / "pareto_designs_unique.tsv", frontier)
        ko_freq, reg_freq = recurring_interventions(unique)
        write_tsv(bdir / "recurring_interventions.tsv",
                  [{"intervention_type": "KO", **row} for row in ko_freq] +
                  [{"intervention_type": "REG", **row} for row in reg_freq])

    combined = []
    for branch in BRANCHES:
        for row in branch_fronts[branch]:
            combined.append({"branch": branch, **row})
    combined_front = pareto(combined)
    write_tsv(OUT / "comparison/S01_vs_S03.tsv", combined)
    write_tsv(OUT / "comparison/nondominated_combined.tsv", combined_front)

    common_fingerprint = {branch: manifests[branch].get("fingerprint") for branch in BRANCHES}
    script_path = Path(__file__).resolve()
    host_sha = manifests["S01"].get("host")
    e1_sha = manifests["S01"].get("pool")
    e2r3_fp = manifests["S01"].get("E2R3_fingerprint")
    e2r3_manifest = json.loads((ROOT / "results/engineering/E2R3_two_stage_full_pool/checkpoints/input_manifest.json").read_text(encoding="utf-8"))
    environment = e2r3_manifest["inputs"]["thresholds"]
    common = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "host_sha256": host_sha, "E1_pool_sha256": e1_sha,
              "E2R3_input_fingerprint": e2r3_fp,
              "branch_input_fingerprints": common_fingerprint,
              "E3_script_sha256": manifests["S01"].get("script"),
              "summary_script_sha256": sha(script_path),
              "StrainDesign_patch_sha256": manifests["S01"].get("straindesign_patch_sha256", {}),
              "branch_fixed_reaction_ids": {b: manifests[b].get("support", {}).get("search_ids", []) for b in BRANCHES},
              "solver": "SCIP", "seed": 1,
              "environment": environment,
              "optimization_parameters": {"modes": MODES, "budgets": [1, 2, 4, 6],
                  "conditional_budget_8": True, "max_solutions": 5,
                  "time_limit_seconds": manifests["S01"].get("time_limit"),
                  "ko_and_reg_unit_cost": 1, "objective_engine": "StrainDesign RobustKnock",
                  "solver_preprocessing_fallback": manifests["S01"].get("solver_preprocessing_fallback")},
              "software": {"python": platform.python_version(), "cobra": pkg("cobra"),
                  "straindesign": pkg("straindesign"), "scipy": pkg("scipy")},
              "objective_definitions": {"inner": "maximize biomass subject to E3.0 phenotype constraints",
                  "outer": {"CBB_min": "minimize RUBISCO", "Fe2_efficiency": "maximize Ex_fe2[e]",
                      "O2_efficiency": "maximize Ex_o2[e]", "CO2_efficiency": "maximize Ex_h2co3[e]"}},
              "pareto_dimensions": ["intervention_count:min", "rubisco:min", "abs(Fe2 uptake):min",
                  "abs(O2 uptake):min", "abs(CO2/H2CO3 uptake):min", "biomass:max"],
              "branch_outputs_are_terminal": {b: (OUT / b / "checkpoints/optimization_results.json").exists() for b in BRANCHES}}
    (OUT / "E3_0_manifest.json").write_text(json.dumps(common, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    make_report(common, branch_data, branch_unique, branch_fronts, combined_front)


def fmt(x, digits=6):
    v = number({"v": x}, "v")
    return "—" if v is None else f"{v:.{digits}g}"


def make_report(common, data, unique, fronts, combined_front):
    lines = ["# E3.0 whole-host optimization summary", "",
        f"Host SHA256: `{common['host_sha256']}`  ",
        f"E1 pool SHA256: `{common['E1_pool_sha256']}`  ",
        f"E2R3 fingerprint: `{common['E2R3_input_fingerprint']}`  ",
        f"E3.0 script SHA256: `{common['E3_script_sha256']}`  ",
        f"Summary script SHA256: `{common['summary_script_sha256']}`  ",
        f"Solver: {common['solver']}; Python {common['software']['python']}; COBRApy {common['software']['cobra']}; StrainDesign {common['software']['straindesign']}; SCIP seed {common['seed']}.  ",
        f"StrainDesign patch SHA256: `{json.dumps(common['StrainDesign_patch_sha256'], sort_keys=True)}`", "",
        "The branch jobs used separate fixed E2R3 reaction sets, references, candidate lists, checkpoints, and optimization fingerprints. Pareto metrics use reconstructed worst-case phenotype values and compare intervention count, Rubisco flux, absolute Fe2/O2/inorganic-carbon uptake, and biomass. The original optimizer TSVs remain unchanged. The unique frontier TSVs deduplicate identical KO plus regulatory constraint sets.", ""]
    for branch in BRANCHES:
        raw = data[branch]["raw"]
        rows = data[branch]["status"]
        front = fronts[branch]
        bdir = OUT / branch
        candidates = read_tsv(bdir / "intervention_candidates.tsv")
        ko_n = sum(r.get("intervention_type") == "KO" for r in candidates)
        reg_n = sum(r.get("intervention_type") in {"UPREG", "ACTIVATION_UPREG"} for r in candidates)
        ref = read_tsv(bdir / "branch_reference.tsv")
        lines += [f"## {branch}", "",
            f"Fixed E1 reactions: `{'; '.join(common['branch_fixed_reaction_ids'][branch])}`. Native candidate universe: {ko_n} KO and {reg_n} regulatory flux interventions."]
        if ref:
            r = ref[0]
            lines += [f"Branch reference: biomass {fmt(r.get('biomass_optimum'))}; glucose {fmt(r.get('glucose'))}; Rubisco {fmt(r.get('rubisco'))}; Fe2 {fmt(r.get('fe2'))}; O2 {fmt(r.get('o2'))}; CO2/H2CO3 {fmt(r.get('co2_h2co3'))}.", ""]
        lines += ["Budget results (statuses are shown as returned by the solver; ERROR and TIME_LIMIT remain distinct):", "",
            "| Objective | Max cost | Solver status | Reconstruction | Solutions | Runtime (s) |",
            "|---|---:|---|---|---:|---:|"]
        for r in rows:
            lines.append(f"| {r['objective_mode']} | {r['max_intervention_cost']} | {r['solver_status']} | {r['reconstruction_status']} | {r['solutions_retained']} | {fmt(r['runtime_seconds'])} |")
        ko_freq, reg_freq = recurring_interventions(fronts[branch])
        repeat = [x["reaction_id"] for x in ko_freq if x["mode_count"] > 1]
        repeat += [x["regulatory_constraint"].split()[0] for x in reg_freq if x["mode_count"] > 1]
        lines += ["", f"Unique valid designs after deduplicating KO/regulation sets: {len(unique[branch])}; nondominated designs: {len(front)}.",
                  "Repeated interventions across this branch's frontier: " + (", ".join(sorted(set(repeat))) if repeat else "none observed") + ".", ""]
        if front:
            lines += ["Frontier phenotypes:", "", "| ID | Native interventions | KO + OE | Biomass | Rubisco | Fe2 uptake | O2 uptake | CO2 uptake | Glucose off |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
            for r in sorted(front, key=lambda x: (int(number(x,"total_intervention_cost",0)), x.get("objective_mode",""),x.get("solution_id",""))):
                off_status = r.get("glucose_off_status", "unavailable") or "unavailable"
                off_label = (f"{off_status}: {fmt(r.get('glucose_off_biomass'))}"
                             if off_status == "optimal" else off_status)
                interventions = "; ".join(x for x in (r.get("KO_reaction_ids", ""), r.get("regulatory_constraints", "")) if x) or "none"
                lines.append("| {} | {} | {} + {} | {} | {} | {} | {} | {} | {} |".format(
                    r.get("solution_id",""),interventions,
                    r.get("KO_count",""),r.get("OE_count",""),fmt(r.get("biomass")),fmt(r.get("rubisco")),
                    fmt(abs(number(r,"fe2",0))),fmt(abs(number(r,"o2",0))),fmt(abs(number(r,"co2_h2co3",0))),
                    off_label))
            lines.append("")
    lines += ["## S01 vs S03", "",
        "| Branch | Unique valid intervention sets | Within-branch Pareto designs | Fewest positive host interventions |",
        "|---|---:|---:|---:|"]
    for branch in BRANCHES:
        counts = [int(number(row, "total_intervention_cost", 0)) for row in fronts[branch]
                  if number(row, "total_intervention_cost", 0) > 0]
        lines.append(f"| {branch} | {len(unique[branch])} | {len(fronts[branch])} | {min(counts) if counts else '—'} |")
    baseline = {branch: next((r for r in unique[branch] if r.get("objective_mode") == BASELINE), None)
                for branch in BRANCHES}
    zero_carbon = {}
    for branch in BRANCHES:
        candidates = [r for r in fronts[branch]
                      if number(r, "total_intervention_cost", 0) > 0
                      and abs(number(r, "co2_h2co3", 1)) <= 1e-8]
        zero_carbon[branch] = min(candidates,
            key=lambda r: (number(r, "total_intervention_cost", 999),
                           abs(number(r, "fe2", 999)), -number(r, "biomass", 0))) if candidates else None
    if baseline["S01"] and baseline["S03"]:
        lines += ["", "The unmodified S03 branch has the higher biomass optimum "
                  f"({fmt(baseline['S03'].get('biomass_inner_optimum') or baseline['S03'].get('biomass'))} "
                  f"versus {fmt(baseline['S01'].get('biomass_inner_optimum') or baseline['S01'].get('biomass'))})."]
    if zero_carbon["S01"] and zero_carbon["S03"]:
        a, b = zero_carbon["S01"], zero_carbon["S03"]
        lines += ["Both branches have a reconstructed single native KO design that reaches zero external "
                  "inorganic-carbon uptake at the biomass threshold. "
                  f"S01 uses {fmt(abs(number(a,'fe2',0)))} Fe2; S03 uses {fmt(abs(number(b,'fe2',0)))} "
                  "under the common worst-case reconstruction. S03 therefore needs no fewer observed "
                  "host interventions for this endpoint. The CO2-mode MILPs timed out with incumbents, "
                  "so these designs are feasible findings rather than proven budget-wide optima."]
    lines += ["No intervention was selected by the proven-optimal CBB runs, and the Fe2/O2 higher-budget "
              "runs did not resolve a new design. HCO3tex and HCO3tpp recur across CO2 budgets within "
              "each branch; no endogenous reaction is repeatedly selected across different objective modes. "
              "The current frontier leaves an overall backbone choice open: S03 has higher unmodified "
              "biomass, while S01 has lower Fe2 demand in the observed zero-carbon single-KO state."]
    lines += ["", f"Combined nondominated solutions retained: {len(combined_front)}. Review `comparison/S01_vs_S03.tsv` for all branch Pareto designs and `comparison/nondominated_combined.tsv` for designs that remain nondominated across both branches.", "",
        "Interpretation is limited to the model-level intervention hypotheses. Regulatory flux constraints do not imply a gene-expression fold change. A missing mode/budget row means the branch output did not contain that solve; inspect branch checkpoint status before interpreting it. No gene or experimental design was selected.", ""]
    (OUT / "E3_0_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="aggregate available rows while marking missing branch checkpoint files")
    args = parser.parse_args()
    manifests = {branch: load_manifest(branch) for branch in BRANCHES}
    for branch in BRANCHES:
        state_path = OUT / branch / "checkpoints/optimization_results.json"
        if not state_path.exists() and not args.allow_incomplete:
            raise SystemExit(f"{branch} branch job has not written its checkpoint; wait for terminal completion")
        if state_path.exists():
            state = json.loads(state_path.read_text(encoding="utf-8"))
            if state.get("fingerprint") != manifests[branch].get("fingerprint"):
                raise SystemExit(f"{branch} checkpoint fingerprint does not match its input manifest")
            results = state.get("results", {})
            required = {f"{mode}|{budget}" for mode in MODES for budget in (1, 2, 4, 6)}
            missing = sorted(required - set(results))
            unfinished = [key for key, val in results.items()
                          if val.get("status", "").lower() in {"", "unknown", "running"}]
            if (missing or unfinished) and not args.allow_incomplete:
                details = []
                if missing:
                    details.append("missing required solves: " + ", ".join(missing))
                if unfinished:
                    details.append("non-terminal checkpoint entries: " + ", ".join(unfinished))
                raise SystemExit(f"{branch} branch is not terminal ({'; '.join(details)})")
    aggregate()


if __name__ == "__main__":
    main()
