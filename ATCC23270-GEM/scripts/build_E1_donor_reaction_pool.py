#!/usr/bin/env python
"""Build the E1 non-native intracellular reaction pool from three donor GEMs.

Run from any directory with the repository's .venv Python. Outputs are written
under results/engineering/E1_donor_reaction_pool; the host model is read-only.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import cobra


ROOT = Path(__file__).resolve().parents[1]
HOST = ROOT.parent / "before final" / "updatedv3.1_verified" / "updatedv3.1.xml"
EXPECTED_HOST_SHA256 = "026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3"
DONORS = [
    ("iML1515", ROOT / "external/reference/donor_models/iML1515.json"),
    ("iJN1463", ROOT / "external/reference/donor_models/iJN1463.xml"),
    ("iCN1361", ROOT / "external/reference/donor_models/iCN1361.xml"),
]
DONOR_URLS = {
    "iML1515": "Reused local file; reference page https://bigg.ucsd.edu/models/iML1515",
    "iJN1463": "Downloaded source https://bigg.ucsd.edu/static/models/iJN1463.xml",
    "iCN1361": "Reused local file; author model repository https://github.com/SBRCNottingham/CnecatorGSM/tree/main/Model",
}
OUT = ROOT / "results/engineering/E1_donor_reaction_pool"
ALIAS_FIELDS = {
    "bigg.metabolite": "bigg",
    "metanetx.chemical": "mnx",
    "mnx": "mnx",
    "kegg.compound": "kegg",
    "kegg": "kegg",
    "biocyc": "biocyc",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def values(v):
    if v is None:
        return []
    return v if isinstance(v, (list, tuple, set)) else [v]


def clean_alias(kind: str, value) -> str:
    s = str(value).strip()
    s = s.rsplit("/", 1)[-1].split(":")[-1]
    s = re.sub(r"\s+", "", s)
    if kind == "bigg":
        s = re.sub(r"_\[[^]]+\]$", "", s)
    return f"{kind}:{s.upper()}"


def metabolite_aliases(met, host_legacy=False):
    result = []
    for field, kind in ALIAS_FIELDS.items():
        for value in values(met.annotation.get(field)):
            alias = clean_alias(kind, value)
            if alias not in result:
                result.append(alias)
    # Explicit BiGG annotation always wins. Older host IDs encode BiGG's
    # double underscore as _DASH_; use the ID only when no BiGG annotation exists.
    if not any(x.startswith("bigg:") for x in result):
        raw_id = str(met.id)
        # The host's legacy BiGG encoding is explicit and safe to normalize.
        legacy = re.match(r"^(?:M_)?(.*?)(?:_DASH_|__)([A-Za-z0-9]+)_([cpe])$", raw_id, re.I)
        if legacy:
            alias = clean_alias("bigg", f"{legacy.group(1)}__{legacy.group(2)}")
        elif host_legacy and re.search(r"_[cpe]$", raw_id, re.I):
            base = re.sub(r"^M_", "", raw_id, flags=re.I)
            base = re.sub(r"_[cpe]$", "", base, flags=re.I)
            alias = clean_alias("bigg", base)
        elif not raw_id.startswith("M_") and re.search(r"_[cpe]$", raw_id, re.I):
            # Generic fallback is limited to known compartment suffixes.
            base = re.sub(r"_[cpe]$", "", raw_id, flags=re.I)
            alias = clean_alias("bigg", base)
        else:
            alias = None
        if alias and alias not in result:
            result.insert(0, alias)
    return result


def compartment_class(met) -> str:
    c = str(met.compartment or "").lower()
    name = str(getattr(met, "name", "")).lower()
    if c in {"c", "cytosol", "cytoplasm", "cytoplasmic"} or any(x in name for x in ("cytosol", "cytoplasm")):
        return "cytosol"
    return c or "unknown"


def is_pseudo_or_biomass(rxn) -> str | None:
    rid = rxn.id.lower()
    name = (rxn.name or "").lower()
    if re.search(r"(?:^|[_-])(?:biomass|bio1)(?:$|[_-])", rid) or "biomass" in name or abs(rxn.objective_coefficient) > 0:
        return "biomass"
    if rxn.boundary or rid.startswith(("ex_", "dm_", "sk_", "sink_", "demand_")):
        return "boundary_pseudo"
    return None


def build_alias_index(models):
    """Alias is usable only if it identifies one metabolite per model.

    Distinct compartments of the same mapped metabolite are allowed; conflicting
    chemical identities under one alias make that alias unusable in that model.
    """
    index = defaultdict(lambda: defaultdict(set))
    for label, model in models:
        for met in model.metabolites:
            for alias in metabolite_aliases(met, host_legacy=(label == "host")):
                index[alias][label].add(met.id)
    return index


def build_alias_canonical(models, alias_index):
    """Join co-annotated xrefs after rejecting model-local alias collisions."""
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    for label, model in models:
        for met in model.metabolites:
            aliases = metabolite_aliases(met, host_legacy=(label == "host"))
            safe = []
            for alias in aliases:
                ids = alias_index[alias].get(label, set())
                bases = {re.sub(r"_[A-Za-z0-9]+$", "", re.sub(r"^M_", "", x)) for x in ids}
                met_base = re.sub(r"_[A-Za-z0-9]+$", "", re.sub(r"^M_", "", met.id))
                if met.id in ids and len(bases) == 1 and met_base in bases:
                    safe.append(alias)
            for alias in safe:
                find(alias)
            for alias in safe[1:]:
                union(safe[0], alias)
    components = defaultdict(list)
    for alias in parent:
        components[find(alias)].append(alias)
    # A transitive xref component is unsafe if it lands on multiple distinct
    # model metabolites in the same model and compartment. Disable the whole
    # component's cross-reference joins in that case.
    component_memberships = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for label, model in models:
        for met in model.metabolites:
            aliases = metabolite_aliases(met, host_legacy=(label == "host"))
            safe_roots = {find(a) for a in aliases if a in parent}
            base_id = re.sub(r"_[A-Za-z0-9]+$", "", re.sub(r"^M_", "", met.id))
            for root in safe_roots:
                component_memberships[root][label][compartment_class(met)].add(base_id)
    conflicting_roots = {
        root for root, by_model in component_memberships.items()
        if any(len(ids) > 1 for by_compartment in by_model.values() for ids in by_compartment.values())
    }
    def rank(alias):
        p = 0 if alias.startswith("bigg:") else 1 if alias.startswith("mnx:") else 2 if alias.startswith("kegg:") else 3
        return p, alias
    canonical = {}
    for root, members in components.items():
        representative = min(members, key=rank)
        for alias in members:
            canonical[alias] = alias if root in conflicting_roots else representative
    return canonical, len(conflicting_roots)


def metabolite_key(met, model_label, alias_index, alias_canonical):
    aliases = metabolite_aliases(met, host_legacy=(model_label == "host"))
    shared = [a for a in aliases if len(alias_index[a]) > 1]
    ordered = sorted(shared, key=lambda a: (0 if a.startswith("bigg:") else 1 if a.startswith("mnx:") else 2, a))
    ordered += [a for a in aliases if a not in ordered]
    for alias in ordered:
        ids = alias_index[alias].get(model_label, set())
        # One model must not assign an alias to distinct internal IDs.
        base_ids = {re.sub(r"_[A-Za-z0-9]+$", "", re.sub(r"^M_", "", x)) for x in ids}
        met_base = re.sub(r"_[A-Za-z0-9]+$", "", re.sub(r"^M_", "", met.id))
        if met.id in ids and len(base_ids) == 1 and met_base in base_ids:
            return alias_canonical.get(alias, alias), "mapped"
    return f"unmapped:{model_label}:{met.id}", "unmapped"


def reaction_signature(rxn, label, alias_index, alias_canonical):
    terms = []
    quality = []
    for met, coeff in rxn.metabolites.items():
        key, status = metabolite_key(met, label, alias_index, alias_canonical)
        quality.append(status)
        # Canonicalize cytosol labels so c/cytosol naming differences compare.
        compartment = compartment_class(met)
        if compartment == "cytosol":
            key += "@cytosol"
        else:
            key += "@" + compartment
        q = Fraction(float(coeff)).limit_denominator(100000)
        terms.append((key, q))
    # Convert coefficients to primitive integer ratios, independent of scale.
    den = math.lcm(*(q.denominator for _, q in terms))
    ints = [(k, int(q * den)) for k, q in terms]
    divisor = math.gcd(*(abs(n) for _, n in ints)) or 1
    ints = [(k, n // divisor) for k, n in ints]
    forward = tuple(sorted(ints))
    reverse = tuple(sorted((k, -n) for k, n in ints))
    if reverse < forward:
        return reverse, "reverse", "unmapped" if "unmapped" in quality else "mapped"
    return forward, "forward", "unmapped" if "unmapped" in quality else "mapped"


def canonical_directions(rxn, orientation):
    directions = set()
    if rxn.upper_bound > 0:
        directions.add(orientation)
    if rxn.lower_bound < 0:
        directions.add("reverse" if orientation == "forward" else "forward")
    return directions


def equation(rxn):
    return rxn.reaction


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    observed_sha = sha256(HOST)
    if observed_sha != EXPECTED_HOST_SHA256:
        raise SystemExit(f"Host SHA256 mismatch: expected {EXPECTED_HOST_SHA256}, observed {observed_sha}")
    # The verified historical SBML has empty optional reactant/product lists in
    # a few reactions. COBRApy's validator still parses the model successfully.
    # libSBML on Windows may fail to open Unicode absolute paths. Use paths
    # relative to the repository root for COBRApy, while retaining absolute
    # paths for the report and SHA verification.
    os.chdir(ROOT)
    host, host_validation = cobra.io.validate_sbml_model("../before final/updatedv3.1_verified/updatedv3.1.xml")
    if host is None:
        raise SystemExit("COBRApy could not parse the verified host model")
    donors = [(label, cobra.io.load_json_model(path.relative_to(ROOT).as_posix()) if path.suffix.lower() == ".json"
              else cobra.io.read_sbml_model(path.relative_to(ROOT).as_posix())) for label, path in DONORS]
    all_models = [("host", host), *donors]
    alias_index = build_alias_index(all_models)
    alias_canonical, alias_component_conflicts = build_alias_canonical(all_models, alias_index)

    host_signatures = set()
    host_directions = defaultdict(set)
    for rxn in host.reactions:
        if is_pseudo_or_biomass(rxn):
            continue
        sig, orientation, _ = reaction_signature(rxn, "host", alias_index, alias_canonical)
        host_signatures.add(sig)
        host_directions[sig].update(canonical_directions(rxn, orientation))

    stats = {}
    donor_pool = {}
    deferred = {}
    mapping_problems = Counter()
    reaction_xref_signatures = defaultdict(set)
    for label, model in donors:
        counts = Counter()
        counts["read"] = len(model.reactions)
        for rxn in model.reactions:
            pseudo = is_pseudo_or_biomass(rxn)
            if pseudo:
                counts["pseudo_biomass_removed"] += 1
                continue
            compartments = {compartment_class(m) for m in rxn.metabolites}
            # Keep transport reactions distinct from reactions wholly in a
            # single non-cytosolic compartment.
            if len(compartments) > 1 or compartments != {"cytosol"}:
                deferred_type = "transport" if len(compartments) > 1 else "non_cytosolic"
                counts[deferred_type] += 1
                sig, orient, quality = reaction_signature(rxn, label, alias_index, alias_canonical)
                for xref in values(rxn.annotation.get("bigg.reaction")):
                    reaction_xref_signatures[clean_alias("bigg.reaction", xref)].add(sig)
                rec = deferred.setdefault(sig, {"signature": repr(sig), "donors": [], "reactions": [], "equations": [], "mapping_quality": quality, "category": deferred_type})
                rec["donors"].append(label); rec["reactions"].append(f"{label}:{rxn.id}"); rec["equations"].append(f"{label}: {equation(rxn)}")
                if quality == "unmapped":
                    rec["mapping_quality"] = "unmapped"
                if quality == "unmapped":
                    mapping_problems["deferred_unmapped"] += 1
                continue
            counts["intracellular_considered"] += 1
            sig, orient, quality = reaction_signature(rxn, label, alias_index, alias_canonical)
            for xref in values(rxn.annotation.get("bigg.reaction")):
                reaction_xref_signatures[clean_alias("bigg.reaction", xref)].add(sig)
            if quality == "unmapped":
                mapping_problems["intracellular_unmapped"] += 1
            rec = donor_pool.setdefault(sig, {"signature": repr(sig), "donors": [], "reactions": [], "equations": [], "mapping_quality": quality, "orientation_to_canonical": []})
            rec["donors"].append(label)
            rec["reactions"].append(f"{label}:{rxn.id}")
            rec["equations"].append(f"{label}: {equation(rxn)}")
            rec["orientation_to_canonical"].append(f"{label}:{rxn.id}={orient}")
            if quality == "unmapped":
                rec["mapping_quality"] = "unmapped"
            if sig in host_signatures:
                counts["host_equivalent"] += 1
                if not canonical_directions(rxn, orient).issubset(host_directions[sig]):
                    counts["host_direction_difference"] += 1
        stats[label] = counts
        classified = counts["pseudo_biomass_removed"] + counts["transport"] + counts["non_cytosolic"] + counts["intracellular_considered"]
        if classified != counts["read"]:
            raise RuntimeError(f"Reaction classification mismatch for {label}: {classified} != {counts['read']}")

    # Deduplicate donor chemistry first, then remove the unique host-equivalent
    # signatures. This makes the requested sequence and counts explicit.
    donor_unique_count = len(donor_pool)
    host_equivalent_unique = set(donor_pool).intersection(host_signatures)
    candidates = {sig: rec for sig, rec in donor_pool.items() if sig not in host_equivalent_unique}
    deferred_unique_before_host = len(deferred)
    deferred_host_equivalent = set(deferred).intersection(host_signatures)
    deferred = {sig: rec for sig, rec in deferred.items() if sig not in deferred_host_equivalent}
    intracellular_source_rows = sum(c["intracellular_considered"] for c in stats.values())
    donor_dedup_removed = intracellular_source_rows - donor_unique_count
    fields = ["candidate_id", "donors", "source_reactions", "source_equations", "mapping_quality", "canonical_signature", "orientation_to_canonical"]
    candidate_ids = {}
    with (OUT / "candidate_reactions.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("\t".join(fields) + "\n")
        for n, (sig, rec) in enumerate(sorted(candidates.items(), key=lambda x: repr(x[0])), 1):
            record_id = f"E1R{n:05d}"
            candidate_ids[sig] = record_id
            row = [record_id, ";".join(sorted(set(rec["donors"]))), ";".join(rec["reactions"]),
                   ";".join(rec["equations"]), rec["mapping_quality"], rec["signature"], ";".join(rec["orientation_to_canonical"])]
            f.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in row) + "\n")
    fields = ["deferred_id", "category", "donors", "source_reactions", "source_equations", "mapping_quality", "canonical_signature"]
    deferred_ids = {}
    with (OUT / "deferred_transport.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("\t".join(fields) + "\n")
        for n, (sig, rec) in enumerate(sorted(deferred.items(), key=lambda x: repr(x[0])), 1):
            record_id = f"E1T{n:05d}"
            deferred_ids[sig] = record_id
            row = [record_id, rec["category"], ";".join(sorted(set(rec["donors"]))), ";".join(rec["reactions"]),
                   ";".join(rec["equations"]), rec["mapping_quality"], rec["signature"]]
            f.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in row) + "\n")

    unresolved_candidates = [(sig, rec) for sig, rec in candidates.items() if rec["mapping_quality"] == "unmapped"]
    unresolved_deferred = [(sig, rec) for sig, rec in deferred.items() if rec["mapping_quality"] == "unmapped"]
    with (OUT / "unresolved_mapping.tsv").open("w", encoding="utf-8", newline="") as f:
        f.write("record_type\trecord_id\tdonors\tsource_reactions\tsource_equations\tcanonical_signature\n")
        for kind, rows, ids in (("candidate", unresolved_candidates, candidate_ids), ("deferred", unresolved_deferred, deferred_ids)):
            for sig, rec in sorted(rows, key=lambda x: repr(x[0])):
                row = [kind, ids[sig], ";".join(sorted(set(rec["donors"]))), ";".join(rec["reactions"]),
                       ";".join(rec["equations"]), rec["signature"]]
                f.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in row) + "\n")

    report = [
        "# E1 donor reaction pool",
        "",
        f"- Host: `{HOST}`",
        f"- Observed host SHA256: `{observed_sha}` (expected hash matched)",
        f"- Host reactions/metabolites: {len(host.reactions)}/{len(host.metabolites)}",
        "- Donor files and sources: " + "; ".join(f"{label}: `{path}` (SHA256 `{sha256(path)}`; {DONOR_URLS[label]})" for label, path in DONORS),
        "",
        "## Counts",
        "",
        "| Donor | Read | Pseudo/biomass removed | Cross-compartment transport | Single non-cytosolic | Intracellular checked | Host equivalent rows | Direction difference source rows |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for label, _ in DONORS:
        c = stats[label]
        report.append(f"| {label} | {c['read']} | {c['pseudo_biomass_removed']} | {c['transport']} | {c['non_cytosolic']} | {c['intracellular_considered']} | {c['host_equivalent']} | {c['host_direction_difference']} |")
    report += [
        "",
        f"- Unique donor reactions after donor-to-donor deduplication and host-equivalence removal: {len(candidates)}",
        f"- Intracellular donor source rows before deduplication: {intracellular_source_rows}",
        f"- Unique donor reactions after donor-to-donor deduplication, before host filtering: {donor_unique_count}",
        f"- Unique donor reactions equivalent to host and removed: {len(host_equivalent_unique)}",
        f"- Host-equivalent intracellular source rows with donor direction capacity beyond host bounds: {sum(c['host_direction_difference'] for c in stats.values())}",
        f"- Shared BiGG reaction xrefs associated with multiple canonical signatures (diagnostic only): {sum(1 for sigs in reaction_xref_signatures.values() if len(sigs) > 1)}",
        f"- Transitive metabolite alias components rejected for same-model, same-compartment conflicts: {alias_component_conflicts}",
        f"- Source rows collapsed by donor-to-donor chemical deduplication: {donor_dedup_removed}",
        f"- Deferred unique donor rows before host filtering: {deferred_unique_before_host}",
        f"- Deferred rows equivalent to host and removed: {len(deferred_host_equivalent)}",
        f"- Deferred unique rows retained (transport and non-cytosolic): {len(deferred)}",
        "",
        "## Mapping and boundaries",
        "",
        "Metabolite identifiers use BiGG, MetaNetX, KEGG, and BioCyc annotations as a conflict-checked alias graph: co-annotations connect aliases only when each alias identifies one chemical per model. A transitive component is disabled if it maps to multiple distinct metabolites in the same model and compartment. The graph uses a deterministic BiGG/MNX/KEGG/BioCyc canonical label. Legacy host `_DASH_` IDs are normalized; generic `M_*` donor IDs are not used as fallback aliases. Names and formulas do not force matches. Reaction stoichiometry is reduced to primitive integer ratios and canonicalized over forward/reverse orientation. Host equivalence is therefore assessed without direction; differing donor/host bounds are counted above for later review. Shared BiGG reaction xrefs are counted as split-signature diagnostics and never force a merge.",
        "",
        f"Unmapped intracellular reaction-source rows: {mapping_problems['intracellular_unmapped']}; after donor deduplication and host filtering, {len(unresolved_candidates)} final candidate rows still contain unresolved metabolite mappings. Unmapped deferred source rows: {mapping_problems['deferred_unmapped']}; these yield {len(unresolved_deferred)} unresolved deferred rows. The exact donor/reaction IDs and equations are listed in `unresolved_mapping.tsv`.",
        "",
        "Exchanges, demands, sinks, boundary reactions, and biomass reactions are excluded. Cross-compartment reactions are categorized as transport; reactions wholly within one non-cytosolic compartment are categorized separately. Host-equivalent deferred signatures are removed. This output is a reaction candidate pool only; no strain design or model edit was run.",
        "",
        "## Reproduction",
        "",
        f"Run with Python and COBRApy 0.32.1 installed (for example `pip install cobra==0.32.1`): `python \"{ROOT / 'scripts/build_E1_donor_reaction_pool.py'}\"`.",
    ]
    (OUT / "E1_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"host_reactions": len(host.reactions), "host_metabolites": len(host.metabolites), "donors": {k: dict(v) for k, v in stats.items()}, "intracellular_source_rows": intracellular_source_rows, "donor_unique_before_host_filter": donor_unique_count, "host_equivalent_unique_removed": len(host_equivalent_unique), "candidate_rows": len(candidates), "deferred_unique_before_host_filter": deferred_unique_before_host, "deferred_host_equivalent_removed": len(deferred_host_equivalent), "deferred_rows": len(deferred), "donor_duplicate_rows": donor_dedup_removed, "unresolved_candidate_rows": len(unresolved_candidates), "unresolved_deferred_rows": len(unresolved_deferred), "alias_component_conflicts": alias_component_conflicts, "mapping_problems": dict(mapping_problems), "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
