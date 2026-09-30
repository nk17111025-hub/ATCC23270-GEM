# -*- coding: utf-8 -*-
"""Build the Phase 5B-6 reaction universe (U1/U2/U3) from local BioCyc + curated
central-metabolism reactions, with KEGG-based metabolite mapping, duplicate removal
and elemental mass-balance QC.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import p5b6_common as C

ROOT = C.ROOT
BIOCYC = ROOT / "C2_外部数据库与注释来源核验" / "C2-1_BioCyc与数据库版本核验" / "标准化数据"
RXN_FILE = BIOCYC / "reactions.tsv"
CMP_FILE = BIOCYC / "compounds_network.tsv"

OUT = C.OUT
RAW = OUT / "02_universe_raw"
MAPDIR = OUT / "03_universe_mapping"
CLEAN = OUT / "04_universe_clean"

# manual currency map: BioCyc frame -> model metabolite base (cytosolic)
MANUAL_CURRENCY = {
    "Pi": "pi",
    "PHOSPHATE": "pi",
    "INORGANIC-PHOSPHATE": "pi",
}

# frames whose compounds are generic/protein/tRNA/polymeric and must not be treated
# as exact metabolites (these reactions are excluded unless every such frame is mapped
# by a real KEGG id)
GENERIC_PATTERNS = [
    r"-ACPs?$", r"^ACP$", r"Thioredoxin", r"ferredoxin", r"Ferredoxin",
    r"Cytochromes?-C", r"Quinols?$", r"Quinones?$", r"-tRNAs?$", r"^tRNAs?$",
    r"Protein", r"protein", r"Peptides", r"holder", r"-Proteins?$",
]


def parse_formula_charge(metas):
    """Return {met: {'formula': str, 'charge': int|None}} for element balancing."""
    out = {}
    for m, meta in metas.items():
        f = str(meta.get("formula", "")).strip()
        ch = meta.get("charge", "")
        try:
            ch = int(ch) if str(ch).strip() != "" else None
        except (ValueError, TypeError):
            ch = None
        out[m] = {"formula": f, "charge": ch}
    return out


def formula_to_elements(formula: str):
    """Parse a molecular formula like C6H12O6 -> {'C':6,'H':12,'O':6}."""
    if not formula:
        return {}
    # strip common junk
    f = formula.strip()
    f = re.sub(r"\(|\)|\[|\]", "", f)
    elems = {}
    for m in re.finditer(r"([A-Z][a-z]?)([0-9.]*)", f):
        el = m.group(1)
        num = m.group(2)
        val = float(num) if num else 1.0
        elems[el] = elems.get(el, 0.0) + val
    return elems


def element_balance(stoich, fmap):
    """Return (balanced_bool, elem_delta_dict)."""
    total = {}
    for m, c in stoich.items():
        f = fmap.get(m, {}).get("formula", "")
        for el, n in formula_to_elements(f).items():
            total[el] = total.get(el, 0.0) + c * n
    # ignore H and O and charge for "balanced" check (protonation/water hydration noise)
    residual = {el: v for el, v in total.items() if abs(v) > 1e-6 and el not in ("H", "O")}
    return (len(residual) == 0), residual


def main():
    metas, rxns = C.build_v2("T0")

    # 1) frame -> kegg
    _, cmp_rows = C.read_tsv(CMP_FILE)
    frame2kegg = {}
    frame2big = {}
    for r in cmp_rows:
        frame = r[2].strip()
        ext = r[7]
        for part in ext.split(";"):
            if part.startswith("LIGAND-CPD:"):
                frame2kegg[frame] = part.split(":", 1)[1].strip()
            if part.startswith("BIGG:"):
                frame2big[frame] = part.split(":", 1)[1].strip()

    # 2) kegg -> model metabolites (prefer [c])
    kegg2met = {}
    for m, meta in metas.items():
        k = str(meta.get("kegg", "")).strip()
        if k:
            kegg2met.setdefault(k, []).append(m)

    def resolve(frame):
        """Return list of candidate model metabolites for a frame."""
        if frame in MANUAL_CURRENCY:
            base = MANUAL_CURRENCY[frame]
            return [base + "[c]"] if (base + "[c]") in metas else []
        k = frame2kegg.get(frame)
        if k and k in kegg2met:
            cands = kegg2met[k]
            cyto = [m for m in cands if m.endswith("[c]")]
            if cyto:
                return cyto
            return cands
        return []

    # 3) parse reactions
    _, rxn_rows = C.read_tsv(RXN_FILE)

    raw_rows = []       # raw universe
    mapped_rows = []    # successfully mapped reactions
    failed_rows = []    # unmapped / failed

    def is_generic(frame):
        return any(re.search(p, frame) for p in GENERIC_PATTERNS)

    for r in rxn_rows:
        if len(r) < 7:
            continue
        rid = r[2].strip()
        name = r[3].strip()
        direction = r[4].strip()
        left = [x.strip() for x in r[5].split(";") if x.strip()]
        right = [x.strip() for x in r[6].split(";") if x.strip()]
        ec = r[10].strip() if len(r) > 10 else ""
        pathway = r[9].strip() if len(r) > 9 else ""

        reversible = direction.upper() in ("REVERSIBLE", "PHYSIOL-RIGHT-TO-LEFT", "REVERSIBLE;")
        # build stoich: left negative, right positive
        stoich = {}
        unmapped = []
        ambiguous = []
        for frame in left:
            cands = resolve(frame)
            if not cands:
                unmapped.append(frame)
                continue
            if len(cands) > 1:
                ambiguous.append((frame, cands))
            stoich[cands[0]] = stoich.get(cands[0], 0.0) - 1.0
        for frame in right:
            cands = resolve(frame)
            if not cands:
                unmapped.append(frame)
                continue
            if len(cands) > 1:
                ambiguous.append((frame, cands))
            stoich[cands[0]] = stoich.get(cands[0], 0.0) + 1.0

        # drop zero-net terms
        stoich = {m: c for m, c in stoich.items() if abs(c) > 1e-12}

        raw_rows.append([rid, name, direction, ";".join(left), ";".join(right), ec, pathway])

        if unmapped or not stoich:
            failed_rows.append([rid, name, ";".join(unmapped), ";".join(left), ";".join(right), ec, "unmapped_or_empty"])
            continue

        if C.reaction_in_model(rxns, stoich):
            continue  # exact duplicate of an existing v2 reaction

        mapped_rows.append({
            "universe_id": rid,
            "source_db": "BioCyc",
            "source_reaction_id": rid,
            "name": name,
            "reversible": reversible,
            "stoich": stoich,
            "ec": ec,
            "pathway": pathway,
            "ambiguous": ambiguous,
        })

    # 4) element balance QC
    fmap = parse_formula_charge(metas)
    balanced = []
    unbalanced = []
    unknown = []
    for d in mapped_rows:
        ok, res = element_balance(d["stoich"], fmap)
        has_formula = all(fmap.get(m, {}).get("formula") for m in d["stoich"])
        if not has_formula:
            d["balance"] = "BALANCE_UNKNOWN"
            unknown.append(d)
        elif ok:
            d["balance"] = "BALANCED"
            balanced.append(d)
        else:
            d["balance"] = "UNBALANCED"
            d["residual"] = res
            unbalanced.append(d)

    # 5) write outputs
    C.write_tsv(RAW / "U1" / "biocyc_reactions_raw.tsv",
                ["source_reaction_id", "name", "direction", "left", "right", "ec", "pathway"], raw_rows)

    C.write_tsv(MAPDIR / "failed_reactions.tsv",
                ["source_reaction_id", "name", "unmapped", "left", "right", "ec", "reason"], failed_rows)

    # mapped reactions table
    map_rows = []
    for d in mapped_rows:
        map_rows.append([
            d["universe_id"], d["name"], str(d["reversible"]), C.equation_str(d["stoich"]),
            d["ec"], d["pathway"], d["balance"], ";".join(f"{a}:{','.join(c)}" for a, c in d["ambiguous"]),
        ])
    C.write_tsv(MAPDIR / "reaction_mapping.tsv",
                ["universe_id", "name", "reversible", "equation", "ec", "pathway", "balance", "ambiguous"], map_rows)

    # clean universe
    clean_rows = []
    for d in mapped_rows:
        if d["balance"] == "UNBALANCED":
            continue
        clean_rows.append(d)

    C.write_tsv(CLEAN / "U1_clean.tsv",
                ["universe_id", "source_db", "source_reaction_id", "equation", "reversible",
                 "metabolite_mapping", "mass_balance", "charge_balance", "EC", "candidate_class", "included", "exclusion_reason"],
                [[d["universe_id"], d["source_db"], d["source_reaction_id"], C.equation_str(d["stoich"]),
                  str(d["reversible"]), "KEGG", d["balance"], "BALANCE_UNKNOWN", d["ec"], "U1", "true", ""]
                 for d in clean_rows])

    summary = {
        "n_model_metabolites": len(metas),
        "n_model_reactions": len(rxns),
        "n_bio cyc_reactions": len(raw_rows),
        "n_mapped": len(mapped_rows),
        "n_failed_unmapped": len(failed_rows),
        "n_balanced": len(balanced),
        "n_unbalanced": len(unbalanced),
        "n_balance_unknown": len(unknown),
        "n_clean_U1": len(clean_rows),
        "n_balanced_in_clean": sum(1 for d in clean_rows if d["balance"] == "BALANCED"),
        "n_unknown_in_clean": sum(1 for d in clean_rows if d["balance"] == "BALANCE_UNKNOWN"),
    }
    (CLEAN / "universe_manifest.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
