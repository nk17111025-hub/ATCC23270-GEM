# -*- coding: utf-8 -*-
"""Universe builder: BioCyc direction parsing, metabolite mapping, balance, duplicates."""
from __future__ import annotations

import re
from collections import defaultdict

import p5b6_common as C
import phase4a_lib as lib
import fw_common as FW

BIOCYC_DIR = C.ROOT / "C2_外部数据库与注释来源核验" / "C2-1_BioCyc与数据库版本核验" / "标准化数据"
RXN_FILE = BIOCYC_DIR / "reactions.tsv"
CMP_FILE = BIOCYC_DIR / "compounds_network.tsv"

DIRECTION_MAP = {
    "LEFT-TO-RIGHT": ("fwd", False),
    "PHYSIOL-LEFT-TO-RIGHT": ("fwd", False),
    "IRREVERSIBLE-LEFT-TO-RIGHT": ("fwd", False),
    "RIGHT-TO-LEFT": ("rev", False),
    "PHYSIOL-RIGHT-TO-LEFT": ("rev", False),
    "IRREVERSIBLE-RIGHT-TO-LEFT": ("rev", False),
    "REVERSIBLE": ("fwd", True),
    "BIDIRECTIONAL": ("fwd", True),
}

GENERIC_PATTERNS = [r"-ACPs?$", r"^ACP$", r"Thioredoxin", r"ferredoxin", r"Ferredoxin",
                    r"Cytochromes?-C", r"Quinols?$", r"Quinones?$", r"-tRNAs?$", r"^tRNAs?$",
                    r"Protein", r"protein", r"Peptides", r"holder", r"-Proteins?$"]

MANUAL_CURRENCY = {"Pi": "pi", "PHOSPHATE": "pi", "INORGANIC-PHOSPHATE": "pi"}


def read_biocyc_directions():
    """Return list of (rid, name, direction_label, left, right, ec)."""
    _, rows = FW.read_tsv(RXN_FILE)
    out = []
    for r in rows:
        if len(r) < 7:
            continue
        out.append((r[2].strip(), r[3].strip(), r[4].strip(),
                    r[5].split(";"), r[6].split(";"), r[10].strip() if len(r) > 10 else ""))
    return out


def enumerate_directions():
    from collections import Counter
    c = Counter()
    for rid, name, d, l, r, ec in read_biocyc_directions():
        c[d] += 1
    return c


def normalize_direction(direction_label):
    """Return (mode, reversible) or None for unknown. mode in {'fwd','rev'}."""
    return DIRECTION_MAP.get(direction_label.strip())


def parse_reaction_side(side, resolve):
    """side = list of frame ids. resolve(frame)->list of model metabolites.
    Returns (stoich_updates, unmapped, ambiguous)."""
    stoich = {}
    unmapped = []
    ambiguous = []
    for frame in side:
        frame = frame.strip()
        if not frame:
            continue
        cands = resolve(frame)
        if not cands:
            unmapped.append(frame)
            continue
        if len(cands) > 1:
            ambiguous.append((frame, cands))
            continue
        stoich[cands[0]] = stoich.get(cands[0], 0.0) + 1.0
    return stoich, unmapped, ambiguous


def build_metabolite_mapper(metas):
    """frame -> list of model metabolites via KEGG C-number (cytosolic preferred)."""
    _, cmp_rows = FW.read_tsv(CMP_FILE)
    frame2kegg = {}
    for r in cmp_rows:
        frame = r[2].strip()
        ext = r[7] if len(r) > 7 else ""
        for part in ext.split(";"):
            if part.startswith("LIGAND-CPD:"):
                frame2kegg[frame] = part.split(":", 1)[1].strip()
    t2 = lib._read_table2()
    kegg2met = defaultdict(list)
    for abbr, meta in t2.items():
        k = str(meta.get("kegg", "")).strip()
        if k and abbr in metas:
            kegg2met[k].append(abbr)

    def resolve(frame):
        if frame in MANUAL_CURRENCY:
            base = MANUAL_CURRENCY[frame]
            return [base + "[c]"] if (base + "[c]") in metas else []
        k = frame2kegg.get(frame)
        if not k or k not in kegg2met:
            return []
        cands = kegg2met[k]
        cyto = [m for m in cands if m.endswith("[c]")]
        return cyto if cyto else cands
    return resolve, frame2kegg, kegg2met


def canonical_signature(stoich):
    """Normalized stoichiometric signature for duplicate detection (forward-normalized)."""
    return tuple(sorted((m, round(c, 6)) for m, c in stoich.items()))


def is_equivalent(a, b):
    """True if a and b are same-direction or scalar-multiple equivalent."""
    if canonical_signature(a) == canonical_signature(b):
        return True
    return False


def is_reverse_equivalent(a, b):
    """True if b is exact reverse of a."""
    return canonical_signature({m: -c for m, c in b.items()}) == canonical_signature(a)


def classify_reaction(stoich):
    """Return (class, note). class in intra/transport/exchange/sink/generic/shortcut."""
    comps = set(m[m.index("[") + 1 : -1] for m in stoich)
    if len(comps) == 1:
        return "intracellular_metabolic", ""
    if len(comps) > 1:
        # detect cross-compartment shortcut (external <-> cytosolic without transport)
        if "e" in comps and "c" in comps:
            return "shortcut", "cross-compartment without transporter"
        return "transport", ""
    return "other", ""
