#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Apply the APPROVED Phase 3.1 HCO3E-bypass patch to a working copy."""
import hashlib
import datetime as dt
from pathlib import Path
import libsbml

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
FROZEN = ROOT / "baseline_v0" / "original" / "mmc3.xml"
WORKDIR = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "working"
FREEZE = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze"

def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    FREEZE.mkdir(parents=True, exist_ok=True)

    frozen_hash = sha256(FROZEN)
    raw = FROZEN.read_bytes()
    doc = libsbml.readSBMLFromString(raw.decode("utf-8", "replace"))
    m = doc.getModel()

    macpd = m.getReaction("R_MACPD")
    macpd.getKineticLaw().getParameter("UPPER_BOUND").setValue(0.0)
    acoata = m.getReaction("R_ACOATA")
    acoata.getKineticLaw().getParameter("LOWER_BOUND").setValue(0.0)

    working = WORKDIR / "minimal_trustworthy_baseline_working.xml"
    working.write_text(libsbml.writeSBMLToString(doc), encoding="utf-8")
    working_hash = sha256(working)

    now = dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds")
    changelog = [
        "# minimal_trustworthy_baseline changelog",
        "",
        "generated_at=" + now,
        "source_model=" + str(FROZEN.relative_to(ROOT)),
        "source_sha256=" + frozen_hash,
        "working_model=" + str(working.relative_to(ROOT)),
        "working_sha256=" + working_hash,
        "",
        "## Approved changes (Phase 3.1 HCO3E bypass patch)",
        "",
        "### 1. MACPD (R_MACPD)",
        "- action: disable",
        "- original_bounds: [0.0, 1000.0]",
        "- final_bounds: [0.0, 0.0]",
        "- original_stoichiometry: h[c] + malACP[c] -> acACP[c] + co2[c] (unchanged)",
        "- GPR: none",
        "- confidence: 1",
        "- reason: malonyl-ACP decarboxylase; EC 2.3.1.38 inconsistent; no PMID; no GCF_049532655.1 gene evidence; closes artificial HCO3E bypass",
        "- rollback: restore UPPER_BOUND=1000.0",
        "",
        "### 2. ACOATA (R_ACOATA)",
        "- action: prevent reverse flux (retain forward as conservative assumption)",
        "- original_bounds: [-1000.0, 1000.0]",
        "- final_bounds: [0.0, 1000.0]",
        "- original_stoichiometry: ACP[c] + accoa[c] <=> acACP[c] + coa[c] (unchanged)",
        "- GPR: none",
        "- confidence: 1",
        "- reason: reverse direction closes artificial HCO3E bypass; forward retained and labeled conservative modeling assumption",
        "- rollback: restore LOWER_BOUND=-1000.0",
        "",
    ]
    (WORKDIR / "MINIMAL_TRUSTWORTHY_BASELINE_CHANGELOG.md").write_text("\n".join(changelog), encoding="utf-8")

    diff = [
        "reaction_id\tchange_type\toriginal_lb\toriginal_ub\tfinal_lb\tfinal_ub\tstoichiometry_changed",
        "MACPD\tbounds\t0.0\t1000.0\t0.0\t0.0\tFalse",
        "ACOATA\tbounds\t-1000.0\t1000.0\t0.0\t1000.0\tFalse",
    ]
    (WORKDIR / "reaction_level_diff.tsv").write_text("\n".join(diff) + "\n", encoding="utf-8")

    print("patch applied")
    print("frozen_sha256", frozen_hash)
    print("working_sha256", working_hash)
    print("working", str(working.relative_to(ROOT)))

if __name__ == "__main__":
    main()
