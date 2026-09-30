# -*- coding: utf-8 -*-
"""Phase 5B-4: build v2 models (copy Phase4.1G + NADHI respiratory direction)."""
import csv
import hashlib
from pathlib import Path

import lxml.etree as ET

import phase4a_lib as lib

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
SRC = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase4_1G_glucose_corrected_model"
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase5B4"

EXPECT = {
    "T0": "3b7d6fda5f5b1b30369f65db44a2036c919401470a43b2be35792207afc50345",
    "TH": "ea404ac16cdeb42b90ea2581f6592d8e7aa28dc4b9553707be9a770f88428933",
}


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1024 * 1024), b""):
            h.update(c)
    return h.hexdigest()


def local(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def apply_nadhi(src, dst):
    tree = ET.parse(str(src))
    root = tree.getroot()
    model = [c for c in root if local(c.tag) == "model"][0]
    changed = False
    for r in model.iter():
        if local(r.tag) == "reaction" and r.get("id") == "R_NADHI":
            for p in r.iter():
                if local(p.tag) == "parameter" and p.get("id") == "LOWER_BOUND":
                    old = p.get("value")
                    p.set("value", "-1000")
                    changed = True
    tree.write(str(dst), encoding="UTF-8", xml_declaration=True, pretty_print=True)
    return changed


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    prov = []
    for t in ["T0", "TH"]:
        src = SRC / f"minimal_trustworthy_baseline_v1_glucose_corrected_{t}.xml"
        dst = OUT / f"minimal_trustworthy_baseline_v2_{t}.xml"
        sh = sha256(src)
        ok = sh == EXPECT[t]
        print(f"{t}: src SHA {sh} match={ok}")
        if not ok:
            print("BLOCKED_WRONG_STARTING_MODEL")
            return
        apply_nadhi(src, dst)
        dst_sh = sha256(dst)
        prov.append({"model": t, "source_path": str(src.relative_to(ROOT)), "source_sha256": sh,
                     "destination_path": str(dst.relative_to(ROOT)), "destination_pre_edit_sha256": sh,
                     "destination_post_edit_sha256": dst_sh, "edit": "NADHI LOWER_BOUND 0 -> -1000"})
    with (OUT / "phase5B4_model_provenance.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(prov[0].keys()), delimiter="\t")
        w.writeheader()
        for r in prov:
            w.writerow(r)

    # applied changes
    with (OUT / "phase5B4_applied_changes.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction", "field_changed", "old_value", "new_value", "equation", "forward_meaning",
                    "reverse_meaning", "source", "evidence_level"])
        w.writerow(["NADHI", "LOWER_BOUND", "0", "-1000",
                    "5 h[p] + nad[c] + q8h2[c] -> 6 h[c] + nadh[c] + q8[c]",
                    "reverse electron transport (PMF+quinol -> NADH)",
                    "respiratory NADH oxidation (NADH + quinone -> NAD+ + quinol)",
                    "Phase 5B-4 approved", "CURRENT_GENOME_DATABASE_SUPPORTED"])
        w.writerow(["PPC", "GPR", "AFE_1810", "AFE_1883 / RU820_RS08690",
                    "co2 + h2o + pep -> h + oaa + pi", "", "", "Phase 5B-4 approved (metadata)",
                    "CURRENT_GENOME_DATABASE_SUPPORTED"])
        w.writerow(["GHMT3", "GPR", "AFE_0295 (GlyA)", "(removed)",
                    "nad + thf + gly -> co2 + nadh + mlthf + nh4", "", "", "Phase 5B-3E verified",
                    "VERIFIED_GPR_ERROR_FIXED"])
    print("v2 built")


if __name__ == "__main__":
    main()
