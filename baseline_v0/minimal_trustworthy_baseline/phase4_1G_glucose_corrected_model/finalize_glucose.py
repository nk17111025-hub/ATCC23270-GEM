# -*- coding: utf-8 -*-
"""Write validation summary JSON + copy glucokinase resolution to evidence dir."""
import json
import datetime as dt
import shutil
from pathlib import Path

from build_glucose_model import sha256

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
OUT = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase4_1G_glucose_corrected_model"
EVID = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "phase4_1G_glucose_evidence"
FROZEN = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"


def main():
    summary = {
        "generated_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds"),
        "frozen_v1_sha256": sha256(FROZEN),
        "frozen_v1_sha256_match": sha256(FROZEN) == "bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c",
        "glucokinase_locus": {
            "legacy_locus": "AFE_2841", "old_refseq": "AFE_RS13045",
            "current_locus": "RU820_RS13140", "protein": "WP_012537424.1",
            "resolution": "GLUCOKINASE_CANDIDATE_ONLY",
            "bdgk_gpr": "GLUCOKINASE_GPR_UNRESOLVED",
        },
        "corrections": {
            "BDGK": {"lb": "0", "reason": "forward-only (block reverse ATP-generating glucose formation)"},
            "G6PDH2": {"lb": "-1000", "reason": "open oxidative PPP direction (G6PDH EC 1.1.1.49)"},
            "GAPD1": {"lb": "-1000", "reason": "open glycolytic direction (type-I GAPDH reversible)"},
            "GAPD2": {"lb": "-1000", "reason": "open glycolytic direction (AFE_3251 NADP variant)"},
            "PGK": {"lb": "-1000", "reason": "open glycolytic direction (reviewed reversible PGK)"},
            "glucose_uptake": "ADD provisional native uptake; NATIVE_GLUCOSE_UPTAKE_GPR_UNRESOLVED (T0/TH)",
        },
        "baseline_regression": {
            "FIM": 0.052076387, "TTM": 0.052076387, "TSM": 0.052076387,
            "pass": True, "free_ATP_NADH_NADPH": 0.0,
        },
        "glucose_result": {
            "biomass_increases": True,
            "FeS_oxidation_retained": True,
            "oxidative_PPP_used": False,
            "lower_EMP_glycolysis_used": False,
            "RUBISCO_dependent": True,
            "RUBISCO_KO": "infeasible",
            "detour_disappears": False,
        },
        "artifact_flags": [
            "NADTRHD large flux up to ~9.3 (confidence 1, no GPR) at high glucose",
            "glucose growth is RUBISCO/CO2-fixation dependent, not heterotrophic",
            "Fe2 uptake saturates -1000 bound at high glucose (energy-limited)",
            "no glyoxylate-shunt activation (MALS=0); no free ATP/NADH/NADPH cycle",
        ],
        "final_status": "PHASE4_1G_PASS_WITH_TRANSPORT_UNCERTAINTY",
        "candidate_model_sha": {
            "T0": sha256(OUT / "minimal_trustworthy_baseline_v1_glucose_corrected_T0.xml"),
            "TH": sha256(OUT / "minimal_trustworthy_baseline_v1_glucose_corrected_TH.xml"),
        },
    }
    (OUT / "phase4_1G_validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    shutil.copy2(OUT / "glucokinase_locus_resolution.md", EVID / "glucokinase_locus_resolution.md")
    print("written validation summary + copied glucokinase resolution")
    print("status:", summary["final_status"])
    print("T0 sha", summary["candidate_model_sha"]["T0"])
    print("TH sha", summary["candidate_model_sha"]["TH"])


if __name__ == "__main__":
    main()
