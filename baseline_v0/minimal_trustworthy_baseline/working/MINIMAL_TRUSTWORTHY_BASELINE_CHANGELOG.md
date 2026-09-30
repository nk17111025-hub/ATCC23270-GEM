# minimal_trustworthy_baseline changelog

generated_at=2026-09-26T07:55:56+08:00
source_model=baseline_v0\original\mmc3.xml
source_sha256=0dbaacfc07aafb892db0731a517160cf46006bc502b0dedd56f950c14b840246
working_model=baseline_v0\minimal_trustworthy_baseline\working\minimal_trustworthy_baseline_working.xml
working_sha256=99c32e73f0025e086b4426ec4e60e5968659ec22f8aae829169fde09f1fc4da7

## Approved changes (Phase 3.1 HCO3E bypass patch)

### 1. MACPD (R_MACPD)
- action: disable
- original_bounds: [0.0, 1000.0]
- final_bounds: [0.0, 0.0]
- original_stoichiometry: h[c] + malACP[c] -> acACP[c] + co2[c] (unchanged)
- GPR: none
- confidence: 1
- reason: malonyl-ACP decarboxylase; EC 2.3.1.38 inconsistent; no PMID; no GCF_049532655.1 gene evidence; closes artificial HCO3E bypass
- rollback: restore UPPER_BOUND=1000.0

### 2. ACOATA (R_ACOATA)
- action: prevent reverse flux (retain forward as conservative assumption)
- original_bounds: [-1000.0, 1000.0]
- final_bounds: [0.0, 1000.0]
- original_stoichiometry: ACP[c] + accoa[c] <=> acACP[c] + coa[c] (unchanged)
- GPR: none
- confidence: 1
- reason: reverse direction closes artificial HCO3E bypass; forward retained and labeled conservative modeling assumption
- rollback: restore LOWER_BOUND=-1000.0
