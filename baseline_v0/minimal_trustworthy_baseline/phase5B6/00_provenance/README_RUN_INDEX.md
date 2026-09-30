# Phase 5B-6 — Run Index

Project: *Acidithiobacillus ferrooxidans* ATCC 23270 (GCF_049532655.1, RU820_RS locus tags)

Authoritative executable model: `minimal_trustworthy_baseline_v3` — now `V3_UPDATED_NATIVE_REDOX`
(Phase 5B-6R + native-redox update, built from v2).

v3 is the **sole executable baseline** for all subsequent Phase 5B-6 / 5B-7 work.
v2, Phase 4.1G and v1 remain unchanged and are retained only for traceability.

## v3 native-redox update (Phase 5B-6R.2) — current state

Two native reactions were added directly to the authoritative v3 (not a side candidate):

| reaction | gene (old -> current) | EC | equation | class |
|---|---|---|---|---|
| NDH2_NATIVE | AFE_1854 -> RU820_RS08555 (WP_012536802.1) | 1.6.5.9 | nadh[c] + h[c] + q8[c] -> nad[c] + q8h2[c] | NATIVE_MODEL_COMPLETION |
| NOX_NATIVE | AFE_1803 -> RU820_RS08315 (WP_009568931.1) | 1.6.3.4 | nadh[c] + h[c] + 0.5 o2[c] -> nad[c] + h2o[c] | NATIVE_MODEL_COMPLETION_PROVISIONAL |

Both are mass/charge balanced, non-proton-pumping, non-PMF-generating, distinct from
Complex I (NADHI). Updated v3 SHA256:

- T0 `5ba88e6dd31d375ea0acbd94442e190b0de9f67d2e61169f8b73380d57becb77`
- TH `4c5c5474c43d215bea6d4b24f999f1f21eab3a75dea10a997e4ffeaf47fbe1b2`

Pre-edit snapshots are in `00_provenance/v3/archive_pre_native_redox/`.
Artifact QC passes (free ATP/NADH/NADPH = 0; WT biomass unchanged; no PMF/quinone/ATP-synthase
loop). Native R0 (Rubisco=0, zero heterologous additions) remains infeasible.

## v3 promotion (Phase 5B-6R) — first

v3 removes the Table-1 bound-precedence defect: historical Table-1 internal-reaction bounds
no longer overwrite the approved direction corrections (BDGK / G6PDH2 / GAPD1 / GAPD2 / PGK /
NADHI). The v3 execution rule is:

`load v3 -> apply environmental condition (exchange-only) -> apply experiment/search override`

Details: `00_provenance/v3/V3_BUILD_REPORT.md`, runtime `scripts/v3_runtime.py`.
Any attempt to apply a legacy internal bound over v3 raises `BOUND_PRECEDENCE_VIOLATION`.

| stage | status | main input | main output | checkpoint | verdict |
|---|---|---|---|---|---|
| 01 inputs / model verify | done | v2 T0/TH XML + provenance | model_provenance.tsv, source_file_hashes.tsv | cp_inputs | v2 SHA256 match |
| 02 universe raw | done | BioCyc reactions.tsv | U1/biocyc_reactions_raw.tsv | cp_raw | 1313 reactions |
| 03 mapping | done | compounds_network.tsv | metabolite_mapping.tsv, reaction_mapping.tsv, failed_reactions.tsv | cp_map | 138 mapped |
| 04 universe clean | done | mapped reactions | U1_clean.tsv, universe_manifest.json | cp_clean | 116 balanced + 27 curated |
| 05 R0 search | done | cleaned universe | phase5B6_search_summary.tsv | cp_milp | INFEASIBLE (no rescue) |
| 06 artifact QC / flux / fva / precursor / cofactor | done (partial) | — | precursor_rescue, cofactor_balance, final_summary | cp_qc | redox bottleneck localised |
| 07 Rubisco restore | skipped | — | — | — | no bypass to restore |
| 08 KO search | skipped | — | — | — | no bypass exists |
| 09 OE search | skipped | — | — | — | no bypass exists |
| 10 robustness | skipped | — | — | — | no accepted solution |
| 11 engineered models | skipped | — | — | — | none |
| 12 pareto | skipped | — | — | — | none |
| 13 final | done | all above | PHASE5B6_FINAL_REPORT.md, final_summary.tsv | cp_final | PARTIAL |

## Key result (one line)

Native v3 R0 (Rubisco=0, zero heterologous additions) is **still infeasible** (biomass 0).
The corrected EMP core is flux-capable, but the amino-acid/cofactor biosynthesis and the
NADH/NADPH redox imbalance diagnosed in Phase 5B-6 remain the bottleneck. All subsequent
heterologous search must start from v3, not v2.

## Where things live

- v3 build report: `00_provenance/v3/V3_BUILD_REPORT.md`
- v3 runtime: `scripts/v3_runtime.py`
- v3 XMLs + provenance + effective bounds: `00_provenance/v3/`
- Final report (Phase 5B-6): `13_final/PHASE5B6_FINAL_REPORT.md`
- Search summary: `05_R0_search/phase5B6_search_summary.tsv`
- Cleaned universe: `04_universe_clean/U1_clean.tsv`, `04_universe_clean/universe_manifest.json`
- Precursor rescue: `06_solution_validation/precursor_rescue/phase5B6_precursor_rescue.tsv`
- Scripts: `scripts/`
- Recovery log: `logs/phase5B6_autonomous_recovery_log.tsv`

## Model repair (do not lose this)

The v2 direction corrections (NADHI/GAPD1/GAPD2/PGK/G6PDH2 opened) were silently overridden
by Table-1 FIM condition bounds. This is now permanently fixed in v3 (execution
architecture correction), so no downstream override hack is needed.

## Legacy reference-only machinery

`phase4a_lib.py` `COND_COL` / `run_fba(with_table1=True)` and `p5b6_common.setup` apply full
Table-1 internal+exchange bounds. These are `legacy_reference_only`. New scripts must use
`v3_runtime.setup_v3` / `v3_runtime.pfba_v3`.
