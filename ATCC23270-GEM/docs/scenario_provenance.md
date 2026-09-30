# Scenario Provenance — ATCC23270-GEM (FIM / TTM / TSM + perturbations)

This document records where the validated physiological scenario definitions and
regression logic were recovered from, so the new repository can reproduce them without
reading the historical baseline directory.

## 1. Historical source files

Priority order applied (executable > machine-readable config > frozen artifact > report):

| role | path |
|---|---|
| exchange condition bounds (authoritative) | `baseline_v0\minimal_trustworthy_baseline\phase5B6\00_provenance\v3\v3_condition_definition.tsv` |
| execution rule (exchange-only conditions) | `baseline_v0\minimal_trustworthy_baseline\phase5B6\scripts\v3_runtime.py` |
| condition API / perturbations | `baseline_v0\minimal_trustworthy_baseline\phase5B6\scripts\conditions.py` |
| build + regression QC | `baseline_v0\minimal_trustworthy_baseline\phase5B6\scripts\build_v3.py` |
| COBRApy regression reference | `baseline_v0\minimal_trustworthy_baseline\regression_qc.py` |

Relevant functions / sections:

- `v3_runtime.py::condition_environmental_bounds()` — applies Table-1 exchange-only
  condition bounds; internal reaction base bounds are authoritative and never overridden.
- `v3_runtime.py::setup_v3()` — the `load -> environmental condition -> override` order.
- `build_v3.py::wt_reference()` / `no_donor_override()` — WT (CO2 fixed -2, glucose
  closed), no-donor (all donor exchanges closed).
- `regression_qc.py::run()` — the COBRApy scenario runner (CO2 fixed at `-hco3`, donor
  open, optional donor/O2 caps).

## 2. Objective

- Objective reaction: biomass exchange `Ex_bio[e]`
  (SBML id `R_Ex_bio_LSQBKT_e_RSQBKT_`), maximized.
- This is the `OBJECTIVE_COEFFICIENT=1` reaction in the SBML.

## 3. Condition exchange bounds (FIM / TTM / TSM)

Recovered verbatim from `v3_condition_definition.tsv`. Only exchange reactions receive
condition-specific bounds; internal reactions keep their SBML base bounds.

The three conditions share most exchange bounds and differ on the electron donor and the
sulfur/hydrogen settings:

| exchange | FIM | TTM | TSM |
|---|---|---|---|
| `Ex_fe2[e]` (Fe2 donor) | -1000 / 1000 | -1000 / 1000 | -1000 / 1000 |
| `Ex_ttton[e]` (tetrathionate) | 0 / 0 | -1000 / 0 | 0 / 0 |
| `Ex_tsul[e]` (thiosulfate) | 0 / 0 | 0 / 0 | -1000 / 0 |
| `Ex_so4[e]` (sulfate) | 0 / 0 | 0 / 1000 | 0 / 1000 |
| `Ex_h2co3[e]` (inorganic carbon) | -2 / -2 | -2 / -2 | -2 / -2 |
| `Ex_o2[e]` (oxygen) | -1000 / 1000 | -1000 / 1000 | -1000 / 1000 |
| `Ex_h2[e]` (hydrogen) | 0 / 1000 | 0 / 0 | 0 / 0 |
| `Ex_eps_AFE[e]` | -1000 / 1000 | 0 / 1000 | 0 / 1000 |
| `Ex_h2s[e]`, `Ex_s[e]` | 0 / 0 | 0 / 0 | 0 / 0 |

The full 28-reaction exchange table is stored in `config/scenarios.json`
(`conditions.<cond>.exchange_bounds`).

Electron donors: FIM = Fe2, TTM = tetrathionate, TSM = thiosulfate.
Carbon source: H2CO3 (inorganic), fixed at uptake 2 mmol/gDW/h in the base condition.

## 4. Perturbations (applied on top of a condition)

| perturbation | exact override |
|---|---|
| `WT` | `Ex_h2co3[e]` = -2 / -2; glucose and donor left open (SBML default) |
| `no_donor` | close the condition donor exchange (`Ex_fe2[e]` / `Ex_ttton[e]` / `Ex_tsul[e]` = 0 / 0); CO2 = -2 / -2 |
| `no_carbon` | `Ex_h2co3[e]` = 0 / 0; glucose and donor left open |
| `glucose_closed` | `Ex_glc-B[e]` = 0 / 0 and `Ex_h2co3[e]` = -2 / -2 |

Note on labels: the historical v3 build report called the glucose-closed + CO2-fixed
state "WT" (objective 0.0520764). In the frozen 12-point regression, that state is named
`glucose_closed`, while `WT` means glucose left open (FIM 0.6735...). The recovered
definitions reproduce all 12 targets to <= 1e-13.

## 5. Reproducibility result

Run `python tests/test_regression.py`. Result against `models/current/updatedv3.1.xml`:
**12/12 PASS** (tolerance 1e-9). Output: `tests/regression_v3.1.tsv`.

## 6. Ambiguity / conflict

- The `WT` label differs between the v3 build report (glucose-closed) and the frozen
  regression table (glucose-open). This is a naming difference, not a bounds conflict;
  both states are recovered and the 12 targets match exactly.
- Solver: the historical scripts used scipy `linprog(method="highs")` (v3 runtime) and
  COBRApy pFBA (older regression). The new `scenarios.py` uses scipy `linprog("highs")`;
  objective values agree to floating-point precision (<= 1e-13).

## 7. Model integrity

The metabolic model was NOT modified. Current model SHA256:
`026352372D0B04D2CC1518B0E37D92075F7AF7F5112B94EB48D6302FC680F3A3`.
