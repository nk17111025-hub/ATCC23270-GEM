# 00_READ_FIRST — before final 入口说明

## 1. Current purpose

`before final/` is the working entry for the current reconstruction stage.

The goal of this stage is:

> integrate the completed G51–G510 evidence audit with the existing v3 model to construct `updatedv3.1`.

This is an evidence-driven update of v3, not a new reconstruction from scratch.

## 2. Source locations

```text
PARENT_V3=baseline_v0/minimal_trustworthy_baseline/phase5B6/00_provenance/v3
G5X_EVIDENCE_ROOT=G5X_2929基因逐基因证据链
WORKSPACE=before final
```

- `PARENT_V3` is the existing parent model/runtime state.
- G51–G510 are the evidence source for the current genome `GCF_049532655.1`.
- current gene coordinates are `RU820_RS...` / `WP_...`.
- old `AFE_xxxx` identifiers are historical mappings for the old model/literature.

## 3. Integration workflow

```text
completed G51-G510 evidence
        ↓
map each model-relevant finding to the corresponding v3 object
        ↓
compare evidence-supported state with current v3 state
        ↓
adjudicate the required action
        ↓
build the complete updatedv3.1 state
        ↓
generate and validate updatedv3.1
```

The next stage compares evidence against v3; it does not blindly import G5X conclusions.

## 4. Adjudication outcomes

- `APPLY` — a supported model change is absent from v3.
- `ALREADY_APPLIED` — the supported correction is already present in v3.
- `ANNOTATION_ONLY` — metadata/evidence/identifier update without changing flux structure.
- `CONFIDENCE_ONLY` — reaction confidence changes without changing reaction structure.
- `HOLD` — evidence is relevant but insufficient to define the exact model change.
- `NO_CHANGE` — retain the current v3 representation.

One key rule:

> Change only the model property actually supported by the evidence.

Examples:

- GPR evidence → GPR change.
- intrinsic direction evidence → base reaction direction/bound change.
- stoichiometric evidence → reaction equation change.
- fully resolved missing native reaction → reaction addition.
- evidence-strength change only → confidence/annotation update only.

Culture-condition-specific constraints are not part of this reconstruction decision layer; they remain reproduction/QC evidence.

## 5. Working products

`01_G5X_to_v3_改动对照集合.xlsx`

- detailed G5X ↔ v3 comparison, adjudication, exact proposed change, and provenance.

`02_updatedv3.1_master.xlsx`

- complete current state of updatedv3.1, not only changed reactions.
- this will be the human-readable authoritative model table.

`03_updatedv3.1.xml`

- executable updatedv3.1 model generated from the accepted master state.

These files are not created in Task 1.

## 6. Traceability

```text
updatedv3.1 object
→ patch/change ID
→ G5X-to-v3 adjudication row
→ exact G5X source record
→ original database/literature evidence
```

The original G5X evidence remains in its existing location rather than being duplicated into `before final/`.

## 7. Current status

```text
PROJECT=ATCC23270_GEM_RECONSTRUCTION
STAGE=BEFORE_FINAL_INITIALIZED
PARENT_V3=baseline_v0/minimal_trustworthy_baseline/phase5B6/00_provenance/v3
G5X_EVIDENCE_ROOT=G5X_2929基因逐基因证据链
TARGET_MODEL=updatedv3.1
UPDATEDV3_1_STATUS=NOT_BUILT
NEXT_TASK=BUILD_G5X_TO_V3_ADJUDICATION_TABLE
```
All paths in this document are relative to the ATCC23270 project root.