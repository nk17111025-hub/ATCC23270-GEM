# PHASE 1.6 FULL CONDITION-BOUND RECONCILIATION REPORT

生成时间：2026-09-25。冻结 SBML 未改动（SHA256 与 Phase 0 一致）。场景定义已更新为 schema 3.0。

## 结论（六个问题）

1. 每个公开条件有多少反应 bound 与原始 SBML 不同？每个条件 19 个（FIM/TTM/TSM 各 19）。
2. 之前漏了哪些 bound？漏了 6 个内部反应 + 1 个特殊反应（见下）。Phase 1 只套用了 exchange + CYT2。
3. Htpp 是否三条件下都固定为 0？是。Table 1 三条件下 Htpp=0/0，FVA 确认 [0,0]。
4. ATP EGC 在实际 2016 公开模拟条件下是否可达？不可达。Htpp=0/0 后，关闭全部 exchange、最大化 ATPM 的目标为 0。
5. 哪些 Phase 1/1.5 结果仍有效、哪些被取代？见下。
6. 完整条件是否足以开始严格 Phase 2？是，但需把 Afe_biomass 的 Table 1 0/0 视为占位符（保留 SBML 0/1000），否则生长为 0。

## 1. 每条件 19 个差异

`phase1_6_bound_reconciliation.tsv` 完整列出。差异构成（每条件）：

- 13 个 exchange 反应（Ex_*）的条件化 bound（已在 Phase 1 套用）。
- CYT2（TTM/TSM 0/0；FIM 与 XML 一致）。
- 内部反应：Htpp、SQRED1、NH4pp、PPTtpp（三条件均 0/0）；SO41tpp、SO42tpp（FIM 0/0，TTM/TSM 0/1000）。
- Afe_biomass_mc507_WT_139p0M（三条件 0/0，但原始 SBML 0/1000）。

## 2. 之前遗漏的 bound

Phase 1/1.5 只套用了 exchange + CYT2，遗漏：

- FIM：Htpp、NH4pp、PPTtpp、SO41tpp、SO42tpp、SQRED1、Afe_biomass（7 个）。
- TTM：Htpp、NH4pp、PPTtpp、SQRED1、Afe_biomass（5 个，CYT2 已套用）。
- TSM：同 TTM（5 个）。

其中 Htpp 是导致 EGC 被错误保留在 Phase 1.5 场景里的关键遗漏。

## 3. Htpp 固定为 0 确认

完整 Table 1 bound 下，Htpp FVA = [0,0]（三条件）。Htpp 在公开最小培养基中确实被禁用。

## 4. ATP EGC 可达性

完整条件（Htpp=0/0）下，关闭全部 exchange、最大化 ATPM：ATPM 最大 = 0（FIM/TTM/TSM 均 0）。EGC 不可达。EGC 只存在于原始 SBML 默认 bound（Htpp 可逆 [-1000,1000]），是冻结文件的默认态缺陷，不是公开模拟条件。

## 5. 结果有效性判定

仍有效：

- 结构/计数、质量/电荷、GPR、dead-end、数值可行性、max|S·v|（Phase 1）。
- 生长/biomass（Phase 1 与 1.6 均为 0.052076，Htpp 不影响生长）。
- Phase 1.5 第 1、2 节（最小 EGC 环 + Htpp 方向）：作为“原始 SBML 默认态缺陷”仍成立。

被取代（场景相关，因 Htpp 未按 0/0 套用）：

- Phase 1.5 第 3 节基线生长里的 Htpp 通量（约 120–155）与 ATPS5rpp FVA 上限（388/953）——完整条件下降为 Htpp=0、ATPS5rpp 上限 188/753。
- Phase 1.5 第 5 节 Table 5 预检的电子供体/O2 摄取（Htpp 开放时 Fe2 ≈ -19，完整条件为 -303.9，后者接近实验 296.67）。
- Phase 1 blocked 计数：FIM 70→71、TTM 66→70、TSM 66→70（新增固定 0 的反应进入 blocked）。

## 6. 完整条件是否足以开始 Phase 2

基本可以，但有一个必须明确的决策点：Afe_biomass 的 Table 1 0/0 是占位符（biomass objective function 不是培养基组分）。按字面套用 0/0 会使三条件生长全为 0（已记录在 phase1_6_literal_full_bounds.tsv）。因此 Phase 2 复现应：

- 套用全部 Table 1 条件 bound（含 Htpp=0/0、SQRED1=0/0 等），
- 但 Afe_biomass 保留 SBML 默认 0/1000，
- 并把这一例外作为“已文档化的 2016 复现场景定义”记录。

该例外不改变冻结 SBML，也不属于修复。

## 产物

- `baseline_v0/qc/phase1_6_bound_reconciliation.tsv`
- `baseline_v0/qc/phase1_6_rediagnosis.tsv`
- `baseline_v0/qc/phase1_6_literal_full_bounds.tsv`
- `baseline_v0/qc/phase1_6_table5_check.tsv`
- `baseline_v0/qc/phase1_6_summary.json`
- `baseline_v0/phase1_6_reconcile.py`
- 场景定义更新：`baseline_v0/scenario_definition/2016_scenario_constraints.json` → schema 3.0

未修复模型，未改冻结 SBML，未 gap-fill，未开始工程优化。
