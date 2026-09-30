# PHASE 3.1 HCO3E BYPASS AUDIT REPORT

生成时间：2026-09-25。诊断专用，未修复模型、未 gap-fill、未改工程。全部 pFBA 作为主通量追踪解。

## 结论（六个问题）

1. 绕过 HCO3E 的确切反应路线：hco3 → ACCOAC（乙酰辅酶 A 羧化酶）→ 丙二酰辅酶 A → 脂肪酸延长脱羧（3OAS 系列 + MACPD）→ 再生 CO2 → RUBISCO。
2. 绕过 HCO3E 会增加供体/O2/ATP 成本：是。
3. 能量受限下 HCO3E 缺失会降低生长：是（FIM cap_donor / cap_both 直接 infeasible）。
4. 零碳供体/O2 消耗完全由 NGAM 解释：是。
5. 修正后的泄漏分类无意外碳摄取或能量源：是。
6. Phase 3 可正式关闭，进入 Phase 4：是（HCO3E 冗余记为动力学盲点，不是失败）。

## 1. 确切绕过路线

正常：h2co3 → hco3（HCO3tpp）→ co2（HCO3E，-1.90）→ RUBISCO 固定。
HCO3E=0：hco3 → ACCOAC（乙酰辅酶 A 羧化酶，0.068 → 1.970，增 +1.90）→ malcoa → 脂肪酸合成脱羧（3OAS/MACPD 产生 co2）→ RUBISCO 固定（2.166，不变）。

CBPS（氨甲酰磷酸合酶）0.030 不变；BIOC2（生物素羧化酶）为 0。绕过点是 ACCOAC：它直接消耗 hco3 生成丙二酰辅酶 A，再由脂肪酸链延长反应的脱羧步骤再生 CO2 供 RUBISCO。这是化学计量冗余，不是碳酸酐酶在生物学上可省。

## 2. 绕过 HCO3E 的能量成本（每单位 biomass）

| 条件 | 供体/biomass 正常→关闭 | O2/biomass 正常→关闭 | ATPS5rpp 正常→关闭 | 总绝对通量 正常→关闭 |
|---|---:|---:|---:|---:|
| FIM | 3159→3342 | 742→788 | 20.75→22.65 | 1320→1401 |
| TTM | 100.2→107.0 | 303→327 | 20.74→22.64 | 380→409 |
| TSM | 175.3→187.3 | 303→327 | 20.74→22.64 | 342→369 |

供体 +6~7%、O2 +6~8%、ATPS5rpp +1.90（对应 ACCOAC 多耗的 1 ATP）、总通量 +6~8%。

## 3. 能量受限下的生长

| 条件 | cap_donor | cap_o2 | cap_both |
|---|---:|---:|---:|
| FIM | infeasible | 0.0422 | infeasible |
| TTM | 0.0423 | 0.0425 | 0.0423 |
| TSM | 0.0423 | 0.0425 | 0.0423 |

基线生长 0.052。能量受限后 HCO3E 缺失使生长降到 ~0.042（约 -19%），FIM 在 cap_donor/cap_both 下 infeasible。证明 HCO3E 表面非必需只是外部能量未受限造成的。

## 4. 零碳维持代谢

| 条件 | ATPM=3.475 供体 | ATPM=0 供体 | ATPM=3.475 ATPS5rpp | ATPM=0 ATPS5rpp |
|---|---:|---:|---:|---:|
| FIM | -17.375 | 0.0 | 3.475 | 0.0 |
| TTM | -0.653 | 0.0 | 3.475 | 0.0 |
| TSM | -1.143 | 0.0 | 3.475 | 0.0 |

ATPM=0 后所有供体/O2/ATPS5rpp 通量归零。零碳下供体/O2 消耗完全由 NGAM（ATPM=3.475）解释，归类 `NGAM_MAINTENANCE_ENERGY`，不是泄漏。

## 5. 修正后的泄漏分类

- 无意外碳摄取（UNEXPECTED_CARBON_UPTAKE 为空）。
- TTM/TSM 的 Fe2/Fe3 摄取 = biomass 矿物需求（biomass 系数 fe2=0.067、fe3=0.067，0.067×0.052=0.00349，与观测吻合），标 `BIOMASS_MINERAL_REQUIREMENT`。
- FIM Fe2 摄取 = INTENDED_DONOR；TTM ttton / TSM tsul = INTENDED_DONOR；O2 = INTENDED_ACCEPTOR；H2CO3 = INTENDED_CARBON。
- H+ = ACID_BASE_BALANCE。
- 4hba 输出 = BIOSYNTHETIC_BYPRODUCT_EXPORT（源自 TTMS 硫胺素合成副产物，非碳泄漏）。

## 6. 关闭判定

Phase 3 核心生理通过；HCO3E 冗余已确认为化学计量冗余/动力学盲点，且能量受限下会降低生长。无隐藏碳源、无隐藏能量源、无自由 ATP。Phase 3 可正式关闭，Phase 4 可开始。

## 产物

- `baseline_v0/qc/physiology/PHASE3_1_HCO3E_BYPASS_AUDIT.md`
- `hco3e_flux_trace.tsv` / `hco3e_flux_delta.tsv` / `hco3e_energy_cost.tsv`
- `zero_carbon_maintenance.tsv` / `exchange_leakage_inventory_corrected.tsv`
- `phase3_1_summary.json`
- `baseline_v0/phase3_1_hco3e_audit.py`

未修复模型、未 gap-fill、未开始工程优化。
