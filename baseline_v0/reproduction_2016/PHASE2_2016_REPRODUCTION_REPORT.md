# PHASE 2 STRICT 2016 REPRODUCTION REPORT

生成时间：2026-09-25。冻结 SBML 未改动。所有模拟使用 Phase 1.6 确立的完整 Table 1 条件 bound（published_operational_bounds，含 Htpp=0/0，Afe_biomass 保留 SBML 0/1000）。目标为最大化 Ex_bio/biomass，GAM=139、NGAM/ATPM=3.475。不约束实验生长率、供体摄取、O2 摄取（这些是验证目标）。

## 最终分类

`PASS_WITH_SOLVER_OR_ALTERNATE_OPTIMUM_DIFFERENCES`

论文平均 R² ≈ 0.92 被复现（本研究 1-SSE/SST 均值 = 0.9225）。差异来源：4 个低生长点的 GLPK 浮点残差略超 1e-8，以及 TSM 硫代硫酸盐摄取系统性偏高（非 alternate optimum）。

## 1. 20 个场景是否全部可行

是，20/20 均 optimal。

## 2. 实验生长复现精度

| 系列 | n | slope | Pearson r | r²(Pearson) | R²(1-SSE/SST) |
|---|---:|---:|---:|---:|---:|
| growth FIM | 9 | 1.187 | 0.9964 | 0.9928 | 0.9184 |
| growth TTM | 5 | 0.937 | 1.0000 | 1.0000 | 0.9688 |
| growth TSM | 6 | 0.937 | 1.0000 | 1.0000 | 0.9804 |

FIM 高生长点略偏高（约 +18%），低生长点略偏低；TTM/TSM 略偏低（约 -6%）。

## 3. 电子供体摄取复现精度

| 系列 | n | slope | R²(1-SSE/SST) |
|---|---:|---:|---:|
| FIM Fe2 | 9 | 1.033 | 0.9842 |
| TTM tetrathionate | 5 | 1.047 | 0.8767 |
| TSM thiosulfate | 6 | 1.301 | 0.7738 |

FIM Fe2 摄取复现非常好（slope 1.033）；TSM 硫代硫酸盐摄取系统性偏高约 30%。

## 4. FIM O2 摄取复现精度

n=9，slope=1.124，R²(1-SSE/SST)=0.9554。O2 摄取略偏高（约 +12%）。

## 5. 七个验证 R² 值与均值

| 系列 | R²(1-SSE/SST) | r²(Pearson) |
|---|---:|---:|
| growth FIM | 0.9184 | 0.9928 |
| growth TTM | 0.9688 | 1.0000 |
| growth TSM | 0.9804 | 1.0000 |
| donor FIM Fe2 | 0.9842 | 0.9928 |
| donor TTM | 0.8767 | 0.9914 |
| donor TSM | 0.7738 | 0.9945 |
| O2 FIM | 0.9554 | 0.9789 |
| 均值 | 0.9225 | 0.9929 |

## 6. 论文 ≈0.92 平均 R² 是否复现

是。本研究用经典决定系数 `1 - SSE/SST` 的均值 = 0.9225，与论文报告的 ≈0.92 一致。论文未明确给出 R² 定义，故同时报告两种口径：Pearson r² 均值 0.9929（相关性口径），决定系数均值 0.9225（绝对误差口径）。

## 7. 哪些差异来自 alternate optima

- 生长（目标函数）在最大生长下是唯一值，无 alternate optimum 问题。
- FIM Fe2 摄取：多数点实验值落在最优 FVA 区间内（EXPERIMENT_WITHIN_OPTIMAL_FVA）；点 1、13 落在区间外。
- TSM 硫代硫酸盐摄取：6/6 点均落在最优 FVA 区间外（EXPERIMENT_OUTSIDE_OPTIMAL_FVA），说明这是模型结构性的偏高，不是 alternate optimum 假象。
- O2（FIM）：多数点在区间内，点 1、4 在区间外。

## 8. EGC 是否在每个公开条件下完全失活

是。所有 20 点 Htpp 通量 = 0，Htpp FVA = [0,0]。原始 SBML 的 ATP 能量生成环在公开条件下不可达。

记录为：`KNOWN_RAW_SBML_DEFECT_NOT_ACTIVE_UNDER_PUBLISHED_CONDITIONS`。

## 9. 公开能量参数是否一致编码

全部一致（按反应式 H+/膜质子净移位核对）：

| 反应 | 公开 | 编码 | 一致 |
|---|---:|---:|---|
| ATPS5rpp | 5 | 5 | 是 |
| CYTAA31 | 0 | 0 | 是 |
| NADHI | 5 | 5 | 是 |
| CYTBC1 | 5 | 5 | 是 |
| CYTAA32 | 1.3 | 1.3 | 是 |
| CYTRED | 0.5 | 0.5 | 是 |
| CYTBO3 | 1.8 | 1.8 | 是 |
| CYTBD | 1.8 | 1.8 | 是 |
| GAM | 139 | 139.1 | 是（差 0.1，属估算取整） |
| NGAM/ATPM | 3.475 | 3.475 | 是 |

两个标签性小问题（不改变化学计量）：CYTBO3/CYTBD 的反应名写“2.6 protons”，但实际化学计量为 1.8（与公开一致）；CYTRED 在 SBML 与 Table 1 中的质子方向相反（幅度均为 0.5）。作为 ISSUE 标记，不修改。

## 10. 冻结基线是否足够可信以进入生理 sanity 测试与后续 2026 重建

基本可信。理由：20/20 可行；生长与 FIM 供体/O2 复现良好；平均决定系数 0.9225 复现论文 ≈0.92；EGC 在公开条件下失活；能量参数一致。

需携带的已知差异：

- max|S·v| = 4.87e-7（4 个低生长点略超 1e-8，GLPK 浮点精度；GLPK-exact 降到 ~1e-15）。16 个点在 ~1e-14。
- TSM 硫代硫酸盐摄取系统性偏高（slope 1.301，R² 0.774，6/6 在 FVA 外），指向 S 氧化/TSQOC/四硫代硫酸盐循环的模型差距，需后续证据审查，不在本阶段修复。
- GAM 编码 139.1（名义 139）、CYTBO3/CYTBD 名称“2.6”与化学计量 1.8 不符、CYTRED 方向歧义，作为 ISSUE 记录。

## 产物

- `baseline_v0/reproduction_2016/PHASE2_2016_REPRODUCTION_REPORT.md`
- `baseline_v0/reproduction_2016/table5_20point_reproduction.tsv`
- `baseline_v0/reproduction_2016/table5_optimal_fva.tsv`
- `baseline_v0/reproduction_2016/table5_pfba_diagnostic.tsv`
- `baseline_v0/reproduction_2016/table5_target_classification.tsv`
- `baseline_v0/reproduction_2016/seven_validation_series.tsv`
- `baseline_v0/reproduction_2016/seven_validation_metrics.tsv`
- `baseline_v0/reproduction_2016/published_energy_parameter_check.tsv`
- `baseline_v0/reproduction_2016/table5_glpk_exact_crosscheck.tsv`
- `baseline_v0/reproduction_2016/scenario_interpretation_notes.md`
- `baseline_v0/reproduction_2016/phase2_summary.json`
- `baseline_v0/phase2_reproduce_2016.py`
- `baseline_v0/logs/`（运行日志）

未修改模型，未重拟合 GAM/ETC，未 gap-fill，未纠正 Htpp，未开始工程优化。
