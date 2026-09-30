# PHASE 2.1 PUBLISHED ARTIFACT CONSISTENCY AUDIT REPORT

生成时间：2026-09-25。冻结 SBML 未改动。仅做 2016 官方 artifacts 的比对与诊断，不修复、不重拟合、不 gap-fill。

## 结论（七个问题）

1. mmc1.xls 与 mmc3.xml 之间有 1 个真实反应式不匹配。
2. CYTRED 的质子方向确实相反。
3. CYTRED 方向不解释 TSM 供体摄取偏差。
4. 两个 CYTRED 表示对已发表验证行为没有影响（CYTRED 在验证条件下通量为 0）。
5. 其他 ETC/硫通路无 stoichiometric 不一致；TSM 偏差是硫氧化通路的真实模型差距。
6. 四个数值弱点在 GLPK-exact 下全部满足 max|S·v| ≤ 1e-8。
7. 最终分类：PARTIAL_REPRODUCTION_WITH_DOCUMENTED_CAUSES。

## 1. 反应式不匹配数量

对全部 615 反应做规范化比较（代谢物集、系数、方向、质子区室），仅 1 个真实差异：

- CYTRED：PROTON_DIRECTION_MISMATCH。

其余 614 反应均为 EXACT_MATCH。没有 DIRECTION_REPRESENTATION_ONLY、COMPARTMENT_MISMATCH、STOICHIOMETRIC_MISMATCH 或 OTHER。

（CYTBO3/CYTBD 的反应名写“2.6 protons”但实际化学计量为 1.8，这是名称标签问题，不是反应式差异；两个 artifact 的化学计量都是 1.8。）

## 2. CYTRED 质子方向确认

mmc3.xml：`0.5 h[c] + q8h2[c] + etpcycAox[p] -> 0.5 h[p] + q8[c] + etpcycArd[p]`（h[c] 在反应物侧）。

mmc1.xls Table 1：`0.5 h[p] + q8h2[c] + etpcycAox[p] -> 0.5 h[c] + q8[c] + etpcycArd[p]`（h[p] 在反应物侧）。

h[c] 与 h[p] 的 0.5 系数被互换，其余（q8h2/q8、etpcycAox/etpcycArd）一致。这是真实的方向相反，不是标签问题。

## 3. CYTRED 方向是否解释 TSM 供体偏差

否。三个临时变体（SBML_ORIGINAL / TABLE1_CYTRED / CYTRED_NO_PROTON）复现 20 个 Table 5 场景，七个验证曲线完全一致：

- growth_FIM slope 1.187 R² 0.9184
- growth_TTM slope 0.937 R² 0.9688
- growth_TSM slope 0.937 R² 0.9804
- donor_FIM_fe2 slope 1.033 R² 0.9842
- donor_TTM slope 1.047 R² 0.8767
- donor_TSM slope 1.301 R² 0.7738
- o2_FIM slope 1.124 R² 0.9554

均值 R²(1-SSE/SST) 三变体均为 0.92252。TSM 供体摄取三个变体都是 slope 1.301、R² 0.774、6/6 实验值落在最优 FVA 区间外。

原因：在公开 TTM/TSM 条件下，硫氧化电子流走 CYTBD 分支（ubiquinol 氧化酶），不走 CYTRED→CYTAA32 分支。实测通量：CYTRED=0、CYTAA32=0、CYTBO3=0、CYTBD=18.07、TSQOC=20.54。因此 CYTRED 的 0.5 H+ 方向（甚至去掉质子）对验证结果零影响。

## 4. 哪个表示更符合论文行为

两者对已发表验证行为无差别（CYTRED 惰性）。CYTRED 的质子方向是真实的 artifact 不一致，需在 2026 重建时裁定，但不是 TSM 偏差的来源。

## 5. 其他 ETC/硫通路不一致

比较结果显示除 CYTRED 外无任何 stoichiometric 不一致。TSM 硫代硫酸盐摄取系统性偏高约 30%（slope 1.301、R² 0.774、6/6 在 FVA 外）不是 artifact 不一致，而是硫氧化通路的真实模型差距，方向指向 TSQOC / 4THASE / CYTBD 分支的化学计量或能量耦合，需后续证据审查，不在本阶段处理。

## 6. GLPK-exact 数值闭合

四个弱点在 glpk_exact 下：

- 点 2（TTM）：max|S·v| = 1.78e-14
- 点 3（TSM）：max|S·v| = 3.55e-14
- 点 19（FIM）：max|S·v| = 1.42e-14
- 点 20（FIM）：max|S·v| = 7.11e-15

全部 ≤ 1e-8，达标。之前 1e-7 级残差是 GLPK 浮点精度所致。

## 7. 最终分类

`PARTIAL_REPRODUCTION_WITH_DOCUMENTED_CAUSES`

理由：20/20 可行、平均 R²≈0.92 复现、能量参数一致、EGC 失活、GLPK-exact 数值闭合；但 TSM 硫代硫酸盐供体摄取系统性偏高约 30%，既不是 solver 效应，也不是 alternate-optimum 效应，也非 CYTRED 方向，而是硫氧化通路（TSQOC/4THASE/CYTBD）的真实模型差距，已记录原因、待后续证据审查。

## 产物

- `baseline_v0/reproduction_2016/PHASE2_1_ARTIFACT_CONSISTENCY_AUDIT.md`
- `baseline_v0/reproduction_2016/mmc1_vs_mmc3_reaction_equation_diff.tsv`
- `baseline_v0/reproduction_2016/cytred_variant_20point_reproduction.tsv`
- `baseline_v0/reproduction_2016/cytred_variant_seven_metrics.tsv`
- `baseline_v0/reproduction_2016/cytred_variant_TTM_TSM_fva.tsv`
- `baseline_v0/reproduction_2016/phase2_exact_four_problem_points.tsv`
- `baseline_v0/reproduction_2016/phase2_1_summary.json`
- `baseline_v0/phase2_1_audit.py`

未修改冻结模型，未重拟合，未 gap-fill，未创建校正模型，未开始工程优化。
