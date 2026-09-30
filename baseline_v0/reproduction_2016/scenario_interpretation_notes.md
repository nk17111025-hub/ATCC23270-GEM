# 2016 复现场景解释说明（Phase 2 前置登记）

## A. Afe_biomass Table 1 例外

两个表示：

1. `literal_table1_bounds`：完全按 mmc1.xls Table 1 第 11–16 列字面套用（含 Afe_biomass=0/0）。三条件生长全为 0（见 baseline_v0/qc/phase1_6_literal_full_bounds.tsv）。
2. `published_operational_bounds`：套用全部已验证的 Table 1 条件 bound，但 Afe_biomass 保留 SBML 默认 0/1000，使 biomass objective 可以工作。

证据：

- mmc3.xml：R_Afe_biomass_mc507_WT_139p0M bounds = 0/1000（不可逆，正向产 bio[c]）。
- mmc1.xls Table 1：Afe_biomass 三条件 lb/ub = 0/0。
- 论文 Methods 2.3：不可逆反应设为 0/1000；Methods/目标描述明确最大化 BOF（biomass objective function）。

判定：Table 1 的 0/0 是“非培养基字段占位符/导出伪影”，不是可运行的生长约束。biomass objective function 不是培养基组分，其 media 列被写成 0/0。故复现用 published_operational_bounds（Afe_biomass 保持 0/1000），literal 仅作演示。

## B. Phase 1.6 计数口径（直接重算）

每条件与原始 SBML 不同的反应数（来自 phase1_6_bound_reconciliation.tsv）：

- FIM：12 exchange + 6 内部（Htpp、NH4pp、PPTtpp、SO41tpp、SO42tpp、SQRED1）+ 1 Afe_biomass = 19。
- TTM：13 exchange + 5 内部（CYT2、Htpp、NH4pp、PPTtpp、SQRED1）+ 1 Afe_biomass = 19。
- TSM：同 TTM = 19。

Phase 1.6 报告里“13 exchange + 1 CYT2 + 6 内部 + 1 Afe_biomass”的文字是把 FIM 与 TTM/TSM 的类别混写了；正确口径见上表。这是记账性修正，不影响任何数值。
