# PHASE 3 PHYSIOLOGICAL SANITY REPORT

生成时间：2026-09-25。诊断专用，未修复模型、未 gap-fill、未改 GPR、未引入有机碳工程、未创建 2026 校正模型。全部用 schema 3.0 published_operational_bounds，GLPK 筛选。

## 结论摘要

冻结 2016 iMC507 在全部核心生理控制上表现合理：需碳、需 RUBISCO、需电子供体、需电子受体、氧化还原方向正确、无自由 ATP、碳响应单调。仅一处需生物学复查：HCO3E（碳酸酐酶）非必需。

## A. 阳性对照好氧表型（H2CO3=2，FBA/pFBA）

| 条件 | growth | O2 | 供体 | RUBISCO | ATPS5rpp | ATPM |
|---|---:|---:|---:|---:|---:|---:|
| FIM | 0.0521 | -38.66 | Fe2 -164.5 | 2.166 | 20.75 | 3.475 |
| TTM | 0.0521 | -15.79 | ttton -5.22 | 2.166 | 20.74 | 3.475 |
| TSM | 0.0521 | -15.79 | tsul -9.13 | 2.166 | 20.74 | 3.475 |

FIM：Fe2 摄取 -164.5、Fe3 输出 +164.5、O2 摄取、CO2 固定；ETC 走 CYT2(77.3)/CYT1(82.3)/CYTAA31(77.3) 下坡 + CYTBC1(4.94)/NADHI(5.21) 上坡（NADH 生成）。TTM/TSM：硫供体 → 硫酸盐（+20.8/+18.2）、O2 摄取、CYTBD(31.6)/TSQOC(5.2)/4THASE1(10.4) 硫氧化，无 Fe 氧化。pFBA 通量与 FBA 一致。

## B. 碳依赖

H2CO3 关闭后三条件生长均为 0。无任何含碳 exchange 被意外摄取（关闭后仅 fe2/o2/h/h2o/ttton/tsul 等非碳维持能摄取，growth=0）。

## C. Calvin 依赖

| 处理 | FIM | TTM | TSM |
|---|---:|---:|---:|
| baseline | 0.0521 | 0.0521 | 0.0521 |
| RUBISCO 关闭 | infeasible 0 | infeasible 0 | infeasible 0 |
| HCO3tpp 关闭 | infeasible 0 | infeasible 0 | infeasible 0 |
| HCO3E 关闭 | 0.0521 | 0.0521 | 0.0521 |

RUBISCO 与 HCO3tpp 必需（PASS）。HCO3E 非必需（生长不变）——见下方 NEEDS BIOLOGICAL REVIEW。

## D. 电子供体依赖

单供体关闭 → 生长 0（三条件）。全部还原态无机供体（fe2/ttton/tsul/h2/h2s/s）关闭 → 生长 0。无替代供体被意外开启。

## E. 电子受体依赖

- E1 仅关 O2：三条件生长 0；FIM 中 Fe3 为输出（+10.3，非摄取），无替代末端受体被使用。
- E2 关 O2 + Fe3：三条件生长 0。

无末端电子受体时无热力学捷径生长。

## F. 氧化还原产物方向

全部正确：FIM Fe2 消耗 / Fe3 输出 / O2 消耗 / 无机碳固定；TTM/TSM 硫供体消耗 / 硫酸盐输出 / O2 消耗 / 无机碳固定。无反向或无解释方向。

## G. 能量 sanity

三条件 Htpp=0；全 exchange 关闭最大化 ATPM = 0（EGC 不可达）。各诊断扰动后未见其它 ATP 生成环。

## H. 碳与能量泄漏清单

三条件只有 h2co3 为含碳摄取（intended carbon），供体/受体/矿物/氮素摄取均正常。唯一非零含碳 exchange 是 4hba 微量输出（~1e-6，可忽略）。FIM 的质子摄取 -154.7 是 Fe2 氧化耗酸（pH 平衡），属生理正常，不是能量泄漏。

## I. 碳响应敏感性（H2CO3 = 0/0.5/1/2/4）

三条件生长单调线性（0 → 0.013 → 0.026 → 0.052 → 0.104），供体/O2 需求同步线性增长。零碳下生长为 0，无跳变、无通路切换。

## 异常表型详情

`SCENARIO`：HCO3E（碳酸酐酶）关闭。
`EXPECTED BIOLOGICAL INTERPRETATION`：论文描述 HCO3E 把 hco3 转为 co2 供 RUBISCO 固定。
`OBSERVED MODEL BEHAVIOR`：三条件下生长不变（0.0521），RUBISCO 仍 2.166。
`EXACT FLUX SOURCE`：HCO3tpp 仍把 h2co3 转为 hco3；RUBISCO 的 co2 由代谢脱羧/其他含碳路径回收供给，模型对 HCO3E 存在冗余。
`IMPACT ON FUTURE ENGINEERING`：未来有机碳工程需注意模型碳固定路径对碳酸酐酶不敏感，避免据此误判碳流瓶颈。
`STATUS`：NEEDS BIOLOGICAL REVIEW（不修复）。

## J. 需携带、不解决的已知问题

- 原始 SBML Htpp EGC（默认态缺陷，公开条件失活）
- CYTRED mmc1/mmc3 质子方向不一致（验证条件惰性）
- CA2tpp 质量/电荷不平衡
- 电荷不平衡反应集（19 条）
- TSM 硫代硫酸盐摄取定量偏高
- 508 vs 507 GPR 计数未解
- ETC 质子系数不确定
- biomass 组成不确定

## 产物

- `baseline_v0/qc/physiology/PHASE3_PHYSIOLOGICAL_SANITY_REPORT.md`
- `positive_control_phenotypes.tsv` / `positive_control_fva.tsv`
- `carbon_dependency.tsv` / `calvin_dependency.tsv`
- `electron_donor_dependency.tsv` / `electron_acceptor_dependency.tsv`
- `redox_product_directions.tsv` / `energy_sanity.tsv`
- `exchange_leakage_inventory.tsv` / `h2co3_response.tsv`
- `phase3_summary.json`
- `baseline_v0/phase3_physiological_sanity.py`

未修复模型，未 gap-fill，未开始工程优化。
