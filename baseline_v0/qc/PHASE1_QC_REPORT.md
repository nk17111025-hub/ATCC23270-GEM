# PHASE 1 STRUCTURAL & NUMERICAL QC REPORT

生成时间：2026-09-25。所有结果来自 `baseline_v0/phase1_qc.py`（可复跑），
未修改权威 SBML，未 gap-fill，未改 GPR / biomass / 持久 bounds，未启动复现或工程优化。

## 1. SBML 解析状态

- libSBML 解析 0 错误，Level 2 Version 1，3 区室（c/e/p），573 物种，615 反应，model id 为空。
- COBRApy 成功导入：615 反应、573 代谢物、0 基因（SBML 无 GPR，符合预期）。
- 目标函数：旧式 `OBJECTIVE_COEFFICIENT=1.0` 指向 `R_Ex_bio_LSQBKT_e_RSQBKT_`，COBRA 识别为最大化 biomass exchange。

## 2. 结构计数是否匹配 Phase 0

全部匹配：总反应 615、代谢物 573、区室 3、exchange 28、非 exchange 587、transport 52、biomass objective 1、其他内部 534。
重复 reaction ID 与 metabolite ID：0。`lb>ub`：0。空化学计量反应：0。固定（零宽）反应：ATPM、SCCR（XML 默认即固定，非异常）。

## 3. 质量/电荷问题数量与性质

| 分组 | 数量 |
|---|---:|
| BALANCED | 529 |
| CHARGE_IMBALANCED | 19 |
| MASS_AND_CHARGE_IMBALANCED | 1 |
| NOT_ASSESSABLE | 66（28 exchange + 38 内部） |

无“仅质量不平衡”反应。19 个电荷不平衡均元素已平衡、仅电荷差 1–4；集中在脂肪酸合成（3OAS4、EAR4/8、BHAH8、CFAS190）、醌/硫/辅因子氧化还原（DUQM、OMBQM、SQRED1、SULR2、PAPSR、TRDR、RNDR4、BTS5、HACT、CYSS）、核苷酸（NDPK6）、糖原（GLBRAN1、GLDBRAN2）以及转运 HCO3tpp。这些多是还原/氧化型辅因子电荷赋值不一致。
`CA2tpp` 为唯一 MASS_AND_CHARGE_IMBALANCED（Na +1、电荷 +1），属于已知不确定项。

## 4. 缺失公式/电荷信息

66 条 NOT_ASSESSABLE：28 条 exchange（边界反应，本不要求内部平衡）＋38 条内部反应。38 条内部反应主要是：氧化磷酸化电子载体（CYT1/CYT2/CYTA2/CYTAA31/CYTAA32/CYTBC1/CYTRED、SCCR 等，rusticyanin、细胞色素无公式/电荷）、LPS/脂质/磷脂/糖原/EPS/生物量大分子、bio/eps 转运与 biomass 目标函数。均标记为 NOT_ASSESSABLE，未当作平衡。

## 5. GPR 参考完整性

非空 GPR：461；空 GPR：154；疑似布尔语法错误：0。唯一基因（AFE_####）：508。与论文“507 genes”相比多 1 个唯一 AFE ID，属需复核的轻微差异，本阶段不处理。

## 6. Dead-end 与 blocked 反应

结构 dead-end 代谢物：13 个（多为核苷酸补救途径中间体，如 din/dad-2、R/RSH、bccp、dcyt/csn、pppi/dgsn、urate、dhptd、thymd、ppt 等）。多为“只有生产无消耗”或“只有消耗无生产”，提示补救/降解途径不完整，但不自动判为错误。

Blocked 反应（FVA，zero_cutoff=1e-7）：FIM 70、TTM 66、TSM 66。其中约 47–48 条标“potentially suspicious”，主要是次级/补救代谢（核苷酸补救、硫氧化、谷胱甘肽、芳香族/色氨酸、部分氧化磷酸化与糖代谢分支）。可能是正常条件性封锁，也可能提示通路缺失，需后续证据审查。

## 7. 可行性与 biomass objective

| 条件 | 状态 | objective | max|S·v| |
|---|---|---|---:|
| baseline_untouched（全 exchange 开放） | optimal | 3.6166 | 3.6e-12 |
| FIM（base medium） | optimal | 0.052076 | 1.4e-14 |
| TTM（base medium） | optimal | 0.052076 | 8.6e-14 |
| TSM（base medium） | optimal | 0.052076 | 4.2e-14 |

三个 base 条件均可行，biomass objective ≈ 0.052076 h^-1（碳源 h2co3 固定 2、电子供体开放，碳限制，故三者相同）。baseline_untouched 的 3.6166 是全部 exchange 开放下的非生理上限，仅作解析可用性检验。

## 8. 最大 |S·v|

所有求解的最大 |S·v| ≤ 3.6e-12，远低于目标 1e-8，数值残差达标。

## 9. GLPK / GLPK-exact 对比

baseline：GLPK 3.6166297847755096 vs GLPK-exact 3.6166297847754376（差 ~7e-14）。FIM：GLPK 0.05207638710274161 vs GLPK-exact 0.05207638710274114（差 ~5e-16）。两者一致，主结果得到交叉验证。

## 10. 能量生成/内部循环

内部循环 POSITIVE：全部 exchange 关闭后仍有约 140 条内部可逆反应可承载非零通量（典型 FBA 闭合回路，多数是合法可逆代谢环），列于 cycle_diagnostics.tsv。自由 ATP 生成 POSITIVE：全部 exchange 关闭、ATPM 放宽到 [0,1000] 并最大化 ATPM，目标 = 200，核心是 ATPS5rpp + ATPM + 可逆质子转运形成的热力学不可行能量环。自由 NADH / NADPH 生成 NEGATIVE（0.0）。

## 11. 需在 Phase 2 前审查的问题

`ISSUE`
自由 ATP 生成环（energy-generating cycle）。
`EVIDENCE`
全 exchange 关闭、最大化 ATPM，目标 200；环由 ATPS5rpp + ATPM + 可逆质子转运（Htpp 型）构成。
`MODEL IMPACT`
预测可能凭空产生 ATP，影响后续能量/生长定量。
`PROPOSED TEST`
在临时副本上对质子转运方向/bounds 做敏感性测试，核对 CYTBD/CYTAA3/ETC 质子系数。
`DO NOT MODIFY YET`

`ISSUE`
CA2tpp 质量+电荷不平衡（Na +1、电荷 +1）。
`EVIDENCE`
`ca2[p] -> ca2[c] + na1[p]`，Na 净 +1、电荷净 +1。
`MODEL IMPACT`
Ca/Na 转运化学计量影响电荷与 Na 平衡，属已知不确定项。
`PROPOSED TEST`
对照原文/数据库核对 CA2tpp 是否应为 Ca2+:Na+ 反向转运及其化学计量。
`DO NOT MODIFY YET`

`ISSUE`
19 条电荷不平衡反应（氧化还原辅因子电荷赋值不一致）。
`EVIDENCE`
SQRED1、SULR2、PAPSR、TRDR、HCO3tpp、NDPK6、GLBRAN1/GLDBRAN2 及脂肪酸/醌类反应电荷差 1–4。
`MODEL IMPACT`
影响电荷平衡校验与质子相关能量计算。
`PROPOSED TEST`
逐条核对还原/氧化型代谢物（NADPH/NADP、trdox/trdrd、q8/q8h2、hco3/h2co3 等）电荷。
`DO NOT MODIFY YET`

`ISSUE`
GPR 唯一基因 508 vs 论文 507。
`EVIDENCE`
461 条反应非空 GPR，抽取到 508 个唯一 AFE_####。
`MODEL IMPACT`
可能多出 1 个基因 ID（错写或未纳入模型），影响 GPR 一致性。
`PROPOSED TEST`
比对 508 与 507 基因集合，定位多出的 ID。
`DO NOT MODIFY YET`

`ISSUE`
38 条内部反应因缺失公式/电荷标 NOT_ASSESSABLE，主要是 ETC 电子载体与 LPS/脂质/糖原/EPS 大分子。
`EVIDENCE`
CYT1/CYT2/CYTA2/CYTAA31/CYTAA32/CYTBC1/CYTRED/SCCR 及 LPS/LIPS/PL-AFE/FA-PL-ACP-AFE/EPS 等无公式。
`MODEL IMPACT`
这部分无法做质量/电荷校验，其中 ETC 正是后续 Fe/S 能量分析关键。
`PROPOSED TEST`
对电子载体与大分子补元数据或按组分拆分，单独做 ETC 方向/质子系数敏感性。
`DO NOT MODIFY YET`

`ISSUE`
约 47–48 条 blocked 反应标“potentially suspicious”，多为补救/次级代谢。
`EVIDENCE`
核苷酸补救、硫氧化、谷胱甘肽、芳香族/色氨酸、部分氧化磷酸化分支在 FIM/TTM/TSM 均被封锁。
`MODEL IMPACT`
可能指示通路缺失或方向/bounds 问题，需区分“条件性正常封锁”与“真实缺陷”。
`PROPOSED TEST`
对 blocked 子集做 FVA 与基因证据对照，判断哪些应激活。
`DO NOT MODIFY YET`

已知不确定项（AFE_0269/SULDO、AFE_0660 malic enzyme、AFE_2841 BDGK、AFE_2024 6PGDH、AFE_0048 DoxD/TSQOC、CA2tpp、NA1tpp、CYTBD 质子、ETC 方向/系数、biomass 糖原系数）在本阶段仅标记，不自动解决；留待后续证据审查与敏感性分析。

## 产物

- `baseline_v0/qc/PHASE1_QC_REPORT.md`
- `baseline_v0/qc/structural_qc.tsv`
- `baseline_v0/qc/reaction_balance_qc.tsv`
- `baseline_v0/qc/gpr_reference_qc.tsv`
- `baseline_v0/qc/dead_end_metabolites.tsv`
- `baseline_v0/qc/blocked_reactions_FIM.tsv` / `_TTM.tsv` / `_TSM.tsv`
- `baseline_v0/qc/numerical_feasibility.tsv`
- `baseline_v0/qc/cycle_diagnostics.tsv`
- `baseline_v0/qc/solver_probe.txt`
- `baseline_v0/qc/phase1_summary.json`
- `baseline_v0/phase1_qc.py`
- `baseline_v0/logs/phase1_qc.log`

冻结 SBML（`baseline_v0/original/mmc3.xml`）保持字节一致，未改动。
