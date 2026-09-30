# MINIMAL_TRUSTWORTHY_BASELINE_V1_FREEZE_REPORT

冻结结论：`MINIMAL_TRUSTWORTHY_BASELINE_V1 — FREEZE PASS`

本报告确认：已批准的三个结构性修复形成稳定的最小工程基线，通过全部验收测试，无新 BLOCKER。

## 冻结前完整性

- 冻结 2016 基线（`baseline_v0/original/mmc3.xml`）SHA256 冻结前后均为
  `0dbaacfc07aafb892db0731a517160cf46006bc502b0dedd56f950c14b840246`，字节未变。
- v1 模型由冻结基线直接生成，仅含 3 个已批准的 bound 改动，无任何化学计量/GPR/代谢物/目标函数/交换/GAM/NGAM 变动。

## 三个已批准改动

| 反应 | 原 bound | 新 bound | 说明 |
|---|---|---|---|
| MACPD | 0..1000 | 0..0 | 禁用 malonyl-ACP 脱羧酶（Confidence=1、无 GPR/PMID、EC 2.3.1.38 矛盾），关闭人工 HCO3E bypass |
| ACOATA | -1000..1000 | 0..1000 | 禁止逆向（逆向闭合 bypass），正向保留为保守假设 |
| Htpp | -1000..1000 | 0..0 | 默认禁用；消除自由 ATP EGC（Htpp+ATPS5rpp+ATPM）；公开 Table 1 本就 0/0；状态 `DISABLED_BY_DEFAULT_DUE_TO_UNCONSTRAINED_EGC; BIOLOGICAL_PROTON_LEAK_DIRECTION_AND_CAPACITY_UNRESOLVED` |

反应级 diff 仅 3 行，且 `stoich_changed=False`、`reversible_changed=False`（纯 bound 改动）。

## 验收测试结果

### A. EGC 消除

| 测试 | 目标值 |
|---|---:|
| 自由 ATP（全 exchange 关闭，最大化 ATPM） | 0.0 |
| 自由 NADH | 0.0 |
| 自由 NADPH | 0.0 |
| 自由 ATP 复检 | 0.0 |

Htpp=0/0 后无任何替代自由 ATP/还原力环出现。

### B. 公开条件回归（FIM/TTM/TSM，schema 3.0 operational bounds）

| 条件 | biomass | Htpp FVA | MACPD FVA | ACOATA FVA |
|---|---|---|---|---|
| FIM | 0.052076 | [0,0] | [0,0] | [0.009, 0.009] |
| TTM | 0.052076 | [0,0] | [0,0] | [0.009, 0.009] |
| TSM | 0.052076 | [0,0] | [0,0] | [0.009, 0.009] |

biomass 保持参考值 0.052076，无回归；Htpp 确认为 [0,0]。

### C. HCO3E 补丁回归

- HCO3E KO：FIM/TTM/TSM 均 infeasible（人工 bypass 已消除，碳酸酐酶恢复必需）。
- MACPD FVA = [0,0]（已禁用）。
- ACOATA 仅正向（FVA 下限 0.009 > 0，无逆向）。
- 脂肪酸支持的 biomass 形成仍可行（biomass=0.052076），ACCOAC/MCOATA/正常 3OAS 保持。

### D. 基础生理 sanity

- 零无机碳 → 生长 0。
- 零必需供体 → infeasible。
- 零末端电子受体（O2+Fe3）→ infeasible（无有氧生长；NGAM=3.475 无法在无受体下满足）。
- 无隐藏有机碳摄取、无隐藏外部能量源、无闭交换自由 ATP。

## 冻结产物

- `baseline_v0/minimal_trustworthy_baseline/freeze/minimal_trustworthy_baseline_v1.xml`
- `.../freeze/MINIMAL_TRUSTWORTHY_BASELINE_V1_CHANGELOG.md`
- `.../freeze/reaction_level_diff_vs_2016.tsv`
- `.../freeze/egc_diagnostic.tsv`
- `.../freeze/regression.tsv`
- `.../freeze/sha256_manifest.json`
- `.../freeze/freeze_validation_summary.json`

## 版本身份

- 源模型：`baseline_v0/original/mmc3.xml`，SHA256 `0dbaacfc…b840246`。
- v1 模型：`baseline_v0/minimal_trustworthy_baseline/freeze/minimal_trustworthy_baseline_v1.xml`，SHA256 `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`。
- 生成时间见 `sha256_manifest.json`。
- 求解器：cobra 0.32.1 + GLPK 5.0（optlang）。

## 后续建议

下一阶段：`glucose local pathway audit`（有机碳局部通路审计）。

## 继续携带的非阻塞问题（不阻断本冻结）

- CYTRED mmc1/mmc3 质子方向差异
- CA2tpp 不平衡
- 剩余电荷不平衡反应
- TSM 硫代硫酸盐定量偏高
- 508 vs 507 GPR 计数
- ETC 质子系数不确定
- biomass 组成不确定
- 未完成的 G5X 全基因组复审
