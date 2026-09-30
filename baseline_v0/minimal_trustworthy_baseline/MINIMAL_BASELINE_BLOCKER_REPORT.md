# MINIMAL_BASELINE_BLOCKER_REPORT

生成时间：2026-09-26。集中问题报告，供主项目规划对话裁决下一轮修复。

## 汇总表

| ID | Severity | Trigger | Core reactions | Engineering impact | Decision required |
|---|---|---|---|---|---|
| MBL-BLOCK-001 | BLOCKER | 原始默认 Htpp 可逆（未施加公开条件 Htpp=0/0） | Htpp + ATPS5rpp + ATPM | 自由 ATP 生成，污染能量平衡/Fe-S 能量分析 | 是否把 Htpp=0/0 烘焙进 v1 默认（或其它已批准的 Htpp 修复） |
| MBL-CARRY-002 | DEFER | 恒定存在 | CYTRED | 质子方向 mmc1/mmc3 相反；验证条件惰性 | 无需（Phase 2.1 已裁定） |
| MBL-CARRY-003 | MEDIUM | 恒定存在 | CA2tpp | Na/Ca 质量+电荷不平衡 | 后续证据审查 |
| MBL-CARRY-004 | MEDIUM | 恒定存在 | 19 条电荷不平衡反应 | 氧化还原辅因子电荷不一致 | 后续证据审查 |
| MBL-CARRY-005 | MEDIUM | TSM 硫氧化 | TSQOC/4THASE/CYTBD | 硫代硫酸盐摄取定量偏高 | 后续证据审查 |
| MBL-CARRY-006 | LOW | 恒定存在 | GPR | 508 vs 507 唯一基因 | 后续机械核对 |

## 已批准补丁状态（Phase 3.1 HCO3E bypass）

已正确应用并验证通过：

- MACPD：禁用（bounds 0..1000 → 0..0）。
- ACOATA：禁止逆向（-1000..1000 → 0..1000，保留正向作为保守假设）。

回归结果：

- FIM/TTM/TSM 生长 = 0.052076（与补丁前一致，无回归）。
- HCO3E KO = infeasible（三条件，人工 bypass 已消除）。
- MACPD FVA = [0,0]（已禁用）。
- ACOATA FVA = [0.009, 0.009]（仅正向，供脂肪酸合成引物）。
- 脂肪酸合成仍正常（biomass 可形成，ACCOAC/MCOATA/3OAS 保持）。
- 零碳生长 0；零供体 infeasible。

标记：PHASE3.1_HCO3E_PATCH_PASS。

---

## ISSUE MBL-BLOCK-001

## Short title
Htpp 可逆质子转运 + ATPS5rpp + ATPM 形成自由 ATP 生成环（EGC）。

## Trigger condition
模型使用原始默认 bound（Htpp = [-1000, 1000] 可逆）且关闭全部外部营养/供体摄取时，最大化 ATPM。

## Observable phenotype
全部 exchange 关闭、最大化 ATPM，目标 = 200 mmol/gDW/h。核心三反应：
- Htpp 通量 +1000（c→p 方向）
- ATPS5rpp 通量 +200
- ATPM 通量 +200

施加公开条件 Htpp=0/0 后，同一测试目标 = 0（EGC 被阻断）。

## Why it matters
有机碳工程的核心是能量平衡（有机碳是否减轻 Fe/S 供能负担）。若模型默认态允许无中生有生成 ATP，会污染后续所有能量/生长定量，无法可信判断“保留 Fe/S 能量代谢”与“有机碳供能”的相对贡献。

## Reproduction steps
1. 载入工作模型（默认 Htpp 可逆）。
2. 将所有 Ex_* exchange 设 0。
3. ATPM 放宽到 [0, 1000] 并设为目标，最大化。
4. 观察目标 = 200。

## Relevant reactions
Htpp（h[c]<=>h[p]，可逆）、ATPS5rpp（ATP 合酶，5 h[p]→4 h[c]）、ATPM（ATP 维持）。

## Relevant metabolites
h[c]、h[p]、atp_c、adp_c、pi_c、h2o_c。

## Key fluxes
Htpp +1000；ATPS5rpp +200；ATPM +200。

## Net cycle / net stoichiometry
ATPS5rpp + ATPM + 可逆 Htpp 净效果：5 h[p] → 5 h[c]（质子梯度跨膜闭合），净 ATP 凭空生成。Htpp ub=1000 / 5 = 200 解释了观测上限。

## GPR / confidence / provenance
Htpp：转运（内质膜质子转运），2016 默认可逆，无方向约束。公开 Table 1（mmc1.xls 第 11–16 列）将 Htpp 设为 0/0。

## Current evidence
- Phase 1.5：确认最小环 Htpp+ATPS5rpp+ATPM，方向诊断确认仅 c→p 方向闭合。
- Phase 1.6：公开最小培养基 Table 1 明确 Htpp=0/0，EGC 在公开条件下失活。
- 本次复测：published Htpp=0 → free ATP=0；raw Htpp 可逆 → free ATP=200。

## Severity
BLOCKER（冻结门槛 #5：不允许存在未解决的 EGC/自由 ATP 问题）。

## Impact area
energy（Fe oxidation / sulfur oxidation / biomass / 未来 glucose engineering 的能量解释均受影响）。

## Interpretation
- 事实：raw 默认 Htpp 可逆允许自由 ATP；公开条件 Htpp=0/0 阻断。
- 可能解释：Htpp 是质子跨膜转运，生物学上不应自由可逆；公开培养基已将其关闭。
- 未决：v1 冻结模型是否应把 Htpp=0/0 烘焙进默认 bound，还是继续依赖“公开条件运行时施加”。

## Candidate repair options（仅提案，不实施）
1. 将 Htpp 默认 bound 改为 0/0（对齐公开 Table 1，最小证据改动）。
2. 将 Htpp 改为不可逆（c→p 或 p→c 单向，需方向证据）。
3. 保持 Htpp 可逆，但在所有工程运行强制施加 Htpp=0/0，并在文档显著标记。

## Recommended decision needed from main project conversation
是否授权把 Htpp 默认 bound 设为 0/0（或其它已批准方向约束）作为 v1 的一部分；否则 v1 默认态仍含自由 ATP 环，不可冻结。

---

## 附：非阻塞待办（DEFER，不阻止冻结但需记录）

- CYTRED mmc1/mmc3 质子方向相反（验证条件惰性，Phase 2.1 已裁定）。
- CA2tpp 质量+电荷不平衡。
- 19 条电荷不平衡反应（氧化还原辅因子电荷）。
- TSM 硫代硫酸盐摄取定量偏高（硫氧化通路）。
- 508 vs 507 GPR 计数。
- ETC 质子系数不确定、biomass 组成不确定。

以上均在 Phase 1–3 中诊断并记录，对“有机碳→biomass 并保留 Fe/S 能量”的当前工程目标非阻塞；除 MBL-BLOCK-001 外无需本阶段裁决。
