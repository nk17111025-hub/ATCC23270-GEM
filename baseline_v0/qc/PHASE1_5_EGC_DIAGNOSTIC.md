# PHASE 1.5 EGC IMPACT DIAGNOSTIC REPORT

生成时间：2026-09-25。诊断仅在临时模型副本上执行，未修改冻结 SBML，未改持久 bounds，未做任何修复/校正，未进入 Phase 2。

## 结论摘要（五个问题）

1. 自由 ATP 环是否被确凿复现？是。
2. 支撑它的最小反应集：Htpp + ATPS5rpp + ATPM（仅 3 个）。
3. 在 FIM/TTM/TSM 最大生长下：可选（alternate optimum），不是必需。
4. 关闭/限制 Htpp 是否实质改变 20 个 Table 5 预测？生长几乎不变，但电子供体摄取和 O2 摄取被实质改变。
5. 是否必须在 Phase 2 前解决？不阻塞先复现原模型；但 Table 5 的电子供体/O2 验证会受影响，需作为已记录缺陷携带，并在复现时采用一致的 flux 分配方法或公开的 Table 1 最小培养基 Htpp=0 约束。

## 1. 最小环复现

全部营养/电子供体 exchange 关闭、最大化 ATPM，目标 = 200。完整非零通量向量仅 3 个反应：

| 反应 | 通量 |
|---|---:|
| Htpp | 1000.0（c→p 方向） |
| ATPS5rpp | 200.0 |
| ATPM | 200.0 |

这与 `Htpp ub=1000 / 5 H+ per ATP = 200` 一致。

逐个中断：

| 被置 0 的反应 | ATPM 最大值 |
|---|---:|
| Htpp=0 | 0.0 |
| ATPS5rpp=0 | 0.0 |

两者任一置 0 均使 ATPM 归零，环被确认。

## 2. Htpp 方向敏感性（全部 exchange 关闭，最大化 ATPM）

| 情况 | ATPM 最大 | ATPS5rpp | Htpp |
|---|---:|---:|---:|
| 原始 [-1000,1000] | 200 | 200 | 1000 |
| 仅 c→p [0,1000] | 200 | 200 | 1000 |
| 仅 p→c [-1000,0] | 0 | 0 | 0 |
| 关闭 [0,0] | 0 | 0 | 0 |

结论：该环只在 Htpp 允许 c→p 方向时存在。本阶段不判定哪个方向生物正确。

## 3. 最大生长下 EGC 是否被使用

用 Phase 1 相同的公开模型约束（exchange + CYT2），FIM/TTM/TSM 最大化 biomass：

| 条件 | biomass | ATPM | ATPS5rpp | Htpp |
|---|---:|---:|---:|---:|
| FIM | 0.052076 | 3.475 | 20.751 | 154.62 |
| TTM | 0.052076 | 3.475 | 20.705 | 119.64 |
| TSM | 0.052076 | 3.475 | 20.705 | 119.64 |

100% 最大生长下 FVA：

| 条件 | Htpp [min,max] | ATPS5rpp [min,max] |
|---|---:|---:|
| FIM | [-835.5, 1000] | [20.75, 388.0] |
| TTM | [-1000, 1000] | [20.70, 953.0] |
| TSM | [-1000, 1000] | [20.70, 953.0] |

判定：ATPS5rpp 下限约 20.7 是真实生长所需（ATP 合酶必需）；Htpp 的 FVA 区间跨 0，说明支撑 EGC 的
c→p 通量是可选的 alternate optimum，不是达到最大生长的必需条件。

`EGC_USED_IN_OPTIMAL_SOLUTION = ALTERNATIVE_OPTIMUM_POSSIBLE`

## 4. 对 biomass 的影响（临时副本）

| 条件 | 处理 | 状态 | biomass 变化 |
|---|---|---:|---:|
| FIM/TTM/TSM | Htpp 关闭 | optimal | ~1e-13 到 1e-9（可忽略） |
| FIM/TTM/TSM | Htpp 仅 c→p | optimal | 0 |
| FIM/TTM/TSM | Htpp 仅 p→c | optimal | ~1e-13 到 1e-9 |
| FIM/TTM/TSM | ATPS5rpp 关闭 | infeasible | -100% |

biomass 对 Htpp 关闭/限制不敏感；ATPS5rpp 是生长必需。

## 5. Table 5 预检（EGC 影响）

每个点固定其实验 H2CO3 输入，比较原始模型 vs Htpp 关闭副本。

- 生长：差异约 1e-13 到 1e-6（相对可忽略）。
- 电子供体摄取：被实质改变。FIM 例：Fe2 摄取原始 -19.25 vs Htpp 关闭 -303.9（实验 296.67），约 16 倍差异。
- O2 摄取（FIM）：-0.015 vs -71.18，被实质改变。
- TTM/TSM 电子供体摄取同样被改变约 9 倍（原始远低于实验，Htpp 关闭后接近实验）。

即：EGC 不改变生长，但实质改变电子供体/O2 摄取这两个 Table 5 验证目标。

## 6. GPR 508 vs 507（机械核查）

508 个唯一 AFE ID 全部格式规范（AFE_####），未发现重复、错写或可疑符号；AFE_0697 不在其中。
仅凭 mmc1.xls 无法机械定位多出的 1 个 ID；最可能是某 GPR 引用了未计入“507 genes”模型基因列表的
一个 ID，或论文的 507 计数口径略不同。保持 unresolved，未改 GPR。

## 关键证据

公开最小培养基（mmc1.xls Table 1 第 11–16 列）将 Htpp 在三条件下均设为 0/0（另把 SQRED1、NH4pp、
PPTtpp、SO41tpp/SO42tpp 设为 0/0）。这说明原始 2016 模拟很可能 Htpp 已关闭，EGC 未污染其验证。
冻结 SBML 默认 Htpp 为 [-1000,1000]，是允许该环的原因。

## 产物

- `baseline_v0/qc/PHASE1_5_EGC_DIAGNOSTIC.md`
- `baseline_v0/qc/phase1_5_minimal_egc_flux.tsv`
- `baseline_v0/qc/phase1_5_egc_interruption.tsv`
- `baseline_v0/qc/phase1_5_htpp_direction.tsv`
- `baseline_v0/qc/phase1_5_baseline_growth_egc.tsv`
- `baseline_v0/qc/phase1_5_egc_fva.tsv`
- `baseline_v0/qc/phase1_5_biomass_impact.tsv`
- `baseline_v0/qc/phase1_5_table5_egc_precheck.tsv`
- `baseline_v0/qc/phase1_5_gpr_count_check.tsv`
- `baseline_v0/qc/phase1_5_summary.json`
- `baseline_v0/phase1_5_egc_diagnostic.py`

未修复模型，未改持久 bounds，未 gap-fill，未开始工程优化。
