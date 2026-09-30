# PHASE 0 FREEZE REPORT

生成时间：2026-09-25

## 1. 输入冻结结果

四个 2016 权威文件已复制到 `baseline_v0/original/`，源文件在 `01_原始数据` 和 `02_论文`
未做任何改动。逐一重新计算 SHA256 后，与 A0 登记值逐项一致，全部 `match=True`。

| 文件 | 大小(byte) | SHA256 比对 |
|---|---:|---|
| mmc3.xml | 903066 | 一致 |
| mmc1.xls | 816640 | 一致 |
| mmc2.docx | 4871639 | 一致 |
| 2016主要模型GEM.pdf | 2711985 | 一致 |

完整清单见 `baseline_v0/input_manifest.tsv`。

## 2. 587 / 615 计数口径

已解决，不是错误。615 = 587（非 exchange）+ 28（exchange）。

587 = 534 代谢反应 + 52 转运反应 + 1 生物量目标函数（`Afe_biomass_mc507_WT_139p0M`）。

论文正文写 `587 metabolic and transport reactions`，SBML 与 mmc1.xls 额外包含 28 个
exchange 边界反应，所以总数为 615。详见 `baseline_v0/reaction_count_convention.md`。
模型不需要修改。

## 3. FIM / TTM / TSM 运行条件

### 3.1 直接 REPORTED 的条件

| 参数 | 值 | 适用条件 |
|---|---|---|
| Ex_h2co3[e]（唯一碳源 CO2/H2CO3 摄取） | -2 ~ -2 mmol/gDW/h | FIM / TTM / TSM |
| Ex_o2[e]（好氧呼吸） | -1000 ~ 1000（自由摄取） | FIM / TTM / TSM |
| Ex_fe2[e]（铁电子供体） | -1000 ~ 1000 | FIM |
| Ex_ttton[e]（四硫代硫酸盐电子供体） | -1000 ~ 0 | TTM |
| Ex_tsul[e]（硫代硫酸盐电子供体） | -1000 ~ 0 | TSM |
| CYT2（关闭铁氧化/rusticyanin 复合体） | 0 ~ 0 | TTM / TSM |
| ATPM（NGAM 维持能） | 3.475 ~ 3.475 mmol/gDW/h | FIM / TTM / TSM |
| GAM（生物量目标函数内） | 139 mmol ATP/gDW | FIM / TTM / TSM |
| Objective | 最大化 Ex_bio[e] | FIM / TTM / TSM |

无机离子（K/Na/Mg/Ca/Mn/Zn/Cu/Pi）在论文 Methods 2.4 中被描述为自由进出，
对应 mmc1.xls 中这些 exchange 的 ±1000，属于 REPORTED。

### 3.2 MODEL_DEFAULT 的条件

以下参数来自 mmc1.xls Table 1 默认边界，论文方法未逐字写明：

- Ex_n2[e] 开放、Ex_nh4[e] 关闭：氮源为 N2 固氮（NITF），NH4 不可用。
- Ex_so4[e]：FIM 下 0/0，TTM/TSM 下 0/1000（硫酸盐仅在 S 氧化条件下允许排出）。
- Ex_h2[e]、Ex_h2s[e]、Ex_s[e]、Ex_polypi[e]、Ex_ppt[e] 等：三条件下均关闭或按默认。
- 电子传递链（ETC）各反应质子化学计量（ATPS5rpp、CYTAA31、NADHI、CYTBC1、CYTAA32、
  CYTRED、CYTBO3/CYTBD）：数值已写入 mmc3.xml/mmc1.xls 反应式，属于模型默认；其
  拟合来源见论文 Supp Note 1/3。

### 3.3 INFERRED

- NGAM 数值 3.475 来自 `GAM 139 × 2.5%`（论文写 NGAM = 2.5% GAM）。该数值本身在
  mmc1.xls ATPM 中固定为 3.475，因此同时是 REPORTED 规则与 MODEL_DEFAULT 数值。

### 3.4 UNKNOWN

- 求解器可行性/最优性/整数容差：GLPK 通过 optlang 未暴露，标记 UNKNOWN，Phase 1 再核实。
- 论文报告的逐反应 FBA 预测通量（除生长/摄取外）没有完整数值化，暂不用于 Phase 0。

## 4. Supplementary Table 5 的 20 个实验点

已完整机器化提取到 `2016_scenario_constraints.json` 的 `experimental_points_table5`。

- FIM（fe2）：9 点（µ 0.086 → 0.01）
- TTM（s4o6）：5 点（µ 0.02 → 0.063）
- TSM（s2o3）：6 点（µ 0.026 → 0.125）
- 合计 20 点。

每点记录 µ、电子供体摄取速率、CO2 消耗；FIM 额外记录 O2 摄取。TTM/TSM 的 O2 摄取
论文表内未给，记为 null，未补猜。

## 5. 是否足以定义可复现 2016 场景

可以定义三个基础场景（FIM/TTM/TSM 最小培养基 + 目标函数 + GAM/NGAM），以及 20 个
实验验证点。基础场景的 exchange 边界全部来自 mmc1.xls Table 1 第 11–16 列，已机器化。

尚缺、且本阶段不自行补齐的：

- 论文 FVA 分析中“维持 100% 最大生长”以外的具体附加约束口径（Phase 2 复现时再对照原文）。
- 论文原求解器是 Insilico Discovery 3.3 + COBRA Toolbox 2.0 + GUROBI 5.5.0，本项目用
  cobra 0.32.1 + GLPK，复现时会有求解器数值差异，需要容差策略，不当作已解决。

## 6. 环境是否就绪

就绪。既有 venv（Python 3.14.6 / cobra 0.32.1 / optlang 1.9.1 / libsbml 5.21.2 / GLPK 5.0）
可以加载 mmc3.xml 并运行 FBA。环境快照和 requirements 锁定文件已写入
`baseline_v0/environment/`。未改动任何包版本。

## 7. 产物清单

- `baseline_v0/PHASE0_FREEZE_REPORT.md`
- `baseline_v0/input_manifest.tsv`
- `baseline_v0/environment/environment_snapshot.txt`
- `baseline_v0/environment/requirements_locked.txt`
- `baseline_v0/environment/rebuild_environment.md`
- `baseline_v0/scenario_definition/2016_scenario_constraints.json`
- `baseline_v0/scenario_definition/gpr_2016_reference.tsv`
- `baseline_v0/reaction_count_convention.md`
- `baseline_v0/qc/reaction_classification.json`
- `baseline_v0/logs/phase0_freeze.log`
- `baseline_v0/phase0_freeze.py`（本阶段生成脚本，可复跑）

## 8. 边界声明

本阶段未修改 SBML，未 gap-fill，未改动 biomass，未改动反应边界（除显式记录 2016
条件下的 exchange/CYT2 约束定义），未运行工程优化。原始 SBML 冻结副本与权威
mmc3.xml 字节一致。
