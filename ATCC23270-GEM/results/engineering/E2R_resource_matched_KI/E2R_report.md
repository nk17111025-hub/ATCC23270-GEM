# E2R 资源匹配 KI 搜索

- 权威宿主 SHA256：026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3
- WT glucose-closed FIM pFBA：B_WT=0.0520763871；Fe2=-163.781529；O2=-38.4815427；Rubisco=2.16569779；CO2=-2。
- 混合无 KI：biomass=0.0911667123；glucose=-0.5；R_REF=0.791344958；CO2=-0.501268707；Fe2=-159.179299；O2=-38.4815427。
- 资源上限：Fe2≥-163.781529、O2≥-38.4815427；CO2/H2CO3为[-2,0]（以当前 FIM 的已定义无机碳供应上限为范围）；glucose 为[-0.5,0]且至少摄取0.05。历史配置 perturbation 中名为 WT 的状态开启了 glucose；本报告的 B_WT 明确使用 glucose-closed FIM。
- FIM 情景还要求至少 0.01 Fe2 uptake，以保留铁氧化条件。
- StrainDesign KI-only：显式 ko_cost={}，单位 ki_cost，max_cost=6，best，SCIP，seed=1；已通过正控。

## 结果概览
- 粗网格：21 个；无 KI 可行 13 个；六 KI 内搜索确认无解 8 个；全池 LP 不可行 0 个；超时 0 个。
- 细化：4 个；无 KI 可行 3 个；六 KI 内无解 1 个；超时 0 个。
- 当前共确认 0 个非零 KI 最小解。9 个场景的 StrainDesign 搜索在六 KI 内确认无解，其中 9 个全池 LP 可行；这些场景的设计若存在，需超过六个 KI。
- 主要边界：biomass ×1: Rubisco 25% → 10%; biomass ×1.25: Rubisco 50% → 25%; biomass ×1.5: Rubisco 50% → 25%; Rubisco 25%: biomass ×1 → ×1.25
- 六 KI 内搜索无解仅表示当前上限内无解；E2R 不据此推断更高 KI 数仍不可行。

## 场景摘要

| 场景 | 阶段 | B倍数 | Rubisco比例 | 状态 | KI数 | B | glucose | glucose-off B | ΔB | Rubisco | CO2 | Fe2 | O2 |
|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| B1_R100 | coarse | 1.0 | 1.0 | zero_KI_feasible | 0 | 0.0911667123679318 | -0.5 | 0.01902868740518603 | 0.07213802496274577 | 0.7913449619096874 | -0.5012687108297501 | -159.17929907582993 | -38.48154267515412 |
| B1_R90 | coarse | 1.0 | 0.9 | zero_KI_feasible | 0 | 0.0892629586223035 | -0.5 | 0.017125818664674632 | 0.07213713995762887 | 0.7122106909675091 | -0.4281936034325292 | -155.0764743623571 | -37.54586317874935 |
| B1_R75 | coarse | 1.0 | 0.75 | zero_KI_feasible | 0 | 0.08640865594275243 | -0.5 | 0.01427188619034121 | 0.07213676975241122 | 0.5935089648007544 | -0.3185737563015416 | -147.0120779685701 | -35.664806828855895 |
| B1_R50 | coarse | 1.0 | 0.5 | zero_KI_feasible | 0 | 0.0816523685596069 | -0.5 | 0.009514590793544771 | 0.07213777776606213 | 0.39567247875600003 | -0.135869176120751 | -134.3608323776732 | -32.72706890603389 |
| B1_R25 | coarse | 1.0 | 0.25 | zero_KI_feasible | 0 | 0.06217702919209398 | -0.39798606528683544 | 0.00475762013697272 | 0.05741940905512126 | 0.19783623937800002 | 0.0 | -103.02119071137345 | -25.201492692381468 |
| B1_R10 | coarse | 1.0 | 0.1 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1_R0 | coarse | 1.0 | 0.0 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p25_R100 | coarse | 1.25 | 1.0 | zero_KI_feasible | 0 | 0.09116582187522045 | -0.5 | 0.01902868740518685 | 0.0721371344700336 | 0.7913450735457185 | -0.5012734082095355 | -159.17929297384453 | -38.48154267515412 |
| B1p25_R90 | coarse | 1.25 | 0.9 | zero_KI_feasible | 0 | 0.0892629586223035 | -0.5 | 0.017125818664674632 | 0.07213713995762887 | 0.7122106909675091 | -0.4281936034325292 | -155.0764743623571 | -37.54586317874935 |
| B1p25_R75 | coarse | 1.25 | 0.75 | zero_KI_feasible | 0 | 0.08640865594275243 | -0.5 | 0.01427188619034121 | 0.07213676975241122 | 0.5935089648007544 | -0.3185737563015416 | -147.0120779685701 | -35.664806828855895 |
| B1p25_R50 | coarse | 1.25 | 0.5 | zero_KI_feasible | 0 | 0.0816523685596069 | -0.5 | 0.009514590793544771 | 0.07213777776606213 | 0.39567247875600003 | -0.13586917612075106 | -134.3608323776732 | -32.72706890603389 |
| B1p25_R25 | coarse | 1.25 | 0.25 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p25_R10 | coarse | 1.25 | 0.1 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p25_R0 | coarse | 1.25 | 0.0 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p5_R100 | coarse | 1.5 | 1.0 | zero_KI_feasible | 0 | 0.09116582187522045 | -0.5 | 0.01902868740518685 | 0.0721371344700336 | 0.7913450735457185 | -0.5012734082095355 | -159.17929297384453 | -38.48154267515412 |
| B1p5_R90 | coarse | 1.5 | 0.9 | zero_KI_feasible | 0 | 0.0892629586223035 | -0.5 | 0.017125818664674632 | 0.07213713995762887 | 0.7122106909675091 | -0.4281936034325292 | -155.0764743623571 | -37.54586317874935 |
| B1p5_R75 | coarse | 1.5 | 0.75 | zero_KI_feasible | 0 | 0.08640865594275243 | -0.5 | 0.01427188619034121 | 0.07213676975241122 | 0.5935089648007544 | -0.3185737563015416 | -147.0120779685701 | -35.664806828855895 |
| B1p5_R50 | coarse | 1.5 | 0.5 | zero_KI_feasible | 0 | 0.0816523685596069 | -0.5 | 0.009514590793544771 | 0.07213777776606213 | 0.39567247875600003 | -0.135869176120751 | -134.3608323776732 | -32.72706890603389 |
| B1p5_R25 | coarse | 1.5 | 0.25 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p5_R10 | coarse | 1.5 | 0.1 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p5_R0 | coarse | 1.5 | 0.0 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1_R17p5 | refined | 1.0 | 0.175 | infeasible_within_6_KI | None | None | None | None | None | None | None | None | None |
| B1p125_R25 | refined | 1.125 | 0.25 | zero_KI_feasible | 0 | 0.06217702919209398 | -0.39798606528683544 | 0.00475762013697272 | 0.05741940905512126 | 0.19783623937800002 | 0.0 | -103.02119071137345 | -25.201492692381468 |
| B1p25_R37p5 | refined | 1.25 | 0.375 | zero_KI_feasible | 0 | 0.07927378263398152 | -0.5 | 0.007136015165955748 | 0.07213776746802578 | 0.296754359067 | -0.044519293458944356 | -128.24812862203896 | -31.311428692228667 |
| B1p5_R37p5 | refined | 1.5 | 0.375 | zero_KI_feasible | 0 | 0.07927378263398152 | -0.5 | 0.007136015165955748 | 0.07213776746802578 | 0.296754359067 | -0.044519293458944356 | -128.24812862203896 | -31.311428692228667 |

## 自适应细化

根据粗网格相邻可行性或 KI 数变化增加 4 个场景（最多 12）。
- B1_R17p5：biomass×1.0，Rubisco 17.5%
- B1p125_R25：biomass×1.125，Rubisco 25%
- B1p25_R37p5：biomass×1.25，Rubisco 37.5%
- B1p5_R37p5：biomass×1.5，Rubisco 37.5%

## 最小 KI 组合
- B1_R100 / 0 KI / S00：候选 ，供体 ；B=0.0911667123679318，glucose=-0.5，glucose-off B=0.01902868740518603，ΔB=0.07213802496274577，贡献验证=True。
- B1_R90 / 0 KI / S00：候选 ，供体 ；B=0.0892629586223035，glucose=-0.5，glucose-off B=0.017125818664674632，ΔB=0.07213713995762887，贡献验证=True。
- B1_R75 / 0 KI / S00：候选 ，供体 ；B=0.08640865594275243，glucose=-0.5，glucose-off B=0.01427188619034121，ΔB=0.07213676975241122，贡献验证=True。
- B1_R50 / 0 KI / S00：候选 ，供体 ；B=0.0816523685596069，glucose=-0.5，glucose-off B=0.009514590793544771，ΔB=0.07213777776606213，贡献验证=True。
- B1_R25 / 0 KI / S00：候选 ，供体 ；B=0.06217702919209398，glucose=-0.39798606528683544，glucose-off B=0.00475762013697272，ΔB=0.05741940905512126，贡献验证=True。
- B1p25_R100 / 0 KI / S00：候选 ，供体 ；B=0.09116582187522045，glucose=-0.5，glucose-off B=0.01902868740518685，ΔB=0.0721371344700336，贡献验证=True。
- B1p25_R90 / 0 KI / S00：候选 ，供体 ；B=0.0892629586223035，glucose=-0.5，glucose-off B=0.017125818664674632，ΔB=0.07213713995762887，贡献验证=True。
- B1p25_R75 / 0 KI / S00：候选 ，供体 ；B=0.08640865594275243，glucose=-0.5，glucose-off B=0.01427188619034121，ΔB=0.07213676975241122，贡献验证=True。
- B1p25_R50 / 0 KI / S00：候选 ，供体 ；B=0.0816523685596069，glucose=-0.5，glucose-off B=0.009514590793544771，ΔB=0.07213777776606213，贡献验证=True。
- B1p5_R100 / 0 KI / S00：候选 ，供体 ；B=0.09116582187522045，glucose=-0.5，glucose-off B=0.01902868740518685，ΔB=0.0721371344700336，贡献验证=True。
- B1p5_R90 / 0 KI / S00：候选 ，供体 ；B=0.0892629586223035，glucose=-0.5，glucose-off B=0.017125818664674632，ΔB=0.07213713995762887，贡献验证=True。
- B1p5_R75 / 0 KI / S00：候选 ，供体 ；B=0.08640865594275243，glucose=-0.5，glucose-off B=0.01427188619034121，ΔB=0.07213676975241122，贡献验证=True。
- B1p5_R50 / 0 KI / S00：候选 ，供体 ；B=0.0816523685596069，glucose=-0.5，glucose-off B=0.009514590793544771，ΔB=0.07213777776606213，贡献验证=True。
- B1p125_R25 / 0 KI / S00：候选 ，供体 ；B=0.06217702919209398，glucose=-0.39798606528683544，glucose-off B=0.00475762013697272，ΔB=0.05741940905512126，贡献验证=True。
- B1p25_R37p5 / 0 KI / S00：候选 ，供体 ；B=0.07927378263398152，glucose=-0.5，glucose-off B=0.007136015165955748，ΔB=0.07213776746802578，贡献验证=True。
- B1p5_R37p5 / 0 KI / S00：候选 ，供体 ；B=0.07927378263398152，glucose=-0.5，glucose-off B=0.007136015165955748，ΔB=0.07213776746802578，贡献验证=True。

## Pareto 前沿与候选拐点
去重后非支配工程点 6 个。
- B1p5_R100：0 KI，B=0.09116582187522045，Rubisco cap=0.7913449575127987，实际 Rubisco=0.7913450735457185，glucose-off B=0.01902868740518685，ΔB=0.0721371344700336，CO2=-0.5012734082095355。
- B1p5_R37p5：0 KI，B=0.07927378263398152，Rubisco cap=0.2967543590672995，实际 Rubisco=0.296754359067，glucose-off B=0.007136015165955748，ΔB=0.07213776746802578，CO2=-0.044519293458944356。
- B1p5_R50：0 KI，B=0.0816523685596069，Rubisco cap=0.3956724787563993，实际 Rubisco=0.39567247875600003，glucose-off B=0.009514590793544771，ΔB=0.07213777776606213，CO2=-0.135869176120751。
- B1p5_R75：0 KI，B=0.08640865594275243，Rubisco cap=0.593508718134599，实际 Rubisco=0.5935089648007544，glucose-off B=0.01427188619034121，ΔB=0.07213676975241122，CO2=-0.3185737563015416。
- B1p5_R90：0 KI，B=0.0892629586223035，Rubisco cap=0.7122104617615188，实际 Rubisco=0.7122106909675091，glucose-off B=0.017125818664674632，ΔB=0.07213713995762887，CO2=-0.4281936034325292。
- B1p125_R25：0 KI，B=0.06217702919209398，Rubisco cap=0.19783623937819966，实际 Rubisco=0.19783623937800002，glucose-off B=0.00475762013697272，ΔB=0.05741940905512126，CO2=0.0。
- 最低复杂度前沿点：B1p125_R25（0 KI）；最强 CBB 限制前沿点：B1p125_R25（Rubisco cap=0.19783623937819966）。
- 高增长拐点：B1p5_R37p5，Rubisco cap=0.2967543590672995，biomass=0.07927378263398152（目标阈值 0.0781146）。
- 再加强 Rubisco 限制至 0.19783623937819966 后，biomass 降至 0.06217702919209398，减少 0.0170968（21.57%）。
- KI 数量拐点：本次有效解均为 0 KI，暂未出现可比较的 KI 复杂度拐点。

## TTM/TSM 固定 KI 验证
抽查高增长点、满足 1.5×B_WT 的 Rubisco 限制拐点和前沿中最强 Rubisco 限制点（相同点合并）。native_Rubisco_unrestricted 使用 TTM/TSM 原条件，不施加 FIM Rubisco cap；strict_FIM_cap 额外施加 FIM 场景 cap；donor_closed_control 关闭对应硫供体。三种模式都使用 Fe2 微量营养范围、FIM O2 资源上限、glucose 上限和固定 CO2。
- B1p5_R100 / TTM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775688668，Rubisco=2.4142444854100282，CO2=-2.0，Ex_ttton[e] flux=-8.897549823265008。
- B1p5_R100 / TTM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p5_R100 / TTM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p5_R100 / TSM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775687497，Rubisco=2.4142444854100544，CO2=-2.0，Ex_tsul[e] flux=-15.57071219071332。
- B1p5_R100 / TSM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。
- B1p5_R100 / TSM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。
- B1p5_R37p5 / TTM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775688668，Rubisco=2.4142444854100282，CO2=-2.0，Ex_ttton[e] flux=-8.897549823265008。
- B1p5_R37p5 / TTM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p5_R37p5 / TTM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p5_R37p5 / TSM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775687497，Rubisco=2.4142444854100544，CO2=-2.0，Ex_tsul[e] flux=-15.57071219071332。
- B1p5_R37p5 / TSM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。
- B1p5_R37p5 / TSM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。
- B1p125_R25 / TTM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775688668，Rubisco=2.4142444854100282，CO2=-2.0，Ex_ttton[e] flux=-8.897549823265008。
- B1p125_R25 / TTM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p125_R25 / TTM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_ttton[e] flux=None。
- B1p125_R25 / TSM / native_Rubisco_unrestricted：optimal，biomass=0.13019096775687497，Rubisco=2.4142444854100544，CO2=-2.0，Ex_tsul[e] flux=-15.57071219071332。
- B1p125_R25 / TSM / strict_FIM_cap：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。
- B1p125_R25 / TSM / donor_closed_control：infeasible，biomass=None，Rubisco=None，CO2=None，Ex_tsul[e] flux=None。

未解决场景 0 个；超时均保留为未解决，未记为不可行。
