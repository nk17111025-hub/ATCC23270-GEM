# PHASE 3.2-R VERIFICATION REPORT

审计对象：`minimal_trustworthy_baseline_v1`（SHA256 `bd9715e6…3689af8c`）。只验证、不修复。关闭 Phase 3.2 的四个方法学缺口。

## 结论

`PHASE3_2R_PASS — PHASE 3.2 MAY BE FORMALLY CLOSED`

未发现 BLOCKER。四个缺口全部闭合。

## 完整性

- v1 SHA256 校验通过，与冻结一致。
- 615 反应 / 573 代谢物。
- 三处已批准 bound 确认：MACPD=0/0、ACOATA=0/1000、Htpp=0/0。
- 参考生长（非 null）：FIM=0.052076、TTM=0.052076、TSM=0.052076。

## Gap 1 闭合：三条件 KO 救援复现

27 个 KO × 3 条件 = 81 次测试。KO 分类在三条件完全一致：

- lethal/infeasible：13（HCO3E、PGI1、PYK、MCOATA、TKT1、RUBISCO、PDH、ENO、TPI、PGK、RPE、RPI、PRUK）。
- zero-growth：2（CS、ACS）。
- partial rescue：2（MDH、FUM）。
- near-complete rescue：10（PFK、FBP、GAPD1、PPC、ACCOAC、G6PDH2、TALA、SBPASE、NADTRHD、BDGK）。

未出现任何条件特异救援。CONFIRMED。

## Gap 2 闭合：全网络救援追踪

对每个保留有意义的 KO，比较扰动 pFBA 与参考 pFBA 在全部 615 反应上的通量差（|delta|>1e-4）。

结果：救援路径未离开中央碳 + 已知能量代谢。外部高 delta 反应均为：

- 外部 exchange（Ex_h/Ex_fe2/Ex_fe3/Ex_h2o/Ex_so4/Ex_o2）；
- ETC（FIM：CYT2/CYTAA31/CYT1；TTM/TSM：CYTBD/TSQOC）；
- 转运（Htex/H2Otex/H2Otpp/O2tpp/O2tex/SO4tex）。

这些都是 KO 造成氧化还原/质子失衡后，能量代谢的正常补偿，非人工跨子系统救援。无重大跨子系统人工救援。CONFIRMED。

## Gap 3 闭合：CO2/HCO3 净化学计量

- ACCOAC(reverse) + BIOC1 + BIOC2 净化学计量 = 零（纯内部 futile 环）。
- HCO3E(reverse)+RUBISCO = 经典 CO2 固定。
- ACCOAC(reverse) = 丙二酰辅酶A脱羧酶样（产 ATP + HCO3）。
- HCO3E KO 在三条件均 infeasible，无任何候选组合恢复 HCO3E 无关生长。

未发现新的无机碳捷径。CONFIRMED。

## Gap 4 闭合：标准 vs loopless FVA

- `find_cyclic_reactions` 识别 30 个参与内部环的反应。
- `loopless_solution`（无环解）显示：ACCOAC=0.068（正向，非极端的 -999.93）、BIOC1=0、BIOC2=0。
- 结论：ACCOAC[-999.93,0.068] 与 BIOC1/BIOC2[0,1000] 的极端标准 FVA 范围是 `LOOP-DRIVEN FVA FREEDOM`，不是生物休眠容量。

方法说明：使用 cobra 的 `loopless_solution`（无环解）+ `find_cyclic_reactions` + 净化学计量作为无环/去循环判据；`loopless_fva_iter` 的逐反应 min/max 未稳定产出，故以无环解为无环参考（已文档化，未静默替代为普通 FVA）。

## 氧化还原扰动（Test E）

GAPD1 KO：GAPD2 + 电子传递链补偿（三条件一致），NADTRHD 参与。NADTRHD KO：通过增加电子传递链/质子交换补偿，生长不变。MDH/PPC 扰动：经能量代谢补偿。无低置信度反应成为主导氧化还原救援。CONFIRMED（措辞：替代在 2016 模型分析中被记载，非独立实验验证）。

## 最终分类

- BLOCKER：0。
- WATCH：5（BDGK、ACCOAC 可逆性、乙醛酸循环、NADTRHD、MDH/FUM 部分必需）。
- LOOP-DRIVEN FVA FREEDOM：ACCOAC/BIOC1/BIOC2、FBP/FBA/FBA3/TALA/SBPASE、PGMT/PGMT2。
- EXPLAINED/ACCEPTABLE：上述净零环、GAPD2 替代、BIOC1/BIOC2 替代（均为可解释化学计量）。
- DEFERRED：CYTRED、CA2tpp 等 Phase 1–3 已记录项。

## 结论

FIM/TTM/TSM 复现基线生长；HCO3E KO 三条件 infeasible；三条件 KO 无隐藏救援；全网络无重大跨子系统人工救援；loopless 消除/解释大循环 FVA；无 loopless 支持的 BDGK/ACCOAC/乙醛酸/氧化还原捷径；无 CO2/HCO3 组合恢复旧人工 bypass；无扰动相关自由/氧化还原捷径。所有剩余发现均归 WATCH / EXPLAINED / LOOP-DRIVEN / DEFERRED。

Phase 3.2 可正式关闭。葡萄糖工程启动时须将 5 项 WATCH 纳入敏感性分析。
