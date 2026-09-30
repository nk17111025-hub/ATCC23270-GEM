# PHASE 3.2 CENTRAL CARBON STRUCTURAL SANITY AUDIT

审计对象：`minimal_trustworthy_baseline_v1`（SHA256 `bd9715e6…3689af8c`）。只诊断、不修复。

## 结论

`PHASE3_2_PASS — CENTRAL CARBON STRUCTURALLY ACCEPTABLE FOR LOCAL GLUCOSE ENGINEERING`

未发现 BLOCKER。中央碳区域结构上可安全进入葡萄糖通路工程，但有 5 项 WATCH 需带入敏感性分析。

## 完整性

- v1 SHA256 校验通过，与冻结一致。
- 615 反应 / 573 代谢物 / 0 基因。
- 三个已批准 bound 确认：MACPD=0/0、ACOATA=0/1000、Htpp=0/0。
- 审计中央碳反应 60 个。

## 参考回归（未变）

FIM/TTM/TSM 生长 = 0.052076，HCO3E KO = infeasible，Phase 3.1 行为保持不变。

## Test A 闭碳循环筛选

三条件内部环均为净零 futile 循环，无净碳/能量/CO2 产出：

- FBP/FBA/FBA3/TALA/SBPASE：CBB 备用最优环（论文 Supp Note 4 已记载）。
- PGMT/PGMT2：磷酸甘油酸变位酶可逆平衡。
- ACCOAC/BIOC1/BIOC2：丙二酰辅酶 A 合成两条途径间的平衡（净零）。

分类：无害内部平衡，EXPLAINED / ACCEPTABLE。

## Test B 单反应必需性救援

必需（KO infeasible）：HCO3E、PGI1、PYK、MCOATA、TKT1、RUBISCO、PDH、ENO、TPI、PGK、RPE、RPI、PRUK。
生长归零：CS、ACS。
非必需（生长不变，均有已知替代）：PFK（自养走糖异生）、FBP/SBPASE（TALA/FBA3 环）、GAPD1（GAPD2+NADTRHD）、PPC（其它回补）、ACCOAC（BIOC1/BIOC2 生物素羧化酶）、G6PDH2、TALA、NADTRHD。
部分必需：MDH、FUM（KO 后生长 0.0098，经低置信度乙醛酸/苹果酸路线部分补偿）—— WATCH。

这些替代大多在 2016 论文 Supp Note 4 中有记载，非人工救援。

## Test C 可逆性 watchlist

中央碳可逆且需关注：
- BDGK（葡萄糖激酶，Confidence=1、无 GPR，可逆）：葡萄糖→G6P 入口，逆向等效 G6P 磷酸酶产 ATP—— WATCH。
- ACCOAC（乙酰辅酶 A 羧化酶，无 PMID，可逆）：逆向等效丙二酰辅酶 A 脱羧酶产 ATP—— WATCH。
- GLXCL / GLYCK / MALS / GLXR（乙醛酸循环，Confidence=1/无 GPR）—— WATCH。

## Test D CO2/HCO3 循环

CO2 消耗：DBTS、HCO3E、PPC、RUBISCO。HCO3 消耗：ACCOAC、BIOC2、CBPS。
人工 HCO3→CO2 bypass（原 MACPD/ACOATA 反向）已消除。合法的脂肪酸合成脱羧（3OAS）保留，但 HCO3E KO 仍 infeasible，故无新的 HCO3E 型救援。PASS。

## Test E 氧化还原捷径

free ATP=0、free NADH=0、free NADPH=0。NADTRHD（转氢酶，Confidence=1）提供 NADH↔NADPH 平衡，是 GAPD1→GAPD2 替代的关键—— WATCH。

## Test F 丙酮酸/乙酰辅酶A/OAA 闭合

PDH、PYK 必需；PPC 可被回补替代；MDH/FUM 部分必需（0.0098）；乙醛酸循环（MALS/GLXCL/GLYCK）是部分替代来源。WATCH。

## Test G FVA 休眠 bypass 容量

零参考通量但有大 FVA 范围的休眠容量：PFK（糖酵解，自养休眠）、GAPD2、BDGK（葡萄糖激酶 0..55.7）、BIOC1/BIOC2（0..1000）、MALS/GLXCL/GLYCK/GLXR（乙醛酸循环）、FDH、RUBISCOX、GLYCH。这些是引入葡萄糖/乙酸时才激活的休眠通路，非当前人工救援。

## Test H 低置信度中央碳反应

BDGK、HCO3tpp、MALS、MACPD（已禁用）、ACOATA（已单向化）、GLXCL、GLYCK、GLXR、NADTRHD 共 9 项 Confidence=1 或低证据。其中与葡萄糖/有机碳工程直接相关的：BDGK、乙醛酸循环（MALS/GLXCL/GLYCK/GLXR）、NADTRHD、ACCOAC。

## 最终分类

BLOCKER：0。
WATCH：5。

| ID | 反应 | 问题 | 影响 |
|---|---|---|---|
| CC-WATCH-001 | BDGK | 葡萄糖激酶 Confidence=1、无 GPR、可逆（逆向 G6P→葡萄糖+ATP） | 葡萄糖入口方向/热力学 |
| CC-WATCH-002 | ACCOAC | 可逆等效丙二酰辅酶 A 脱羧酶产 ATP | 乙酰辅酶A/丙二酰辅酶A界面 |
| CC-WATCH-003 | 乙醛酸循环 (MALS/GLXCL/GLYCK/GLXR) | 低证据，MDH/FUM 部分替代来源 | 乙酸/有机碳同化与 TCA 回补 |
| CC-WATCH-004 | NADTRHD | 转氢酶 Confidence=1、无 GPR | NADH↔NADPH 平衡 |
| CC-WATCH-005 | MDH/FUM 部分必需 | KO 后经低置信度路线仍 0.0098 生长 | 不完整 TCA 的可替代性 |

DEFERRED：CYTRED、CA2tpp、剩余电荷不平衡、TSM 定量偏高、508vs507、ETC 质子系数、biomass 组成（Phase 1–3 已记录，非本审计范围）。

## 是否可开始 Phase 4

是（NO BLOCKER）。但葡萄糖工程启动时需把 5 项 WATCH 纳入敏感性分析，尤其是 BDGK 的方向/证据与 ACCOAC 的可逆性。
