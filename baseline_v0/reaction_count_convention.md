# 2016 iMC507 反应计数口径

## 结论

615（SBML / mmc1.xls 总反应数）= 587（非 exchange 反应）+ 28（exchange 反应）。

587 的构成：534 代谢反应 + 52 转运反应 + 1 生物量目标函数。

## 精确计数（来源：mmc1.xls Table 1，共 615 数据行）

| 类别 | 数量 | 说明 |
|---|---:|---|
| exchange（Ex_*） | 28 | 形如 `metabolite[e] <=>` 的边界反应 |
| transport | 52 | Subsystem 含 `Transport`（含 tex/tpp） |
| biomass objective | 1 | `Afe_biomass_mc507_WT_139p0M`（Subsystem 标为 Exchange） |
| metabolic（其余酶促反应） | 534 | 非 exchange、非 transport、非 biomass |
| 合计 | 615 | |

非 exchange 反应数 = 52 + 1 + 534 = 587。

## 28 个 exchange 反应

Ex_4hba[e], Ex_bio[e], Ex_ca2[e], Ex_cu2[e], Ex_eps_AFE[e], Ex_fe2[e], Ex_fe3[e],
Ex_h2[e], Ex_h2co3[e], Ex_h2o[e], Ex_h2s[e], Ex_h[e], Ex_k[e], Ex_mg2[e], Ex_mn2[e],
Ex_n2[e], Ex_na1[e], Ex_nh4[e], Ex_o2[e], Ex_pi[e], Ex_polypi[e], Ex_ppt[e], Ex_s[e],
Ex_so4[e], Ex_so4aa[e], Ex_tsul[e], Ex_ttton[e], Ex_zn2[e]

## 论文原文措辞

摘要：
> A total of 587 metabolic and transport/exchange reactions, 507 genes and 573
> metabolites organized in over 42 subsystems were incorporated into the model.

方法（第 2.1 节附近）：
> contained 507 genes, 587 metabolic and transport reactions, and 573 nonunique
> metabolites, which were distributed over 42 subsystems and three different
> cellular compartments: extracellular, periplasm and cytoplasm.

## 说明

- 论文正文明确写 `587 metabolic and transport reactions`，即 587 个非 exchange 反应。
- SBML（mmc3.xml）和 mmc1.xls Table 1 额外包含 28 个 exchange 边界反应，因此总数为 615。
- mmc1.xls Table 1 中 `Afe_biomass_mc507_WT_139p0M` 被归入 `Exchange` 子系统，但它不是
  `Ex_` 边界反应，只是生物量目标函数；因此在按子系统统计时 `Exchange` 显示 29 条
  （28 个 Ex_ + 1 个 biomass），而真正的 exchange 反应是 28 个。
- 此差异只是计数口径，模型无需修改。
