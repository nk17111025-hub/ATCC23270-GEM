# G53 专项证据复核：587–592、616、646、653–654

日期：2026-09-25。判断仅依据下列权威数据库、原始论文和本地 2016 原表；表达证据不等同精确催化证据。

## 裁决

| 基因 | 2016 表示及证据 | 动作 |
|---|---|---|
| 587 RU820_RS03040 / AFE_0630 | Amidotransferase domain-containing protein；未找出精确反应或 GPR | 暂缓，不猜反应 |
| 588–591 AFE_0631–0634 | iMC CYTBO3 已有公式、AND GPR（AFE_0631–0634）、confidence 3 | 不新增、不改 GPR；AFE_0632 注释冲突另记 |
| 592 RU820_RS03065 / AFE_0635? | PHFT 已有公式；GPR=(AFE_0635 or AFE_3143)。WP_012536251.1 同时映射 AFE_0635/3143 | 不改反应/GPR；映射 QC |
| 616 AFE_0660 | KEGG NADP malic enzyme 与 BioCyc NAD-malate dehydrogenase 冲突；iMC MDH GPR 是 AFE_3000 | 暂不新增、不加进 MDH GPR |
| 646 AFE_0692 | iMC FDH: for[c]+nad[c]→co2[c]+nadh[c]；GPR=(AFE_0692 or AFE_1652)，confidence 2 | 保留旧假设；无依据重写或增反应 |
| 653–654 AFE_0701/0702 | uptake hydrogenase 小/大亚基候选；iMC HYD1pp/HYD3pp 由其他基因关联 | 不并入旧 GPR；电子受体与耦联未定 |

## 证据细节

### 587–592 呼吸酶

587 的当前产品为 type 1 glutamine amidotransferase domain-containing protein；未见确定底物/产物或 2016 反应链接。

588–591 的 B1 完全 WP 映射将 RU820_RS03045/050/055/060 对应 AFE_0631–0634；旧注释为 cyoA/B/C/D，KEGG K02297/K02298 支持前两亚基。当前 RU820_RS03050 / AFE_0632 / WP_012536248.1 的 NCBI product 标为 cbb3-type cytochrome c oxidase subunit I，与旧 CyoB 和 K02298 不一致，作为注释冲突记录。

2016 CYTBO3 公式为 1.8 h[c] + q8h2[c] + 0.5 o2[c] -> 1.8 h[p] + h2o[c] + q8[c]，GPR=(AFE_0631 and AFE_0632 and AFE_0633 and AFE_0634)，confidence 3，bounds 0–1000。Campodonico 2016 用遗传算法拟合呼吸参数；1.8 H+/2e− 是模型拟合值，不是本株纯化 bo3 直接实测。Rhea:30251 为一般质子转位式；Rhea:30255 为非电生形式，不决定本株系数。

Quatrini 2009 虽讨论 AFE_0631–0634 cyoABCD，但对应 bo3 实验证据来自 ATCC 19859。Bellenberg 2019 的 DSM 14882T (= ATCC 23270) 蛋白组检出 AFE_0631/WP_012536247.1（pyrite/Fe log2=5.03,q=0.00）、AFE_0632/WP_012536248.1（2.56,0.01）、AFE_0633/WP_012536249.1（0.73,0.48）。前两者丰度显著上升，第三者检出但变化不显著；只支持表达，不证明复合体组装/独立催化/1.8 系数；AFE_0634 不在该表。这组反应在 iMC 已表示，现无可执行变更。

592 的 heme O synthase 反应 PHFT 公式为 h2o[c]+pheme[c]+frdp[c]→ppi[c]+hemeO[c]，GPR=(AFE_0635 or AFE_3143)，bounds 0–1000。B1 将 WP_012536251.1 同时对应 RU820_RS03065 与 RU820_RS14610，故不能专将 RU820_RS03065 接入某一个旧位点；记映射 QC。

### 616 AFE_0660 苹果酸酶候选

RefSeq 当前为 malic enzyme-like NAD(P)-binding protein；旧 AFE/BioCyc 标为 malate dehydrogenase。KEGG 对 afr:AFE_0660 关联 K00029、EC 1.1.1.40、R00216（NADP 苹果酸酶）；D3 的 BioCyc 字段却继承到 1.1.1.39-RXN/MALATE-DEH-RXN 和 NAD/NADP 多种 Rhea 反应，并标注 EC 冲突。苹果酸酶式为 (S)-malate + NAD(P)+ -> pyruvate + CO2 + NAD(P)H；MDH 式为 (S)-malate + NAD+ <=> oxaloacetate + NADH + H+。Rhea:12653、18253、21432 分别显示这些化学差别，不能确定 AFE_0660 的辅酶或产物。

2016 MDH 的式为 nad[c]+mal-L[c]<=>h[c]+nadh[c]+oaa[c]，GPR=AFE_3000，不是 AFE_0660。Bellenberg 同株补表检测到 AFE_0660/WP_009560818.1，log2=1.44、q=0.09；未达 q<0.05，也未测具体反应或酶活。因此不新增反应、不接入 MDH GPR。

### 646 AFE_0692

WP_012536289.1 注释为 FdhF/YdeP family oxidoreductase / molybdopterin oxidoreductase alpha。Valdés 2008 将 AFE_0690–0692 作为潜在 FDH cluster，并推测可能与 group 4 hydrogenase 成复合体，属基因组推论。Bellenberg 2019 同株补表列 AFE_0692（log2=-0.35,q=0.29），仅为表达/丰度记录。

2016 FDH 已有 for+NAD→CO2+NADH 与 OR GPR，score 2；无本株位点专属实验证明 NAD 为受体或产品。保留旧建模假设，不升证据、不另造反应；需专属电子受体和产物测定后重审。

### 653–654 AFE_0701/0702

WP 完全匹配对应 RU820_RS03380/WP_012536297.1→AFE_0701、RU820_RS03385/WP_012536298.1→AFE_0702。KEGG K23548/K23549 为 uptake hydrogenase 小/大亚基，EC 1.12.99.6。D3 的 R08034 链接不是氢化酶反应；KEGG 该号实际为二氨基硝基甲苯转化。Rhea:12116 只有 H2 + acceptor -> reduced acceptor，未指定受体。

Drobner 1990 直接证明 ATCC 23270 可用 H2 生长且 H2 诱导 hydrogenase，但没归因 AFE_0701/0702。Valdés 2008 按基因组把它们列作 group 2 uptake [NiFe] hydrogenase。Bellenberg 2019 同株补表列 AFE_0701（log2=-0.88,q=0.05）和 AFE_0702（-1.01,0.00），仅支持蛋白检出/丰度变化。Kucera 2020 用的是 CCM 4253（异株）；报告 group 2 Hup 及 H2 氧化多组学，Hup 到醌池的 Fe-S 接力仍为推测。

iMC 的 HYD1pp/HYD3pp 已表达 H2→Q8H2，并各有 2 H+ 胞质至周质转位；GPR 分别为 AFE_3283–3286 和 AFE_0937–0940，未含 AFE_0701/0702。HYD4 是 AFE_2149–2154 关联的甲酸→CO2+H2。总体 H2 通量有表示，但 Hup 的具体受体、Fe-S 伙伴和质子耦联未证，不能直接并入旧 HYD GPR或另造化学式。暂不改模型。

## DOI 与原始链接

- Campodonico 2016 iMC507：<https://doi.org/10.1016/j.meteno.2016.03.003>。本地原始 mmc1.xls、mmc3.xml；对照 A0_复现2016与2024模型/A0_主任务/2016_到2024_真实差异.tsv。
- Valdés 2008：<https://doi.org/10.1186/1471-2164-9-597>；<https://pmc.ncbi.nlm.nih.gov/articles/2621215/>。
- Quatrini 2009：<https://doi.org/10.1186/1471-2164-10-394>。
- Brasseur 2004：<https://doi.org/10.1016/j.bbabio.2004.02.008>。
- Bellenberg 2019，DSM 14882T (= ATCC 23270)：<https://doi.org/10.3389/fmicb.2019.00592>；补表 <https://www.frontiersin.org/api/v4/articles/434788/file/Data_Sheet_2.pdf/434788_supplementary-materials_datasheets_2_pdf/2>。
- Drobner 1990：<https://doi.org/10.1128/aem.56.9.2922-2923.1990>；<https://pmc.ncbi.nlm.nih.gov/articles/PMC184866/>。
- Kucera 2020，CCM 4253：<https://doi.org/10.3389/fmicb.2020.610836>；<https://pmc.ncbi.nlm.nih.gov/articles/PMC7735108/>。
- KEGG AFE_0660 <https://www.kegg.jp/entry/afr:AFE_0660>；Rhea <https://www.rhea-db.org/rhea/12653>, <https://www.rhea-db.org/rhea/18253>, <https://www.rhea-db.org/rhea/21432>, <https://www.rhea-db.org/rhea/30251>, <https://www.rhea-db.org/rhea/30255>, <https://www.rhea-db.org/rhea/12116>；KEGG R08034 <https://www.kegg.jp/entry/R08034>。
- 当前组装 <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/>。
