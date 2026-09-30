# ATCC23270 混合营养工程项目交接文档

> 用途：把当前聊天中的项目背景、固定文献编号、已形成的方法学判断和下一步技术路线交给新的聊天继续讨论。  
> 研究对象：**Acidithiobacillus ferrooxidans ATCC 23270**。  
> 核心目标：在尽量保持 **Fe²⁺/RISC 氧化与生物浸矿能力** 的前提下，让外源有机碳逐渐承担更多生物量碳来源，降低对 CO₂ 固定的依赖，构建稳定的混合营养状态。

## 1. 研究问题的准确表述

不要再把项目简单定义成“gap filling”。更准确的问题是：

**在 ATCC23270 的基因组尺度代谢网络中，寻找最少的反应添加、过表达（OE）、下调（KD）或敲除（KO），使有机碳能够稳定贡献生物量，同时保持 Fe²⁺/RISC 氧化、生长和能量代谢。**

更接近：
- phenotype-constrained network augmentation
- metabolic network rewiring
- strain design for trophic rewiring
- mixotrophy engineering

理想模型问题：
> 在 ATCC23270 的 iMC507/更新模型中，当 Fe²⁺/RISC oxidation 必须保持、CO₂供给逐渐降低、有机碳摄取逐渐提高时，哪一组 transport / central-carbon / redox / respiratory changes 能维持 biomass 与 bioleaching？

输出最终应分类为：
1. 异源反应添加
2. 本地反应 OE
3. 本地反应 KD
4. KO
5. 调控瓶颈
6. 酶容量瓶颈

关键警告：
**基因存在 ≠ 反应真实存在 ≠ 酶容量足够 ≠ 生理上真的有通量。**
**有机碳 uptake > 0 ≠ 有机碳直接进入 biomass。**它可能先被氧化为 CO₂，再通过 CBB 重新固定。

## 2. 当前最大的科学问题

ATCC23270 已有一些有机碳利用迹象，但远没有形成强异养代谢。真正的问题可能同时包括：
- 有机碳转运能力不足；
- PFK 等关键酶容量弱；
- PPP/EMP/pyruvate/TCA 连接不足；
- TCA 不是典型完整异养型循环；
- CBB 与 PPP/中央碳共享部分反应；
- NADH/NADPH/醌池/质子动力势重新分配；
- Fe/S 电子传递与有机碳氧化之间可能竞争；
- 有机碳加入后可能发生调控互斥；
- 普通 FBA 会把“存在但很弱的反应”误当作可以承担高通量。

因此，不能只做传统 gap filling。

## 3. 固定文献编号

### A 系列：ATCC23270 本体生理、CO₂固定、Fe/S 能量代谢
- A1 Barreto et al. 2003 — partial genome / predicted physiology
- A2 Appia-Ayme et al. 2006 — Fe²⁺ vs S⁰ 条件下碳代谢转录组
- A3 Quatrini et al. 2006 — Fe/S 电子传递模型
- A4 Valdés et al. 2008 — 全基因组代谢蓝图
- A5 Quatrini et al. 2009 — 扩展 Fe/S 氧化模型
- A6 Esparza et al. 2010 — CO₂ fixation / cbb 操纵子
- A7 Wang et al. 2012 — markerless gene replacement + pfkB mutant
- A8 Osorio et al. 2013 — anaerobic sulfur metabolism + Fe(III) reduction
- A9 Esparza et al. 2019 — CO₂浓度响应
- A10 Wang et al. 2024 — electroautotrophy vs chemoautotrophy transcriptome/metabolome

A 系列关键点：
- Fe²⁺氧化有 downhill / uphill 电子流；downhill 主要供 ATP，uphill 代价高、主要供还原力。
- Fe 与 S 使用不同但耦联的 ETC 模块。
- CBB 不是简单开关，多个 cbb operons、Rubisco form I/II、carboxysome 有条件依赖。
- A7 很关键：`pfkB/AFE1807` 有 PFK 活性但很弱；ΔpfkB 影响 S⁰ 生长而对 Fe²⁺氧化影响较小；加葡萄糖能提高 WT 生物量但葡萄糖消耗有限。
- A2/A3/A10 共同说明：碳代谢和电子传递具有强状态依赖调控，未来必须考虑“有机碳加入后是否把 Fe/S 能量模块压下去”。

### B 系列：ATCC23270 代谢模型

**B1 Hold et al. 2009**  
小型 stoichiometric/MFA 模型，不是 GEM。价值是估计 Fe²⁺电子分配：绝大多数 downhill，小部分 uphill，但 uphill 能量成本很大。

**B2 Campodonico et al. 2016 — iMC507**
- 507 genes
- 587 reactions
- 573 metabolites
- 3 compartments
- FBA + chemiosmotic constraints + GA 校准
- 对 Fe²⁺、thiosulfate、tetrathionate 等已有生理数据拟合
- 约 7 条拟合曲线平均 R²≈0.92

注意：这不是“模型准确率92%”。它主要是校准条件下宏观表型拟合较好，不代表内部 flux 正确，也不代表有机碳条件可靠。

iMC507 强项：
- Fe/S 生物能量学
- 极低 pH 场景
- 电子传递
- chemolithoautotrophy

弱项：
- 有机碳转运
- 中央碳高通量
- 酶容量
- 调控
- 混合营养

**B3 Khaleque et al. 2024**  
基于更新模型做 compatible-solute/ectoine 工程。价值是“在已有 GEM 上加合成路径并比较多个环境场景”的范式；工程预测仍主要是 in silico。

## 4. C 系列：外部菌株给项目的方法学启发

### C1 Hommes et al. 2003 — Nitrosomonas europaea
传统上是专性化能无机自养：NH₃供能、CO₂供碳。用标记果糖证明在低/无外源无机碳时，新生物量碳可以大量来自果糖，但抑制氨氧化后即使有果糖也不能正常生长。  
**借鉴点：**证明“无机物继续供能、有机物开始供碳”的目标表型本身可行。最终 ATCC23270 必须同时验证 Fe/S 氧化仍在、有机碳进入 biomass、CO₂固定需求下降。

### C2 Kappler & Nouwens 2013 — Starkeya novella
天然就是灵活的异养 + 甲基营养 + 硫化能无机自养/混合营养菌。具有 PPP+ED、完整 TCA、glyoxylate shunt、Sox、CBB。多营养条件下 Sox、甲醇代谢和 CBB 蛋白可以长期保留。  
**借鉴点：**不要预设 CBB 必须降到 0，更合理的是找 Fe/S oxidation + organic carbon assimilation + residual CBB 的合适通量区域。

### C3 Quasem et al. 2017 — Hydrogenovibrio crunogenus
严格硫化能无机自养。多种有机物不能支持有效异养；关键中央代谢酶活和基因注释不完全一致，TCA 更偏 biosynthetic。  
**借鉴点：**这是“为什么只看 GEM 会被骗”的反例。每个关键候选反应都要检查 substrate specificity、direction、cofactor、enzyme capacity、pathway redundancy。

### C4 Jahn et al. 2021 — Cupriavidus necator H16
天然兼性化能无机自养 + 异养。作者做多条件 chemostat、多生长速率、LC-MS/MS 定量蛋白质组、RBA、约 6 万条码转座子突变体竞争实验。

RBA 比 FBA 多的核心约束：
`v_j <= k_app,j × E_j`

即反应通量必须由足够的酶量支持，还加入总蛋白池、膜蛋白池、核糖体、RNA polymerase、DNA polymerase、chaperone 等资源成本。

`k_app` 大致来自：
`reaction flux / measured enzyme abundance`

但并非所有参数都自己测；部分机器参数、总蛋白比例等来自文献，测不到的 k_app 也有缺省处理。

**借鉴点：**iMC507 最大的问题之一是“反应存在就可以跑很大通量”。项目早期不必复制完整 RBA，可以先给少数关键反应加局部酶容量约束：
- organic-carbon transporter
- PFK
- PPP
- pyruvate node
- CBB
- Complex I
- bc₁
- ATP synthase

### C5 Lawson et al. 2021 — Nitrospira moscoviensis
核心是 GEM + ^13C tracer。模型预测甲酸可经还原型甘氨酸途径直接同化，但 ^13C-formate 实验显示真实细胞主要是 formate → CO₂ → rTCA 重新固定。  
**借鉴点：**模型可行路线 ≠ 真实生理路线。未来 ATCC23270 必须区分：
1. organic carbon 直接进入 biomass；
2. organic carbon 先被氧化成 CO₂，再通过 CBB 固定。
理想验证：^13C-organic carbon、胞内代谢物标记、biomass labeling、^13CO₂。

### C6 Canto-Encalada et al. 2022 — Nitrosomonas europaea
这个菌本来就有实验报道可利用一定果糖/丙酮酸。C6 不是预测“怎样把不会吃有机碳的菌改出来”，而是建立 GEM 后扫描有机碳增加时，网络何时从自养型转向有机碳主导。

iGC535：
- 535 genes
- 1149 reactions
- 1114 metabolites

模型构建：
本菌蛋白 → 与 iML1515 / iHN637 / iPC815 做同源匹配 → 迁移候选 reaction/GPR → 得到含模板基因的初稿 → 人工清洗 GPR → 用本菌基因替代或删除模板基因 → 按 N. europaea 文献补氨氧化/CBB/ETC 等物种特异反应 → gap filling / consistency → iGC535。

**最大借鉴点：**做 `organic-carbon uptake × CO₂ availability` 二维扫描，每个格子记录：
- max biomass
- minimum required Fe²⁺/RISC oxidation
- CBB flux
- PPP/EMP/TCA
- NAD(P)H source
- Complex I
- bc₁
- ATP synthase
- CO₂ secretion/fixation
- FVA range

目标是划分：
- 纯自养区
- 混合营养区
- organic-carbon dominant 区
- 模型外推危险区

C6 是目前最适合直接移植的“干实验骨架”。

### C7 Pu et al. 2026 — Cupriavidus necator
本来既能 fructose heterotrophy，也能 H₂+CO₂ autotrophy，但 WT 不能很好同时运行。作者发现有机碳条件会压低 hydrogenase 模块；通过 HoxA 等调控让 H₂利用系统在有机碳存在时保持开启，实现更真正的 mixotrophy。  
**借鉴点：**即使 GEM 说两套通路可以同时运行，真实细胞也可能因转录调控而互斥。ATCC23270 后续需要检查有机碳加入以后是否压低：
- Fe/S ETC
- CBB
- reverse electron transport
- respiratory modules
如果存在类似 HoxA 的上游调控节点，可能需要 regulatory engineering，而不是只加代谢酶。

### C 系列一句话分工
- C1：证明目标生理状态可以存在。
- C2：说明 mixotrophy 不一定需要关掉 CBB。
- C3：警告 reaction annotation ≠ 真实酶功能。
- C4：给 FBA 增加酶容量和蛋白资源现实性。
- C5：用 ^13C 检验模型预测碳流真假。
- C6：用 GEM + FBA/FVA 扫描自养→混合营养→有机碳主导。
- C7：解决“模型允许，但真实调控不允许”的问题。

最有用的组合：
**C6 找状态 → C3 审反应 → C4 查容量 → C7 查调控 → C5 做 ^13C 强验证。**

## 5. D 系列：ATCC23270 工程工具

D 系列主要回答“模型给出靶点以后，现实里怎么改 ATCC23270”。

**D1 Kernan et al. 2016**  
证明 ATCC23270 可以表达异源代谢酶、加入人工代谢通路。  
适用：宿主缺失关键反应时做异源表达。

**D2 Inaba et al. 2018**  
Tn5 介导染色体整合。  
适用：需要稳定加入外源 transporter/enzyme/多基因模块。

**D3 Jung et al. 2024**  
CRISPR/dCas12a knock-down bc₁ complex。证明 ATCC23270 可以做 CRISPRi/KD，并且能量代谢支路下调会造成真实 Fe/S 和浸矿表型变化。  
适用：模型给出 MUST↓ / KD 靶点。

**D4 Jung et al. 2025**  
SQR overexpression。证明本地能量代谢反应可 OE，并能改变硫氧化和矿物浸出。  
适用：模型给出 MUST↑ / OE 靶点。

映射关系：
- 缺失反应 → D1 异源表达
- 稳定加入外源功能 → D2 染色体整合
- 本地反应需降低 → D3 CRISPRi/KD
- 本地反应需提高 → D4 OE
- Fe/S 电子流重分配 → D3/D4

## 6. 目前最重要的新判断：不要局限于 gap filling

用户数据很少，不适合依赖大样本监督机器学习，但这不等于没有算法。这个课题反而非常适合：
- FBA
- FVA
- MILP
- pathway addition
- OptStrain
- OptForce
- constrained Minimal Cut Sets（cMCS）
- OptCouple / StrainDesign
- TFA
- flux sampling
- 局部 enzyme-capacity constraint

这些方法本来就不依赖大量训练数据。

## 7. 推荐的核心算法框架

### 第一层：宿主模型整理
以 iMC507/更新模型为底盘，重点人工重审：
- organic carbon transport
- glycolysis/EMP
- PPP
- pyruvate node
- acetyl-CoA
- TCA
- anaplerosis
- CBB ↔ PPP shared reactions
- NADH/NADPH
- quinone pool
- Complex I
- bc₁
- terminal oxidase
- ATP synthase

检查：
- GPR
- reaction direction
- cofactor
- compartment
- substrate specificity
- evidence level

### 第二层：构建候选 reaction universe
来源可考虑：
- Rhea
- MetaCyc
- KEGG
- BiGG
- 近缘酸硫杆菌/化能无机营养菌
- 已知能利用目标有机碳的细菌

优先筛：
- transporter
- first-step activation enzymes
- PPP/EMP bridge
- pyruvate/acetyl-CoA bridge
- TCA bypass
- cofactor-balancing reactions

### 第三层：OptStrain 式最小反应添加
目标是寻找满足 phenotype 所需的最少新增反应。

约束示例：
- `growth >= 0.3~0.5 × WT`
- `Fe/RISC oxidation >= beta × WT`
- `organic carbon uptake > 0`
- `CO2 fixation <= delta × WT`

最小化：
- `number of added reactions`

最好做 multi-scenario robust design，要求同一套改造在：
- Fe²⁺ background
- S / tetrathionate / thiosulfate background
- low/medium/high CO₂
- low/medium organic carbon
- 30/50/70/90% growth requirement
中都成立。

### 第四层：C6 式状态扫描
对每个候选 engineered network 做：
**organic-carbon uptake × CO₂ availability 二维扫描**

输出：
- biomass
- Fe/S oxidation
- CBB
- PPP/EMP
- TCA
- NADH/NADPH
- Complex I
- bc₁
- ATP synthase
- CO₂ secretion

配合 FVA 区分：
- 必须通量
- 可选通量
- 多解状态

### 第五层：OptForce 找 OE/KD/KO
比较 WT 与目标 mixed state 的 flux ranges：
- target minimum > WT maximum → MUST↑ → OE
- target maximum < WT minimum → MUST↓ → KD
- target flux forced to zero → KO

### 第六层：cMCS 排除不想要的表型
最怕：
**organic carbon 让模型长，但 Fe/S oxidation 掉到接近 0。**

因此定义：
undesired region，例如：
- `growth > 0`
- `organic carbon uptake > 0`
- `Fe/S oxidation < 0.2 × WT`

desired region，例如：
- `growth >= 0.3 × WT`
- `Fe/S oxidation >= 0.7 × WT`
- `organic carbon assimilation > 0`

用 constrained Minimal Cut Sets 找最少 intervention，使坏状态不可行，同时保留好状态。

### 第七层：TFA 热力学过滤
对最终少量候选检查：
- reaction thermodynamic directionality
- ΔG feasibility
- metabolite concentration constraints（若有）

ATCC23270 极低 pH，因此要特别注意：
- pH
- proton balance
- ionic strength
- membrane potential

### 第八层：C4 式局部 enzyme-capacity constraints
目前不建议直接做 full RBA，因为缺：
- 多条件定量 proteomics
- 本菌 k_app
- 全局 protein allocation data

但可以先给关键反应加：
`v_j <= k_app,j × E_j`

优先：
- sugar transporter
- PFK
- PPP enzymes
- pyruvate node
- CBB
- Complex I
- bc₁
- ATP synthase

A7 的 PFK 酶活数据尤其值得利用。

### 第九层：C7 式调控审查
模型筛出候选后，检查有机碳加入是否转录下调：
- petI/petII
- cyc2
- rus
- bc₁
- nuo
- terminal oxidases
- cbb genes
- organic carbon transporter
- PPP/EMP regulators

如果存在上游 master regulator，应考虑：
**regulatory engineering > 单个酶 OE**

### 第十层：C5 式湿实验验证
最终不能只测 growth。

理想最小验证：
1. growth
2. organic carbon consumption
3. Fe²⁺/RISC oxidation
4. CO₂ uptake/release
5. ^13C-organic-carbon tracing

最好进一步：
- intracellular metabolomics
- biomass labeling
- ^13CO₂

真正回答：
**organic carbon 是直接进入 biomass，还是先变 CO₂ 再被 CBB 固定。**

## 8. 当前最推荐的“论文级”技术路线

**模型人工审查 → 候选反应库 → OptStrain式最小反应添加 → 多场景鲁棒设计 → C6式 organic-C × CO₂ 状态扫描 → OptForce 找 MUST↑/↓/0 → cMCS 排除丢失 Fe/S 氧化的坏状态 → FVA/flux sampling 检查鲁棒性 → TFA 过滤 → 局部酶容量约束 → 调控层检查 → 湿实验 ^13C 验证。**

核心创新点不是“又做一次 gap filling”，而是：

> **在保证 Fe/S 能量代谢持续存在的条件下，寻找使 ATCC23270 获得稳定有机碳参与生物量合成能力的最小网络重布线方案。**

## 9. 下一步优先深入的算法

下一聊天优先系统比较：
1. OptStrain
2. OptForce
3. constrained Minimal Cut Sets（cMCS）
4. OptCouple
5. StrainDesign
6. TFA / pyTFA
7. RobustKnock
8. FastKnock
9. Growth Coupling Suite
10. FluxRETAP
11. gDel_minRN
12. DBgDel
13. 是否存在专门针对 trophic switch / mixotrophy engineering / network augmentation 的新算法
14. reaction addition + regulation intervention 的联合 MILP
15. atom mapping / carbon-source tracing 的计算方法

重点比较：
- 是否支持 reaction insertion + KO/OE/KD 联合优化
- 是否可强制保留 Fe/S oxidation
- 是否支持 multi-scenario
- 是否支持 Python
- 是否需要 Gurobi/CPLEX
- 是否能接 COBRApy/iMC507
- 是否依赖实验 omics
- 是否有活跃代码仓库

## 10. 当前研究策略优先级

第一优先：
- GEM curation
- reaction-universe construction
- OptStrain / reaction addition
- FBA/FVA
- C6-style condition scanning
- OptForce
- cMCS

第二优先：
- TFA
- flux sampling
- local enzyme-capacity constraints

第三优先：
- full RBA
- ML

原因：
**当前瓶颈不是训练数据不足，而是物理/生化约束和网络结构。**

## 11. 用户偏好

新聊天请注意：
- 不喜欢每两句话就分段；
- 尽量使用连续、紧凑的中文；
- 英文只在首次解释术语时括号补充；
- 技术路线导向，不要泛泛重复论文背景；
- 对论文明确区分 prediction / calibration / independent validation / wet-lab confirmation；
- 2018年前老论文一般抓主结论，不必逐图细读；
- 近年方法论文可以细讲；
- 最关心“这篇对我的模型能直接借什么”；
- 不要把“模型能跑”当成“生理真实”；
- 对算法优先讲输入、输出、数据要求、优化类型、软件生态和能否落地。

## 12. 下一聊天可直接从这里继续

建议直接提问：

**请以代谢工程/合成生物学技术路线设计为目标，系统比较 OptStrain、OptForce、cMCS、OptCouple、StrainDesign、TFA 这几类方法，重点判断它们是否能支持 ATCC23270 的 reaction addition + OE/KD/KO + Fe/S retention + multi-scenario design，并给出一个真正能落地的算法组合。**

也可以进一步要求：

**不要只总结算法，要把每个算法对应到我的变量、约束、目标函数和输入输出。**
