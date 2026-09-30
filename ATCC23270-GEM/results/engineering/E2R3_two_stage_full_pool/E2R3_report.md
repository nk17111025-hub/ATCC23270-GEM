# E2R3 两阶段全池筛选

- 权威宿主 SHA256：026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3
- E1 pool SHA256：23b45be4333e32f2b6535813d9bb017ef019afa417f3f494267145ad1b6df0d6；候选行：2424；构建后方向变体：2427。
- WT glucose-closed FIM：B_WT=0.0520763871027；pFBA Fe2=163.781528744、O2=38.4815426752。
- 资源匹配混合零 KI：biomass=0.0911667122963；glucose=-0.5；R_REF=0.791344957513。
- 两阶段共同 biomass 下限=0.0781145806541；Rubisco cap：Stage 1=0.237403487254（30%），Stage 2=0.0791344957513（10%）。
- Fe2 bounds=[-163.781528744,0] 且 uptake≥0.01；O2 bounds=[-38.4815426752,0]；glucose bounds=[-0.5,0] 且 uptake≥0.05；CO2/H2CO3=[-2,0]。
- Stage 1 全池 FBA：optimal，biomass=0.13435120673715378；pFBA：optimal，biomass=0.07811458065411052。
- FVA 分类（阈值容差 1e-08）：active_in_pfba=114, inactive_stage1=1218, optional_stage1=1095
- Stage-1 pFBA 保留 biomass 下限和全部其余表型约束，以零主目标最小化总通量；它是一个简约通量代表。FVA mandatory 表示区间不含零。对 FVA 可活动范围最大的 32 个 optional variants，逐个强制方向一致的非零 flux 后再做阈值条件 pFBA，以发现需要多个共同反应的备选路线。
- 每个支持集经 Stage 1 独立 FBA 重建，并逐反应删除验证不可再单删；该性质不代表全局 KI 数最小。

- reaction_count 与优先排序按 donor-supported search reaction variants 计数；E1_candidate_count 按去重后的 E1 candidate_id 计数，两个数值均保存在支持集表。
## Stage-1 validated irreducible support sets

- **S01**（3 search variants；low_pfba_flux_first）：E1R01016; E1R01493; E1R01642。B=0.084562227，glucose=-0.5，Rubisco=0.23740349，Fe2=-137.81877，O2=-33.453882；glucose-off B=0.006201070634220194，ΔB=0.078361157。
- **S02**（4 search variants；high_pfba_flux_first）：E1R00037; E1R00038; E1R01343; E1R01344。B=0.079719114，glucose=-0.5，Rubisco=0.23740349，Fe2=-125.1971，O2=-30.603128；glucose-off B=0.00584591807757693，ΔB=0.073873196。
- **S03**（2 search variants；candidate_id_descending）：E1R00732; E1R00741。B=0.08519834，glucose=-0.5，Rubisco=0.23740349，Fe2=-142.57092，O2=-34.611824；glucose-off B=0.006247717662623617，ΔB=0.078950622。

## Stage 2: same-set 10% challenge

- S01：PASS_BOTH；status=optimal；B=0.07950623572120567，glucose=-0.5，Rubisco=0.07913449575127987，glucose-off B=0.002067023544739748，ΔB=0.07743921217646592。Stage 2 使用同一反应集合，未增加 KI。
- S02：PASS_30_ONLY；status=infeasible；B=None，glucose=None，Rubisco=None，glucose-off B=None，ΔB=None。Stage 2 使用同一反应集合，未增加 KI。
- S03：PASS_BOTH；status=optimal；B=0.080089852603804，glucose=-0.5，Rubisco=0.07913449575127987，glucose-off B=0.0020825725542079846，ΔB=0.078007280049596。Stage 2 使用同一反应集合，未增加 KI。

## 优先级

- A_PASS_BOTH / S03（2 reactions）：Stage1 B=0.08519834 / Rubisco=0.23740349；Stage2 B=0.080089852603804 / Rubisco=0.07913449575127987；Stage1 glucose contribution=True。
- A_PASS_BOTH / S01（3 reactions）：Stage1 B=0.084562227 / Rubisco=0.23740349；Stage2 B=0.07950623572120567 / Rubisco=0.07913449575127987；Stage1 glucose contribution=True。
- B_PASS_30_ONLY / S02（4 reactions）：Stage1 B=0.079719114 / Rubisco=0.23740349；Stage2 B=None / Rubisco=None；Stage1 glucose contribution=True。

## Reproducibility
- Input fingerprint：`6a47044c69e9408352200ba87b45cf8a4e526822a2ecd4ab68c83f3f7940e52c`
- Donor model SHA256：`{'iML1515': 'b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5', 'iJN1463': 'd573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1', 'iCN1361': 'c4e418d831f428a93702f2d1b6c62d6fe3538545474f12397640028169d4ce0c'}`
- Package versions：`{"cobra": "0.32.1", "optlang": "1.9.1", "pyscipopt": "6.2.1", "scipy": "1.18.1", "straindesign": "1.19.1", "swiglpk": "5.0.13"}`
- pFBA/FVA 设置：WT: maximize biomass then COBRApy pFBA at exact optimum. Stage1: retain explicit phenotype constraints including biomass >=1.5xB_WT; objective=0 and COBRApy add_pfba(fraction_of_optimum=0) minimize sum of forward/reverse flux without forcing maximum biomass.；COBRApy flux_variability_analysis; fraction_of_optimum=0; exact explicit Stage1 phenotype constraints; loopless=None。
- FVA 使用精确 Stage-1 表型约束（biomass 下限、Rubisco 上限、glucose 最低摄取、Fe2 trace、Fe2/O2 caps、CO2 bounds），fraction_of_optimum=0。
- 结果、输入指纹检查点及单独 E2R3 日志均保存在本目录；未触碰 E1/E2/E2R 结果。
