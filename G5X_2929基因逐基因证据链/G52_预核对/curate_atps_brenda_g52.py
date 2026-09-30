import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE / name).read_text(encoding='utf-8'))

def write(name, value):
    (HERE / name).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')

manual = read('g52_manual_adjudications.json')
manual['RU820_RS02615'] = {
    'issue': 'Jaramillo et al. 2012 克隆 ATCC 23270 atpS，在大肠杆菌表达物中测得 APS + PPi 生成 ATP；FM177944.1 与现行 RU820_RS02615 CDS 同为 1674 nt，差 7 nt/4 aa，蛋白一致率 99.2819%。实测为同株历史等位序列，现行 WP 的催化身份需保留版本差异。APS kinase 结构域仅预测，未测其活性；SFAT2 未测。',
    'recommendation': 'SFAT1 的 APS + PPi -> ATP + sulfate 核心反应有同株近同序列重组表达实验证据；2016 原式保留历史分数，2026 对核心反应可评 4，质子系数、细胞区室和网络方向待独立核对。ADSK 保留注释候选；Hold SFAT2 的该基因 GPR。',
    'chemistry': 'SFAT1: APS + PPi -> ATP + sulfate，2016 完整式 ppi[c] + aps[c] -> h[c] + atp[c] + so4[c] 已按旧补表配平；实验仅测 ATP 形成。ADSK: ATP + APS -> ADP + PAPS 未实测。SFAT2: APS + Pi -> ADP + sulfate 未实测。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S027;G52-S028'
}
manual['RU820_RS01775'] = {
    'issue': 'BRENDA 的 EC 1.8.1.7 条目指向 Sugio et al. 1995 的 glutathione reductase 纯化，实验菌株为 AP19-3；现行 G52 RU820_RS01775 是 ATCC 23270，不能跨菌株转为此 WP 的直接催化证据。',
    'recommendation': '保留当前注释和数据库支持的 glutathione reductase 候选；暂不提升为同株直接实验分数。',
    'chemistry': 'GSSG + NADPH + H+ -> 2 GSH + NADP+；当前 WP 的精确化学和区室按既有数据库核对。',
    'evidence': 'G52-S001;G52-S006;G52-S007;G52-S027'
}
manual['RU820_RS02045'] = {
    'issue': '2016 GLUTRS 使用 tRNA(Glu)；KEGG 将历史 AFE_0422 标作 gltX-2。Salazar et al. 2003 及 Núñez et al. 2004 的同物种重组 GluRS2 优先谷氨酰化 tRNA(Gln)-UUG，GluRS1 可作用于 tRNA(Glu)。原文早于完整基因注释，尚未找到可与现行 WP_012536128.1 精确比对的克隆序列或引物，因而该实验与本 WP 的身份关联仍须复核。',
    'recommendation': 'Hold GLUTRS 的 AFE_0422/RU820_RS02045 OR 分支；先核对旧 GluRS2 与当前 WP 的序列，再决定 tRNA(Glu) 或 tRNA(Gln) 的反应归属。旧反应的另两条 GPR 分支由所在组独立审查。',
    'chemistry': '2016 GLUTRS: L-glutamate + ATP + tRNA(Glu) -> L-glutamyl-tRNA(Glu) + AMP + PPi；GluRS2 实验主要指向 tRNA(Gln)-UUG，同一底物类型不能直接替代。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S027;G52-S029'
}
manual['RU820_RS02255'] = {
    'issue': 'BRENDA EC 6.3.1.2 引用 Barros et al. 1986 早期 Thiobacillus ferrooxidans glutamine synthetase 克隆纯化实验；原文未提供能与现行 RU820_RS02255/WP 精确核对的序列身份或 ATCC 23270 菌株证明。',
    'recommendation': '保留 glutamine synthetase 化学候选；早期同物种纯化论文仅作功能线索，暂不作为本 WP 的直接实验分数。',
    'chemistry': 'L-glutamate + ammonium + ATP -> L-glutamine + ADP + phosphate；精确质子和区室以模型代谢物式核对。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S027'
}
manual['RU820_RS02030'] = {
    'issue': 'BRENDA EC 2.2.1.2 所引 He et al. 2005 是磷限制条件下差异蛋白组，报告一个 transaldolase 蛋白点；原文该点标识 STY3758，未与现行 RU820_RS02030/WP 建立序列对应。Osorio et al. 2013 同株蛋白组明确记录 AFE_0419 表达，仍未测该 WP 的催化底物。',
    'recommendation': '保留 2016 transaldolase 关联候选和同株表达信息；精确反应分数不因差异蛋白组记录上调。',
    'chemistry': 'transaldolase 反应按 2016 式子与 KEGG/Rhea 逐式核对；蛋白丰度变化不决定底物或方向。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S018;G52-S027'
}
manual['RU820_RS02195'] = {
    'issue': '2016 BILG 把 biotin + ATP + BCCP 写为整体反应；现行 WP 注释为 biotin--[acetyl-CoA-carboxylase] ligase，KEGG 分别给 biotinyl-AMP 生成与蛋白赖氨酸生物素化步骤。2016 代谢物 bccp[c] 缺元素式，旧式无法完整配平。',
    'recommendation': '保留整体生物素连接反应候选；Hold 2016 精确计量式和 2026 分数，先补 BCCP 及生物素化蛋白形式，再核对合并两步时 AMP/PPi 系数。',
    'chemistry': 'biotin + ATP -> biotinyl-AMP + PPi；biotinyl-AMP + apo-BCCP lysine -> biotinyl-BCCP + AMP；整体式取决于蛋白载体定义。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S025;G52-S026'
}
manual['RU820_RS02370'] = {
    'issue': '2016 TDSK 使用 lipidA_AFE[c]、lipid4A_AFE[c] 两个模型特定脂质；旧补表没有足够元素式/电荷，KEGG R04657 使用通用 lipid A disaccharide 与 lipid IVA。当前注释支持 LpxK 家族，但酰链结构未与模型代谢物逐一对应。',
    'recommendation': '保留 tetraacyldisaccharide 4-prime-kinase 功能候选；Hold 2016 精确式，先确定两个模型脂质的结构与质子平衡。',
    'chemistry': 'lipid A disaccharide + ATP -> lipid IVA + ADP；具体 lipidA_AFE/lipid4A_AFE 的酰链及磷酸化结构待定。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016'
}
manual['RU820_RS02870'] = {
    'issue': '现行 WP 仅注释为 class I SAM-dependent methyltransferase；KEGG KO 的磷脂甲基化条目连接三种底物。2016 PMEAS 的 pe_AFE[c] 与 ptdmeeta_AFE[c] 元素式缺失，旧式无法完整配平，也没有当前 WP 特定磷脂底物实验。',
    'recommendation': 'Hold PMEAS 该基因 GPR 和精确式子；先确认具体磷脂底物、甲基化步骤以及两个模型脂质定义。',
    'chemistry': 'SAM + phosphatidylethanolamine -> SAH + phosphatidyl-N-methylethanolamine 是一种可能步骤；不可由 class I 甲基转移酶家族名直接指定。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S025;G52-S026'
}
manual['RU820_RS02245'] = {
    'issue': 'Acosta et al. 2005 对 ATCC 23270 重组 rhodanese P16.2 测得 thiosulfate:cyanide sulfurtransferase 活性。论文同名序列 AY863108.1 与现行 RU820_RS02245 CDS 450/451 nt 对齐后仅多一个 G（历史序列第 361–363 位 G 连续区），导致蛋白 C 端读码不同；现行 WP_009564831.1 并非被测的完全相同序列。P15 的序列 AY863107 对应本组外 WP，不得混同。',
    'recommendation': '记录 P16.2 的体外 thiosulfate + cyanide -> sulfite + thiocyanate 直接实验，反应核心可评 4；当前 WP 的 GPR 仍需序列版本和活性核验。cyanide 是体外受体，生理受体与代谢网络连通性未证实，暂不直接纳入 2026 模型。',
    'chemistry': 'thiosulfate + cyanide -> sulfite + thiocyanate；体外 TST 测定。细胞区室和生理硫受体未定。',
    'evidence': 'G52-S001;G52-S005;G52-S030'
}
manual['RU820_RS02155'] = {
    'issue': 'Barreto et al. 2005 用 ATCC 23270 galU 克隆互补 E. coli galU 缺失株，并在带克隆的细胞提取物中测得 UDP-glucose pyrophosphorylase 活性：45±7 对阴性 0.5±0.1 nmol/min/mg。论文提交 AY789511.1 与现行 RU820_RS02155/WP_009567323.1 的 897-nt CDS 完全相同。',
    'recommendation': '确认本 WP 的 GalU 催化和 GPR；2016 GALU 历史分数 2 保留原列，2026 对 UTP + glucose-1-phosphate -> UDP-glucose + PPi 的核心反应评 4。旧式中的 H+ 系数已按补表配平，但未被实验单独测定；细胞区室及模型通量方向另以网络证据核对。',
    'chemistry': 'alpha-D-glucose 1-phosphate + UTP -> UDP-alpha-D-glucose + diphosphate；2016 完整式 h[c] + utp[c] + g1p[c] -> ppi[c] + udpg[c] 配平。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S031'
}
for loc in ['RU820_RS02555','RU820_RS02995']:
    note='Barreto et al. 2005 克隆并实测的 phosphoglucomutase 提交序列 AY789512.1 精确对应当前 WP_012537123.1，位于 G52 范围之外，不能作为此 G52 PGM1 分支的直接实验依据。'
    if note not in manual[loc]['issue']:
        manual[loc]['issue'] += ' '+note
    manual[loc]['evidence'] = ';'.join(dict.fromkeys(manual[loc]['evidence'].split(';')+['G52-S031']))
write('g52_manual_adjudications.json', manual)

overrides = read('g52_reaction_overrides.json')
overrides['SFAT1'] = {
    'action': '保留并核正证据等级',
    'scope': '2012 同株历史 atpS 等位序列重组表达实测 APS+PPi 生成 ATP；与现行 CDS 差 7 nt/4 aa。核心反应可评 4；质子、区室与现行 WP 版本差异单列。'
}
overrides['ADSK']['scope'] = '当前 WP 注释 APS kinase 反应；2012 atpS 论文仅预测该结构域，未测 kinase 活性。'
overrides['SFAT2']['scope'] = '2012 atpS 论文未测 APS+Pi -> ADP+sulfate，当前 WP 亦无此精确反应对应。'
overrides['GLUTRS'] = {
    'action': 'Hold该基因OR分支',
    'scope': '历史 AFE_0422 在 KEGG 标作 gltX-2；同物种 GluRS2 实验偏好 tRNA(Gln)-UUG，而旧 GLUTRS 使用 tRNA(Glu)；原文蛋白序列与当前 WP 尚未精确对齐。'
}
overrides['BILG'] = {'action': 'Hold旧式子', 'scope': 'BCCP 与生物素化 BCCP 缺完整分子式；KEGG 对应两步反应，需核对整体计量。'}
overrides['TDSK'] = {'action': 'Hold旧式子', 'scope': '两个 AFE 专用 lipid A 代谢物缺结构/式子，无法核对酰链与电荷。'}
overrides['PMEAS'] = {'action': 'Hold该基因GPR及旧式子', 'scope': '当前 WP 仅为通用 SAM methyltransferase；模型脂质缺式子且具体底物未确认。'}
write('g52_reaction_overrides.json', overrides)
print({'manual': len(manual), 'atps_current_sequence': read('g52_atps_2012_sequence_comparison.json')['current_wp']})
