import csv
import json
from collections import defaultdict
from pathlib import Path

out = Path(__file__).resolve().parent / 'G57_审查输出'
scan = json.loads((out / 'EuropePMC_gene_scan_20260927.json').read_text(encoding='utf-8'))
current_scan = json.loads((out / 'EuropePMC_current_ids_20260927.json').read_text(encoding='utf-8')) if (out / 'EuropePMC_current_ids_20260927.json').exists() else {}

# 文献结论限于当前已核对的标题、原文实验段落及基因对象；“待核”不可用作反应评分。
decision = {
 '10.1007/s12223-013-0244-8': ('A. ferrooxidans DC', '异株金属胁迫表达', 'AFE_2126 在 DC 株克隆/表达；不可视为 ATCC 23270 的同株反应实测'),
 '10.1016/j.jbc.2024.107703': ('ATCC 23270', '同株 CRISPRi；命中基因仅背景', '敲低对象是 petA/petB；AFE_2321/2323/2324 为背景中的 EPS 候选，不能给这些基因赋实验反应'),
 '10.1038/ncomms9883': ('ATCC 23270 基因来源；异源蛋白', '重组酶与活化因子实验', 'CbbM 正向羧化及 CbbQ2/CbbO2 活化有直接实验；逆向和体内区室未证'),
 '10.1074/mcp.m700042-mcp200': ('ATCC 23270', '周质蛋白质组', '原文 AFE_2266 蛋白质组有检出；旧文功能名与当前 GNAT 注释不一致，需肽段/WP 复核；不证明催化反应'),
 '10.1128/aem.03057-12': ('ATCC 23270T', '厌氧硫/铁培养与表达', 'SreABCD 表达与生理可用作间接证据；其他命中仅组学记录；不能整体继承硫还原或氧化的精确式'),
 '10.1128/aem.07230-11': ('ATCC 23270', 'pfkB 敲除及异源酶活', '直接对象为 AFE_1807（本组以外）；本组 AFE_2025/2155/2324 命中仅背景'),
 '10.1186/1471-2164-10-394': ('ATCC 23270 基因组', '铁硫代谢模型', 'AFE_2157/2324 为模型或背景提及；无本基因催化实测'),
 '10.1186/1471-2180-11-259': ('实验为 A. ferrooxidans LR', '旧编号假阳性；已用引物排除', '文中 Afe_2172 两引物落当前 796313/796390 的 Hsp20 区域，不落 G57 RU820_RS10045'),
 '10.1371/journal.pone.0112226': ('ATCC 23270 与 53993', '毒素抗毒素系统研究', 'AFE_2129/2130 为毒素抗毒素基因；不形成代谢反应；具体实验对象与蛋白活性需按原文细核'),
 '10.3389/fbioe.2018.00157': ('综述；多菌种', '二次文献', 'AFE_2018 仅线索，不能作为 ATCC 23270 单基因实验'),
 '10.3389/fbioe.2020.543807': ('ATCC 23270 基因来源；E. coli 宿主', '异源 CbbM/CbbQO 表达', '支持 CbbM 羧化及 CbbQO 活化；宿主是 E. coli，不能据此指定本株体内方向/区室'),
 '10.3389/fmicb.2016.01365': ('ATCC 23270T', '群体感应扰动与表达组学', '命中基因有表达变化；不等于对应蛋白催化精确反应'),
 '10.3389/fmicb.2019.00592': ('DSM 14882T＝ATCC 23270T', '黄铁矿生物膜蛋白质组', '蛋白丰度及全细胞生理为同株间接证据；不等于单酶反应测定'),
 '10.3389/fmicb.2019.00603': ('A. ferrooxidans；菌株需逐篇确认', 'CO2 条件转录研究', 'CbbM/Q/O 等转录线索；不等于这些蛋白的催化活性测定'),
 '10.3390/genes11080844': ('A. ferrivorans ACH', '异种铜抗性研究', 'AFE_2021 作为 ferrooxidans 同源对照；不可给本株赋直接实验'),
 '10.3390/ijms23073580': ('A. ferrooxidans；菌株需逐篇确认', '辉锑矿条件表达研究', 'AFE_2312 表达/同源性线索；未确定被转运的糖、耦联离子或跨膜方向'),
 '10.1039/d4cb00195h': ('ATCC 23270 来源 WP_012536942.1；E. coli 表达', '重组酶纯化、PLP 晶体结构、半胱氨酸/Cd2+ 体系纳米粒子生成', '直接支持该 WP 的蛋白身份和体系活性；未完整测定半胱氨酸分解化学计量，也未测旧 HACT 的 O-乙酰高丝氨酸反应'),
}

by_doi = defaultdict(lambda: {'afe': set(), 'paper': None})
for afe, item in scan.items():
    for paper in item['papers']:
        doi = paper.get('doi', '').lower()
        if not doi:
            continue
        by_doi[doi]['afe'].add(afe)
        by_doi[doi]['paper'] = paper
for key, item in current_scan.items():
    for paper in item.get('papers', []):
        doi = paper.get('doi', '').lower()
        if not doi:
            continue
        by_doi[doi]['afe'].add(key)
        by_doi[doi]['paper'] = paper
assert set(by_doi) == set(decision), (set(by_doi) - set(decision), set(decision) - set(by_doi))

rows = []
for doi in sorted(by_doi):
    item = by_doi[doi]
    organism, typ, boundary = decision[doi]
    paper = item['paper']
    rows.append({'DOI':doi, '标题':paper['title'], '年份':paper.get('year',''), '命中编号':';'.join(sorted(item['afe'])), '命中编号数':len(item['afe']), '论文实验对象':organism, '证据类型':typ, '对 G57 的使用边界':boundary, '原文':f'https://doi.org/{doi}'})
with (out / 'G57_文献命中逐篇判读.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print(f'{len(rows)} DOI classified')
