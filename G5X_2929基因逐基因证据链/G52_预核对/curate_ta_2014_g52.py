"""Record paper-specific limits for four historical TA labels in G52."""
import json
from pathlib import Path

DIR = Path(__file__).resolve().parent
path = DIR / 'g52_manual_adjudications.json'
manual = json.loads(path.read_text(encoding='utf-8'))
details = {
    'RU820_RS02005': ('AFE_0413', 'VapC 毒素', 'RU820_RS02010'),
    'RU820_RS02010': ('AFE_0414', 'Phd 抗毒素', 'RU820_RS02005'),
    'RU820_RS02310': ('AFE_0477', 'VapB 抗毒素', 'RU820_RS02315'),
    'RU820_RS02315': ('AFE_0478', 'VapC 毒素', 'RU820_RS02310'),
}
for loc, (afe, role, partner) in details.items():
    manual[loc] = {
        'issue': f'Bustamante 等 2014 原文 Table 2 将 {afe} 列为 ATCC 23270 的 {role}，对应伴侣 {partner}；该文直接抑菌与 RNase 实验针对 ICEAfe1 上其他 TA 对，本位点仅属基因组/家族归类。',
        'recommendation': '保留现行 TA 功能注释线索；本轮无可写入代谢模型的独立小分子反应或 GPR 修改。',
        'chemistry': '无可确认的小分子代谢反应；毒素/抗毒素作用属于 RNA 或蛋白调控过程。',
        'evidence': 'G52-S001;G52-S002;G52-S024',
    }
path.write_text(json.dumps(manual, ensure_ascii=False, indent=2), encoding='utf-8')
print('TA paper boundaries updated:', ', '.join(details))
