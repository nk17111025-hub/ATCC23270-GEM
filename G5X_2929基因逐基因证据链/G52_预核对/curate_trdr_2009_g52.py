"""Add the primary TrxR assay after validating its primers against current CDS."""
import json
from pathlib import Path

DIR = Path(__file__).resolve().parent
primer = json.loads((DIR / 'g52_trxr_2009_exact_identity.json').read_text(encoding='utf-8'))
assert primer['both_primer_ends_match_current_cds']

path = DIR / 'g52_manual_adjudications.json'
manual = json.loads(path.read_text(encoding='utf-8'))
manual['RU820_RS01820'] = {
    'issue': 'Wang et al. 2009 直接测定 ATCC 23270 重组 TrxR 对氧化型 thioredoxin 的 NADPH 依赖还原；论文两端基因特异引物均与现行 RU820_RS01820 CDS 精确匹配。2016 TRDR 式子把还原型 thioredoxin 作为反应物，并按旧补表电荷计算净差 +4。',
    'recommendation': '确认本 WP 的 thioredoxin-disulfide reductase 催化功能；2016 TRDR 原式应替换为已配平的 h[c] + nadph[c] + trdox[c] -> nadp[c] + trdrd[c] 候选，完成全模型网络核验后落表。',
    'chemistry': '原始实验证明 NADPH + oxidized thioredoxin -> NADP+ + reduced thioredoxin；按 2016 代谢物式与电荷，候选 h[c] + nadph[c] + trdox[c] -> nadp[c] + trdrd[c] 原子与电荷均平衡。',
    'evidence': 'G52-S001;G52-S004;G52-S006;G52-S007;G52-S016;G52-S022',
}
path.write_text(json.dumps(manual, ensure_ascii=False, indent=2), encoding='utf-8')

path = DIR / 'g52_reaction_overrides.json'
overrides = json.loads(path.read_text(encoding='utf-8'))
overrides['TRDR'] = {
    'action': '替换2016错误式子候选',
    'scope': '同株重组蛋白催化实验和论文引物与当前 CDS 双端精确匹配；原式方向及电荷错误；候选 h[c] + nadph[c] + trdox[c] -> nadp[c] + trdrd[c] 配平，待网络核验。',
}
path.write_text(json.dumps(overrides, ensure_ascii=False, indent=2), encoding='utf-8')
print('TRDR manual evidence and override updated')
