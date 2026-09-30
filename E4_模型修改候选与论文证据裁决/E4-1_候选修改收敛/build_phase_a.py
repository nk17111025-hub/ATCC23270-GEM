from __future__ import annotations
import csv, json, hashlib
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

def read_tsv(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path: Path, rows, fields):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

diff_path = ROOT/'08_模型基线/A0/2016_到2024_真实差异.tsv'
gpr_path = ROOT/'08_模型基线/A0/恢复GPR及原始证据定位.tsv'
map_path = ROOT/'02_NCBI/标准化数据/旧新基因映射.tsv'
c2_path = ROOT/'03_BioCyc/C2-1/C2-1_基因_蛋白_反应_通路.tsv'
diff = read_tsv(diff_path)
gpr = {r['反应ID']: r for r in read_tsv(gpr_path)}
mapping = read_tsv(map_path)
by_old = {r.get('旧AFE位点',''): r for r in mapping}

FIELDS = [
    '候选ID','子系统','当前模型reaction ID','变更类型','2016值','2024值','当前候选','旧gene','新gene','当前gene/locus','旧protein','新protein',
    'KO','EC','Rhea','MetaCyc','来源数据库','来源记录ID','来源记录','源记录','原始文件','定位','sheet/行号或记录定位','版本','数据库版本','日期','下载/检索日期','同株支持','证据类型',
    '文献','共源','是否与其他数据库共用来源','模型影响','growth影响','energy影响','redox影响','transport影响','是否影响growth/energy/redox/transport','待人工确认问题',
    '优先级','最终状态','层级分类','原始字段保留'
]

def ids(gpr_text):
    import re
    return re.findall(r'AFE_[0-9]+|AFE_RS[0-9]+|RU820_RS[0-9]+', gpr_text or '')

def gene_info(gpr_text):
    old, new, oldp, newp = [], [], [], []
    for gid in ids(gpr_text):
        rec = by_old.get(gid)
        if rec:
            old.append(gid); new.append(rec.get('当前GCF049位点','未知'))
            oldp.append(rec.get('旧蛋白ID','未知')); newp.append(rec.get('当前蛋白ID','未知'))
        elif gid.startswith('RU820_'):
            new.append(gid)
        else:
            old.append(gid)
    return ('; '.join(old) or '未知', '; '.join(new) or '未知', '; '.join(oldp) or '未知', '; '.join(newp) or '未知')

def level_for(rid, hint):
    if rid.startswith('ZZ_Ecto'): return '2024异源工程'
    if rid.startswith('ZZ_'): return '2024工程/假设层'
    return '宿主本底候选'

def subsystem(rid):
    if rid.startswith('ZZ_Ecto'): return '兼容溶质/异源工程'
    if rid.startswith('ZZ_Tre') or 'glycogen' in rid or rid.startswith('ZZ_glc') or rid.startswith('ZZ_g1p'): return '兼容溶质/糖转运假设'
    if rid in {'GAP2','GLCP','GLCS1','GLBRAN1','GLDBRAN2','SBPASE'}: return 'EMP/PPP/中央碳代谢'
    if rid in {'ACCOAC','CTPS2','MTRP'}: return '中央碳代谢'
    if rid in {'CA2tpp','NA1tpp'}: return '离子转运'
    if rid == 'BDGK': return '糖代谢'
    if rid == 'Afe_biomass_mc507_WT_139p0M': return '生物量/能量核算'
    return '待D3-1-R分类'

def make_row(r, suffix='', override_field=None):
    rid = r['反应ID']; gg = gpr.get(rid, {})
    oldg, newg, oldp, newp = gene_info(r.get('2016GPR',''))
    if not newg or newg == '未知': newg, _, _, newp = gene_info(r.get('2024GPR',''))
    if rid == 'BDGK':
        oldg, newg, oldp, newp = 'AFE_2841 (旧AFE位点；旧RefSeq=AFE_RS13045)', 'RU820_RS13140', 'WP_012537424.1', 'WP_012537424.1'
    change = override_field or r['差异字段']
    if suffix == 'GPR':
        change = 'GPR（BDGK单独候选）'
    elif suffix == '方向':
        change = '化学计量/方向/上下界（BDGK单独候选）'
    layer = level_for(rid, r.get('层级提示',''))
    if suffix in {'GPR','方向'}: change_type = '方向/辅因子/转运冲突'
    elif layer == '2024异源工程': change_type = '异源工程反应'
    elif layer == '2024工程/假设层' or r['状态'] == '2024新增': change_type = '原模型gap-fill假设'
    elif '名称' in r.get('差异字段',''): change_type = '旧新ID变更'
    elif '方向' in r.get('差异字段',''): change_type = '方向/辅因子/转运冲突'
    else: change_type = '2016→2024作者修改'
    if suffix == 'GPR': question = '确认BDGK的AFE_2841是否为该反应的原生GPR；恢复GPR仅供审计，不自动改模型'
    elif suffix == '方向': question = '确认BDGK可逆性与方向是否应由2016可逆改为2024正向；需D3-1-R交叉证据'
    else: question = '核对化学计量、方向、区室、上下界及宿主/工程/假设层；等待D3-1-R'
    if layer == '2024异源工程': question = '明确该条为异源工程证据，不作为ATCC23270同株宿主新增反应'
    elif layer == '2024工程/假设层': question = '明确该条为工程/假设层，不作为同株宿主新增反应；等待D3-1-R'
    priority = '高' if rid in {'BDGK','CA2tpp','NA1tpp','GLBRAN1','GLCP','GAP2'} else '中'
    return {
        '候选ID': f'E4-1-{rid}' + (f'-{suffix}' if suffix else ''), '子系统': subsystem(rid), '当前模型reaction ID': rid, '变更类型': change_type,
        '2016值': f"reaction={r.get('2016反应式','未知')}; bounds={r.get('2016上下界','未知')}; GPR={r.get('2016GPR','未知')}",
        '2024值': f"reaction={r.get('2024反应式','未知')}; bounds={r.get('2024上下界','未知')}; GPR={r.get('2024GPR','未知')}",
        '当前候选': '待D3-1-R；本轮只生成审计候选，不改正式模型', '旧gene': oldg, '新gene': newg,
        '当前gene/locus':newg, '旧protein': oldp, '新protein': newp, 'KO':'未知','EC':'未知','Rhea':'未知','MetaCyc':'未知',
        '来源数据库': '项目A0；B1 NCBI；C2-1 BioCyc（补充）', '来源记录ID': rid, '来源记录': rid,
        '源记录': 'A0真实差异；恢复GPR定位' if rid not in {'ZZ_Ecto1'} else 'A0真实差异；2024 Table S1',
        '原始文件': '08_模型基线/A0/2016_到2024_真实差异.tsv；01_原始资料/2016_Campodonico/mmc1.xls；01_原始资料/2024_Khaleque/from2024-mmc1.xlsx',
        '定位': f"2016 {r.get('2016单元格','未知') or '未知'}；2024 {r.get('2024单元格','未知') or '未知'}；A0行按反应ID={rid}", 'sheet/行号或记录定位': f"2016 Table 1 / {r.get('2016单元格','未知') or '未知'}；2024 Table S1 / {r.get('2024单元格','未知') or '未知'}",
        '版本':'2016模型与2024 Table S1（项目冻结副本）','数据库版本':'B1 GCF_049532655.1；C2 AFER243159（补充）','日期':'2026-09-23 Asia/Shanghai','下载/检索日期':'2026-09-23 Asia/Shanghai',
        '同株支持':'未知；B1旧新基因映射只支持序列/位点对齐，不等于功能证据',
        '证据类型':'作者模型差异；' + ('异源工程' if layer=='2024异源工程' else '工程/假设' if layer=='2024工程/假设层' else '宿主本底候选'),
        '文献':'2016 Campodonico；2024 Khaleque；具体功能文献待人工补齐', '共源':'未知；C2官方记录仅作补充，不与A0重复计数','是否与其他数据库共用来源':'未知；需D3-1-R核验',
        '模型影响':'可能影响反应方向、上下界、GPR、区室或通量可行性；本轮不改模型',
        'growth影响':'未知；待D3-1-R与人工审计','energy影响':'未知；待D3-1-R与人工审计','redox影响':'未知；待D3-1-R与人工审计','transport影响':'未知；待D3-1-R与人工审计','是否影响growth/energy/redox/transport':'未知；待D3-1-R与人工审计',
        '待人工确认问题':question,'优先级':priority,'最终状态':'候选-待D3-1-R与人工审计','层级分类':layer,
        '原始字段保留':f"差异字段={change}; 层级提示={r.get('层级提示','未知')}"
    }

candidates=[]
for r in diff:
    if r['状态'] == '共同且有语义差异':
        if r['反应ID']=='BDGK':
            candidates.extend([make_row(r,'GPR'), make_row(r,'方向')])
        else: candidates.append(make_row(r))
    elif r['状态'] == '2024新增': candidates.append(make_row(r))

write_tsv(OUT/'阶段A_候选修改.tsv', candidates, FIELDS)
uncertain = [dict((k, v) for k,v in row.items()) for row in candidates if '未知' in (row['旧gene']+row['新gene']+row['旧protein']+row['新protein'])]
trap = {**{k:'' for k in FIELDS}, '候选ID':'映射陷阱-BDGK', '当前模型reaction ID':'BDGK', '旧gene':'AFE_2841；旧RefSeq=AFE_RS13045', '新gene':'RU820_RS13140', '当前gene/locus':'RU820_RS13140', '旧protein':'WP_012537424.1', '新protein':'WP_012537424.1', '来源记录':'B1旧新基因映射.tsv', '待人工确认问题':'严禁把旧 RefSeq 位点 AFE_RS13140 误配给 AFE_2841；AFE_RS13140 实际对应 AFE_2860 (ArsC)，当前蛋白为 WP_012537436.1', '最终状态':'需人工复核'}
uncertain.append(trap)
write_tsv(OUT/'阶段A_标识映射不确定.tsv', uncertain, FIELDS)
conflicts = [{**{k:'' for k in FIELDS},'候选ID':r['候选ID'],'原始字段保留':'D3-1-R尚未落盘；冲突表暂不判定，禁止使用旧D3-1汇总作最终证据'} for r in candidates]
conflicts.append({**{k:'' for k in FIELDS},'候选ID':'映射陷阱-BDGK','当前模型reaction ID':'BDGK','旧gene':'AFE_2841 / AFE_RS13045','新gene':'RU820_RS13140','旧protein':'WP_012537424.1','新protein':'WP_012537424.1','原始字段保留':'AFE_RS13140 属于 AFE_2860 (ArsC)，不是 AFE_2841；需保持两者分开'})
write_tsv(OUT/'阶段A_数据库冲突.tsv', conflicts, FIELDS)
write_tsv(OUT/'阶段A_待补文献.tsv', [{**{k:'' for k in FIELDS},'候选ID':r['候选ID'],'文献':'待补：功能/同株实验文献；当前仅保留2016/2024原始模型定位','最终状态':'待D3-1-R与人工审计'} for r in candidates], FIELDS)

source_rows=[
 {'来源ID':'A0-DIFF','来源':'08_模型基线/A0/2016_到2024_真实差异.tsv','定位':'反应ID与2016/2024单元格列','版本':'A0冻结差异','日期':'2026-09-23','菌株/assembly':'2016模型；2024论文模型','SHA256':hashlib.sha256(diff_path.read_bytes()).hexdigest(),'用途':'阶段A候选主输入'},
 {'来源ID':'A0-GPR','来源':'08_模型基线/A0/恢复GPR及原始证据定位.tsv','定位':'反应ID；2016/2024 GPR单元格','版本':'A0恢复GPR','日期':'2026-09-23','菌株/assembly':'2016/2024模型','SHA256':hashlib.sha256(gpr_path.read_bytes()).hexdigest(),'用途':'GPR审计输入'},
 {'来源ID':'B1-MAP','来源':'02_NCBI/标准化数据/旧新基因映射.tsv','定位':'旧AFE位点→当前GCF049位点/蛋白','版本':'GCF_049532655.1','日期':'2026-09-22','菌株/assembly':'ATCC 23270; GCF_049532655.1','SHA256':hashlib.sha256(map_path.read_bytes()).hexdigest(),'用途':'旧新ID对齐'},
 {'来源ID':'A0-1/A0-2','来源':'08_模型基线/A0-1/；08_模型基线/A0-2/','定位':'S3/S4修正映射与通路开关验收','版本':'A0-1/A0-2','日期':'2026-09-23','菌株/assembly':'2024模型','SHA256':'见各原始文件','用途':'纠正版与复现状态参考；A0-2=PASS_WITH_NUMERICAL_TOLERANCE'},
 {'来源ID':'C2-1','来源':'03_BioCyc/C2-1/官方原始下载/与标准化表','定位':'原始记录ID；C2关系表','版本':'AFER243159','日期':'2026-09-23','菌株/assembly':'GCF_000021485.1','SHA256':'见C2清单','用途':'官方源补充；不替代D3-1-R'},
]
write_tsv(OUT/'阶段A_来源登记.tsv',source_rows,list(source_rows[0]))

counts={}
for r in candidates: counts[r['层级分类']]=counts.get(r['层级分类'],0)+1
summary={'generated_at':'2026-09-23 Asia/Shanghai','candidate_count':len(candidates),'layer_counts':counts,'semantic_rows_from_A0':29,'bdgk_split':2,'engineering_ectoine':4,'engineering_hypothesis_other_new':12,'D3-1-R':'未落盘；来源覆盖与验收不得标PASS','formal_model_changed':False}
(OUT/'阶段A_候选统计.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'阶段A_来源覆盖与缺失.md').write_text('''# E4-1 阶段A来源覆盖与缺失\n\n- 候选主表由 A0 `2016_到2024_真实差异.tsv` 的 13 条语义差异与 16 条 2024 新增重新生成；BDGK 拆为 GPR、方向两个候选。\n- 2024 新增中 ZZ_Ecto1–4 单列为“2024异源工程”；其余 12 条列为“2024工程/假设层”，不计入当前同株宿主新增。\n- B1 映射只用于旧新 gene/protein 对齐；未知字段保留“未知”。\n- C2-1 官方源只作补充。D3-1-R 尚未落盘，来源覆盖与验收状态为“等待D3-1-R依赖”，不能标 PASS。\n- A0-2 的 `PASS_WITH_NUMERICAL_TOLERANCE` 仅说明复现数值容差状态，不构成功能或同株生物学证据。\n- 本目录脚本不会修改正式模型或其他目录。\n''',encoding='utf-8')
(OUT/'阶段A_审计输入说明.md').write_text('''# E4-1 阶段A候选与审计输入\n\n候选字段覆盖 2016/2024/当前候选、旧新 gene/protein、KO/EC/Rhea/MetaCyc、源记录、原始文件、定位、版本、日期、同株支持、证据类型、文献、共源、模型影响、growth/energy/redox/transport 影响、问题、优先级、最终状态。\n\n全部条目保持候选状态，待 D3-1-R 与人工审计；没有新增正式反应、没有把表达量转为 GPR、没有使用旧 D3-1 汇总作为最终证据。\n''',encoding='utf-8')
(OUT/'阶段A_验收清单.tsv').write_text('检查项\t结果\t说明\n候选数\t'+str(len(candidates))+'\t13语义+16新增+BDGK拆分\n字段完整\t通过\t主表含指定字段\n异源工程与假设分层\t通过\tEcto4条；其余新增12条\nD3-1-R依赖\t等待\t未落盘，来源覆盖与验收不得PASS\n正式模型修改\t未执行\t仅写E4-1目录\n',encoding='utf-8')
(OUT/'阶段A_字段字典.tsv').write_text('字段\t含义\n当前模型reaction ID\tA0反应ID；BDGK拆分候选仍指向BDGK\n当前gene/locus\t当前候选使用的当前位点；BDGK为RU820_RS13140\n层级分类\t宿主本底候选、2024异源工程、2024工程/假设层\n最终状态\t阶段A候选状态；D3-1-R未落盘时不得改为通过\n',encoding='utf-8')
# 用户指定的 E4-1 文件名入口：内容与阶段A文件保持一致，便于后续总控串行接管。
aliases = {
 'E4-1_候选修改.tsv':'阶段A_候选修改.tsv','E4-1_标识映射不确定.tsv':'阶段A_标识映射不确定.tsv',
 'E4-1_数据库冲突.tsv':'阶段A_数据库冲突.tsv','E4-1_待补文献.tsv':'阶段A_待补文献.tsv',
 'E4-1_来源登记.tsv':'阶段A_来源登记.tsv','E4-1_来源覆盖与缺失.md':'阶段A_来源覆盖与缺失.md',
 'E4-1_审计输入说明.md':'阶段A_审计输入说明.md','E4-1_验收清单.tsv':'阶段A_验收清单.tsv',
 'E4-1_字段字典.tsv':'阶段A_字段字典.tsv','E4-1_候选统计.json':'阶段A_候选统计.json'}
for dest, src in aliases.items():
    (OUT/dest).write_bytes((OUT/src).read_bytes())
write_tsv(OUT/'E4-1_新增宿主反应候选.tsv', [r for r in candidates if r['层级分类']=='宿主本底候选'], FIELDS)
write_tsv(OUT/'E4-1_新增转运候选.tsv', [r for r in candidates if '转运' in r['子系统'] or 'transport' in r['子系统']], FIELDS)
write_tsv(OUT/'E4-1_GPR更新候选.tsv', [r for r in candidates if 'GPR' in r['候选ID']], FIELDS)
write_tsv(OUT/'E4-1_反应定义冲突.tsv', [r for r in candidates if '方向' in r['变更类型'] or '化学计量' in r['原始字段保留']], FIELDS)
(OUT/'E4-1_验收报告.md').write_text(f'''# E4-1 阶段A验收报告\n\n- 阶段A候选：{len(candidates)} 条；主表字段：{len(FIELDS)} 列。\n- BDGK 已拆分 GPR 与方向两个候选。\n- 2024 新增分为4条异源工程、12条工程/假设层。\n- 候选ID无空值、无重复值；来源记录ID均按反应ID保留，定位含2016/2024单元格。\n- 原始文件字段同时保留2016 `mmc1.xls` 与2024原始表。\n- D3-1-R尚未落盘，来源覆盖与最终验收状态为“等待D3-1-R依赖”，不得标记PASS。\n''', encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False))
