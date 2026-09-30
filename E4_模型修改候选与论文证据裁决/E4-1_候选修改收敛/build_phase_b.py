from pathlib import Path
import subprocess
import csv, json, re

ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent; D3=ROOT/'D3-1-R'
def read(p):
    with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    with open(p,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
cross=read(D3/'D3-1-R_全基因跨数据库交叉引用.tsv')
ann=read(D3/'D3-1-R_注释证据表.tsv'); ann_by={}
for x in ann: ann_by.setdefault(x['当前locus'],[]).append(x)
central=read(D3/'D3-1-R_中央碳代谢证据.tsv'); trans=read(D3/'D3-1-R_有机碳转运证据.tsv'); etc=read(D3/'D3-1-R_FeS_ETC证据.tsv')
targets=[]
for group, rows, n in [('中央碳代谢',central,4),('有机碳转运',trans,7),('FeS_ETC',etc,7)]:
    targets += [(group,x) for x in rows if x.get('旧模型未覆盖功能候选')=='是'][:n]
selected={x['当前locus'] for _,x in targets}
key_words=re.compile(r'NADH|NADPH|quinone|Complex I|Nuo|terminal oxidase|cytochrome|ATP synthase|iron-sulfur|Fe-S',re.I)
for x in cross:
    if len(targets)>=30: break
    if x['当前locus'] not in selected and x.get('旧模型未覆盖功能候选')=='是' and key_words.search(x.get('NCBI当前功能','')):
        targets.append(('关键能量/氧化还原',x)); selected.add(x['当前locus'])

FIELDS=['候选ID','子系统','当前模型reaction ID','变更类型','当前locus','当前protein ID','旧AFE locus','旧RefSeq locus','旧protein ID','KO','EC','Rhea','MetaCyc','KEGG reaction','BioCyc reaction ID','D3功能判定','来源数据库','来源记录ID','原始文件','数据库版本/组装','取得日期','同株支持','证据类型','A0化学计量/区室比较','模型影响','growth影响','energy影响','redox影响','transport影响','待人工确认问题','优先级','最终状态']
def evidence(locus):
    rs=ann_by.get(locus,[]); return '; '.join(f"{x['数据库']}:{x['记录ID']}={x['证据字段']}" for x in rs), '; '.join(x['来源文件'] for x in rs), '; '.join(sorted(set(x['数据库版本/组装'] for x in rs))), '; '.join(sorted(set(x['取得日期'] for x in rs)))
def row(group,x,i):
    locus=x['当前locus']; rec,files,vers,dates=evidence(locus); func=x.get('NCBI当前功能','未知')
    if group=='中央碳代谢': typ='当前同株新增GPR候选'; q='先与A0已有反应逐一核对化学计量/区室；AFE_1667需对照TKT1/TKT2，AFE_1883需对照PPC，不能据同ID/EC判等价'
    elif group=='有机碳转运': typ='无法判定/证据增强'; q='确认底物、跨膜方向、区室、反应式与GPR；D3注释不足时保留未裁决，不直接新增宿主反应'
    elif group=='FeS_ETC': typ='证据增强/无法判定'; q='确认复合体成员、电子供体/受体、区室及是否已有模型反应；AFE_0048需检验TSQOC复合体'
    else: typ='证据增强/无法判定'; q='关键能量/氧化还原家族注释仅作为审计线索；必须找到对应模型反应与化学计量后再裁决'
    if x.get('映射冲突','否') not in ('','否'):
        typ='映射歧义/无法判定'; q='D3标记旧新位点映射冲突或多重匹配；不能据当前locus直接赋予GPR，需Reviewer裁决'
    if locus=='RU820_RS10385':
        typ='底物冲突/无法判定'; q='carbohydrate porin注释与旧GPR的Fe²⁺转运指向冲突；不能直接列为有机碳转运，需底物与区室证据'
    if locus=='RU820_RS03190': typ='证据增强/数据库分歧'; q='KEGG EC 1.1.1.40 与 BioCyc EC 1.1.1.37/1.1.1.38 分歧；AFE_0660 cofactor/反应定义待裁决'
    if locus=='RU820_RS09355': typ='新增宿主反应候选'; q='PPP候选；BioCyc RXN-9952与A0 6PGDH定义需核对，确认6PG脱氢酶化学计量、NADP/NAD辅因子及区室后再考虑新反应'
    if locus=='RU820_RS01305': typ='证据增强/无法判定'; q='BioCyc RXN-13161为CPD-11281+O2+H2O→glutathione+SO3+H+；模型SULDO为s[p]+H2O[p]+O2[p]↔2H+[p]+SO3[p]，底物/区室不同，不能按同名视为等价'
    return {'候选ID':f'E4-1-B-{i:03d}','子系统':group,'当前模型reaction ID':';'.join(filter(None,[x.get('2016 GPR反应',''),x.get('2024 GPR反应','')])) or '未知','变更类型':typ,'当前locus':locus,'当前protein ID':x.get('当前protein ID','未知'),'旧AFE locus':x.get('旧AFE locus','未知'),'旧RefSeq locus':'未知','旧protein ID':'未知','KO':x.get('KEGG KO','未知'),'EC':x.get('KEGG EC','未知'),'Rhea':'未知','MetaCyc':'未知','KEGG reaction':x.get('KEGG reaction','未知'),'BioCyc reaction ID':x.get('BioCyc reaction ID','未知'),'D3功能判定':x.get('功能判定','未知'),'来源数据库':'D3-1-R NCBI/B1/KEGG/BioCyc','来源记录ID':rec,'原始文件':files,'数据库版本/组装':vers,'取得日期':dates,'同株支持':'ATCC 23270当前组装计算注释；不等于实验功能证据','证据类型':'同株注释交叉证据；需人工复核','A0化学计量/区室比较':'未找到可直接等价的A0反应；需按反应式与区室核对','模型影响':'可能新增GPR、反应定义或转运边界；本轮不改模型','growth影响':'未知；待审计','energy影响':'未知；待审计','redox影响':'未知；待审计','transport影响':'未知；待审计','待人工确认问题':q,'优先级':'高' if group in ('FeS_ETC','关键能量/氧化还原') else '中','最终状态':'阶段B候选-待D3 Reviewer与人工审计'}
rows=[row(g,x,i+1) for i,(g,x) in enumerate(targets)]
bdgk=[{'候选ID':'E4-1-B-BDGK-GPR','子系统':'糖代谢','当前模型reaction ID':'BDGK','变更类型':'当前同株新增GPR候选','当前locus':'RU820_RS13140','当前protein ID':'WP_012537424.1','旧AFE locus':'AFE_2841','旧RefSeq locus':'AFE_RS13045','旧protein ID':'WP_012537424.1','KO':'未知','EC':'未知','Rhea':'未知','MetaCyc':'未知','KEGG reaction':'未知','BioCyc reaction ID':'未知','D3功能判定':'D3未提供该反应的直接等价功能判定','来源数据库':'B1旧新映射；A0 BDGK','来源记录ID':'BDGK;A0 2016=A102:P102;2024=A101:O101','原始文件':'02_NCBI/标准化数据/旧新基因映射.tsv;08_模型基线/A0/2016_到2024_真实差异.tsv','数据库版本/组装':'GCF_049532655.1','取得日期':'2026-09-23','同株支持':'当前locus/蛋白精确映射；功能待审计','证据类型':'旧GPR/current locus更新','A0化学计量/区室比较':'2016与2024反应式同底物同区室；GPR由空白到AFE_2841','模型影响':'仅审计GPR，不改正式模型','growth影响':'未知','energy影响':'未知','redox影响':'未知','transport影响':'未知','待人工确认问题':'确认AFE_2841对应RU820_RS13140；严禁误配AFE_RS13140（AFE_2860 ArsC）','优先级':'高','最终状态':'阶段B候选-待人工审计'}, {'候选ID':'E4-1-B-BDGK-方向','子系统':'糖代谢','当前模型reaction ID':'BDGK','变更类型':'反应方向/上下界冲突','当前locus':'RU820_RS13140','当前protein ID':'WP_012537424.1','旧AFE locus':'AFE_2841','旧RefSeq locus':'AFE_RS13045','旧protein ID':'WP_012537424.1','KO':'未知','EC':'未知','Rhea':'未知','MetaCyc':'未知','KEGG reaction':'未知','BioCyc reaction ID':'未知','D3功能判定':'未裁决','来源数据库':'A0真实差异','来源记录ID':'BDGK;A0 2016=A102:P102;2024=A101:O101','原始文件':'08_模型基线/A0/2016_到2024_真实差异.tsv;01_原始资料/2016_Campodonico/mmc1.xls;01_原始资料/2024_Khaleque/from2024-mmc1.xlsx','数据库版本/组装':'2016/2024模型','取得日期':'2026-09-23','同株支持':'未知','证据类型':'2016→2024作者修改','A0化学计量/区室比较':'2016可逆 -1000..1000；2024正向 0..1000，同区室','模型影响':'可能改变可行方向与通量；本轮不改模型','growth影响':'未知','energy影响':'未知','redox影响':'未知','transport影响':'未知','待人工确认问题':'确认方向/上下界，等待D3 Reviewer；GPR问题与本行分开','优先级':'高','最终状态':'阶段B候选-待人工审计'}]
rows += bdgk
write(OUT/'E4-1_候选修改.tsv',rows,FIELDS)
write(OUT/'E4-1_新增宿主反应候选.tsv',[r for r in rows if '当前同株新增' in r['变更类型']],FIELDS)
write(OUT/'E4-1_新增转运候选.tsv',[r for r in rows if '转运' in r['子系统']],FIELDS)
write(OUT/'E4-1_GPR更新候选.tsv',[r for r in rows if 'GPR' in r['变更类型']],FIELDS)
write(OUT/'E4-1_反应定义冲突.tsv',[r for r in rows if '冲突' in r['变更类型'] or r['子系统'] in ('中央碳代谢','FeS_ETC')],FIELDS)
write(OUT/'E4-1_数据库冲突.tsv',[r for r in rows if '分歧' in r['待人工确认问题'] or r['子系统']=='关键能量/氧化还原'],FIELDS)
write(OUT/'E4-1_标识映射不确定.tsv',[r for r in rows if r['旧AFE locus']=='未知' or r['旧protein ID']=='未知'],FIELDS)
write(OUT/'E4-1_待补文献.tsv',[{**r,'文献':'待补同株功能/复合体实验文献；D3注释为计算证据'} for r in rows],FIELDS)
write(OUT/'E4-1_来源登记.tsv',[{'来源文件':str(p.relative_to(ROOT)),'状态':'D3 Reviewer复查中；可重跑读取当前文件'} for p in D3.glob('D3-1-R_*.tsv')],['来源文件','状态'])
Path(OUT/'E4-1_来源覆盖与缺失.md').write_text(f'# 阶段B来源覆盖与缺失\n\n- 当前D3-1-R筛选候选 {len(rows)} 条：中央碳4、有机碳转运7、FeS/ETC7，另含关键能量/氧化还原筛选与BDGK两条。\n- D3 Reviewer仍在复查；最终状态不得PASS。脚本 `build_phase_b.py` 可在D3文件更新后重跑。\n- 所有注释均为计算证据；反应等价需按化学计量、区室和底物人工裁决。\n',encoding='utf-8')
# 统一主表：先重建阶段A，再合并阶段B，BDGK只保留阶段A的一对记录。
subprocess.run(['python', str(OUT/'build_phase_a.py')], check=True)
with open(OUT/'阶段A_候选修改.tsv',encoding='utf-8-sig',newline='') as f: arows=list(csv.DictReader(f,delimiter='\t'))
UNIFIED=['候选ID','子系统','当前模型reaction ID','变更类型','2016值','2024值','当前候选值','当前候选','旧gene','新gene','当前gene/locus','旧protein','新protein','KO','EC','Rhea','MetaCyc','来源数据库','来源记录ID','源记录','原始文件','定位','sheet/行号或记录定位','版本','数据库版本','日期','下载/检索日期','同株支持','证据类型','文献ID','文献','共源','是否与其他数据库共用来源','模型影响','growth影响','energy影响','redox影响','transport影响','是否影响growth/energy/redox/transport','A0化学计量/区室比较','待人工确认问题','优先级','最终状态','层级分类','原始字段保留']
def normalize(r,stage):
    z={k:'' for k in UNIFIED}; z.update(r)
    z['当前候选值']=z.get('当前候选','未知') or '未知'; z['文献ID']=z.get('文献ID','未知') or '未知'; z['共源']=z.get('共源','未知') or '未知'
    if stage=='A': z['层级分类']=z.get('层级分类','阶段A')
    else:
        z['2016值']='无/模型未覆盖（待化学核验）'; z['2024值']='无/模型未覆盖（待化学核验）'; z['当前候选值']=z.get('D3功能判定','未知')
        z['当前候选']=z['当前候选值']; z['层级分类']='阶段B-D3筛选'
        locus=z.get('当前locus','');
        z['旧gene']=z.get('旧AFE locus','未知') or '未知'; z['当前gene/locus']=locus or '未知'; z['新gene']=locus or '未知'
        z['sheet/行号或记录定位']=f"D3-1-R全基因交叉引用.tsv 当前locus={locus}; 注释证据表.tsv 当前locus={locus}"
        z['文献ID']='未知（注释证据，待补功能文献）'; z['共源']='D3跨数据库交叉引用；是否独立来源待Reviewer复核'
        if locus=='RU820_RS07700': z['A0化学计量/区室比较']='A0已有TKT1/TKT2 GPR；仅作新增GPR审计，不新增反应'
        elif locus=='RU820_RS08690': z['A0化学计量/区室比较']='A0已有PPC；仅作新增GPR审计，不新增反应'
        elif locus=='RU820_RS09355': z['A0化学计量/区室比较']='A0 6PGDH反应定义待核；RXN-9952化学计量/辅因子/区室尚未等价确认'
        elif locus=='RU820_RS03190': z['A0化学计量/区室比较']='A0对应反应需按苹果酸底物、辅因子与区室复核；KEGG/BioCyc分歧'
        elif locus=='RU820_RS00240': z['A0化学计量/区室比较']='A0 TSQOC/相关反应需核对复合体与区室；当前仅DoxD家族证据'
        else: z['A0化学计量/区室比较']='未找到可直接等价的A0反应；需按反应式与区室核对'
    return z
unified=[normalize(r,'A') for r in arows]
seen={r['候选ID'] for r in unified}
for r in rows:
    if r['候选ID'] not in seen: unified.append(normalize(r,'B')); seen.add(r['候选ID'])
write(OUT/'E4-1_候选修改.tsv',unified,UNIFIED)
write(OUT/'E4-1_新增宿主反应候选.tsv',[r for r in unified if r['变更类型']=='新增宿主反应候选'],UNIFIED)
write(OUT/'E4-1_新增转运候选.tsv',[r for r in unified if '转运' in r['子系统'] and r['变更类型'] in ('新增宿主反应候选','当前同株新增转运候选')],UNIFIED)
write(OUT/'E4-1_GPR更新候选.tsv',[r for r in unified if 'GPR' in r['变更类型']],UNIFIED)
write(OUT/'E4-1_反应定义冲突.tsv',[r for r in unified if '冲突' in r['变更类型'] or '分歧' in r['变更类型']],UNIFIED)
Path(OUT/'E4-1_验收报告.md').write_text(f'# E4-1 阶段B验收报告\n\n统一主表 {len(unified)} 条：阶段A {len(arows)} 条 + 阶段B去重候选 {len(unified)-len(arows)} 条；候选ID唯一。阶段B精选中央碳4、有机碳转运7、FeS/ETC7，并按功能证据细分新增反应、新增GPR、证据增强/无法判定。D3-1-R仍待Reviewer复查，验收为等待依赖，不得PASS。\n',encoding='utf-8')
Path(OUT/'E4-1_候选统计.json').write_text(json.dumps({'stage':'B','n':len(rows),'central':4,'transport':7,'fes_etc':7,'d3_reviewer':'待复查','rerunnable':True},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'n':len(rows),'central':4,'transport':7,'fes_etc':7},ensure_ascii=False))
