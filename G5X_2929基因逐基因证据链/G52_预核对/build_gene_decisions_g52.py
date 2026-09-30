import json
import re
from pathlib import Path

DIR=Path(__file__).resolve().parent
data=json.loads((DIR/'g52_precheck.json').read_text(encoding='utf-8'))
manual=json.loads((DIR/'g52_manual_adjudications.json').read_text(encoding='utf-8'))
balance={r['reaction']:r for r in json.loads((DIR/'g52_2016_balance.json').read_text(encoding='utf-8'))}
literature=json.loads((DIR/'g52_literature_search.json').read_text(encoding='utf-8'))['by_afe']
wp_literature=json.loads((DIR/'g52_wp_literature_search.json').read_text(encoding='utf-8'))['by_wp']
function_literature=json.loads((DIR/'g52_function_literature_search.json').read_text(encoding='utf-8'))['by_locus']

structural=re.compile(r'ribosomal|transcription|translation|tRNA|rRNA|RNA polymerase|DNA polymerase|DNA repair|DNA translocase|recombinase|integrase|antitoxin|toxin|sigma factor|ribosome|chaperon|elongation factor|initiation factor|nuclease|helicase|releasing factor|release factor|maturation|Sec[YE]|BamD|LolA|LptE|FtsK|GroE|preprotein|peptid|protein adenylyl|protein-L-isoaspartate|glycosylase|polyribonucleotide|primosomal|Clp protease|apolipoprotein|metalloprotease',re.I)
transport=re.compile(r'transporter|permease|antiporter|channel|porin|outer membrane|export|uptake|efflux|TolC|ExbD|MotA|MurJ|YhdP',re.I)
catalytic=re.compile(r'ase\b|reductase|oxidase|transferase|synthase|hydrolase|kinase|phosphatase|peptidase|glyoxalase|methyltransferase|rhodanese|lipase|dehydrogenase|mutase|lyase|epimerase',re.I)

rows=[]
for r in data['records']:
    loc=r['当前locus']; oldids=[x for x in r['2016反应ID'].split(';') if x]
    b=[balance[id] for id in oldids if id in balance]
    m=manual.get(loc)
    afe_hit=literature.get(r['旧AFE'],{}).get('hitCount',0)
    wp_hit=wp_literature.get(r['WP'],{}).get('hitCount',0)
    if m:
        decision=m['recommendation']
        kind='人工重点裁决'
        chemistry=m['chemistry']
        issue=m['issue']
        sources=m['evidence']
    elif oldids:
        problems=[x for x in b if x['status']!='balanced']
        if problems:
            decision='Hold 2016 精确式子；保留基因功能线索，先解决代谢物式子/电荷和反应对应。'
            kind='旧式子化学未决'
            issue='; '.join(x['reaction']+':'+x['status'] for x in problems)
        elif r['KEGG live reaction'] or r['UniProt reaction Rhea']:
            decision='维持 2016 历史 GPR 候选；跨库注释支持反应家族，当前 WP 的精确底物、辅酶及方向缺少独立实验证实。本轮不提升分数、不改式。'
            kind='历史关联有跨库支持'
            issue='2016 式子已按补表配平；跨库反应支持家族功能，但未独立证实本 WP 的精确底物、辅酶、方向及区室。'
        else:
            decision='Hold 2016 基因关联；当前 WP 尚无独立的精确催化反应交叉引用。'
            kind='历史关联缺少现行反应支持'
            issue='只有 2016 模型关联和当前蛋白家族注释。'
        chemistry='2016: '+r['2016反应ID']+'；KEGG: '+r['KEGG live reaction']+'；Rhea: '+r['UniProt reaction Rhea']
        sources='G52-S001;G52-S004;G52-S006;G52-S007;G52-S013;G52-S015;G52-S016;G52-S020'
    elif r['KEGG live reaction'] or r['UniProt reaction Rhea'] or r['D3 BioCyc reaction ID']:
        if structural.search(r['当前product']):
            decision='暂不新增小分子代谢反应；该链接涉及 DNA/RNA/蛋白质加工，需另按生物量与聚合反应范围评估。'
            kind='大分子过程'
        elif transport.search(r['当前product']):
            decision='Hold 转运反应；数据库关联尚未确定本 WP 的转运底物、方向与膜两侧代谢物。'
            kind='转运底物未定'
        else:
            decision='Hold 新增反应候选；已有数据库反应链接，尚缺本 WP 精确底物、区室和网络代谢物核查。'
            kind='无旧GPR的反应候选'
        chemistry='KEGG: '+r['KEGG live reaction']+'；Rhea: '+r['UniProt reaction Rhea']+'；BioCyc: '+r['D3 BioCyc reaction ID']
        issue='旧模型无该 AFE GPR；数据库链接不能直接作为 2026 催化证明。'
        sources='G52-S001;G52-S003;G52-S004;G52-S006;G52-S007;G52-S013;G52-S015'
    elif structural.search(r['当前product']):
        decision='本轮不新增独立小分子代谢反应；当前功能属于转录、翻译、修复或蛋白加工，需按生物量与聚合反应范围另核。'
        kind='大分子过程'
        chemistry='无独立小分子反应'
        issue='2016 无对应 GPR；现行资料未确定可直接落模的小分子化学式。'
        sources='G52-S001;G52-S002;G52-S003;G52-S004;G52-S006;G52-S007;G52-S013;G52-S015'
    elif transport.search(r['当前product']):
        decision='Hold 转运反应；当前注释未确定可落模的具体底物、方向和跨膜区室。'
        kind='转运底物未定'
        chemistry='未确认精确转运式子'
        issue='当前 WP/产品名为转运相关，2016 无对应 GPR，未检出 KEGG/Rhea 精确反应。'
        sources='G52-S001;G52-S002;G52-S003;G52-S006;G52-S007;G52-S013;G52-S015'
    elif catalytic.search(r['当前product']):
        decision='Hold 催化功能候选；家族/产品名不足以指定底物与化学式，暂不新增反应。'
        kind='酶家族底物未定'
        chemistry='未确认精确反应'
        issue='当前蛋白注释提示催化家族，2016 无对应 GPR，未检出 KEGG/Rhea 精确反应。'
        sources='G52-S001;G52-S002;G52-S003;G52-S006;G52-S007;G52-S013;G52-S015'
    else:
        decision='本轮不新增化学计量反应；当前资料未提供此 WP 可确认的小分子转化。'
        kind='无可确认小分子反应'
        chemistry='无已确认反应'
        issue='2016 无对应 GPR；现行 KEGG/Rhea/BioCyc 未检出精确反应。'
        sources='G52-S001;G52-S002;G52-S003;G52-S004;G52-S006;G52-S007;G52-S013;G52-S015'
    func_hit=function_literature[loc]['hit_count']
    extra=['G52-S023']
    if r['UniProt reaction Rhea'] or r['D3 Rhea交叉引用']:
        extra.extend(['G52-S025','G52-S026'])
    source_ids=';'.join(dict.fromkeys(sources.split(';')+extra))
    rows.append({'总序号':r['总序号'],'当前locus':loc,'WP':r['WP'],'旧AFE':r['旧AFE'],'基因功能':r['当前product'],'决策类别':kind,'2026建议':decision,'化学边界':chemistry,'证据问题':issue,'2016反应ID':r['2016反应ID'],'KEGG reaction':r['KEGG live reaction'],'Rhea':r['UniProt reaction Rhea'],'同WP蛋白组检出':r['2019同WP命中'],'旧AFE文献索引命中':afe_hit,'当前WP文献索引命中':wp_hit,'功能词文献索引命中':func_hit,'来源ID':source_ids,'身份依据':r['映射依据'] or '现行 RefSeq RU820/WP 直接身份','备注':'文献索引阴性不可证明无研究；Hold 是本次资料条件下的审查结论'})

assert len(rows)==293
assert [x['总序号'] for x in rows]==list(range(294,587))
assert len({x['当前locus'] for x in rows})==293
(DIR/'g52_gene_decisions_draft.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
from collections import Counter
print(json.dumps({'rows':len(rows),'classes':dict(Counter(x['决策类别'] for x in rows))},ensure_ascii=False))
