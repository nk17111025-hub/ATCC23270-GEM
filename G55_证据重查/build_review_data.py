import json
import re
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

import openpyxl

HERE = Path(__file__).parent
ROOT = HERE.parent
genes = json.loads((HERE/'g55_data.json').read_text(encoding='utf-8'))
uniprot = json.loads((HERE/'uniprot_g55.json').read_text(encoding='utf-8'))['by_wp']
literature = json.loads((HERE/'literature_hits.json').read_text(encoding='utf-8'))['results']
wb = openpyxl.load_workbook(ROOT/'00_项目导航与最终成果/ATCC23270-2026_最终模型证据表.xlsx', read_only=True, data_only=True)
alpha = {r[0]:list(r[:16]) for r in list(wb.worksheets[0].values)[2:] if r[0] and r[16]=='2016 original'}

dbsource = {
  'NCBI':'G55-S01','B1':'G55-S02','UniProt':'G55-S03','KEGG':'G55-S04',
  'BioCyc':'G55-S05','Rhea':'G55-S06','iMC507':'G55-S07','Europe PMC':'G55-S08',
}
papers = {
 '10.1074/jbc.M400596200':'G55-S09',
 '10.1021/acschembio.0c00791':'G55-S10',
 '10.1074/jbc.M400597200':'G55-S11',
 '10.1007/s00284-018-1453-9':'G55-S12',
 '10.1128/AEM.03057-12':'G55-S13',
 '10.3389/fmicb.2019.00592':'G55-S14',
 '10.1007/s10529-007-9488-1':'G55-S15',
 '10.1186/1471-2164-9-597':'G55-S16',
 '10.1016/j.meteno.2016.03.003':'G55-S17',
}

nonmetabolic = re.compile(r'(transposase|integrase|recombinase|toxin|antitoxin|restriction|DNA[- ]|RNA[- ]|ribosomal|transcriptional regulator|response regulator|sigma factor|helicase|nuclease|phage|conjugative|secretion system|VirB|Trb[A-Z]|single-stranded DNA|HU family DNA|integration host|methyltransferase|methylase|Hsp20|chaperone|protease|peptidase|protein export|signal peptidase|tRNA|translation|replication|cell division|pilus|pili|flagell|cytoskeletal|ribosome|ATP-dependent proteolysis|heat shock|mercury resistance regulator)',re.I)
hold_override = {
  1206:'SAM依赖型甲基转移酶的受体底物未明确，不能仅凭家族名称排除小分子代谢作用。',
  1249:'PhoB为调控蛋白，旧Pitpp把其列为必需转运成员；同簇PhoR/PhoU/Ppx的角色也须独立复核。',
  1251:'Pitpp的PstS/PstC/PstA/PstB分支需与另一分支分别验证；旧GPR混入调控蛋白。',
  1252:'Pitpp的PstS/PstC/PstA/PstB分支需与另一分支分别验证；旧GPR混入调控蛋白。',
  1253:'Pitpp的PstS/PstC/PstA/PstB分支需与另一分支分别验证；旧GPR混入调控蛋白。',
  1254:'Pitpp的PstS/PstC/PstA/PstB分支需与另一分支分别验证；旧GPR混入调控蛋白。',
  1256:'Ppx/GppA是磷酸酶家族，旧Pitpp把其列为必需转运成员。',
  1276:'ACCOAC旧GPR以OR连接AccA/B/C/D，BIOC1以AND连接同组；须核实复合体与两反应是否重复。',
  1278:'DHORTS旧GPR以AND连接3个二氢乳清酸酶样蛋白；各蛋白是否为必需亚基未证。',
  1303:'FE3tpp的两套ABC候选底物和必需亚基未证，第二分支含钼酸盐结合蛋白。',
  1304:'FE3tpp的两套ABC候选底物和必需亚基未证，第二分支含钼酸盐结合蛋白。',
  1305:'FE3tpp的两套ABC候选底物和必需亚基未证，第二分支含钼酸盐结合蛋白。',
  1307:'当前注释为钼酸盐结合蛋白；旧FE3tpp称铁(III)转运，底物冲突。',
  1309:'FE3tpp第二分支缺明确ATP结合亚基且铁(III)底物未证。',
  1321:'CYSDSS旧分数4所引PMID 17660945研究IscU，不是此NifS精确反应。',
  1329:'NITF旧GPR将NifX/NifEN和铁硫结合蛋白列为催化必需成员；应区分成熟过程与催化复合体。',
  1330:'NITF旧GPR将NifX/NifEN和铁硫结合蛋白列为催化必需成员；应区分成熟过程与催化复合体。',
  1331:'NITF旧GPR将NifX/NifEN和铁硫结合蛋白列为催化必需成员；应区分成熟过程与催化复合体。',
  1332:'NITF旧GPR将铁硫结合蛋白列为催化必需成员；缺少实体复合物证据。',
  1333:'NITF旧GPR将铁硫结合蛋白列为催化必需成员；缺少实体复合物证据。',
  1334:'NITF结构亚基较明确，但旧反应使用NADH供电子，数据库标准式使用还原型铁氧还蛋白。',
  1335:'NITF结构亚基较明确，但旧反应使用NADH供电子，数据库标准式使用还原型铁氧还蛋白。',
  1336:'NITF结构亚基较明确，但旧反应使用NADH供电子，数据库标准式使用还原型铁氧还蛋白。',
  1349:'CYSDSS旧分数4所引PMID 17660945研究IscU，不能证明本蛋白反应。',
  1437:'仅有胰蛋白酶样肽酶结构域，底物、位置及是否属于模型蛋白降解反应均未确定。',
  1442:'FdhF/YdeP家族不证明旧FDH的NAD+受体、方向或辅助亚基。',
  1452:'旧GLYCH以O2产H2O2，UniProt/Rhea仅给未指定受体A；旧名称与EC也不一致。',
  1453:'旧GLYCH以O2产H2O2，UniProt/Rhea仅给未指定受体A；旧名称与EC也不一致。',
  1454:'旧GLYCH以O2产H2O2，UniProt/Rhea仅给未指定受体A；旧名称与EC也不一致。',
  1461:'DHORTS旧GPR以AND连接3个二氢乳清酸酶样蛋白；各蛋白是否为必需亚基未证。',
  1465:'FBA3旧GPR将两种醛缩酶写成AND；必需复合体关系未证。',
}
formal = {
  1228:('修改或补充 GPR','GMHEPAT: AFE_1406→空；GMHEPK: 空→AFE_1406；GMHEPK EC 5.3.1.28→2.7.1.167。原反应式和培养边界保持2016原值。','K21344/RHEA:27473支持7-磷酸激酶；K21345是独立腺苷酸转移酶。AFE_1406缺后者证据。'),
  1235:('修改或补充 GPR','ARGDC: AFE_1417→AFE_1471；原反应式保留。','AFE_1417为GNAT/PF13508，无脱羧酶结构域；AFE_1471为K01585/EC4.1.1.19且UniProt有RHEA:17641自动推注。'),
  1284:('修改或补充 GPR','ARGDC: AFE_1417→AFE_1471；原反应式保留。','本基因为AFE_1471，属生物合成型精氨酸脱羧酶；尚无同株纯化酶实验。'),
  1308:('修改或补充 GPR','DB4PS: (AFE_1494 or AFE_0299)→AFE_0299；原反应式保留。','AFE_1494为ABC底物结合蛋白/PF01497；AFE_0299经B1精确映射为RibB。FE3tpp关联另列待核。'),
}
primary = {
  1243:('10.1007/s00284-018-1453-9','G55-S12'),
  1244:('10.1007/s00284-018-1453-9','G55-S12'),
  1270:('10.1074/jbc.M400597200','G55-S11'),
  1271:('10.1074/jbc.M400596200','G55-S09;G55-S10'),
  1272:('10.1074/jbc.M400596200','G55-S09'),
  1364:('10.3389/fmicb.2019.00592','G55-S14'),
  1456:('10.1128/AEM.03057-12','G55-S13;G55-S14'),
}

candidate_equations={
  1239:'RHEA:25213 UDP-2-N,3-O-bis[(3R)-3-hydroxytetradecanoyl]-alpha-D-glucosamine + H2O = lipid X + UMP + 2 H(+); RHEA:67824 为泛化酰基式',
  1242:'RHEA:47612 a 1,2-diacyl-sn-glycerol + UDP-alpha-D-glucose = a 1,2-diacyl-3-O-(alpha-D-glucopyranosyl)-sn-glycerol + UDP + H(+)',
  1286:'RHEA:26168 3-fumarylpyruvate + H2O = fumarate + pyruvate + H(+)',
  1288:'RHEA:16437 RX + glutathione = an S-substituted glutathione + a halide anion + H(+); RX 为未指定卤代底物',
  1320:'RHEA:12929 acetyl-CoA + 2-oxoglutarate + H2O = (2R)-homocitrate + CoA + H(+)',
  1358:'RHEA:13213 succinate semialdehyde + NADP(+) + H2O = succinate + NADPH + 2 H(+); RHEA:13217 则使用 NAD(+)',
  1362:'RHEA:16125 D-fructose + ATP = D-fructose 6-phosphate + ADP + H(+)',
  1363:'RHEA:22172 beta-D-fructose 6-phosphate + UDP-alpha-D-glucose = sucrose 6(F)-phosphate + UDP + H(+)',
  1364:"RHEA:55092 D-fructose + UDP-alpha-D-glucose = sucrose + UDP + H(+); RHEA:16241 为泛化 NDP-glucose 式",
  1447:'RHEA:22052 A + NH4(+) + H2O = hydroxylamine + AH2 + H(+); A/AH2 为未指定电子载体，方向未定',
  1463:'RHEA:36407 L-threonine + hydrogencarbonate + ATP = L-threonylcarbamoyladenylate + diphosphate + H2O',
}
candidate_boundaries={
  1239:'LpxH蛋白家族与脂质A水解步骤相符；旧模型缺对应特定酰基载体/产物，需核酰基链长。',
  1242:'糖基转移酶家族不足以证明具体二酰甘油底物和产物。',
  1286:'水解酶家族与3-fumarylpyruvate专一性仍需验证。',
  1288:'谷胱甘肽转移酶受体底物RX未定，无法写具体模型反应。',
  1320:'NifV/同源家族支持同化学类型；旧模型缺(2R)-homocitrate池和对应辅因子去向。',
  1358:'NAD(+)与NADP(+)受体候选并存，不能同时作本株精确反应。',
  1362:'PfkB家族有多种糖底物；仅Rhea/KO推注不足以确认果糖专一性。',
  1363:'NCBI产品名为HAD-IIB hydrolase，UniProt/KEGG推注为蔗糖-6-磷酸合酶；存在功能注释冲突。',
  1364:'同株蛋白组检出并报告丰度变化；论文未测UDP-glucose+果糖的酶活。',
  1447:'标准式电子载体A/AH2未指定，模型不能据此选定供电子体和方向。',
  1463:'tRNA修饰前体合成步骤，旧模型未表示tRNA物种与下游t6A途径；UniProt反应为自动推注。',
}

gene_headers=['总序号','Current locus','WP protein_id','坐标','链向','NCBI product','旧AFE','映射状态','映射依据','同WP拷贝数','UniProt accession','UniProt审校','UniProt蛋白证据','Pfam','InterPro','KEGG KO','KEGG EC','KEGG reaction','BioCyc reaction','Rhea交叉引用','2016 linked Reaction ID','2016原分数','是否找到精确reaction','2026 action','具体差异或结论','证据边界','未决点','关联行ID','Source ID','Primary evidence DOI','NCBI链接','UniProt链接','EuropePMC检索命中数','EuropePMC链接','审查状态']
gene_rows=[]
review_rows=[]
candidate_rows=[]
wp_counts=Counter(z['wp'] for z in genes)
for z in genes:
    n=z['num']; x=z['cross']; u=uniprot.get(z['wp'],[]); up=u[0] if u else {}
    old=z['old_reactions']; product=unquote(z['product']); oldids=[r['id'] for r in old]
    doi, paper_ids=primary.get(n,('',''))
    hit_count=sum(literature.get(a,{}).get('hit_count',0) for a in z['afe'].split(';'))
    paper_url='; '.join(literature[a]['url'] for a in z['afe'].split(';') if a in literature)
    sources=['G55-S01','G55-S08']
    if z['afe']: sources.append('G55-S02')
    if up: sources.append('G55-S03')
    if x['KEGG KO'] or x['KEGG EC'] or x['KEGG reaction']: sources.append('G55-S04')
    if x['BioCyc gene ID'] or x['BioCyc reaction ID']: sources.append('G55-S05')
    if x['Rhea精确交叉引用']: sources.append('G55-S06')
    if old: sources.append('G55-S07')
    if paper_ids: sources.extend(paper_ids.split(';'))
    sources=list(dict.fromkeys(sources))
    status='证据链已建；反应结论保守复核'
    if n in formal:
        action,diff,reason=formal[n]
        boundary='直接核实了当前CDS、旧新映射和旧GPR；反应基因归属由相容的结构域/KO推定，尚无同株纯化酶确认。'
        unresolved=reason+'；实施前核查整体模型质量守恒与相关基因跨组关系。'
    elif n in hold_override:
        action='证据不足暂缓'; diff=('旧反应保留为历史建模假设；本轮不提交未证化学式或GPR。' if old else '功能可能涉及代谢，但目前无法确定精确反应；本轮不新增反应。')
        boundary=('当前CDS、旧GPR和数据库注释已核；精确底物、辅因子、复合体或旧文献对应关系存在冲突。' if old else '当前CDS和数据库注释已核；缺少基因专属底物、产物和反应实验。')
        unresolved=hold_override[n]
    elif n==1270:
        action='仅更新证据或编号'; diff='补入Afe LpxA对UDP-GlcNAc3N的体外底物选择性证据；旧UAGAAT式暂不改。'
        boundary='重组Afe LpxA糖核苷酸选择性直接测得；旧式指定的酰基供体及精确化学等价仍未核定。'
        unresolved='核对旧式3hmrsACP与原始试验酰基供体的链长、定位和方向。'
    elif n==1271 or n==1272:
        action='仅更新证据或编号'; diff='AGAD/UACAA反应本身已在2016表；补入同蛋白生化文献，反应分数3→4。'
        boundary='原始论文对重组GnnA/GnnB联合反应及产物作体外测量；区室、体内通量和模型边界未直接测得。'
        unresolved='保留2016方向/区室为模型假设；核对原文单酶中间体定量与反应式质子约定。'
    elif n==1364:
        action='证据不足暂缓'; diff='补入同株AFE_1552蛋白组检出及丰度变化；蔗糖合酶反应暂不新增。'
        boundary='原始论文测得蛋白丰度；其文中蔗糖合酶催化叙述基于功能注释，未测该蛋白底物和产物。'
        unresolved='需测定UDP-glucose/果糖反应，并定义蔗糖代谢物池与下游去向。'
    elif n in (1243,1244):
        action='仅更新证据或编号'; diff='记录tce簇Fe²⁺条件转录证据；不向rus反应加OR/AND。'
        boundary='同株论文直接测得Fe²⁺与硫培养转录差异；未测得纯化蛋白电子传递。'
        unresolved='Fe²⁺供体、受体、膜侧、复合体须实验确认。'
    elif n==1456:
        action='仅更新证据或编号'; diff='AFE_1667蛋白检出证据并入；TKT1/TKT2旧GPR不变。'
        boundary='同株蛋白组直接检出AFE_1667；未测定其TKT1/TKT2底物或独立催化。'
        unresolved='需纯化蛋白、底物实验和复合体鉴定。'
    elif old:
        action='仅更新证据或编号'; diff='登记当前RU820/WP与2016 AFE-GPR对应；旧反应式与边界保留。'
        boundary='官方CDS、映射及2016 GPR直接可复核；UniProt/KEGG/BioCyc多为计算注释，不能视作同株酶实测。'
        unresolved='精确反应的直接实验、方向、区室和GPR组装程度尚未逐项验证。'
    elif x['旧模型未覆盖功能候选']=='是' or x['转运候选']=='是' or x['映射冲突']=='是':
        action='证据不足暂缓'; diff='数据库功能候选；无可直接执行的正式反应或GPR修改。'
        boundary='当前蛋白与数据库命中已记录；未证实该基因的精确底物、产物、偶联或完整GPR。'
        unresolved='需同株催化/转运实验；数据库间继承与旧模型无GPR反应须逐项确认。'
    elif product.lower().startswith('hypothetical') or not z['afe']:
        action='证据不足暂缓'; diff='未定义可评分的精确代谢反应。'
        boundary='只核实当前CDS及可得序列/数据库注释；无同株反应级证据。'
        unresolved='需功能测定或高可信结构/邻域证据。'
    elif nonmetabolic.search(product):
        action='确认无需修改'; diff='现有证据指向DNA/RNA维护、调控、移动元件或分泌功能；未提出代谢反应。'
        boundary='基于当前产物和可得家族注释判断模型范围；未穷尽该蛋白的所有非代谢作用。'
        unresolved='若将来发现独立代谢/转运功能，重新评估。'
    else:
        action='证据不足暂缓'; diff='当前功能注释不能独自确定精确反应。'
        boundary='当前CDS与数据库交叉引用已记录；基因专属催化或转运未直接证实。'
        unresolved='需确定底物、产物、辅因子、方向及模型2016化学等价。'
    exact='2016模型精确定义' if old else ('候选标准反应（基因未实证）' if x['Rhea精确交叉引用'] or (up and up['catalytic']) else '未找到可确认的基因专属精确反应')
    rowid=f'G55-{n}'
    if old: rowid+=';'+ ';'.join(f'G55-{n}-{r["id"]}' for r in old)
    if n==1228: rowid+=';G55-1228-GMHEPK'
    if n==1284: rowid+=';G55-1284-ARGDC'
    accession=up.get('accession','')
    gen_row=[n,z['locus'],z['wp'],f'NZ_CP136162.1:{z["start"]}-{z["end"]}',z['strand'],product,z['afe'] or '未确认',x['B1映射状态'] or '未确认',x['B1映射方法'],wp_counts[z['wp']],accession,up.get('reviewed','未检索到同株条目'),up.get('protein_existence',''),';'.join(up.get('pfam',[])),';'.join(up.get('interpro',[])),x['KEGG KO'],x['KEGG EC'],x['KEGG reaction'],x['BioCyc reaction ID'],x['Rhea精确交叉引用'],';'.join(oldids),';'.join(str(r['score']) for r in old),exact,action,diff,boundary,unresolved,rowid,';'.join(sources),doi,f'https://www.ncbi.nlm.nih.gov/protein/{z["wp"]}',up.get('url',''),hit_count,paper_url,status]
    gene_rows.append(gen_row)
    # A-P stay blank for the gene-level record, preserving Alpha schema.
    ext=['gene review',f'G55-{n}','G55',';'.join(oldids),'',z['afe'] or '未确认',z['locus'],doi,';'.join(sources),boundary,action,'不评分',';'.join(sources),'gene-level，不把旧反应分数转给基因']
    review_rows.append([None]*16+ext)
    for r in old:
        rid=r['id']; score=r['score']; newscore='2'; score_reason='仅家族/KO/EC和旧模型支持；同株精确反应无直接酶实验。'
        if rid in ('AGAD','UACAA'):
            newscore='4'; score_reason='同源于旧AFE的重组GnnA/GnnB在体外测得相应底物与产物，DOI 10.1074/jbc.M400596200。'
        elif rid in ('NITF','FDH','GLYCH','FE3tpp','Pitpp','BIOC1'):
            newscore='1'; score_reason='精确供电子体/受体、底物或复合体归属没有基因专属证据；保留旧建模反应。'
        elif rid=='CYSDSS':
            newscore='2'; score_reason='旧PMID 17660945研究IscU而非本反应；仅NifS家族间接支持。'
        elif rid=='GMHEPAT' and n==1228:
            newscore='1'; score_reason='移除错误AFE_1406关联后，尚无本株腺苷酸转移酶候选；反应作为旧模型缺基因步骤保留。'
        elif rid in ('ARGDC','DB4PS') and n in (1235,1308):
            newscore='2'; score_reason='反应酶类间接支持，但旧GPR指向错误蛋白；反应评分不等于旧基因归属评分。'
        elif rid=='UAGAAT':
            newscore='2'; score_reason='原始LpxA论文测得糖核苷酸选择性，但旧式具体酰基供体及化学等价仍需核。'
        # Alpha old A-P: exact frozen values, separate review row.
        o=alpha[rid]
        ext2=['reaction review',f'G55-{n}-{rid}','G55',rid,rid,z['afe'] or '未确认',z['locus'],doi,';'.join(sources),boundary,action,newscore,';'.join(sources),f'2016={score}; 2026={newscore}; 评分仅限反应。{score_reason}']
        review_rows.append(o+ext2)
    if n==1228:
        o=alpha['GMHEPK']; review_rows.append(o+['reaction review','G55-1228-GMHEPK','G55','GMHEPK','GMHEPK',z['afe'],z['locus'],'',';'.join(sources),'K21344/RHEA:27473与旧GMHEPK式一致；同株直接酶实验未取得。','修改或补充 GPR','2',';'.join(sources),'2016=1; 2026=2；仅反应间接证据，GPR修改单独记录。'])
    if n==1284:
        o=alpha['ARGDC']; review_rows.append(o+['reaction review','G55-1284-ARGDC','G55','ARGDC','ARGDC',z['afe'],z['locus'],'',';'.join(sources),'K01585、RHEA:17641和脱羧酶结构域支持；缺同株直接催化实验。','修改或补充 GPR','2',';'.join(sources),'2016=2; 2026=2；反应级同源支持，不将其视为实验GPR验证。'])
    if x['旧模型未覆盖功能候选']=='是' and n not in (1284,1456):
        candidate_rows.append([n,z['locus'],z['wp'],product,x['KEGG KO'],x['KEGG EC'],x['KEGG reaction'],x['BioCyc reaction ID'],x['Rhea精确交叉引用'],up.get('catalytic',[]),candidate_equations.get(n,'无可确认的基因专属标准式'),candidate_boundaries.get(n,'转运底物、偶联机制或精确反应仍未确定。'),action,unresolved])

sources=[
 ['G55-S01','NCBI RefSeq GCF_049532655.1 genomic.gff; SHA-256 758a8b72423a824840a52e1e4184571474500c7fce155944a28e05fef6577091; https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/','2026','','official annotation','ATCC 23270','CDS/坐标/WP/product','G55-1173:G55-1465','CDS身份、坐标与当前注释','2 for annotation','product是计算注释，非反应实验'],
 ['G55-S02','B1 旧新基因映射.tsv；RefSeq protein_id exact和邻域证据；本地标准化数据','2026','','mapping','ATCC 23270','旧AFE→RU820映射','G55-1173:G55-1465','269条单值映射，2条重复/歧义','not scored','AFE_1360/1364双拷贝不得单值映射'],
 ['G55-S03','UniProtKB RefSeq_Protein ID mapping job pJOiC7Ywlr; UniProtKB entry JSON 2026-09-24; https://rest.uniprot.org/','2026','','protein database','ATCC 23270 taxonomy 243159','216条同株UniProt、Pfam/InterPro、证据代码','G55-1173:G55-1465','2条条目有ECO:0000269反应文献；其余多为自动注释','2 or 4 by reaction','Swiss-Prot审校本身不等于湿实验；需看ECO'],
 ['G55-S04','KEGG afr gene/KO/EC/reaction 本地D3-1-R提取；https://www.kegg.jp/kegg-bin/show_organism?org=afr','2026','','metabolic database','ATCC 23270 old assembly','KO/EC/反应交叉引用','G55-1173:G55-1465','功能和标准反应的间接关联','2','afr采用旧组装，gene→现行RU820须经B1核'],
 ['G55-S05','BioCyc GCF_000021485 PGDB 30.0 本地D3-1-R提取；https://biocyc.org/','2026','','Tier-3 computational PGDB','ATCC 23270 old assembly','gene/protein/reaction交叉引用','G55-1173:G55-1465','计算预测的反应对象','2','与MetaCyc共享Pathway Tools，未独立验证'],
 ['G55-S06','Rhea/ChEBI 本地D3-1-R EC交叉引用；https://www.rhea-db.org/','2026','','reaction database','not strain-specific','标准化反应、化学物种','G55-1173:G55-1465','标准反应式及方向信息','2','EC命中不证明本株基因催化'],
 ['G55-S07','Campodonico et al. iMC507 original Supplementary Table mmc1.xls and Alpha frozen A-P; https://doi.org/10.1016/j.meteno.2016.03.003','2016','10.1016/j.meteno.2016.03.003','model baseline','ATCC 23270','2016反应/GPR/上下界/原分数','G55-1173:G55-1465','旧模型记录，不等于每项反应实验','1-4 legacy','原表错误和假设保留作对照'],
 ['G55-S08','Europe PMC 271个AFE关键词检索，2026-09-24；https://www.ebi.ac.uk/europepmc/','2026','','literature search','mixed','各AFE标签全文索引','G55-1173:G55-1465','49个标签有检索命中','not scored','词检索无命中不等于无文献；命中须分辨直接实验与背景提及'],
 ['G55-S09','Sweet CR, Ribeiro AA, Raetz CRH. J Biol Chem 279:25400-25410 (2004). https://doi.org/10.1074/jbc.M400596200','2004','10.1074/jbc.M400596200','original biochemical experiment','A. ferrooxidans clone corresponding to AFE_1457/1458','重组GnnA/B+标记底物+NMR','AGAD;UACAA;G55-1271;G55-1272','UDP-GlcNAc→UDP-GlcNAc3N；NAD+与谷氨酸参与','4 for assayed reactions','未直接测体内通量、区室或模型边界'],
 ['G55-S10','Manissorn et al. ACS Chem Biol 15:3235-3243 (2020). https://doi.org/10.1021/acschembio.0c00791','2020','10.1021/acschembio.0c00791','original biochemical/structural','A. ferrooxidans GnnA','GnnA生化与结构','AGAD;G55-1271','GnnA反应和结构','4 for GnnA reaction','具体培养株编号需回原文核'],
 ['G55-S11','Sweet et al. J Biol Chem 279:25411-25419 (2004). https://doi.org/10.1074/jbc.M400597200','2004','10.1074/jbc.M400597200','original biochemical experiment','A. ferrooxidans LpxA clone','糖核苷酸底物选择性测定','UAGAAT;G55-1270','Afe LpxA偏好UDP-GlcNAc3N，UDP-GlcNAc低活性','2 for exact legacy model equation','未确认旧式所用酰基供体与实测完全相同'],
 ['G55-S12','Ai et al. Current Microbiology 75:818-826 (2018). https://doi.org/10.1007/s00284-018-1453-9','2018','10.1007/s00284-018-1453-9','original transcriptomics','ATCC 23270','Fe2+对S0转录、序列预测','G55-1243;G55-1244','AFE_1428/1429条件表达','2','未测纯化电子传递、方向或复合体'],
 ['G55-S13','Osorio et al. Appl Environ Microbiol 79:2172-2181 (2013). https://doi.org/10.1128/AEM.03057-12','2013','10.1128/AEM.03057-12','original proteomics','ATCC 23270','硫/氧与硫/Fe3+比较蛋白组','G55-1456','AFE_1667蛋白斑点变化','2 expression only','未测TKT催化和伙伴'],
 ['G55-S14','Bellenberg et al. Front Microbiol 10:592 (2019). https://doi.org/10.3389/fmicb.2019.00592','2019','10.3389/fmicb.2019.00592','original proteomics','DSM 14882 = ATCC 23270','铁培养/黄铁矿生物膜蛋白组','AFE_1667;AFE_1552','相应蛋白检出及丰度','2 expression only','蛋白检出不能证明具体反应'],
 ['G55-S15','Zeng et al. Biotechnol Lett 29:1965-1972 (2007). https://doi.org/10.1007/s10529-007-9488-1','2007','10.1007/s10529-007-9488-1','original biochemical experiment','A. ferrooxidans IscU','重组IscU表达纯化与表征','CYSDSS QC','IscU支架蛋白，不是NifS/IscS的CYSDSS精确反应','not applicable','2016 CYSDSS引用此PMID不能得到分数4'],
 ['G55-S16','Valdes et al. BMC Genomics 9:597 (2008). https://doi.org/10.1186/1471-2164-9-597','2008','10.1186/1471-2164-9-597','genome/pathway reconstruction','ATCC 23270','序列/通路分析','NITF;G55-1318:G55-1350','预测nif基因簇','2','非本株氮酶单酶生化表征'],
 ['G55-S17','Campodonico et al. Metab Eng Commun 3:84-96 (2016). https://doi.org/10.1016/j.meteno.2016.03.003','2016','10.1016/j.meteno.2016.03.003','model reconstruction','ATCC 23270','iMC507构建与培养场景拟合','all 2016 links','模型行为与旧信度','not independent','不能作为新基因反应级湿实验'],
]

issues=[
 ['G55-Q01','GMHEPAT/GMHEPK','RU820_RS06520','正式GPR/EC候选','旧AFE_1406被放在GMHEPAT；本蛋白是单功能7P激酶；GMHEPK缺GPR且EC误写5.3.1.28','GMHEPAT去AFE_1406；GMHEPK加AFE_1406；GMHEPK EC 2.7.1.167','G55-S01;G55-S04;G55-S06;G55-S07','https://www.rhea-db.org/rhea/27473'],
 ['G55-Q02','ARGDC','RU820_RS06560;RU820_RS06805','正式GPR候选','旧AFE_1417属GNAT，AFE_1471属脱羧酶家族','GPR AFE_1417→AFE_1471；保留2016化学式','G55-S01;G55-S03;G55-S04;G55-S07','https://www.rhea-db.org/rhea/17641'],
 ['G55-Q03','DB4PS','RU820_RS06925','正式GPR候选','旧AFE_1494是ABC结合蛋白；另一旧GPR AFE_0299是RibB','(AFE_1494 or AFE_0299)→AFE_0299','G55-S01;G55-S02;G55-S03;G55-S07','https://www.ncbi.nlm.nih.gov/protein/WP_012536073.1'],
 ['G55-Q04','CYSDSS','RU820_RS06990;RU820_RS07130','评分更正','2016原分4引PMID17660945，原文对象IscU，不是脱硫酶','精确反应本轮分2；暂不改化学式和GPR','G55-S07;G55-S15','https://doi.org/10.1007/s10529-007-9488-1'],
 ['G55-Q05','AGAD/UACAA','RU820_RS06740;RU820_RS06745','证据/评分更新','旧分3；同株GnnA/B重组酶反应和产物直接测得','反应分数4；保持原模型行只做扩展审查','G55-S09;G55-S10;G55-S07','https://doi.org/10.1074/jbc.M400596200'],
 ['G55-Q06','Pitpp','RU820_RS06630;RU820_RS06640:RU820_RS06665','暂缓/QC','旧AND含PhoB/PhoR/PhoU/Ppx；第二分支需复核完整性','保留旧反应；转运复合体候选PstS/C/A/B待验证','G55-S01;G55-S03;G55-S04;G55-S07',''],
 ['G55-Q07','NITF','RU820_RS07030:RU820_RS07065','暂缓/QC','旧名称sulfide quinone reductase；NADH受体/供体与标准铁氧还蛋白式不一致；GPR混入成熟蛋白','不直接改式；核实电子供体、成熟成员及菌株实测后再改','G55-S03;G55-S07;G55-S16','https://www.rhea-db.org/rhea/21448'],
 ['G55-Q08','GLYCH','RU820_RS07680:RU820_RS07690','暂缓/QC','旧名D-lactate dehydrogenase、EC1.1.2.4与式中的glycolate+O2不一致；Rhea受体A未定','保留旧低可信假设；实测电子受体后确定反应','G55-S03;G55-S04;G55-S07','https://www.rhea-db.org/rhea/21264'],
 ['G55-Q09','FE3tpp','RU820_RS06900:RU820_RS06930','暂缓/QC','第二ABC分支含钼酸盐结合蛋白且未明确ATP亚基；Fe3+特异性未证','转运底物/复合体实验前不写新GPR','G55-S01;G55-S03;G55-S07',''],
 ['G55-Q10','FDH','RU820_RS07625','暂缓/QC','FdhF/YdeP家族不能确定旧NAD+依赖式','查MoCo、电子受体及蛋白复合体','G55-S01;G55-S03;G55-S07',''],
 ['G55-Q11','DHORTS/FBA3/ACCOAC','RU820_RS06765;RU820_RS06775;RU820_RS07725;RU820_RS07745','暂缓/QC','旧GPR的AND/OR与当前亚基/同工酶注释可能冲突','需生化复合体和序列位点核查，不凭邻近改GPR','G55-S01;G55-S03;G55-S07',''],
 ['G55-Q12','AFE_1360/AFE_1364','RU820_RS06290;RU820_RS06310','映射歧义','两当前HEPN基因共享WP_012607030.1且B1一对多','旧AFE保持未确认，不并入确定GPR','G55-S01;G55-S02',''],
]

candidate_headers=['总序号','Current locus','WP','NCBI product','KEGG KO','KEGG EC','KEGG reaction','BioCyc reaction','Rhea','UniProt catalytic records','标准反应式（数据库）','反应级证据边界','2026 action','未决点']
related_by_source={
 'G55-S01':'G55-1173:G55-1465', 'G55-S02':'G55-1173:G55-1465',
 'G55-S03':'G55-1173:G55-1465', 'G55-S04':'G55-1173:G55-1465',
 'G55-S05':'G55-1173:G55-1465', 'G55-S06':'G55-1173:G55-1465',
 'G55-S07':'2016 linked reactions', 'G55-S08':'G55-1173:G55-1465',
 'G55-S09':'AGAD;UACAA', 'G55-S10':'AGAD', 'G55-S11':'UAGAAT',
 'G55-S12':'G55-1243;G55-1244', 'G55-S13':'G55-1456',
 'G55-S14':'G55-1456;G55-1364', 'G55-S15':'CYSDSS QC',
 'G55-S16':'NITF', 'G55-S17':'all 2016 links',
}
sources=[r[:8]+[related_by_source[r[0]]]+r[8:] for r in sources]
assert all(len(r)==12 for r in sources)
out={'gene_headers':gene_headers,'gene_rows':gene_rows,'alpha_headers':['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe²⁺ lb','Fe²⁺ ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub','Record class','Candidate ID','Candidate group','2016 linked Reaction ID','2026 reaction object','Legacy AFE locus','Current locus','Primary evidence DOI','Evidence type','Evidence boundary','2026 action','2026 reaction score','Source ID','Score scope'],'alpha_rows':review_rows,'source_headers':['Source ID','Full citation','Year','DOI','Source type','Strain','Experimental type','Genes / proteins studied','Related Reaction IDs / candidate IDs','What was directly measured','Applicable confidence level','Evidence limitation'],'source_rows':sources,'issue_headers':['Issue ID','Reaction/subject','Current locus','Class','Finding','Proposed action / experiment','Source IDs','Primary link'],'issue_rows':issues,'candidate_headers':candidate_headers,'candidate_rows':candidate_rows}
(HERE/'review_data.json').write_text(json.dumps(out,ensure_ascii=False,default=str),encoding='utf-8')
assert len(gene_rows)==293 and len(set(r[1] for r in gene_rows))==293
assert all(len(r)==30 for r in review_rows)
print('genes',len(gene_rows),'alpha review rows',len(review_rows),'sources',len(sources),'issues',len(issues),'candidates',len(candidate_rows))
print('actions',Counter(r[23] for r in gene_rows))
