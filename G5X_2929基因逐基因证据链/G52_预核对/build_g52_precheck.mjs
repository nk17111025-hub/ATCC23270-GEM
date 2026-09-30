import fs from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const dir = new URL('.', import.meta.url);
const data = JSON.parse(await fs.readFile(new URL('g52_precheck.json', dir), 'utf8'));
const keggEntries = JSON.parse(await fs.readFile(new URL('g52_kegg_reaction_entries.json', dir), 'utf8')).reactions;
const rheaEntries = JSON.parse(await fs.readFile(new URL('g52_rhea_reaction_entries.json', dir), 'utf8')).reactions;
const rheaKeggRows = JSON.parse(await fs.readFile(new URL('g52_rhea_kegg_xref_comparison.json', dir), 'utf8')).rows;
const rheaKeggByLocus = new Map(rheaKeggRows.map(r=>[r.locus,r]));
const manual = JSON.parse(await fs.readFile(new URL('g52_manual_adjudications.json', dir), 'utf8'));
const literature = JSON.parse(await fs.readFile(new URL('g52_literature_search.json', dir), 'utf8')).by_afe;
const wpLiterature = JSON.parse(await fs.readFile(new URL('g52_wp_literature_search.json', dir), 'utf8')).by_wp;
const functionLiterature = JSON.parse(await fs.readFile(new URL('g52_function_literature_search.json', dir), 'utf8')).by_locus;
const balanceRows = JSON.parse(await fs.readFile(new URL('g52_2016_balance.json', dir), 'utf8'));
const balance = new Map(balanceRows.map(r=>[r.reaction,r]));
const reactionOverrides = JSON.parse(await fs.readFile(new URL('g52_reaction_overrides.json', dir), 'utf8'));
const geneDecisions = JSON.parse(await fs.readFile(new URL('g52_gene_decisions_draft.json', dir), 'utf8'));
const decisionsByLocus = new Map(geneDecisions.map(r=>[r['当前locus'],r]));
const chemistryComparisons = JSON.parse(await fs.readFile(new URL('g52_2016_kegg_chemistry.json', dir), 'utf8'));
const chemistryByPair = new Map(chemistryComparisons.map(r=>[`${r.locus}|${r.old_reaction}`,r]));
const candidateDraft = JSON.parse(await fs.readFile(new URL('g52_candidate_rows_draft.json', dir), 'utf8'));
const tcdbAll = JSON.parse(await fs.readFile(new URL('g52_tcdb_homology.json', dir), 'utf8')).by_locus;
const tcdbFamily = JSON.parse(await fs.readFile(new URL('g52_tcdb_family_validation.json', dir), 'utf8')).by_locus;
const brenda = JSON.parse(await fs.readFile(new URL('g52_brenda_ec_references.json', dir), 'utf8')).by_ec;
const candidateByLocus = new Map();
for(const c of candidateDraft) candidateByLocus.set(c.locus,[...(candidateByLocus.get(c.locus)||[]),c.candidate_id]);
const book = Workbook.create();
const index = book.worksheets.add('G52基因索引');
const alpha = book.worksheets.add('Alpha最终证据表');
const audit = book.worksheets.add('逐基因证据链');
const sources = book.worksheets.add('来源与限制');
for (const sheet of [index, alpha, audit, sources]) sheet.showGridLines = false;

const indexCols = ['总序号','当前locus','WP','旧AFE','映射状态','当前product','染色体','起点','终点','链向','2016反应ID','2016反应数','2016对照状态','D3 KEGG KO','D3 KEGG EC','D3 KEGG reaction','D3 BioCyc reaction ID','D3 Rhea交叉引用','KEGG live KO','KEGG live EC','KEGG live reaction','KEGG live pathways','KEGG AFE URL','UniProt accession','UniProt reviewed','UniProt protein name','UniProt EC','InterPro ID','Pfam ID','UniProt active site','UniProt binding site','UniProt URL','UniProt catalytic activity','UniProt reaction Rhea','UniProt reaction evidence','2019蛋白组命中','2019同WP命中','2019序列追溯','2019补表原文','2019来源DOI','精确reaction确认','2026动作','审查状态','关联行ID','未决点','映射依据','D3 B1状态','官方来源','WP来源','任务来源'];
index.getRangeByIndexes(0,0,1,indexCols.length).values = [indexCols];
index.getRangeByIndexes(1,0,data.records.length,indexCols.length).values = data.records.map(r=>indexCols.map(c=>{
  const d=decisionsByLocus.get(r['当前locus']);
  if(c==='2026动作') return d?.['2026建议']||r[c]||'';
  if(c==='精确reaction确认' && r['当前locus']==='RU820_RS01820') return '同株重组蛋白催化已确认；旧式子需替换';
  if(c==='审查状态') return 'G52证据重查完成；未证实项按当前建议Hold';
  if(c==='未决点') return d?.['证据问题']||r[c]||'';
  if(c==='关联行ID') return [r[c],...(candidateByLocus.get(r['当前locus'])||[])].filter(Boolean).join(';');
  return r[c] ?? '';
}));
index.getRange('A1:AX1').format = {fill:'#17324D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:32};
index.getRange('A1:AX1').format.verticalAlignment='center';
index.getRange('A2:AX294').format.font = {name:'Arial',size:10,color:'#1B2733'};
index.getRange('A2:A294').setNumberFormat('0');
index.getRange('H2:I294').setNumberFormat('0');
index.getRange('L2:L294').setNumberFormat('0');
index.getRange('A:A').format.columnWidth=10;
index.getRange('B:E').format.columnWidth=19;
index.getRange('F:F').format.columnWidth=48;
index.getRange('G:J').format.columnWidth=17;
index.getRange('K:M').format.columnWidth=28;
index.getRange('N:U').format.columnWidth=24;
index.getRange('V:AX').format.columnWidth=40;
index.freezePanes.freezeRows(1);

const core = ['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe²⁺ lb','Fe²⁺ ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub'];
const ext = ['Record class','Candidate ID','Candidate group','2016 linked Reaction ID','2026 reaction object','Legacy AFE locus','Current locus','Primary evidence DOI','Evidence type','Evidence boundary','2026 action','2026 reaction score','Source ID','Score scope'];
alpha.getRangeByIndexes(0,0,1,30).values=[core.concat(ext)];
const byLocus=new Map(data.records.map(r=>[r['当前locus'],r]));
const rows=data.reaction_rows.map(r=>{
  const gene=byLocus.get(r['Current locus']);
  const proteinDetected=gene?.['2019同WP命中']==='是';
  const b=balance.get(r['Reaction ID']);
  let verdict=reactionOverrides[r['Reaction ID']]||{};
  if(r['Reaction ID']==='PGM1' && r['Current locus']==='RU820_RS02555') verdict={action:'保留候选并核正EC',scope:'当前 WP 为 BPG-independent iPGM，EC 5.4.2.12；2016 表填 5.4.2.1'};
  if(r['Reaction ID']==='PGM1' && r['Current locus']==='RU820_RS02995') verdict={action:'Hold该基因OR分支',scope:'当前 WP 仅 histidine phosphatase family，无精确 PGAM 化学注释'};
  if(r['Reaction ID']==='G5SD') verdict={action:'保留候选',scope:'同号他株 CopD 文献不可转移；当前 WP 为 proA'};
  if(!verdict.action && b?.status!=='balanced') verdict={action:'Hold旧式子',scope:`2016 补表配平状态 ${b?.status||'未计算'}；先解决代谢物式子或电荷`};
  if(!verdict.action) verdict={action:'维持2016历史候选；Hold 2026提升或改式',scope:'当前 WP 与旧 AFE 已映射；跨库注释和化学比对仅支持候选，精确底物、方向及区室未获本 WP 独立实验证实'};
  const cmp=chemistryByPair.get(`${r['Current locus']}|${r['Reaction ID']}`);
  const cBest=cmp?.best?.[0];
  const chemistryScope=`2016 配平：${b?.status||'未计算'}；KEGG 最接近：${cBest?`${cBest.kegg_reaction}（化合物集合重叠 ${cBest.similarity}）`:'未链接'}${cmp?.model_missing_kegg_ids?.length?`；旧代谢物缺 KEGG ID：${cmp.model_missing_kegg_ids.join(',')}`:''}；Rhea：${gene?.['UniProt reaction Rhea']||'未链接'}；旧 [c] 区室及方向未获独立实验确认。`;
  const directTrdr=r['Reaction ID']==='TRDR' && r['Current locus']==='RU820_RS01820';
  const directAtps=r['Reaction ID']==='SFAT1' && r['Current locus']==='RU820_RS02615';
  const directGalu=r['Reaction ID']==='GALU' && r['Current locus']==='RU820_RS02155';
  const primarySources=directTrdr?'G52-S004;G52-S016;G52-S022':directAtps?'G52-S004;G52-S016;G52-S027;G52-S028':directGalu?'G52-S004;G52-S016;G52-S031':proteinDetected?'G52-S004;G52-S008;G52-S016':'G52-S004;G52-S016';
  const fullSources=[...new Set([...primarySources.split(';'),...(decisionsByLocus.get(r['Current locus'])?.['来源ID']||'').split(';')].filter(Boolean))].join(';');
  return [
  r['Reaction ID'],r['Reaction Name'],r['Reaction Formula'],r['Confidence Level'],r['EC Number'],r['PMID'],r['Subsystem'],r['Gene-Reaction Association'],r['Gene-Protein-Reaction Association'],r['Protein-Reaction-Association'],r['Fe2 lb'],r['Fe2 ub'],r['tetrathionate lb'],r['tetrathionate ub'],r['sulfur lb'],r['sulfur ub'],
  '2016 GPR 证据重查',r['Candidate ID'],'G52',r['Reaction ID'],directTrdr?'h[c] + nadph[c] + trdox[c] -> nadp[c] + trdrd[c]':directAtps?'aps[c] + ppi[c] -> atp[c] + so4[c]':directGalu?'utp[c] + g1p[c] -> udpg[c] + ppi[c]':'',r['Legacy AFE locus'],r['Current locus'],directTrdr?'10.1007/s00284-009-9390-2':directAtps?'10.6026/97320630008695':directGalu?'10.1128/AEM.71.6.2902-2909.2005':proteinDetected?'10.3389/fmicb.2019.00592':'',directTrdr?'同株重组蛋白底物测定＋现行CDS双端引物核验':directAtps?'同株历史等位序列重组表达物实测ATP形成＋序列版本核验':directGalu?'同株 galU 缺失互补＋细胞提取物酶活＋提交 CDS 与现行完全相同':proteinDetected?'2016模型＋同株蛋白组检出':'2016模型＋历史映射',`${directGalu?'AY789511 与当前 CDS 完全一致；同株缺失互补及提取物活性证实核心催化；':verdict.scope+'；'}${chemistryScope}`,directGalu?'保留 GALU GPR；核心反应证据升级':verdict.action,directTrdr||directAtps||directGalu?4:'',fullSources,directTrdr?'4 仅对 NADPH 还原氧化型 thioredoxin 的精确反应；不转给旧错误式子、区室或参数':directAtps?'4 仅对 APS+PPi 生成 ATP+sulfate 的核心反应；历史等位序列与现行 WP 差 4 aa；质子、区室及网络方向待核':directGalu?'4 仅对 UTP+glucose-1-phosphate 生成 UDP-glucose+PPi 的核心反应；质子系数和区室由模型核对':'仅保留2016原分数；2026不评分'
  ];
});
for(const c of candidateDraft){
  const core=[...c.core_2016];
  const p162=c.locus==='RU820_RS02245';
  if(!c.old_reaction){
    core[0]=c.candidate_id;
    core[1]=byLocus.get(c.locus)?.['当前product']||'';
  }
  rows.push([...core,c.class,c.candidate_id,'G52',c.old_reaction,c.chemistry,c.afe,c.locus,p162?'10.1089/omi.2005.9.13':'',p162?'同株 P16.2 重组蛋白体外 TST 实验＋AY863108 与当前 CDS 序列核验':c.old_reaction?'2016空GPR＋当前数据库反应':'当前数据库与官方蛋白注释',c.scope,c.action,p162?4:'',c.sources,p162?c.scope:'新候选，不继承2016分数']);
}
for(const d of geneDecisions){
  const loc=d['当前locus'];
  const directDoi=loc==='RU820_RS01820'?'10.1007/s00284-009-9390-2':loc==='RU820_RS02615'?'10.6026/97320630008695':loc==='RU820_RS02245'?'10.1089/omi.2005.9.13':loc==='RU820_RS02155'?'10.1128/AEM.71.6.2902-2909.2005':'';
  rows.push([...Array(16).fill(''),'逐基因证据链',`G52-GENE-${loc}`,'G52',d['2016反应ID'],'',d['旧AFE'],loc,directDoi,d['决策类别'],d['证据问题'],d['2026建议'],'',d['来源ID'],'基因层结论；反应评分只见具体反应行，Hold 项不赋分']);
}
alpha.getRangeByIndexes(1,0,rows.length,30).values=rows;
alpha.getRange('A1:AD1').format={fill:'#17324D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:32};
alpha.getRange('A1:AD1').format.verticalAlignment='center';
alpha.getRange(`A2:AD${rows.length+1}`).format.font={name:'Arial',size:10,color:'#1B2733'};
alpha.getRange('A:B').format.columnWidth=26;
alpha.getRange('C:C').format.columnWidth=70;
alpha.getRange('D:P').format.columnWidth=20;
alpha.getRange('Q:AD').format.columnWidth=27;
alpha.freezePanes.freezeRows(1);

const auditCols=['总序号','当前locus','当前WP','旧AFE','当前product','2016反应ID','2016旧式子配平状态','2016配平差或缺少式子','2016/KEGG化合物ID重叠','KEGG KO','KEGG候选反应ID','KEGG候选反应定义','UniProt Rhea','UniProt催化反应','2019同WP蛋白检出','Europe PMC旧AFE索引命中数','旧AFE索引 DOI','现行WP索引命中数','现行WP索引 DOI','文献索引说明','逐基因决策类别','需重点核查的问题','当前建议','化学边界','证据来源','完成状态','TCDB全库近邻（未判定）','TCDB同家族比对及边界','功能词索引命中数','功能词索引 DOI（未核验）','Rhea官方反应式','Rhea ChEBI ID','Rhea/KEGG交叉引用判读','BRENDA同物种EC文献线索','BRENDA适用边界'];
audit.getRangeByIndexes(0,0,1,auditCols.length).values=[auditCols];
const auditRows=data.records.map(r=>{
  const ids=(r['KEGG live reaction']||'').split(';').filter(Boolean);
  const defs=ids.map(id=>`${id}: ${keggEntries[id]?.DEFINITION||'KEGG原始式未读取'}`).join(' | ');
  const m=manual[r['当前locus']]||{};
  const d=decisionsByLocus.get(r['当前locus']);
  const lit=literature[r['旧AFE']]||{};
  const dois=[...new Set((lit.hits||[]).map(h=>h.doi).filter(Boolean))].join(';');
  const wplit=wpLiterature[r['WP']]||{};
  const wpDois=[...new Set((wplit.hits||[]).map(h=>h.doi).filter(Boolean))].join(';');
  const flit=functionLiterature[r['当前locus']]||{};
  const fDois=[...new Set((flit.hits||[]).map(h=>h.doi).filter(Boolean))].join(';');
  const br=(r['2016反应ID']||'').split(';').filter(Boolean).map(id=>balance.get(id)).filter(Boolean);
  const bs=br.map(x=>`${x.reaction}:${x.status}`).join(';');
  const bd=br.filter(x=>x.status!=='balanced').map(x=>`${x.reaction}:差${JSON.stringify(x.difference_products_minus_reactants)};缺${(x.missing||[]).join(',')}`).join(' | ');
  const cm=(r['2016反应ID']||'').split(';').filter(Boolean).map(id=>{
    const x=chemistryByPair.get(`${r['当前locus']}|${id}`);
    return x?`${id}: ${x.best?.[0]?.kegg_reaction||'无'} ${x.best?.[0]?.similarity??'NA'}${x.model_missing_kegg_ids.length?'（旧代谢物无KEGG ID）':''}`:'';
  }).filter(Boolean).join('; ');
  const issue=d?.['证据问题']||'';
  const action=d?.['2026建议']||'';
  const t=tcdbAll[r['当前locus']]?.hits?.[0];
  const f=tcdbFamily[r['当前locus']];
  const tcAll=t?`${t.tcid}；identity ${t.identity}；query coverage ${t.query_coverage}；全库最近邻未经家族核验`:'未对该 WP 进行 TCDB 序列比对';
  const tcFamily=f?`${f.expected_family}；${f.best?.tcid||'无'}；identity ${f.best?.identity??''}；query coverage ${f.best?.query_coverage??''}；${f.family_homology_supported?'支持家族同源':'家族身份未获序列比对支持'}；底物、方向、细胞区室未定`:'未指定可核验家族';
  const rc=rheaKeggByLocus.get(r['当前locus']);
  const rheaDef=(rc?.rhea_ids||[]).map(id=>`${id}: ${rheaEntries[id]?.Equation||''}`).join(' | ');
  const rheaChebi=(rc?.rhea_ids||[]).map(id=>`${id}: ${rheaEntries[id]?.['ChEBI identifier']||''}`).join(' | ');
  const compare=rc?.status==='xref-overlap'?`共享 ${rc.overlap.join(',')}；同一交叉引用链，非两项独立实验证据`:rc?.status==='different-granularity-or-chemistry-to-review'?'无同号交叉引用；可能为总体反应/分步反应或通用/特定底物差异，逐式复核':'无双边链接';
  const ecList=[r['UniProt EC'],r['D3 KEGG EC'],r['KEGG live EC']].flatMap(x=>String(x||'').split(/[;, ]+/)).filter(x=>/^\d+\.\d+\.\d+\.\d+$/.test(x));
  const brendaHits=[...new Set(ecList)].filter(ec=>brenda[ec]?.count).map(ec=>`${ec}: ${(brenda[ec].matching_reference_rows||[]).join(' | ')}`).join(' || ');
  const brendaBoundary=brendaHits?(r['当前locus']==='RU820_RS01775'?'EC 1.8.1.7 纯化研究为 AP19-3，菌株不同；不可提升本 WP 分数':r['当前locus']==='RU820_RS02615'?'2012 原文和 FM177944 已独立核验；见 G52-S028':'仅物种/EC 文献线索，原文基因身份、菌株及实验类型待核'):'该 EC 未找到同物种参考记录或 EC 未指定';
  return [r['总序号'],r['当前locus'],r['WP'],r['旧AFE'],r['当前product'],r['2016反应ID'],bs,bd,cm,r['KEGG live KO'],r['KEGG live reaction'],defs,r['UniProt reaction Rhea'],r['UniProt catalytic activity'],r['2019同WP命中'],lit.hitCount??'',dois,wplit.hitCount??'',wpDois,(lit.hitCount||wplit.hitCount||flit.hit_count)?'索引线索已按当前证据边界处理；不能据此转移实验分数':'索引未命中；不推断无研究',d?.['决策类别']||'',issue,action,d?.['化学边界']||'',[d?.['来源ID']||'','G52-S023',t?'G52-S021':'',rheaDef?'G52-S025;G52-S026':'',brendaHits?'G52-S027':''].filter(Boolean).join(';'), 'G52证据重查完成；未证实项Hold',tcAll,tcFamily,flit.hit_count??'',fDois,rheaDef,rheaChebi,compare,brendaHits,brendaBoundary];
});
audit.getRangeByIndexes(1,0,auditRows.length,auditCols.length).values=auditRows;
audit.getRange('A1:AI1').format={fill:'#17324D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:32};
audit.getRange('A2:AI294').format.font={name:'Arial',size:10,color:'#1B2733'};
audit.getRange('A:A').format.columnWidth=10;
audit.getRange('B:D').format.columnWidth=20;
audit.getRange('E:H').format.columnWidth=28;
audit.getRange('I:AI').format.columnWidth=48;
audit.freezePanes.freezeRows(1);

const sourceCols=['Source ID','Full citation','Year','DOI','Source type','Strain','Experimental type','Genes / proteins studied','Related Reaction IDs / candidate IDs','What was directly measured','Applicable confidence level','Evidence limitation','Original URL'];
const sourceRows=[
 ['G52-S001','NCBI RefSeq assembly GCF_049532655.1, genomic.gff and CDS',2025,'','official genome annotation','ATCC 23270','sequence annotation','G52 293 RU820/WP','G52 index','Genome coordinates and protein IDs','identity only','Product and inferred function are computational annotation','https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/'],
 ['G52-S002','B1 旧新基因映射.tsv, exact WP and qualified mappings',2026,'','project verified mapping','ATCC 23270','sequence/neighbor comparison','mapped AFE/RU820','G52 index','Historical identity mapping','identity only','Ambiguous mappings remain unresolved',''],
 ['G52-S003','D3-1-R 全基因跨数据库交叉引用.tsv',2026,'','project database compilation','ATCC 23270 / older PGDB','secondary database cross references','G52 293 RU820/WP','G52 index','No direct experiment in compilation','2 at most, not assigned','KEGG/BioCyc/Rhea links may share annotation chain; reaction chemistry unverified',''],
 ['G52-S004','Campodonico et al., Acidithiobacillus ferrooxidans comprehensive model driven analysis of the electron transfer metabolism and synthetic strain design for biomining applications',2016,'10.1016/j.meteno.2016.03.003','2016 iMC507 original paper and Supplementary Table 1','ATCC 23270','model reconstruction','2016 AFE GPR','Alpha reaction precheck','Published reaction formulations, confidence and bounds','historical only','Model entries are not gene-specific wet-lab validation','https://doi.org/10.1016/j.meteno.2016.03.003'],
 ['G52-S005','G5 全新基因身份核验及辅助原始AFE核验',2026,'','project mapping review','ATCC 23270','sequence/neighbor comparison','G52 ambiguous B1 mappings','G52 index','Historical locus correspondence','identity only','Some mappings remain unresolved',''],
 ['G52-S006','UniProtKB taxonomy 243159 query, RefSeq WP exact cross references; InterPro/Pfam IDs reproduced from UniProtKB',2026,'','database accession and feature annotation','ATCC 23270 / DSM 14882','database record; catalytic evidence codes extracted','239 matched G52 WP','G52 index','No gene-specific wet experiment in UniProt citations','not scored','54 current WP unmatched; catalytic activity uses ECO:0000255/0000256 automated predictions; reviewed status is not wet-lab proof','https://rest.uniprot.org/uniprotkb/search?query=organism_id%3A243159'],
 ['G52-S007','KEGG REST afr ATCC 23270 KO, EC and pathway links; KO-to-reaction links',2026,'','database cross references','ATCC 23270, old AFE gene catalog','database assignment','279 confirmed AFE mapped G52 loci','G52 index','No direct experiment','not scored','KO-to-reaction link indicates annotation, not exact current-WP catalysis or model-compatible chemistry','https://rest.kegg.jp/info/afr']
 ,['G52-S008','Bellenberg et al. 2019. Proteomics Reveal Enhanced Oxidative Stress Responses and Metabolic Adaptation in Acidithiobacillus ferrooxidans Biofilm Cells on Pyrite',2019,'10.3389/fmicb.2019.00592','original proteomics study','DSM 14882T = ATCC 23270','pyrite biofilm versus Fe(II) proteomics','32 G52 AFE loci in Supplementary Tables 1–5','G52 index and linked model reactions','Protein IDs and relative abundance log2 pyrite/Fe with q-values','expression only','Protein detection and abundance do not establish catalytic substrate, compartment or direction; 2 historic WP differ from current WP','https://doi.org/10.3389/fmicb.2019.00592'],
 ['G52-S009','Primary studies and original Form IV Rubisco-like literature: A. ferrooxidans Form IV assigned to AFE_0434 in uranium-response study; Form IV superfamily studies describe loss of canonical RuBP carboxylation',2016,'','original study plus enzyme-family comparison','ATCC 23270 for locus expression; other strains for catalysis','transcription and comparative biochemistry','AFE_0434 / RU820_RS02105','RUBISCO;DKMPE','Same-strain expression and Form IV assignment; no same-WP catalytic assay','identity/expression only','The 2016 uranium paper does not establish DKMPE activity or direct loss of carboxylation in this WP','https://discovery.ucl.ac.uk/1474235/1/1-s2.0-S0923250816000127-main.pdf'],
 ['G52-S010','Alpha final evidence workbook: ATCC23270-2026_最终模型证据表.xlsx, historical candidate rows',2026,'','project prior evidence table','ATCC 23270','historical curation','RU820_RS02225 and RU820_RS02560','B-020;MTRP','Earlier Hold/Keep decisions and model comparison','historical only','Must re-evaluate against current RU820/WP and original evidence',''],
 ['G52-S011','Valdes et al. 2008, Acidithiobacillus ferrooxidans metabolism: from genome sequence to industrial applications; glgP2 AFE0527 context',2008,'10.1186/1471-2164-9-597','original genome metabolism analysis','ATCC 23270','genome based pathway inference','AFE_0527','MTRP','glgP2 is identified as glucan phosphorylase candidate','annotation/pathway support','No purified AFE_0527 catalysis or maltotetraose specific assay','https://doi.org/10.1186/1471-2164-9-597'],
 ['G52-S012','Chi et al. 2007, Periplasmic Proteins of the Extremophile Acidithiobacillus ferrooxidans: A High Throughput Proteomics Analysis',2007,'10.1074/mcp.M700042-MCP200','original periplasmic proteomics','ATCC 23270','periplasmic fraction LC-MS','paper labels CycA-2 as AFE_0378','CYTAA32 identity warning','Periplasmic cytochromes detected in older numbering','identity conflict only','Paper AFE_0378 text conflicts with current exact AFE_0378 to WP_009566069.1 FadL; protein sequence/accession needs reconciliation before transferring evidence','https://doi.org/10.1074/mcp.M700042-MCP200'],
 ['G52-S013','Europe PMC REST exact historical AFE text search, 279 G52 confirmed AFE labels',2026,'','literature index','multiple strains and article types','title/fulltext index search','279 G52 AFE labels','G52 gene index','35 AFE labels returned at least one hit','discovery only','Index hit is neither exact WP identity nor functional/experimental proof; absent hit does not prove absent literature','https://www.ebi.ac.uk/europepmc/webservices/rest/search'],
 ['G52-S014','Esparza et al. 2019, Effect of CO2 Concentration on Uptake and Assimilation of Inorganic Carbon in the Extreme Acidophile Acidithiobacillus ferrooxidans',2019,'10.3389/fmicb.2019.00603','original expression study','ATCC 23270','RT-qPCR and selected Western blots','AFE_0536/cbbP in G52','PRUK','cbbP RNA and CbbP protein abundance across CO2 conditions','expression only','No purified current-WP phosphoribulokinase reaction measurement','https://doi.org/10.3389/fmicb.2019.00603'],
 ['G52-S015','Europe PMC REST exact current WP protein identifier text search, 293 G52 proteins',2026,'','literature index','multiple strains and article types','title/fulltext index search','293 current WP identifiers','G52 gene index','0 exact WP text hits','discovery only','WP can be shared across strains; index absence does not prove no publication','https://www.ebi.ac.uk/europepmc/webservices/rest/search'],
 ['G52-S016','2016 iMC507 Supplementary Table 1 reaction formulas and Table 2 charged metabolite formulas; G52 reaction mass/charge audit',2026,'','derived balance check','ATCC 23270 model','atom and charge accounting','87 unique 2016 reactions linked to G52 AFE GPR','G52 2016 linked reactions','78 balanced; 4 charge-unbalanced; 5 unknown because metabolite formula or charge missing','chemistry check only','Uses published charged formulas; imbalance may reflect metabolite annotation error and must be resolved before 2026 reaction editing','https://doi.org/10.1016/j.meteno.2016.03.003'],
 ['G52-S017','Wu et al. 2010, Differential gene expression in response to copper in Acidithiobacillus ferrooxidans strains possessing dissimilar copper resistance',2010,'10.2323/jgam.56.491','original expression study, other strains','26# and DC mine isolates','RT-qPCR under copper stress','paper labels a CopD-like gene afe_0454','G5SD identity warning only','Copper-stress expression of another-strain gene called afe_0454','none for current WP','Cannot transfer other-strain early-locus CopD annotation to current ATCC 23270 RU820_RS02200/WP','https://doi.org/10.2323/jgam.56.491'],
 ['G52-S018','Osorio et al. 2013, Anaerobic sulfur metabolism coupled to dissimilatory iron reduction in the extremophile Acidithiobacillus ferrooxidans',2013,'10.1128/AEM.03057-12','original transcript/protein study','ATCC 23270','aerobic sulfur versus anaerobic Fe(III) transcriptomics/proteomics','AFE_0423 and AFE_0424 among G52 loci','ACONT1;ACONT2;ICDHyr;ICITRED','Relative transcript and protein changes','expression only','No purified enzyme chemistry; AFE_0424 historical protein accession must be reconciled with current WP','https://doi.org/10.1128/AEM.03057-12'],
 ['G52-S019','IUBMB enzyme nomenclature EC 3.13.2.1, adenosylhomocysteinase; transferred from EC 3.3.1.1',2026,'','official enzyme nomenclature','general','EC identifier update','RU820_RS02590 / AFE_0534','ADNSE','Current EC number and historical transfer','nomenclature only','EC transfer does not prove this WP enzymatic activity','https://iubmb.qmul.ac.uk/enzyme/EC3/13/2/1.html'],
 ['G52-S020','2016 model metabolite KEGG IDs compared with KEGG REST G52 KO-linked reaction compound IDs',2026,'','derived cross-database chemistry check','ATCC 23270 annotation','compound ID set comparison','90 historical G52 reaction-gene rows','G52 2016 linked reactions','40 rows have >=0.95 two-side set similarity with no missing model KEGG IDs','cross-reference only','Low similarity can result from generic KEGG compounds, proton conventions, missing model KEGG IDs or true substrate mismatch; not a catalyst assay','https://rest.kegg.jp/info/reaction'],
 ['G52-S021','TCDB public sequence set, retrieved 2026-09-25; local current-WP alignment with family-constrained check',2026,'','transporter sequence family database','multiple organisms','local amino-acid sequence comparison','25 G52 transport-related WPs; 10 specified family checks','G52 transporter audit','9 family-consistent sequence homology leads; SecE nearest family reference is SecY and is not accepted','family only','No direct substrate, transport direction, cellular compartment or same-WP transport assay follows from TCDB homology; all-database nearest hits can be unrelated families','https://www.tcdb.org/public/'],
 ['G52-S022','Wang et al. Expression, purification and molecular structure modeling of thioredoxin and thioredoxin reductase from Acidithiobacillus ferrooxidans',2009,'10.1007/s00284-009-9390-2','original biochemical study','ATCC 23270','recombinant enzyme NADPH/oxidized thioredoxin assay and Cys mutants','AFE_0375 / RU820_RS01820 / WP_009566072.1','TRDR','NADPH-dependent reduction of oxidized thioredoxin; both published gene-specific primer ends exactly match current 972-nt CDS','direct same-CDS catalysis','2016 reaction direction/charge remains wrong; candidate corrected equation balances but full-model network pending','https://doi.org/10.1007/s00284-009-9390-2'],
 ['G52-S023','Europe PMC REST phrase search of Acidithiobacillus ferrooxidans and each of the 293 current RefSeq product names',2026,'','literature discovery index','multiple strains and article types','fulltext phrase index search','all 293 G52 current product names','G52 gene audit','127 product phrase queries returned at least one indexed hit','discovery only','Generic products produce many unrelated hits; a positive hit does not identify the current WP or establish an experiment; zero hits does not prove no study','https://www.ebi.ac.uk/europepmc/webservices/rest/search'],
 ['G52-S024','Bustamante et al. Toxin-antitoxin systems in the mobile genome of Acidithiobacillus ferrooxidans',2014,'10.1371/journal.pone.0112226','original genome comparison and ICEAfe1 assay','ATCC 23270 and ATCC 53993','genome family classification; recombinant assays of other ICEAfe1 toxin-antitoxin pairs','AFE_0413, AFE_0414, AFE_0477, AFE_0478','G52 TA annotations','Table 2 lists four G52 loci as chromosomal putative TA pairs; direct functional assays were on ICEAfe1 genes outside G52','annotation only for G52 four','No gene-specific RNase, tRNA target, toxin-antitoxin pairing assay or metabolic reaction demonstrated for these four WPs','https://doi.org/10.1371/journal.pone.0112226'],
 ['G52-S025','Rhea official REST query, 144 G52 Rhea reaction identifiers; equation, ChEBI, EC and curated xrefs',2026,'','curated reaction database','general','reaction definition and curated reference lookup','G52 144 Rhea IDs','G52 gene audit','144 official reaction definitions recovered, none missing','chemistry definition only','Rhea reaction existence is not gene-specific catalysis or strain-specific direction/compartment evidence','https://www.rhea-db.org/help/rest-api'],
 ['G52-S026','G52 Rhea/KEGG cross-reference comparison from official Rhea reaction xrefs and KEGG REST',2026,'','derived annotation-chain check','G52 ATCC 23270 mapping','cross-reference comparison','74 G52 WPs with both Rhea and KEGG links','G52 gene audit','70 share an exact reaction cross-reference; 4 have differing granularity or chemistry requiring review','xref only','Shared Rhea/KEGG cross-reference counts as one inherited annotation chain; lack of same ID is not itself a contradictory experiment','https://www.rhea-db.org/help/rest-api'],
 ['G52-S027','BRENDA EC reference table query for 119 G52 numeric EC labels',2026,'','enzyme literature index','multiple A. ferrooxidans strains','EC-linked reference discovery','7 EC labels with species-name references','G52 gene audit','Seven EC labels returned species-reference rows; individual original studies require gene and strain matching','discovery only','Species/EC record cannot establish current WP catalysis; EC 1.8.1.7 purified glutathione reductase was from AP19-3','https://www.brenda-enzymes.org/'],
 ['G52-S028','Jaramillo et al. Cloning, expression and bioinformatics analysis of ATP sulfurylase from Acidithiobacillus ferrooxidans ATCC 23270 in Escherichia coli',2012,'10.6026/97320630008695','original recombinant expression and ATP formation study','ATCC 23270','APS + PPi to ATP bioluminescence assay of expressing E. coli extracts','historical FM177944.1 atpS; RU820_RS02615 / WP_012536192.1 near-identical current allele','SFAT1;ADSK;SFAT2','ATP formation from APS and PPi in induced extracts versus controls; published cloned allele differs 7 nt and 4 aa from current CDS','4 for SFAT1 core only','APS kinase domain was predicted but not assayed; SFAT2 unassayed; proton coefficient, compartment and current WP allele activity not directly measured','https://pmc.ncbi.nlm.nih.gov/articles/PMC3449377/'],
 ['G52-S029','Salazar et al. 2003 and Nunez et al. 2004, complementary A. ferrooxidans GluRS1/GluRS2 tRNA specificity',2003,'10.1073/pnas.1936123100;10.1016/S0014-5793(03)01460-1','original recombinant tRNA specificity studies','A. ferrooxidans, original study cloned genes before complete annotation','in vitro aminoacylation and E. coli heterologous assay','GluRS1 and GluRS2; historical AFE_0422 later labeled gltX-2 by KEGG','GLUTRS','GluRS2 preferentially charges tRNA(Gln)-UUG; GluRS1 charges tRNA(Glu)','not scored for current WP','No published clone sequence/primers yet aligned to current WP_012536128.1; do not transfer directly to its exact reaction score','https://pmc.ncbi.nlm.nih.gov/articles/PMC283512/'],
 ['G52-S030','Acosta et al. Identification of putative sulfurtransferase genes in the extremophilic Acidithiobacillus ferrooxidans ATCC 23270 genome: structural and functional characterization of the proteins',2005,'10.1089/omi.2005.9.13','original recombinant rhodanese activity study','ATCC 23270','P16.2 thiosulfate:cyanide sulfurtransferase assay','published P16.2 AY863108.1; near current RU820_RS02245 / WP_009564831.1','G52-NEW-RU820_RS02245','P16.2 recombinant enzyme shows in vitro thiosulfate:cyanide sulfurtransferase activity; AY863108 has one extra G versus current 450-nt CDS','4 for historical P16.2 assay only','Single-base insertion causes C-terminal reading-frame difference from current WP; cyanide is assay acceptor, physiological acceptor and compartment not established','https://doi.org/10.1089/omi.2005.9.13'],
 ['G52-S031','Barreto et al. Identification of a gene cluster for the formation of extracellular polysaccharide precursors in the chemolithoautotroph Acidithiobacillus ferrooxidans',2005,'10.1128/AEM.71.6.2902-2909.2005','original same-strain gene complementation and enzyme assays','ATCC 23270','E. coli galU deletion complementation and coupled UDP-glucose pyrophosphorylase assay','AY789511.1 exact current RU820_RS02155 / WP_009567323.1 CDS; AY789512.1 phosphoglucomutase maps outside G52','GALU;PGM1 identity warning','GalU complemented mutant and extract activity 45±7 versus 0.5±0.1 nmol/min/mg negative control; AY789511 exact 897-nt current CDS','4 for GALU core reaction','Coupled assay was in extracts and did not measure proton coefficient or compartment; PGM experiment belongs to a WP outside G52','https://pmc.ncbi.nlm.nih.gov/articles/PMC1151869/']
];
sources.getRangeByIndexes(0,0,1,sourceCols.length).values=[sourceCols];
sources.getRangeByIndexes(1,0,sourceRows.length,sourceCols.length).values=sourceRows;
sources.getRange('A1:M1').format={fill:'#17324D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:32};
sources.getRange('A1:M1').format.verticalAlignment='center';
sources.getRange(`A2:M${sourceRows.length+1}`).format.font={name:'Arial',size:10,color:'#1B2733'};
sources.getRange('A:A').format.columnWidth=14;
sources.getRange('B:B').format.columnWidth=75;
sources.getRange('C:D').format.columnWidth=24;
sources.getRange('E:M').format.columnWidth=32;
sources.freezePanes.freezeRows(1);

book.recalculate();
const preview=await book.render({sheetName:'G52基因索引',range:'A1:F8',scale:1.3,format:'png'});
await fs.writeFile(new URL('g52_precheck_preview.png',dir),new Uint8Array(await preview.arrayBuffer()));
const auditPreview=await book.render({sheetName:'逐基因证据链',range:'A1:H8',scale:1.1,format:'png'});
await fs.writeFile(new URL('g52_audit_preview.png',dir),new Uint8Array(await auditPreview.arrayBuffer()));
const output=await SpreadsheetFile.exportXlsx(book);
await output.save(fileURLToPath(new URL('G52_全量证据重查_Alpha最终格式.xlsx',dir)));
console.log(JSON.stringify({index:data.records.length,alpha:rows.length,audit:auditRows.length,manual:Object.keys(manual).length,sources:sourceRows.length}));
