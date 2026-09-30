import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const here = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(here, 'g57_data.json'), 'utf8'));
const outDir = path.resolve(here, '../../outputs/G57');
await fs.mkdir(outDir, {recursive:true});
const wb = Workbook.create();

function createTableSheet(name, rows, color, widthOverrides={}) {
  const sh = wb.worksheets.add(name);
  sh.showGridLines = false;
  sh.tabColor = color;
  const fields = Object.keys(rows[0] || {});
  const matrix = [fields, ...rows.map(r => fields.map(f => {
    const v = r[f];
    return v === undefined || v === null ? '' : v;
  }))];
  sh.getRangeByIndexes(0,0,matrix.length,fields.length).values = matrix;
  sh.getRangeByIndexes(0,0,matrix.length,fields.length).format.font = {name:'Arial',size:10,color:'#1F2937'};
  const head = sh.getRangeByIndexes(0,0,1,fields.length);
  head.format = {fill:'#17324D',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},rowHeight:28,verticalAlignment:'center'};
  sh.freezePanes.freezeRows(1);
  if (fields.length > 10) sh.freezePanes.freezeColumns(3);
  for (let j=0;j<fields.length;j++) {
    const label = fields[j];
    const w = widthOverrides[label] ?? (label.includes('DOI') || label.includes('URL') ? 38 : label.includes('product') || label.includes('证据') || label.includes('待验证') || label.includes('未决') ? 48 : 22);
    sh.getRangeByIndexes(0,j,matrix.length,1).format.columnWidth = w;
  }
  sh.getRangeByIndexes(1,0,Math.max(1,rows.length),fields.length).format.rowHeight = 20;
  return sh;
}

const summary = wb.worksheets.add('验收摘要');
summary.showGridLines = false;
summary.tabColor = '#12324B';
summary.getRange('A2').values = [['G57 基因证据重查 · 阶段结果']];
summary.getRange('A2').format.font = {name:'Arial',size:15,bold:true,color:'#17324D'};
summary.getRange('A4:B13').values = [
  ['冻结范围','1759—2051'],['当前基因数',293],['同版身份核对通过',293],
  ['UniProt 精确 WP 命中',213],['Europe PMC 编号检索命中',34],
  ['2016 反应逐项对照行',data.reaction_review.length],
  ['无需修改（暂定）',data.counts['无需修改（暂定）'] || 0],
  ['证据不足暂缓',data.counts['证据不足暂缓'] || 0],
  ['仅更新证据或编号',data.counts['仅更新证据或编号'] || 0],
  ['全量证据链验收','未通过'],
];
summary.getRange('A4:A13').format.font = {name:'Arial',size:10,bold:true,color:'#17324D'};
summary.getRange('B4:B13').format.font = {name:'Arial',size:10,color:'#1F2937'};
summary.getRange('A4:B13').format.rowHeight = 25;
summary.getRange('A4:A13').format.columnWidth = 28;
summary.getRange('B4:B13').format.columnWidth = 32;
summary.getRange('A15').values = [['已核实的重点']];
summary.getRange('A15').format.font = {name:'Arial',size:11,bold:true,color:'#17324D'};
summary.getRange('A16:B19').values = [
  ['RUBISCO','AFE_2155 是 CbbM；重组酶正向羧化有实验，2016 原分 3，建议正向反应分 4。'],
  ['CbbQ2/CbbO2','AFE_2156/2157 为 Rubisco 活化因子；不并入羧化催化 GPR。'],
  ['PPTtpp','phnG—phnM 有同株生长和转录支持；旧转运 GPR 与零上下界需独立复核。'],
  ['SQRED1','SreABCD 同株厌氧表达支持硫还原候选；旧硫化物:醌氧化关联需独立复核。'],
];
summary.getRange('A16:A19').format.font = {name:'Arial',size:10,bold:true,color:'#17324D'};
summary.getRange('B16:B19').format.font = {name:'Arial',size:10,color:'#1F2937'};
summary.getRange('B16:B19').format.columnWidth = 95;
summary.getRange('A16:B19').format.rowHeight = 31;
summary.getRange('A21').values = [['状态说明']];
summary.getRange('A21').format.font = {name:'Arial',size:11,bold:true,color:'#17324D'};
summary.getRange('A22').values = [['每个基因均有身份和模型索引记录。绝大多数基因的精确反应、原始论文与化学定位仍需逐项核实；暂定动作不是正式模型修改。']];
summary.getRange('A22').format.font = {name:'Arial',size:10,color:'#374151'};
summary.getRange('A22').format.columnWidth = 105;

createTableSheet('逐基因索引',data.index,'#29607D',{'总序号':12,'RU820':20,'WP':20,'当前 product':52,'未决点':95,'证据边界':75});
createTableSheet('九项证据链',data.chain,'#47748B',{'总序号':12,'RU820':20,'WP':20,'1 身份':56,'2 历史映射':55,'3 蛋白功能':95,'4 反应通路':85,'5 原始论文':65,'6 化学与定位':88,'7 2016 对照':76,'9 复核':65});
createTableSheet('Alpha逐基因',data.gene_review,'#47748B',{'Candidate ID':18,'Current locus':20,'Evidence boundary':76,'2026 reaction object':63});
createTableSheet('2016反应对照',data.reaction_review,'#47748B',{'Reaction ID':22,'Reaction Formula':90,'Gene-Reaction Association':80,'Evidence boundary':75,'Score scope':72});
createTableSheet('来源',data.sources,'#7A8590',{'Full citation':95,'What was directly measured':95,'Evidence limitation':92,'Original URL':68});
createTableSheet('文献逐篇判读',data.literature_triage,'#7A8590',{'DOI':38,'标题':90,'命中编号':90,'论文实验对象':48,'证据类型':48,'对 G57 的使用边界':110,'原文':70});
createTableSheet('论文引物身份',data.primer_audit,'#7A8590',{'论文表1引物':26,'引物序列5to3':35,'当前染色体精确命中起点':28,'匹配方向':18,'当前RU820':24,'证据边界':92});
createTableSheet('同EC旧模型候选',data.ec_model_candidates,'#7A8590',{'总序号':12,'RU820':20,'WP':20,'当前product':55,'基因侧EC':18,'2016同EC反应':24,'2016反应式':95,'2016原GPR':75,'证据边界':100});
createTableSheet('当前编号文献检索',data.current_id_search,'#7A8590',{'总序号':12,'RU820':20,'WP':20,'RU820命中':17,'WP命中':17,'RU820 DOI':44,'WP DOI':44,'RU820查询URL':85,'WP查询URL':85,'检索状态':42});
createTableSheet('TCDB家族核验',data.tcdb_review,'#7A8590',{'总序号':12,'RU820':20,'WP':20,'当前product':55,'查询TCDB大类':25,'最佳TCDB ID':22,'最佳TCDB描述':85,'同一性':15,'当前蛋白覆盖率':20,'家族同源支持':20,'证据边界':90});
const gaps = data.index.map(r => ({总序号:r['总序号'],RU820:r['RU820'],WP:r['WP'],旧AFE:r['旧 AFE'],'2016反应':r['2016 GPR 反应'],初步动作:r['动作'],待验证:r['未决点']}));
createTableSheet('缺证与待验证',gaps,'#7A8590',{'待验证':110,'RU820':20,'WP':20});

wb.recalculate();
for (const [name,range,filename] of [['验收摘要','A1:B23','摘要预览.png'],['逐基因索引','A1:H11','索引预览.png'],['2016反应对照','A1:H9','反应预览.png'],['来源','A1:F8','来源预览.png']]) {
  const blob = await wb.render({sheetName:name,range,scale:1.5,format:'png'});
  await fs.writeFile(path.join(outDir,filename),new Uint8Array(await blob.arrayBuffer()));
}
const check = await wb.inspect({kind:'table',range:'验收摘要!A4:B13',tableMaxRows:12,tableMaxCols:2,maxChars:2000});
await fs.writeFile(path.join(outDir,'验收摘要.inspect.ndjson'),check.ndjson,'utf8');
const output = await SpreadsheetFile.exportXlsx(wb);
const finalPath = path.join(outDir,'G57_阶段证据审查.xlsx');
await output.save(finalPath);
console.log(finalPath);
