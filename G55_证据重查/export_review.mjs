import fs from 'node:fs/promises';
import { SpreadsheetFile, Workbook } from '@oai/artifact-tool';

const folder = 'D:/嗜酸氧化亚铁硫杆菌/G55_证据重查';
const data = JSON.parse(await fs.readFile(`${folder}/review_data.json`, 'utf8'));
const wb = Workbook.create();

function col(n) {
  let s = '';
  for (let x = n; x > 0; x = Math.floor((x - 1) / 26)) s = String.fromCharCode(65 + (x - 1) % 26) + s;
  return s;
}

function makeSheet(name, headers, records, color, widths = {}) {
  const sh = wb.worksheets.add(name);
  const rows = [headers, ...records.map(row => row.map(value => {
    if (value === null || value === undefined) return '';
    if (Array.isArray(value) || typeof value === 'object') return JSON.stringify(value);
    return value;
  }))];
  if (rows.some(r => r.length !== headers.length)) throw new Error(`${name}: column mismatch`);
  const last = col(headers.length);
  sh.getRange(`A1:${last}${rows.length}`).values = rows;
  sh.getRange(`A1:${last}${rows.length}`).format.font = { name: 'Aptos', size: 10 };
  sh.getRange(`A1:${last}1`).format = { fill: color, font: { name: 'Aptos', size: 10, bold: true, color: '#FFFFFF' } };
  sh.getRange(`A1:${last}1`).format.rowHeight = 34;
  sh.getRange(`A1:${last}1`).format.wrapText = true;
  sh.freezePanes.freezeRows(1);
  sh.showGridLines = false;
  sh.tabColor = color;
  for (let i = 1; i <= headers.length; i++) {
    sh.getRange(`${col(i)}:${col(i)}`).format.columnWidth = widths[i] ?? 19;
  }
  return sh;
}

const index = makeSheet('基因索引_293', data.gene_headers, data.gene_rows, '#17365D', {
  1: 10, 2: 20, 3: 22, 4: 29, 5: 9, 6: 43, 7: 21, 8: 20, 9: 25,
  10: 12, 11: 17, 12: 14, 13: 19, 14: 23, 15: 25, 16: 18, 17: 16,
  18: 25, 19: 25, 20: 25, 21: 25, 22: 15, 23: 31, 24: 20, 25: 70,
  26: 70, 27: 70, 28: 36, 29: 29, 30: 27, 31: 45, 32: 45, 33: 15,
  34: 45, 35: 28,
});
index.getRange('A2:A294').format.numberFormat = '0';
index.getRange('J2:J294').format.numberFormat = '0';
index.getRange('AG2:AG294').format.numberFormat = '0';
const actionColor = {
  '修改或补充 GPR': '#E2F0D9',
  '仅更新证据或编号': '#DDEBF7',
  '确认无需修改': '#E7E6E6',
  '证据不足暂缓': '#FFF2CC',
};
for (let i = 0; i < data.gene_rows.length; i++) {
  const color = actionColor[data.gene_rows[i][23]];
  if (color) index.getRange(`X${i + 2}`).format.fill = color;
}

makeSheet('Alpha审查记录', data.alpha_headers, data.alpha_rows, '#375623', {
  1: 18, 2: 41, 3: 77, 4: 17, 5: 22, 6: 32, 7: 34, 8: 85, 9: 72,
  10: 65, 11: 16, 12: 16, 13: 16, 14: 16, 15: 16, 16: 16,
  17: 20, 18: 25, 19: 15, 20: 25, 21: 25, 22: 21, 23: 22,
  24: 27, 25: 27, 26: 75, 27: 23, 28: 18, 29: 31, 30: 65,
});
makeSheet('证据来源', data.source_headers, data.source_rows, '#5B2C6F', {
  1: 15, 2: 115, 3: 11, 4: 34, 5: 33, 6: 32, 7: 47, 8: 56,
  9: 47, 10: 75, 11: 30, 12: 75,
});
makeSheet('冲突与动作', data.issue_headers, data.issue_rows, '#9E480E', {
  1: 15, 2: 29, 3: 53, 4: 22, 5: 103, 6: 92, 7: 45, 8: 66,
});
makeSheet('候选暂缓', data.candidate_headers, data.candidate_rows, '#7F6000', {
  1: 10, 2: 21, 3: 22, 4: 45, 5: 20, 6: 18, 7: 33, 8: 33,
  9: 33, 10: 80, 11: 110, 12: 85, 13: 20, 14: 90,
});
const counts = new Map();
for (const row of data.gene_rows) counts.set(row[23], (counts.get(row[23]) ?? 0) + 1);
const overview = [
  ['项目', '数值', '说明'],
  ['G55索引范围', '1173–1465', 'G55清单全部293条'],
  ['基因总数', data.gene_rows.length, '当前RU820 locus各一行'],
  ['Alpha审查行', data.alpha_rows.length, '293个基因行；72个2016基因-反应行；2个补充反应行'],
  ['2016已关联基因', data.gene_rows.filter(r => r[20]).length, '关联51个不同反应'],
  ['修改或补充 GPR', counts.get('修改或补充 GPR') ?? 0, '涉及GMHEPAT/GMHEPK、ARGDC、DB4PS'],
  ['仅更新证据或编号', counts.get('仅更新证据或编号') ?? 0, '包括AGAD/UACAA直接酶学证据'],
  ['确认无需修改', counts.get('确认无需修改') ?? 0, '依当前功能和iMC507范围判断'],
  ['证据不足暂缓', counts.get('证据不足暂缓') ?? 0, '不据自动注释写入新反应'],
  ['来源记录', data.source_rows.length, '与Alpha来源表相同的12列表头'],
  ['冲突清单', data.issue_rows.length, '待实验或跨组模型整合核查'],
  ['候选暂缓清单', data.candidate_rows.length, '未作正式reaction新增'],
  ['底本规则', 'Alpha A–P', '2016原行A–P逐格冻结；本组判断写在Q–AD'],
];
const check = makeSheet('验收摘要', overview[0], overview.slice(1), '#1F4E78', {1:27, 2:32, 3:95});
check.getRange('A1:C13').format.rowHeight = 24;

const file = await SpreadsheetFile.exportXlsx(wb);
await file.save(`${folder}/G55_逐基因证据重查.xlsx`);
console.log(`${folder}/G55_逐基因证据重查.xlsx`);
