import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const here = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(here, 'G56_审查数据.json'), 'utf8'));
const wb = Workbook.create();

function letter(n) {
  let out = '';
  for (let k = n + 1; k > 0; k = Math.floor((k - 1) / 26)) out = String.fromCharCode(65 + (k - 1) % 26) + out;
  return out;
}

function addSheet(name, headers, body, headerRow = 1) {
  const sh = wb.worksheets.add(name);
  sh.showGridLines = false;
  const n = headers.length;
  if (name === '2026-extended from 2016') {
    const top = Array(n).fill(null);
    top[10] = 'Minimal Media for fe2';
    top[12] = 'Minimal Media for ttton';
    top[14] = 'Minimal Media for tsul';
    top[16] = '2026 evidence extension';
    sh.getRangeByIndexes(0, 0, 1, n).values = [top];
  }
  sh.getRangeByIndexes(headerRow - 1, 0, 1, n).values = [headers];
  const chunk = 80;
  for (let i = 0; i < body.length; i += chunk) {
    const part = body.slice(i, i + chunk).map(r => r.map(v => v === undefined ? null : v));
    sh.getRangeByIndexes(headerRow + i, 0, part.length, n).values = part;
  }
  const end = letter(n - 1);
  const h = sh.getRange(`A${headerRow}:${end}${headerRow}`);
  h.format.fill = '#25364D';
  h.format.font = { name: 'Aptos', size: 10, bold: true, color: '#FFFFFF' };
  h.format.rowHeight = 32;
  h.format.wrapText = true;
  h.format.verticalAlignment = 'center';
  const full = sh.getRange(`A${headerRow + 1}:${end}${headerRow + body.length}`);
  full.format.font = { name: 'Aptos', size: 10, color: '#1E293B' };
  full.format.rowHeight = 21;
  full.format.verticalAlignment = 'center';
  sh.getRange(`A:${end}`).format.columnWidth = 17;
  sh.getRange('A:A').format.columnWidth = 16;
  sh.getRange('B:B').format.columnWidth = 28;
  sh.getRange('C:C').format.columnWidth = 34;
  sh.freezePanes.freezeRows(headerRow);
  return sh;
}

const alpha = addSheet('2026-extended from 2016', data.alpha_headers, data.review, 2);
alpha.getRange('Q:Q').format.columnWidth = 24;
alpha.getRange('R:R').format.columnWidth = 28;
alpha.getRange('Z:Z').format.columnWidth = 56;
alpha.getRange('AA:AA').format.columnWidth = 20;
alpha.getRange('AD:AD').format.columnWidth = 45;
alpha.freezePanes.freezeColumns(3);

const sourceSheet = addSheet('2026-61 papers sources', data.source_headers, data.sources);
sourceSheet.getRange('B:B').format.columnWidth = 60;
sourceSheet.getRange('H:J').format.columnWidth = 38;
sourceSheet.getRange('L:L').format.columnWidth = 52;

const indexHeaders = Object.keys(data.index[0]);
const indexBody = data.index.map(obj => indexHeaders.map(k => obj[k] ?? null));
const indexSheet = addSheet('G56 gene index', indexHeaders, indexBody);
indexSheet.getRange('E:E').format.columnWidth = 48;
indexSheet.getRange('AD:AD').format.columnWidth = 54;
indexSheet.freezePanes.freezeColumns(3);

wb.recalculate();
for (const [name, range, key] of [
  ['2026-extended from 2016', 'Q2:AD6', 'alpha'],
  ['2026-61 papers sources', 'A1:L5', 'sources'],
  ['G56 gene index', 'A1:L5', 'index'],
]) {
  const info = await wb.inspect({ kind: 'table', range: `'${name}'!${range}`, include: 'values', tableMaxRows: 6, tableMaxCols: 14, maxChars: 5000 });
  await fs.writeFile(path.join(here, `G56_${key}_inspect.ndjson`), info.ndjson ?? '', 'utf8');
  const blob = await wb.render({ sheetName: name, range, scale: 1.4, format: 'png' });
  await fs.writeFile(path.join(here, `G56_${key}_preview.png`), new Uint8Array(await blob.arrayBuffer()));
}

const output = await SpreadsheetFile.exportXlsx(wb);
await output.save(path.join(here, 'G56_Alpha格式证据表.xlsx'));
console.log(JSON.stringify({ reviewRows: data.review.length, geneRows: data.index.length, sourceRows: data.sources.length, output: path.join(here, 'G56_Alpha格式证据表.xlsx') }));
