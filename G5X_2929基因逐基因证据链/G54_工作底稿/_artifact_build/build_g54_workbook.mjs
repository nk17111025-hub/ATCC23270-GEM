import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const here = path.dirname(fileURLToPath(import.meta.url));
const column = (n) => { let s = ''; for (; n; n = Math.floor((n - 1) / 26)) s = String.fromCharCode(65 + (n - 1) % 26) + s; return s; };
const data = JSON.parse(await fs.readFile(path.join(here, 'g54_workbook_data.json'), 'utf8'));
const output = process.argv[2];
if (!output) throw new Error('output path required');
const wb = Workbook.create();
const colors = ['#0F766E', '#1967A3', '#4263A9', '#5A57A6', '#6F42A1', '#697681', '#48707D', '#44775B', '#A2672F', '#8A573D', '#AA4F55'];
for (let i = 0; i < data.length; i++) {
  const { name, rows } = data[i];
  const sheet = wb.worksheets.add(name);
  const n = rows.length;
  const m = rows[0].length;
  sheet.showGridLines = false;
  sheet.tabColor = colors[i];
  sheet.getRangeByIndexes(0, 0, n, m).values = rows;
  sheet.getRangeByIndexes(0, 0, n, m).format.font = { name: 'Aptos', size: 10 };
  sheet.getRangeByIndexes(0, 0, 1, m).format = {
    fill: colors[i],
    font: { name: 'Aptos', size: 10, bold: true, color: '#FFFFFF' },
    wrapText: true,
    rowHeight: 42,
  };
  sheet.getRangeByIndexes(0, 0, n, m).format.columnWidth = 23;
  sheet.freezePanes.freezeRows(1);
  sheet.tables.add(`A1:${column(m)}${n}`, true, `G54Table${i + 1}`);
  console.log(name, n - 1, m);
}
await wb.recalculate();
for (const sheet of data) {
  const preview = await wb.render({ sheetName: sheet.name, range: 'A1:D3', scale: 1, format: 'png' });
  await fs.writeFile(path.join(here, `preview_${sheet.name}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const blob = await SpreadsheetFile.exportXlsx(wb);
await fs.mkdir(path.dirname(output), { recursive: true });
await blob.save(output);
console.log('saved', output);
