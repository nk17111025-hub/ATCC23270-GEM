import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const folder = path.dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await fs.readFile(path.join(folder, "G58_中间数据.json"), "utf8"));
const wb = Workbook.create();
const tabs = [
  ["逐基因索引", data.index_header, data.index],
  ["Alpha审查行", data.review_header, data.review],
  ["来源", data.sources_header, data.sources],
  ["D3交叉引用", data.d3_header, data.d3],
  ["UniProt反应模板", data.rhea_header, data.rhea_templates],
  ["TCDB家族线索", data.tcdb_header, data.tcdb],
  ["BioCyc证据链", data.biocyc_header, data.biocyc],
  ["BRENDA酶类线索", data.brenda_header, data.brenda],
];
for (const [name, header, rows] of tabs) {
  const sheet = wb.worksheets.add(name);
  const matrix = [header, ...rows].map(row => row.map(v => v === undefined ? null : v));
  sheet.getRange("A1").write(matrix);
  const last = String.fromCharCode(64 + Math.min(header.length, 26));
  const endCol = header.length <= 26 ? last : "A" + String.fromCharCode(64 + header.length - 26);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill: "#17365D",
    font: { name: "Microsoft YaHei", bold: true, color: "#FFFFFF", size: 10 },
    rowHeight: 32,
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${matrix.length}`).format.font = {
    name: "Microsoft YaHei", size: 10,
  };
  sheet.getRange(`A1:${endCol}1`).format.font = {
    name: "Microsoft YaHei", bold: true, color: "#FFFFFF", size: 10,
  };
  sheet.getRange(`A:${endCol}`).format.columnWidth = 20;
  sheet.freezePanes.freezeRows(1);
}
wb.worksheets.getItem("逐基因索引").getRange("E:E").format.columnWidth = 38;
wb.worksheets.getItem("逐基因索引").getRange("G:G").format.columnWidth = 60;
wb.worksheets.getItem("逐基因索引").getRange("U:W").format.columnWidth = 48;
wb.worksheets.getItem("Alpha审查行").getRange("C:C").format.columnWidth = 64;
wb.worksheets.getItem("Alpha审查行").getRange("H:J").format.columnWidth = 45;
wb.worksheets.getItem("Alpha审查行").getRange("Z:Z").format.columnWidth = 65;
wb.worksheets.getItem("来源").getRange("B:B").format.columnWidth = 75;
wb.worksheets.getItem("来源").getRange("L:M").format.columnWidth = 65;
await wb.recalculate();
const preview = await wb.inspect({ kind:"sheet", include:"id,name", maxChars:1000 });
console.log(preview.ndjson || preview);
const blob = await SpreadsheetFile.exportXlsx(wb);
await blob.save(path.join(folder, "G58_证据审查_阶段版.xlsx"));
console.log("saved", data.counts);
