import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const copyPath = "D:/嗜酸氧化亚铁硫杆菌/08_模型基线/A0-1/2024补充表_S3S4纠正版.xlsx";
const previewPath = "D:/嗜酸氧化亚铁硫杆菌/.tmp.a0_1_work/S3S4_corrected_preview.png";

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(copyPath));
const s1 = workbook.worksheets.getItem("Table S1");
const s1Rows = s1.getRange("A2:C632").values;
const definitions = new Map(s1Rows.map((row) => [row[0], row[2]]));

const correctedIds = [
  "ZZ_glc_transport",
  "ZZ_glc_EX",
  "ZZ_glycogen_EX",
  "ZZ_g1p_transport",
  "ZZ_g1p_EX",
  "ZZ_Tre_amy",
  "ZZ_Tre_cga",
  "ZZ_Tre_EX",
  "ZZ_Tre_transport",
  "ZZ_TreT1_ADP",
  "ZZ_TreT2_UDP",
  "ZZ_TreYZ",
];
const correctedRows = correctedIds.map((id) => {
  if (!definitions.has(id)) throw new Error(`S1 中缺少 reaction definition: ${id}`);
  return [id, definitions.get(id)];
});

for (const sheetName of ["Table S3", "Table S4"]) {
  const sheet = workbook.worksheets.getItem(sheetName);
  sheet.getRange("A622:B633").values = correctedRows;
}

workbook.recalculate();
for (const sheetName of ["Table S3", "Table S4"]) {
  const check = await workbook.inspect({
    kind: "table",
    sheetId: sheetName,
    range: "A620:V633",
    include: "values,formulas",
    tableMaxRows: 20,
    tableMaxCols: 22,
    maxChars: 24000,
  });
  console.log(`===== ${sheetName} AFTER FIX =====`);
  console.log(check.ndjson);
}

const preview = await workbook.render({sheetName: "Table S3", range: "A618:V633", scale: 1.2, format: "png"});
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(copyPath);
console.log(`EXPORTED ${copyPath}`);
