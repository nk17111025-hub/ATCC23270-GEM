import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const sourcePath = "D:/嗜酸氧化亚铁硫杆菌/01_原始资料/2024_Khaleque/from2024-mmc1.xlsx";
const copyPath = "D:/嗜酸氧化亚铁硫杆菌/08_模型基线/A0-1/2024补充表_S3S4纠正版.xlsx";

for (const [label, path] of [["SOURCE", sourcePath], ["COPY", copyPath]]) {
  const input = await FileBlob.load(path);
  const wb = await SpreadsheetFile.importXlsx(input);
  console.log(`===== ${label} =====`);
  const summary = await wb.inspect({
    kind: "workbook,sheet,table,definedName",
    maxChars: 20000,
    tableMaxRows: 8,
    tableMaxCols: 12,
    tableMaxCellChars: 120,
  });
  console.log(summary.ndjson);
  const sheets = await wb.inspect({kind: "sheet", include: "id,name"});
  console.log("SHEETS", sheets.ndjson);
}
