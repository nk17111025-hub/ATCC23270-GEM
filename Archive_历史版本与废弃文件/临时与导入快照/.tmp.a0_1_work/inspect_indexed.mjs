import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";
const path = "D:/嗜酸氧化亚铁硫杆菌/01_原始资料/2024_Khaleque/from2024-mmc1.xlsx";
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(path));
for (const name of ["Table S1", "Table S3", "Table S4"]) {
  const sh = wb.worksheets.getItem(name);
  const range = name === "Table S1" ? "A616:O632" : "A616:V633";
  const rows = sh.getRange(range).values;
  const start = 616;
  console.log(`===== ${name} =====`);
  for (let i = 0; i < rows.length; i++) {
    console.log(JSON.stringify({row:start+i, values:rows[i]}));
  }
}
