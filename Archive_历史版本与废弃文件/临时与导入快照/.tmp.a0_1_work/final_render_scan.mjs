import fs from "node:fs/promises";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";
const path = "D:/嗜酸氧化亚铁硫杆菌/08_模型基线/A0-1/2024补充表_S3S4纠正版.xlsx";
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(path));
wb.recalculate();
for (const sheetName of ["Table S3", "Table S4"]) {
  const check = await wb.inspect({kind:"table", sheetId:sheetName, range:"A622:V633", include:"values,formulas", tableMaxRows:20, tableMaxCols:22, maxChars:22000});
  console.log(`===== ${sheetName} FINAL =====`);
  console.log(check.ndjson);
  const image = await wb.render({sheetName, range:"A622:V633", scale:1.2, format:"png"});
  await fs.writeFile(`D:/嗜酸氧化亚铁硫杆菌/.tmp.a0_1_work/${sheetName.replaceAll(" ","_")}_final.png`, new Uint8Array(await image.arrayBuffer()));
}
const errors = await wb.inspect({kind:"match", searchTerm:"#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options:{useRegex:true,maxResults:300}, summary:"final formula error scan"});
console.log("ERROR_SCAN", errors.ndjson);
