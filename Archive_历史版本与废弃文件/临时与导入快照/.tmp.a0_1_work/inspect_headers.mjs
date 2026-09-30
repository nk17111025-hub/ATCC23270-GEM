import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";
const path = "D:/嗜酸氧化亚铁硫杆菌/01_原始资料/2024_Khaleque/from2024-mmc1.xlsx";
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(path));
for (const name of ["Table S3", "Table S4"]) {
  const sh = wb.worksheets.getItem(name);
  console.log(name, JSON.stringify(sh.getRange("A1:V2").values));
}
