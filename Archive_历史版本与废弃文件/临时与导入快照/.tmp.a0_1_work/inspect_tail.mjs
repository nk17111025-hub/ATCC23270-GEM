import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const path = "D:/嗜酸氧化亚铁硫杆菌/01_原始资料/2024_Khaleque/from2024-mmc1.xlsx";
const wb = await SpreadsheetFile.importXlsx(await FileBlob.load(path));
for (const spec of [
  ["Table S1", "A570:O632"],
  ["Table S2", "A560:K581"],
  ["Table S3", "A600:V633"],
  ["Table S4", "A600:V633"],
]) {
  const [name, range] = spec;
  const sheet = wb.worksheets.getItem(name);
  console.log(`===== ${name} ${range} =====`);
  console.log(JSON.stringify(sheet.getRange(range).values));
}
