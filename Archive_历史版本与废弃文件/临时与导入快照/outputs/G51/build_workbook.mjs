import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const folder = path.dirname(fileURLToPath(import.meta.url));
const load = async name => JSON.parse(await fs.readFile(path.join(folder, name), "utf8"));
const data = await load("G51_working_tables.json");
const sources = await load("G51_sources.json");
const inventory = await load("G51_evidence_inventory.json");
const balance = await load("2016_G51_balance.json");

const wb = Workbook.create();
const blue = "#163A5F";
const pale = "#DCEAF5";
const font = "Aptos";

function writeSheet(name, records, widths = {}) {
  const sh = wb.worksheets.add(name);
  sh.showGridLines = false;
  const headers = Object.keys(records[0]);
  const matrix = [headers, ...records.map(row => headers.map(key => {
    const value = row[key];
    return value === null || value === undefined ? "" : value;
  }))];
  sh.getRangeByIndexes(0, 0, matrix.length, headers.length).values = matrix;
  const all = sh.getRangeByIndexes(0, 0, matrix.length, headers.length);
  all.format.font = {name: font, size: 10, color: "#1D2935"};
  const head = sh.getRangeByIndexes(0, 0, 1, headers.length);
  head.format.fill = blue;
  head.format.font = {name: font, size: 10, bold: true, color: "#FFFFFF"};
  head.format.rowHeight = 30;
  sh.freezePanes.freezeRows(1);
  for (let col = 0; col < headers.length; col++) {
    const width = widths[headers[col]] ?? 18;
    sh.getRangeByIndexes(0, col, matrix.length, 1).format.columnWidth = width;
  }
  sh.tables.add(sh.getRangeByIndexes(0, 0, matrix.length, headers.length), true, "T" + name.replace(/[^A-Za-z0-9]/g, ""));
  return sh;
}

const indexWidths = {"总序号": 9, "当前locus": 19, "当前WP": 19, "坐标": 31, "当前product": 48,
  "旧AFE": 18, "UniProt": 18, "InterPro": 24, "Pfam": 20, "2016直接关联反应": 25,
  "精确反应状态": 44, "原始论文直接事实": 70, "证据边界": 90, "本轮动作": 55,
  "文献候选DOI": 50, "文献检索链接": 80, "未决点": 90};
writeSheet("基因索引", data.index, indexWidths);
writeSheet("Alpha结构审查", data.alpha, {"Reaction ID": 18, "Reaction Name": 44, "Reaction Formula": 75,
  "Gene-Reaction Association": 52, "Gene-Protein-Reaction Association": 58, "Record class": 35,
  "Candidate ID": 21, "2026 reaction object": 75, "Evidence boundary": 95, "2026 action": 63, "Score scope": 68});

const rheaRecords = inventory.flatMap(g => g.rhea_definitions.map(r => ({
  "总序号": g.ordinal, "当前locus": g.locus, "当前WP": g.wp, "当前product": g.product,
  "Rhea ID": r["Reaction identifier"], "Equation": r.Equation, "EC": r["EC number"],
  "ChEBI": r["ChEBI identifier"], "KEGG": r["Cross-reference (KEGG)"],
  "MetaCyc": r["Cross-reference (MetaCyc)"], "2016同EC比较对象": g.model_2016_same_ec.map(x => x["Reaction ID"]).join("; "),
  "证据边界": "反应定义/EC 来源；未证明当前 WP 催化本式，也未证明底物、区室或方向。",
  "原始链接": "https://www.rhea-db.org/rhea/" + r["Reaction identifier"].replace("RHEA:", ""),
  "本地源文件": r.source_file,
})));
writeSheet("Rhea反应定义", rheaRecords, {"当前locus": 19, "当前WP": 19, "当前product": 47, "Equation": 100,
  "ChEBI": 58, "2016同EC比较对象": 29, "证据边界": 85, "原始链接": 45, "本地源文件": 70});

writeSheet("证据来源", sources, {"Source ID": 16, "Full citation": 105, "DOI": 40, "Source type": 32,
  "Strain": 38, "Experimental type": 70, "Genes / proteins studied": 43, "Related Reaction IDs / candidate IDs": 48,
  "What was directly measured": 85, "Applicable confidence level": 58, "Evidence limitation": 90, "URL": 75});

const issueRecords = inventory.map(g => {
  const idx = data.index[g.ordinal - 1];
  const product = g.product.toLowerCase();
  let experiment = "验证具体底物/产物和此 WP 的催化活性；核对区室与辅因子。";
  if (/transport|permease|efflux|porin|channel|sulp|aquaporin/.test(product))
    experiment = "测底物特异性、摄取/外排方向、膜侧、离子或 ATP 耦联；查 TCDB 同源分类。";
  if (/hypothetical|domain-containing|family protein/.test(product))
    experiment = "先确定蛋白家族/结构域和底物候选，再做直接活性或遗传实验。";
  return {"总序号": g.ordinal, "当前locus": g.locus, "当前product": g.product, "未决点": idx["未决点"],
    "推荐验证": experiment, "目前动作": idx["本轮动作"], "审查状态": idx["审查状态"]};
});
writeSheet("冲突缺证", issueRecords, {"当前locus": 19, "当前product": 50, "未决点": 105,
  "推荐验证": 85, "目前动作": 72, "审查状态": 40});
writeSheet("2016平衡检查", balance.map(r => ({"Reaction ID": r.reaction_id, "2016 Reaction Formula": r.formula,
  "检查结果": r.status, "元素差": JSON.stringify(r.element_delta), "电荷差": r.charge_delta,
  "缺失代谢物": r.missing_metabolites.join("; "), "无式代谢物": r.unknown_formula.join("; ")})),
  {"Reaction ID": 20, "2016 Reaction Formula": 100, "检查结果": 23, "元素差": 37, "缺失代谢物": 40, "无式代谢物": 70});

wb.recalculate();
const report = await wb.inspect({kind: "sheet", include: "id,name", maxChars: 3000});
console.log(report.ndjson);
const output = await SpreadsheetFile.exportXlsx(wb);
const target = path.join(folder, "G51_初核工作簿.xlsx");
await output.save(target);
console.log(target);
