import fs from "node:fs/promises";
import path from "node:path";
import crypto from "node:crypto";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "D:/嗜酸氧化亚铁硫杆菌";
const outputPath = path.join(root, "10_人工审计", "人工审计表.xlsx");
const previewPath = path.join(root, ".tmp.artifact_builder_8b93883e86f7405d97b3dadc9c9aebbc", "验收汇总预览.png");

function parseTsv(text) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/).filter(line => line.length > 0);
  if (!lines.length) return [];
  const split = (line) => {
    const cells = [];
    let cur = "";
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        if (quoted && line[i + 1] === '"') { cur += '"'; i++; }
        else { quoted = !quoted; }
      } else if (ch === "\t" && !quoted) { cells.push(cur); cur = ""; }
      else { cur += ch; }
    }
    cells.push(cur);
    return cells;
  };
  const headers = split(lines[0]);
  return lines.slice(1).map(line => {
    const vals = split(line);
    return Object.fromEntries(headers.map((h, i) => [h, vals[i] ?? ""]));
  });
}

function matrix(rows, headers) { return [headers, ...rows.map(r => headers.map(h => r[h] ?? ""))]; }
function colLetter(n) {
  let s = ""; let x = n;
  while (x > 0) { const r = (x - 1) % 26; s = String.fromCharCode(65 + r) + s; x = Math.floor((x - 1) / 26); }
  return s;
}
function applyTableStyle(sheet, rowCount, colCount, freezeRows = 2) {
  const used = sheet.getRange(`A1:${colLetter(colCount)}${rowCount}`);
  used.format.font = { name: "Arial", size: 10, color: "#1F2937" };
  used.format.verticalAlignment = "center";
  used.format.wrapText = true;
  const header = sheet.getRange(`A1:${colLetter(colCount)}1`);
  header.format = { fill: "#1F4E78", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", verticalAlignment: "center", wrapText: true, borders: { preset: "all", style: "thin", color: "#FFFFFF" } };
  used.format.borders = { insideHorizontal: { style: "thin", color: "#D9E2F3" }, bottom: { style: "thin", color: "#9FBAD0" } };
  sheet.showGridLines = false;
  sheet.freezePanes.freezeRows(freezeRows);
  used.format.autofitColumns();
  used.format.autofitRows();
}

async function walk(dir, rel = "") {
  const entries = await fs.readdir(dir, { withFileTypes: true });
  const out = [];
  for (const entry of entries) {
    const name = entry.name;
    if (name.startsWith(".tmp.") || name === ".venv" || name === "node_modules" || name === "__pycache__") continue;
    const abs = path.join(dir, name);
    const relPath = path.join(rel, name).replaceAll("\\", "/");
    if (entry.isDirectory()) {
      out.push({ path: relPath, type: "文件夹" });
      out.push(...await walk(abs, relPath));
    } else if (entry.isFile()) {
      const stat = await fs.stat(abs);
      const hash = crypto.createHash("sha256").update(await fs.readFile(abs)).digest("hex");
      out.push({ path: relPath, type: "文件", size: stat.size, hash });
    }
  }
  return out;
}

const candidates = parseTsv(await fs.readFile(path.join(root, "09_差异比较", "候选修改.tsv"), "utf8"));
const uncertain = parseTsv(await fs.readFile(path.join(root, "10_人工审计", "标识映射不确定.tsv"), "utf8"));
const conflicts = parseTsv(await fs.readFile(path.join(root, "10_人工审计", "数据库冲突.tsv"), "utf8"));
const pending = parseTsv(await fs.readFile(path.join(root, "10_人工审计", "待补文献.tsv"), "utf8"));
const usageSeed = parseTsv(await fs.readFile(path.join(root, "00_项目总控", "文件夹与文件用途.tsv"), "utf8"));
const usageByPath = new Map(usageSeed.map(r => [r["相对路径"], r]));

const tasks = [
  { task: "A0", type: "Codex", status: "部分完成", detail: "原始资料冻结通过；2024场景复现18/24；模型副本已生成并解析，仍有未通过场景。", evidence: "08_模型基线/A0；12_验收报告/A0" },
  { task: "B1", type: "Work", status: "已完成", detail: "ATCC 23270身份、当前/旧组装、下载解压、SHA256和旧新位点映射已落盘。", evidence: "02_NCBI；12_验收报告/B1" },
  { task: "C2", type: "Work", status: "访问失败", detail: "BioCyc官方页面连续超时，PGDB身份/版本/导出未核实；未使用MetaCyc替代。", evidence: "03_BioCyc；12_验收报告/C2" },
  { task: "D3", type: "Work", status: "部分完成", detail: "Rhea样例成功；KEGG返回非目标菌株；BRENDA具体条目未核实；MetaCyc失败；BioCyc依赖阻塞。", evidence: "04_KEGG—07_MetaCyc；12_验收报告/D3" },
  { task: "E4", type: "Work", status: "基础阶段完成", detail: `候选修改${candidates.length}条；标识不确定${uncertain.length}条；数据库冲突${conflicts.length}条；待补文献${pending.length}条；最终差异等待D3收束。`, evidence: "09_差异比较；10_人工审计；12_验收报告/E4" },
  { task: "总控", type: "本机会话", status: "部分完成", detail: "已生成审计工作簿；最终裁决受A0未通过场景、C2访问阻塞和D3部分完成状态限制。", evidence: "10_人工审计/人工审计表.xlsx；12_验收报告/最终验收.md" },
];

const workbook = Workbook.create();
const summary = workbook.worksheets.add("验收汇总");
const audit = workbook.worksheets.add("人工审计表");
const usage = workbook.worksheets.add("文件用途表");
const sources = workbook.worksheets.add("来源与限制");

summary.getRange("A1:F1").values = [["ATCC23270 项目验收汇总", "", "", "", "", ""]];
summary.mergeCells("A1:F1");
summary.getRange("A1:F1").format = { font: { name: "Arial", size: 14, bold: true, color: "#1F2937" }, horizontalAlignment: "left" };
summary.getRange("A2:F2").values = [["截至 2026-09-22 Asia/Shanghai；状态按实际落盘、访问结果和验收记录登记。", "", "", "", "", ""]];
summary.mergeCells("A2:F2");
summary.getRange("A3:F3").values = [["任务", "执行类型", "状态", "实际结果", "证据位置", "人工后续"]];
summary.getRange("A4:F9").values = tasks.map(t => [t.task, t.type, t.status, t.detail, t.evidence, t.task === "总控" ? "继续处理阻塞后再做最终裁决" : t.status.includes("完成") ? "按证据继续人工审计" : "等待官方访问或依赖恢复"]);
summary.getRange("A11:B15").values = [["关键计数", "数量"], ["E4候选修改", candidates.length], ["标识映射不确定", uncertain.length], ["数据库冲突", conflicts.length], ["待补文献", pending.length]];
summary.getRange("D11:F15").values = [["主要限制", "状态", "记录"], ["BioCyc同株PGDB", "访问失败", "官方页面连续超时"], ["MetaCyc相关数据", "访问失败", "请求失败，未作同株替代"], ["KEGG afe", "非目标菌株", "返回ATCC 53993，未计入同株证据"], ["A0场景复现", "部分通过", "18/24通过，未宣告完全复现"]];
summary.getRange("A3:F3").format = { fill: "#1F4E78", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center", verticalAlignment: "center", wrapText: true };
summary.getRange("A11:B11").format = { fill: "#5B9BD5", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center" };
summary.getRange("D11:F11").format = { fill: "#5B9BD5", font: { name: "Arial", size: 10, bold: true, color: "#FFFFFF" }, horizontalAlignment: "center" };
summary.getRange("A1:F15").format.wrapText = true;
summary.getRange("A1:F15").format.verticalAlignment = "center";
summary.getRange("A1:F15").format.autofitColumns(); summary.getRange("A1:F15").format.autofitRows();
summary.getRange("A1:F15").format.borders = { insideHorizontal: { style: "thin", color: "#D9E2F3" }, bottom: { style: "thin", color: "#9FBAD0" } };
summary.showGridLines = false;

const auditHeaders = [...Object.keys(candidates[0] ?? {}), "人工判定", "审计备注", "模型处理"];
const auditRows = candidates.map(r => [...auditHeaders.slice(0, -3).map(h => r[h] ?? ""), "待人工审计", "", "本轮不修改模型"]);
audit.getRange(`A1:${colLetter(auditHeaders.length)}${auditRows.length + 1}`).values = [auditHeaders, ...auditRows];
applyTableStyle(audit, auditRows.length + 1, auditHeaders.length, 1);
audit.getRange(`A1:${colLetter(auditHeaders.length)}1`).format.rowHeight = 34;

const fileEntries = await walk(root);
const usageRows = fileEntries.map(e => {
  const seed = usageByPath.get(e.path);
  let purpose = seed?.["用途"] ?? "";
  if (!purpose) {
    if (e.path.startsWith("02_NCBI/")) purpose = "NCBI同株下载、标准化数据、基因映射与验收记录";
    else if (e.path.startsWith("03_BioCyc/")) purpose = "BioCyc同株PGDB核验与访问阻塞记录";
    else if (e.path.startsWith("04_KEGG/")) purpose = "KEGG来源核验与有限查询结果";
    else if (e.path.startsWith("05_BRENDA/")) purpose = "BRENDA来源核验与访问状态";
    else if (e.path.startsWith("06_Rhea/")) purpose = "Rhea反应标准化、查询结果与交叉引用";
    else if (e.path.startsWith("07_MetaCyc/")) purpose = "MetaCyc来源核验与访问失败记录";
    else if (e.path.startsWith("08_模型基线/")) purpose = "2016基线、2024论文复原、场景配置与复现诊断";
    else if (e.path.startsWith("09_差异比较/")) purpose = "候选变化与差异结构，不直接修改模型";
    else if (e.path.startsWith("10_人工审计/")) purpose = "人工审计输入、冲突、映射不确定和待补文献";
    else if (e.path.startsWith("11_脚本与环境/")) purpose = "项目内脚本、依赖说明与运行环境";
    else if (e.path.startsWith("12_验收报告/")) purpose = "任务验收、失败、缺失和总控记录";
    else if (e.path.startsWith("00_项目总控/")) purpose = "任务书、依赖、配置、状态和文件索引";
    else if (e.path === "启动任务.md") purpose = "本机总控启动任务与执行边界";
    else purpose = "项目根目录原始资料或启动包";
  }
  const state = e.type === "文件夹" ? "已落盘" : (seed?.["真实状态"]?.includes("未写D盘") ? "已落盘，状态字段已过时" : "已落盘");
  return [e.path, e.type, purpose, state, e.size ?? "", e.hash ?? ""];
});
const usageHeaders = ["相对路径", "类型", "用途", "真实状态", "字节数", "SHA256"];
usage.getRange(`A1:${colLetter(usageHeaders.length)}${usageRows.length + 1}`).values = [usageHeaders, ...usageRows];
applyTableStyle(usage, usageRows.length + 1, usageHeaders.length, 1);
usage.getRange("A1:F1").format.rowHeight = 30;

const sourceRows = [
  ["NCBI", "https://www.ncbi.nlm.nih.gov/", "已核实并下载", "ATCC 23270；GCA/GCF_049532655.1与旧版均已保存", "B1"],
  ["BioCyc", "https://biocyc.org/", "访问失败", "官方页面连续超时，PGDB版本和导出未核实", "C2"],
  ["KEGG", "https://www.kegg.jp/", "非目标菌株", "查询返回 ATCC 53993，未计入同株证据", "D3"],
  ["BRENDA", "https://www.brenda-enzymes.org/", "入口可访问，条目未核实", "没有保存项目相关具体官方条目", "D3"],
  ["Rhea", "https://www.rhea-db.org/", "有限成功", "保存5条样例和SHA256，仅用于反应标准化/交叉引用", "D3"],
  ["MetaCyc", "https://metacyc.org/", "访问失败", "请求失败，未用全物种页面替代同株PGDB", "D3"],
];
sources.getRange("A1:E7").values = [["来源", "官方入口", "状态", "可用内容或限制", "任务"], ...sourceRows];
applyTableStyle(sources, 7, 5, 1);
sources.getRange("A1:E7").format.autofitColumns(); sources.getRange("A1:E7").format.autofitRows();

workbook.recalculate();
const inspect = await workbook.inspect({ kind: "sheet,table", maxChars: 5000, tableMaxRows: 5, tableMaxCols: 8, tableMaxCellChars: 100 });
console.log(inspect.ndjson);
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
console.log(errors.ndjson);
const preview = await workbook.render({ sheetName: "验收汇总", autoCrop: "all", scale: 1, format: "png" });
await fs.writeFile(previewPath, new Uint8Array(await preview.arrayBuffer()));
await fs.mkdir(path.dirname(outputPath), { recursive: true });
try { await fs.access(outputPath); throw new Error("目标工作簿已存在，停止以避免覆盖"); } catch (err) { if (err.code !== "ENOENT") throw err; }
const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(JSON.stringify({ outputPath, candidateCount: candidates.length, usageCount: usageRows.length, previewPath }, null, 2));
