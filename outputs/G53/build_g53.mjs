import fs from "node:fs/promises";
import {fileURLToPath} from "node:url";
import {Workbook,SpreadsheetFile} from "@oai/artifact-tool";

const dir = new URL(".", import.meta.url);
const data=JSON.parse(await fs.readFile(new URL("g53_data.json",dir),"utf8"));
const files=["agent_587_684.json","agent_685_782.json","agent_783_879.json"];
let agent=[];
for(const fn of files){
  const raw=JSON.parse(await fs.readFile(new URL(fn,dir),"utf8"));
  agent.push(...(Array.isArray(raw)?raw:raw.records));
}
agent.sort((a,b)=>Number(a.seq)-Number(b.seq));
const byLocus=new Map(agent.map(a=>[a.locus,a]));
const indexHeaders=["总序号","Current locus","WP protein ID","start","end","strand","Current product","Legacy AFE locus","AFE mapping status","GFF identity check","KEGG KO","KEGG EC","KEGG reaction","BioCyc reaction","Rhea exact crossref","Rhea EC-only candidate","2016 GPR reactions","Exact reaction status","2026 action","关联审查行 ID","2016 linked reaction IDs","Unresolved points","Source IDs","UniProt accession","UniProt match basis","UniProt reviewed","UniProt protein name","InterPro IDs","Pfam IDs","TCDB","UniProt EC","Same-strain proteomics","Europe PMC hit count","DOI hits—see paper audit","Europe PMC query URL","D3 BRENDA EC crossref"];
const reviewHeaders=["Reaction ID","Reaction Name","Reaction Formula","Confidence Level","EC Number","PMID","Subsystem","Gene-Reaction Association","Gene-Protein-Reaction Association","Protein-Reaction-Association","Fe²⁺ lb","Fe²⁺ ub","tetrathionate lb","tetrathionate ub","sulfur lb","sulfur ub","Record class","Candidate ID","Candidate group","2016 linked Reaction ID","2026 reaction object","Legacy AFE locus","Current locus","Primary evidence DOI","Evidence type","Evidence boundary","2026 action","2026 reaction score","Source ID","Score scope"];
const sourceHeaders=["Source ID","Full citation","Year","DOI","Source type","Strain","Experimental type","Genes / proteins studied","Related Reaction IDs / candidate IDs","What was directly measured","Applicable confidence level","Evidence limitation"];
const auditHeaders=["总序号","Current locus","分段结论","结论依据","待核事项","来源链接"];
const audit=agent.map(a=>{
  const seq=Number(a.seq);
  const reason=seq===654?
    "同株 DSM 14882 的 AFE_0702 qRT-PCR 显示转录，培养物气相色谱显示低浓度 H2 摄取和 H2 促进生长；该株其他氢化酶也表达。基因专属电子受体和模型反应仍未确定。":a.reason||"";
  const questions=seq===654?
    "需 AFE_0702 敲除/互补、纯化蛋白受体及膜侧/H+ 计量。":
    (a.open_questions||[]).map(x=>typeof x==="string"?x:JSON.stringify(x)).join("；");
  const links=(a.evidence_links||[]).map(x=>typeof x==="string"?x:x.url||JSON.stringify(x));
  if(seq===654)links.push("https://pmc.ncbi.nlm.nih.gov/articles/PMC7784904/");
  return [seq,a.locus,a.decision||"",reason,questions,links.join("；")];
});
const wb=Workbook.create();
const sSummary=wb.worksheets.add("G53验收");
const sIndex=wb.worksheets.add("G53基因索引");
const sReview=wb.worksheets.add("2026-extended from 2016");
const sSources=wb.worksheets.add("2026-61 papers sources");
const sAudit=wb.worksheets.add("G53分段核查");
const sLiterature=wb.worksheets.add("G53论文命中核查");
for(const s of [sSummary,sIndex,sReview,sSources,sAudit,sLiterature])s.showGridLines=false;
function write(s,header,rows){
  s.getRangeByIndexes(0,0,rows.length+1,header.length).values=[header,...rows];
  const head=s.getRangeByIndexes(0,0,1,header.length);
  head.format={fill:"#263B5B",font:{name:"Arial",size:10,bold:true,color:"#FFFFFF"},rowHeight:30,verticalAlignment:"center"};
  s.getRangeByIndexes(1,0,rows.length,header.length).format.font={name:"Arial",size:10,color:"#17243A"};
  s.getRangeByIndexes(1,0,rows.length,header.length).format.rowHeight=23;
  s.getRangeByIndexes(0,0,rows.length+1,header.length).format.verticalAlignment="center";
  s.freezePanes.freezeRows(1);
}
write(sIndex,indexHeaders,data.index);
write(sReview,reviewHeaders,data.review);
write(sSources,sourceHeaders,data.sources);
write(sAudit,auditHeaders,audit);
write(sLiterature,["DOI","关联基因数","G53 序号与位点","论文对象/实验","本组证据边界","原文链接"],data.literature);
sIndex.freezePanes.freezeColumns(2);
sReview.freezePanes.freezeColumns(1);
sAudit.freezePanes.freezeColumns(2);
for(const s of [sIndex,sReview,sSources,sAudit,sLiterature]){
  const n=s===sIndex?36:s===sReview?30:s===sSources?12:6;
  const rowCount=s===sIndex?data.index.length:s===sReview?data.review.length:s===sSources?data.sources.length:s===sAudit?audit.length:data.literature.length;
  s.getRangeByIndexes(0,0,1,n).format.wrapText=true;
  s.getRangeByIndexes(0,0,rowCount+1,n).format.columnWidth=17;
}
sIndex.getRange("B:B").format.columnWidth=22;
sIndex.getRange("C:C").format.columnWidth=20;
sIndex.getRange("G:G").format.columnWidth=44;
sIndex.getRange("S:S").format.columnWidth=28;
sIndex.getRange("V:V").format.columnWidth=70;
sIndex.getRange("AA:AA").format.columnWidth=50;
sIndex.getRange("AF:AF").format.columnWidth=30;
sIndex.getRange("AG:AG").format.columnWidth=21;
sIndex.getRange("AH:AH").format.columnWidth=60;
sIndex.getRange("AI:AI").format.columnWidth=72;
sIndex.getRange("AJ:AJ").format.columnWidth=32;
sReview.getRange("B:B").format.columnWidth=40;
sReview.getRange("C:C").format.columnWidth=62;
sReview.getRange("H:J").format.columnWidth=34;
sReview.getRange("Z:Z").format.columnWidth=66;
sReview.getRange("AC:AD").format.columnWidth=48;
sSources.getRange("B:B").format.columnWidth=72;
sSources.getRange("G:L").format.columnWidth=48;
sAudit.getRange("A:B").format.columnWidth=22;
sAudit.getRange("C:C").format.columnWidth=40;
sAudit.getRange("D:F").format.columnWidth=95;
sLiterature.getRange("A:A").format.columnWidth=38;
sLiterature.getRange("C:C").format.columnWidth=85;
sLiterature.getRange("D:E").format.columnWidth=70;
sLiterature.getRange("F:F").format.columnWidth=75;
const n=data.index.length;
const counts={};
for(const r of data.index)counts[r[18]]=(counts[r[18]]||0)+1;
const exact=data.index.filter(r=>r[24]==="exact WP").length;
const oldOnly=data.index.filter(r=>r[24]==="legacy AFE only; WP mismatch").length;
const seen2016=data.index.filter(r=>r[16]).length;
const pw=data.index.filter(r=>r[31]).length;
const paperHits=data.index.filter(r=>Number(r[32])>0).length;
const rheaPairs=data.review.filter(r=>r[16]==="database-reaction-candidate").length;
const summary=[
  ["G53｜ATCC 23270 第3组",""],
  ["责任范围","总序号 587—879；293 个蛋白编码基因"],
  ["清单核对",`${n} 条；序号连续；RU820 唯一`],
  ["同版 GFF/CDS",`${data.index.filter(r=>r[9]==="一致").length}/293 身份、坐标、链向、WP、product 一致`],
  ["历史 AFE",`${data.index.filter(r=>r[8]==="已核验").length} 条 B1 核验；余者保留未确认/歧义`],
  ["2016 原表",`${seen2016} 个基因有已核验 GPR；${data.review.filter(r=>r[16]==="reaction-review").length} 条逐基因反应关联；PHFT 另有一处歧义映射；A—P 原行与 mmc1 Table 1 全表 615 行一致`],
  ["Rhea 数据库候选",`${rheaPairs} 条基因–反应交叉引用；通用化学式尚未确认为本株功能`],
  ["UniProt",`${exact} 条精确 WP；${oldOnly} 条仅旧 AFE 线索；20 条已审校`],
  ["同株蛋白组",`${pw} 条 AFE/WP 精确配对；仅支持蛋白检出`],
  ["论文检索",`${paperHits} 条 Europe PMC 有 DOI 命中；${data.literature.length} 篇已按菌株/实验/证据边界核查`],
  ["证据不足暂缓",counts["证据不足暂缓"]||0],
  ["仅更新证据或编号（待确认）",counts["仅更新证据或编号（待确认）"]||0],
  ["无需修改（本轮未见代谢反应）",counts["无需修改（本轮未见代谢反应）"]||0],
  ["正式新增/修改 reaction 或 GPR",0],
  ["验收状态","293 条均有本轮判断；精确反应、基因专属实验及模型 GPR 的未决部分保留暂缓"],
  ["不入模型的未决项","CYTBO3 质子数；CYSDSS GPR/原引文；GLUSy 原引文；B-001/MDH；B-021 ATP 合酶 ε；B-022 Fe-S 电子转移；HucL 电子受体；CTPS2 方向；NRAMP 底物；HYD3pp 化学计量"],
  ["来源时间","本地 D3 2026-09-23；UniProt 2026-09-24/25；NCBI GCF_049532655.1"],
  ["交付边界","仅证据索引和候选裁决；未改 SBML，未采纳未经证实的反应"]
];
sSummary.getRangeByIndexes(0,0,summary.length,2).values=summary;
sSummary.getRange("A1").format.font={name:"Arial",size:15,bold:true,color:"#17243A"};
sSummary.getRange("A1:B1").format.rowHeight=34;
sSummary.getRange("A2:A18").format.font={name:"Arial",size:10,bold:true,color:"#263B5B"};
sSummary.getRange("A1:A18").format.columnWidth=35;
sSummary.getRange("B1:B18").format.columnWidth=105;
sSummary.getRange("A2:B18").format.rowHeight=26;
sSummary.getRange("A15:B15").format.fill="#FFF2CF";
sSummary.getRange("A16:B16").format.fill="#FFF2CF";
sSummary.getRange("A15:B16").format.wrapText=true;
sSummary.getRange("A15:B16").format.rowHeight=52;
wb.recalculate();
for(const [tag,sheet,range] of [["summary","G53验收","A1:B18"],["index","G53基因索引","A1:H8"],["review","2026-extended from 2016","Q1:AD5"],["sources","2026-61 papers sources","A1:D5"],["agent","G53分段核查","A1:D6"]]){
  const preview=await wb.render({sheetName:sheet,range,scale:1.3,format:"png"});
  await fs.writeFile(new URL(`preview_${tag}.png`,dir),new Uint8Array(await preview.arrayBuffer()));
}
const blob=await SpreadsheetFile.exportXlsx(wb);
await blob.save(fileURLToPath(new URL("G53_逐基因证据审查_交付稿.xlsx",dir)));
console.log("saved G53 workbook",data.index.length,data.review.length,data.sources.length,agent.length);
