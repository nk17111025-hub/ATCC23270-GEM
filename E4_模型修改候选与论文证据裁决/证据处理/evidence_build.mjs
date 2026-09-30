import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const data=JSON.parse(await fs.readFile('evidence_data.json','utf8'));
const outDir=path.resolve('outputs/01a0d100-c9f0-7122-903b-f8175e453668');
await fs.mkdir(outDir,{recursive:true});
const wb=Workbook.create();
const sh=wb.worksheets.add('2026-extended from 2016');
const src=wb.worksheets.add('2026-61 papers sources');
sh.showGridLines=false; src.showGridLines=false;

const extra=['Record class','Candidate ID','Candidate group','2016 linked Reaction ID','2026 reaction object','Legacy AFE locus','Current locus','Primary evidence DOI','Evidence type','Evidence boundary','2026 action','2026 reaction score','Source ID','Score scope'];
const allHeaders=[...data.headers,...extra];
const row1=Array(30).fill(null);
row1[10]=data.group_titles[0]; row1[12]=data.group_titles[1]; row1[14]=data.group_titles[2]; row1[16]='2026 evidence extension';
sh.getRange('A1:AD1').values=[row1];
sh.getRange('A2:AD2').values=[allHeaders];
const sourceByDoi=new Map(data.sources.map(x=>[x.doi.toLowerCase(),x.id]));
const original=data.original.map(r=>[...r,'2016 original',...Array(13).fill(null)]);
const candidates=data.candidates.map(c=>[...c.core,'2026 candidate',c.id,c.group,c.linked_reaction,c.object,c.legacy_locus,c.locus,c.doi?'https://doi.org/'+c.doi:null,c.evidence_type,c.boundary,c.action,c.score,sourceByDoi.get(c.doi.toLowerCase())??null,c.score_text]);
sh.getRange(`A3:AD${2+original.length}`).values=original;
const candidateStart=3+original.length;
sh.getRange(`A${candidateStart}:AD${candidateStart+candidates.length-1}`).values=candidates;

// The first 16 columns reproduce the original 2016 supplementary table order and values.
sh.mergeCells('K1:L1'); sh.mergeCells('M1:N1'); sh.mergeCells('O1:P1'); sh.mergeCells('Q1:AD1');
sh.getRange('A1:AD1').format={fill:'#E6E6E6',font:{name:'Arial',size:10,bold:true,color:'#252525'},rowHeight:24};
sh.getRange('Q1:AD1').format={fill:'#DCE9F3',font:{name:'Arial',size:10,bold:true,color:'#17365D'},rowHeight:24};
sh.getRange('A2:P2').format={fill:'#F2F2F2',font:{name:'Arial',size:10,bold:true,color:'#1F1F1F'},rowHeight:32,wrapText:true};
sh.getRange('Q2:AD2').format={fill:'#DCE9F3',font:{name:'Arial',size:10,bold:true,color:'#17365D'},rowHeight:32,wrapText:true};
sh.getRange(`A3:AD${candidateStart+candidates.length-1}`).format.font={name:'Arial',size:10,color:'#242424'};
sh.getRange(`A${candidateStart}:AD${candidateStart+candidates.length-1}`).format.fill='#F6FAFD';
sh.getRange(`A${candidateStart}:AD${candidateStart+candidates.length-1}`).format.rowHeight=68;
sh.getRange(`U${candidateStart}:AD${candidateStart+candidates.length-1}`).format.wrapText=true;
sh.getRange(`R${candidateStart}:R${candidateStart+candidates.length-1}`).format.font={name:'Arial',size:10,bold:true,color:'#17365D'};
sh.getRange(`X${candidateStart}:X${candidateStart+candidates.length-1}`).format.font={name:'Arial',size:10,color:'#0563C1'};
sh.getRange('A:A').format.columnWidth=23;
sh.getRange('B:B').format.columnWidth=42;
sh.getRange('C:C').format.columnWidth=79;
sh.getRange('D:F').format.columnWidth=17;
sh.getRange('G:G').format.columnWidth=43;
sh.getRange('H:J').format.columnWidth=46;
sh.getRange('K:P').format.columnWidth=12;
sh.getRange('Q:R').format.columnWidth=19;
sh.getRange('S:T').format.columnWidth=20;
sh.getRange('U:U').format.columnWidth=85;
sh.getRange('V:W').format.columnWidth=31;
sh.getRange('X:X').format.columnWidth=48;
sh.getRange('Y:Y').format.columnWidth=40;
sh.getRange('Z:Z').format.columnWidth=100;
sh.getRange('AA:AD').format.columnWidth=26;
sh.freezePanes.freezeRows(2);

const sourceHeaders=['Source ID','Full citation','Year','DOI','Source type','Strain','Experimental type','Genes / proteins studied','Related Reaction IDs / candidate IDs','What was directly measured','Applicable confidence level','Evidence limitation'];
src.getRange('A1:L1').values=[sourceHeaders];
const srows=data.sources.map(s=>[s.id,s.citation,s.year,'https://doi.org/'+s.doi,s.source_type,s.strain,s.experimental_type,s.genes,s.related,s.measured,s.applicable_confidence,s.limitation]);
src.getRange(`A2:L${srows.length+1}`).values=srows;
src.getRange('A1:L1').format={fill:'#E6E6E6',font:{name:'Arial',size:10,bold:true,color:'#1F1F1F'},rowHeight:34,wrapText:true};
src.getRange(`A2:L${srows.length+1}`).format={font:{name:'Arial',size:10,color:'#242424'},wrapText:true,rowHeight:72};
src.getRange(`D2:D${srows.length+1}`).format.font={name:'Arial',size:10,color:'#0563C1'};
src.getRange('A:A').format.columnWidth=12;
src.getRange('B:B').format.columnWidth=95;
src.getRange('C:C').format.columnWidth=12;
src.getRange('D:D').format.columnWidth=50;
src.getRange('E:E').format.columnWidth=36;
src.getRange('F:F').format.columnWidth=42;
src.getRange('G:G').format.columnWidth=47;
src.getRange('H:H').format.columnWidth=55;
src.getRange('I:I').format.columnWidth=64;
src.getRange('J:J').format.columnWidth=55;
src.getRange('K:K').format.columnWidth=34;
src.getRange('L:L').format.columnWidth=75;
src.freezePanes.freezeRows(1);

wb.recalculate();
const a=await wb.inspect({kind:'table',range:`'2026-extended from 2016'!Q${candidateStart}:AD${candidateStart+3}`,include:'values,formulas',tableMaxRows:4,tableMaxCols:14,maxChars:3000});
console.log('main check',a.ndjson);
const b=await wb.inspect({kind:'table',range:'\'2026-61 papers sources\'!A1:L3',include:'values,formulas',tableMaxRows:3,tableMaxCols:12,maxChars:3000});
console.log('sources check',b.ndjson);
for(const [name,range,file] of [['2026-extended from 2016','A1:P7','2016-header.png'],['2026-extended from 2016',`Q${candidateStart-1}:AD${candidateStart+4}`,'candidate-preview.png'],['2026-61 papers sources','A1:L4','sources-preview.png']]){
  const preview=await wb.render({sheetName:name,range,scale:1.2,format:'png'});
  await fs.writeFile(path.join(outDir,file),new Uint8Array(await preview.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
const dest=path.join(outDir,'ATCC23270-2026_最终模型证据表.xlsx');
await xlsx.save(dest);
console.log('EXPORTED',dest,'rows',original.length,candidates.length,srows.length);
