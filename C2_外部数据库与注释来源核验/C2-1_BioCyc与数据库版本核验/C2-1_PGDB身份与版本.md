# C2-1 PGDB身份与版本

验收对象：Acidithiobacillus ferrooxidans ATCC 23270

## 当前官方身份

- 用户预登记标识：`AFER243159`
- 当前 BioCyc 官方数据库列表和摘要页中的 PGDB ID：`GCF_000021485`
- 当前官方列表：菌种 `Acidithiobacillus ferrooxidans`，菌株 `ATCC 23270`，序列来源 `GCF_000021485.1`
- NCBI Taxonomy ID：`243159`
- 当前 PGDB 版本：`30.0`
- 数据库等级：`Tier 3 Uncurated Organism Database`
- 摘要页：`https://biocyc.org/GCF_000021485/organism-summary`
- 数据库列表原始 XML：`官方原始下载\AFER243159_数据库列表_20260923_082931.txt`
- 版本接口原始 JSON：`官方原始下载\GCF_000021485_kb-version_20260923_083031.json`

`AFER243159` 在当前官方数据库列表中未出现；对该标识调用 `kb-version` 得到 HTTP 404 和 “Organism not found” 页面。为避免把旧标识写成当前官方标识，本快照以官方当前页面和接口实际返回的 `GCF_000021485` 为准，保留 `AFER243159` 失败证据。

## 生成和注释元数据

- PGDB 生成日期：`2021-11-18`
- 生成工具：PathoLogic，生成时 Pathway Tools `25.5`
- 生成时 MetaCyc 版本：`25.5`
- 当前摘要页生成工具：Pathway Tools `30.0`
- 当前摘要页生成时间：`2026-09-22`
- 摘要页发行标记：`BIOCYC14A`
- 序列来源：`GCF_000021485.1`
- 复制子：RefSeq `NC_011761`
- 基因组长度：`2,982,397 bp`
- NCBI BioSample：`SAMN02603974`
- NCBI BioProject：`PRJNA224116`
- NCBI Taxonomy：`243159`
- 注释提供方：`NCBI RefSeq`
- 注释日期：`2021-01-24`
- 注释流水线：`NCBI Prokaryotic Genome Annotation Pipeline (PGAP)`
- 注释流水线版本：`5.0`
- 注释备注：`Best-placed reference protein set; GeneMarkS-2+`

## 官方摘要统计

| 项目 | 数量 |
|---|---:|
| Total Genes | 3,075 |
| Protein Genes | 2,893 |
| RNA Genes | 92 |
| Other Genes | 90 |
| Pathways | 217 |
| Enzymatic Reactions | 1,135 |
| Transport Reactions | 31 |
| Polypeptides | 2,897 |
| Protein Complexes | 7 |
| Enzymes | 704 |
| Transporters | 267 |
| Compounds | 778 |
| Transcription Units | 1,545 |
| tRNAs | 82 |
| Protein Features | 4,185 |

## 批量 XML 实际结果统计

批量 XML 结果保留原始查询口径，不与摘要页统计强行对齐：`All-Genes` 3,084、`Genes` 2,989、`Proteins` 3,235、`Reactions` 1,313、`Pathways` 953、`Enzymatic-Reactions` 1,378；反应左右侧去重后的网络化合物 1,113。摘要页统计和 XML 类别统计存在类别/层级口径差异，均已保留。