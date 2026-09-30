# MetaCyc D3-1 方法记录

本轮没有对 MetaCyc 名称做模糊搜索。按任务要求复用 C2-1 BioCyc PGDB 已有的 Pathway Tools/MetaCyc reaction IDs、反应方向、左右底物、EC、通路、证据代码和文献字段，形成 `MetaCyc_BioCyc_ID精确映射.tsv`。

PGDB 为 `GCF_000021485`，version 30.0、Tier 3 uncurated，依赖旧 assembly `GCF_000021485.1`。因此这些 ID/反应定义可用于标准化候选和通路背景，不能独立证明最新 ATCC 23270 genome 拥有该反应。主表只把 exact reaction ID 计为 ID 覆盖，未将 gene-level enzyme association 当成反应等价。

当前 MetaCyc 定义、文献和许可状态未直接访问 MetaCyc 官方条目，因此覆盖率仅报告 C2-1 已有 ID 的连接，不报告独立 MetaCyc 页面覆盖。
