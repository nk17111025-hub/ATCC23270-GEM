# ATCC23270-GEM 项目资料

本仓库保存 *Acidithiobacillus ferrooxidans* ATCC 23270 代谢模型项目的原始资料、证据链、历史记录、计算脚本和结果。`ATCC23270-GEM/` 是计算工作目录；权威宿主模型位于 `before final/updatedv3.1_verified/updatedv3.1.xml`。

主要计算结果位于 `ATCC23270-GEM/results/engineering/`。第三方工具以 Git 子模块记录，克隆时使用 `git clone --recurse-submodules`。E3.0 所用 StrainDesign 本地兼容性补丁保存在 `patches/straindesign_E3_solver_status.patch`。

本仓库公开。虚拟环境、`node_modules` 和临时同步缓存可从依赖重新生成，未纳入版本记录。
