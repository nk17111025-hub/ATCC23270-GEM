# 环境重建说明

## 当前实际使用的环境

本阶段使用既有环境，未改动任何包版本：

- Python：3.14.6
- 解释器：`Archive_历史版本与废弃文件\A0\运行环境\.venv\Scripts\python.exe`
- cobra：0.32.1
- optlang：1.9.1
- python-libsbml：5.21.2（libsbml 字符串版本 52102）
- xlrd：2.0.2
- openpyxl：3.1.5
- solver：GLPK 5.0（swiglpk 后端，optlang.glpk_interface；另有 glpk_exact）
- cobra Configuration：solver=glpk，tolerance=1e-07，bounds=±1000，processes=13

## 锁定文件

- `requirements_locked.txt`：`pip freeze` 全量输出，用于日后复现同版本。
- `environment_snapshot.txt`：机器可读的运行环境快照。

## 重建命令（参考）

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements_locked.txt
```

## 注意事项

- cobra 0.32.1 在 Python 3.14 上已可用，但属于较新组合；重建时若安装失败，需要先确认
  该包版本与 Python 3.14 的兼容性，不要擅自升级到不同版本。
- GLPK 通过 optlang 的 swiglpk 提供；swiglpk 版本未直接暴露，GLPK 库版本为 5.0。
- 可行性/最优性/整数容差没有通过 cobra.Configuration 暴露，已标记 UNKNOWN，
  需在 Phase 1 用直接 GLPK 探针或另行核实，不要自行填写。
