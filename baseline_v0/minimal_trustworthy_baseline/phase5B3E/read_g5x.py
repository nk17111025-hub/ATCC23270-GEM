# -*- coding: utf-8 -*-
"""Read G5X G56 review data for the target central-carbon genes."""
import json
import pathlib

ROOT = pathlib.Path(r"D:\嗜酸氧化亚铁硫杆菌")


def main():
    d = json.loads((ROOT / "G5X_2929基因逐基因证据链/G56_结果/G56_审查数据.json").read_text(encoding="utf-8", errors="replace"))
    targets = ["AFE_1406", "AFE_1417", "AFE_1471", "AFE_1494", "AFE_0299", "AFE_2025", "AFE_3251", "AFE_3250"]
    keys = ["总序号", "当前locus", "旧AFE", "WP", "当前NCBI product", "2026动作", "2016关联反应",
            "未决点", "证据边界", "KEGG EC", "KEGG reaction", "Rhea交叉引用反应式"]
    for rec in d["index"]:
        afe = str(rec.get("旧AFE", ""))
        if afe in targets:
            out = {k: rec.get(k) for k in keys}
            print(json.dumps(out, ensure_ascii=False)[:1800])
            print("=" * 100)


if __name__ == "__main__":
    main()
