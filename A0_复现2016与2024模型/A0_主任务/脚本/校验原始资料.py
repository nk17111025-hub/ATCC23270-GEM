from pathlib import Path
import argparse, csv, hashlib, sys

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1048576), b""):
            h.update(data)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="核验启动包原始资料；不下载、不修改、不启动代理。")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1], help="项目根目录")
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = root / "00_项目总控" / "原始资料SHA256清单.tsv"
    failures = []
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = (root / row["相对路径"]).resolve()
            if not path.is_relative_to(root):
                failures.append("目录越界：" + row["相对路径"])
                continue
            if not path.is_file():
                failures.append("缺失：" + row["相对路径"])
            elif digest(path) != row["SHA256"]:
                failures.append("校验不一致：" + row["相对路径"])
            else:
                print("通过：" + row["相对路径"])
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("全部原始副本校验通过。本结果不代表模型复原或代理任务已经完成。")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
