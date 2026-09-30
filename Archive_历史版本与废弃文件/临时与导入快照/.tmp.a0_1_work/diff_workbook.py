from pathlib import Path
import sys, openpyxl
sys.stdout.reconfigure(encoding="utf-8")
src=Path(r"D:\嗜酸氧化亚铁硫杆菌\01_原始资料\2024_Khaleque\from2024-mmc1.xlsx")
dst=Path(r"D:\嗜酸氧化亚铁硫杆菌\08_模型基线\A0-1\2024补充表_S3S4纠正版.xlsx")
a=openpyxl.load_workbook(src,data_only=False,read_only=True)
b=openpyxl.load_workbook(dst,data_only=False,read_only=True)
for name in a.sheetnames:
    wa,wb=a[name],b[name]
    ar=[[c.value for c in row] for row in wa.iter_rows()]
    br=[[c.value for c in row] for row in wb.iter_rows()]
    print(name, len(ar), max(len(x) for x in ar), len(br), max([len(x) for x in br] or [0]))
    count=0
    for r in range(max(len(ar),len(br))):
        for c in range(max(len(ar[r]) if r<len(ar) else 0,len(br[r]) if r<len(br) else 0)):
            va=ar[r][c] if r<len(ar) and c<len(ar[r]) else None
            vb=br[r][c] if r<len(br) and c<len(br[r]) else None
            if va!=vb:
                print("DIFF",name,r+1,c+1,repr(va),type(va).__name__,repr(vb),type(vb).__name__)
                count+=1
                if count>=20: break
        if count>=20: break
    print("count_first",count)
