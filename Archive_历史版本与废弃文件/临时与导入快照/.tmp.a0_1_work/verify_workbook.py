from pathlib import Path
import sys, hashlib, openpyxl
sys.stdout.reconfigure(encoding="utf-8")

src=Path(r"D:\嗜酸氧化亚铁硫杆菌\01_原始资料\2024_Khaleque\from2024-mmc1.xlsx")
dst=Path(r"D:\嗜酸氧化亚铁硫杆菌\08_模型基线\A0-1\2024补充表_S3S4纠正版.xlsx")
expected_hash="88735E363ECE058F483F526355ED2325B7B2D4450E0B9AAC381522F278DE86C9"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
    return h.hexdigest().upper()

assert sha(src)==expected_hash, sha(src)
assert dst.exists() and dst.stat().st_size>0
srcwb=openpyxl.load_workbook(src,data_only=False,read_only=True)
dstwb=openpyxl.load_workbook(dst,data_only=False,read_only=True)
assert srcwb.sheetnames==dstwb.sheetnames==["Table S1","Table S2","Table S3","Table S4"]

def rows(ws):
    return [[c.value for c in row] for row in ws.iter_rows()]

def assert_same_matrix(a, b, label):
    assert len(a)==len(b), f"row count differs: {label}"
    for ri,(ra,rb) in enumerate(zip(a,b),1):
        width=max(len(ra),len(rb))
        for ci in range(width):
            va=ra[ci] if ci<len(ra) else None
            vb=rb[ci] if ci<len(rb) else None
            assert va==vb, f"cell differs: {label}!{ri},{ci}: {va!r} != {vb!r}"

for name in ("Table S1","Table S2"):
    a,b=rows(srcwb[name]),rows(dstwb[name])
    assert_same_matrix(a,b,name)

mapping=[
    (622,"ZZ_TreT2_UDP","ZZ_glc_transport"),
    (623,"ZZ_TreT1_ADP","ZZ_glc_EX"),
    (624,"ZZ_TreYZ","ZZ_glycogen_EX"),
    (625,"ZZ_Tre_amy","ZZ_g1p_transport"),
    (626,"ZZ_Tre_cga","ZZ_g1p_EX"),
    (627,"ZZ_Tre_EX","ZZ_Tre_amy"),
    (628,"ZZ_glc_transport","ZZ_Tre_cga"),
    (629,"ZZ_glc_EX","ZZ_Tre_EX"),
    (630,"ZZ_glycogen_EX","ZZ_Tre_transport"),
    (631,"ZZ_g1p_transport","ZZ_TreT1_ADP"),
    (632,"ZZ_g1p_EX","ZZ_TreT2_UDP"),
    (633,"ZZ_Tre_transport","ZZ_TreYZ"),
]
s1_defs={r[0]:r[2] for r in rows(srcwb["Table S1"])[1:] if r and r[0]}
for sheet_name in ("Table S3","Table S4"):
    s=rows(srcwb[sheet_name]); d=rows(dstwb[sheet_name])
    assert len(s)==len(d)==633
    for ri in range(633):
        width=max(len(s[ri]),len(d[ri]))
        for ci in range(width):
            if 622<=ri+1<=633 and ci in (0,1): continue
            sv=s[ri][ci] if ci<len(s[ri]) else None
            dv=d[ri][ci] if ci<len(d[ri]) else None
            assert sv==dv, f"unexpected change {sheet_name}!{ri+1},{ci+1}"
    for row,old,new in mapping:
        assert s[row-1][0]==old
        assert d[row-1][0]==new
        assert d[row-1][1]==s1_defs[new]

print("source_sha256",sha(src))
print("dest_size",dst.stat().st_size)
print("S1_S2_UNCHANGED PASS")
print("S3_S4_ONLY_A_B_ROWS_622_633_CHANGED PASS")
print("MAPPING_DEFINITIONS_FROM_S1 PASS")
srcwb.close(); dstwb.close()
