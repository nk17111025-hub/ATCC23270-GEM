from pathlib import Path
import re, sys
from pypdf import PdfReader
sys.stdout.reconfigure(encoding="utf-8")
pdf = Path(r"D:\嗜酸氧化亚铁硫杆菌\Simulating compatible solute biosynthesis using a metabolic flux model of the biomining acidophile, Acidithiobacillus ferrooxidans ATCC 23270.pdf")
reader = PdfReader(str(pdf))
pat = re.compile(r"TreYZ|TreT|trehalose|Figure\s*3|Fig\.\s*3|30%|maximum biomass|pathway", re.I)
for page_no in range(3, 8):
    text = reader.pages[page_no-1].extract_text() or ""
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines) if pat.search(line)]
    if not hits:
        continue
    print(f"\n===== PAGE {page_no} =====")
    seen = set()
    for i in hits:
        lo, hi = max(0, i-3), min(len(lines), i+4)
        if any(j in seen for j in range(lo, hi)):
            continue
        seen.update(range(lo, hi))
        print("\n".join(lines[lo:hi]))
