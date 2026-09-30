from pathlib import Path
import sys
import re
from pypdf import PdfReader

sys.stdout.reconfigure(encoding="utf-8")

pdf = Path(r"D:\嗜酸氧化亚铁硫杆菌\Simulating compatible solute biosynthesis using a metabolic flux model of the biomining acidophile, Acidithiobacillus ferrooxidans ATCC 23270.pdf")
reader = PdfReader(str(pdf))
terms = re.compile(r"treyz|tret|trehalose|figure\s*3|fig\.\s*3|30%|methods", re.I)
print(f"pages={len(reader.pages)}")
for i, page in enumerate(reader.pages, 1):
    text = page.extract_text() or ""
    if terms.search(text):
        print(f"\n===== PAGE {i} =====\n{text}")
