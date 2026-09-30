from pathlib import Path
import sys
from pypdf import PdfReader
sys.stdout.reconfigure(encoding="utf-8")
p=Path(r"D:\嗜酸氧化亚铁硫杆菌\Simulating compatible solute biosynthesis using a metabolic flux model of the biomining acidophile, Acidithiobacillus ferrooxidans ATCC 23270.pdf")
t=(PdfReader(str(p)).pages[5].extract_text() or "").splitlines()
for i,line in enumerate(t):
    print(f"{i+1:03d}: {line}")
