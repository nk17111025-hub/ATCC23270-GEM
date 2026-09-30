from pathlib import Path
import csv
from urllib.parse import unquote

root = Path(__file__).resolve().parent.parent
fna = root / '01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/GCF_049532655.1_AFEATCC23270_v5.2_genomic.fna'
seq = ''.join(x.strip() for x in fna.read_text().splitlines() if not x.startswith('>')).upper()
gff = root / 'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/genomic.gff'
genes = []
for line in gff.read_text(encoding='utf-8').splitlines():
    fields = line.split('\t')
    if len(fields) < 9 or fields[2] != 'gene':
        continue
    attrs = dict(part.split('=', 1) for part in fields[8].split(';') if '=' in part)
    genes.append((int(fields[3]), int(fields[4]), unquote(attrs.get('locus_tag', ''))))
primers = {
 'phnG_N':'GAAATTCGATGAAGCGCC', 'phnG_R1':'TCATGACCGCATGGTTTC',
 'phnH_N':'TGATGGCATTGGATTGCC', 'phnI_N':'ATGTCCTACGTGGCAGTC',
 'phnJ_N':'ATAATTTCGCCTATCTCG', 'phnK_N':'GAGTTACAGTTCCGTTCC',
 'phnL_N':'ATGAACCTGCTCGAAGTG', 'phnM_N':'ATGACGCAGGGTTTCAGC',
 'paper_Afe2172_F':'AGGTAATCTTCAGCGGCAAC', 'paper_Afe2172_R':'TAGGGGATCTCCAGACGATG',
}
complement = str.maketrans('ACGT','TGCA')
rows = []
for name, primer in primers.items():
    positions = []
    for strand, query in [('+',primer),('-',primer.translate(complement)[::-1])]:
        start = 0
        while (pos := seq.find(query,start)) != -1:
            positions.append((strand,pos+1))
            start = pos+1
    for strand, position in positions:
        locus = next((locus for start, end, locus in genes if start <= position and position + len(primer) - 1 <= end), '基因外')
        rows.append({'论文表1引物':name,'引物序列5to3':primer,'当前染色体精确命中起点':position,'匹配方向':strand,'当前RU820':locus,'证据边界':'精确短序列命中；应与配对引物、基因排列和论文实验对象共同解读'})
    print(name, positions)
out = root / 'G5X_2929基因逐基因证据链/G57_审查输出/G57_论文引物身份核验.tsv'
with out.open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
