import collections
import pathlib
import urllib.request
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from io import StringIO

root = pathlib.Path(__file__).parents[2]
pdb_fasta = urllib.request.urlopen('https://www.rcsb.org/fasta/entry/3KPK/display', timeout=30).read().decode()
pdb = next(SeqIO.parse(StringIO(pdb_fasta), 'fasta'))
print('PDB 3KPK', len(pdb.seq), pdb.description)

files = [
    root / 'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_049532655.1/protein.faa',
    root / 'B1_最新基因组与旧AFE编号映射/标准化数据/GCF_000021485.1/protein.faa',
]
seq = str(pdb.seq).replace('X','')
kmers = {seq[i:i+8] for i in range(len(seq)-7)}
aligner = PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -6
aligner.extend_gap_score = -0.5
for f in files:
    hits = []
    for record in SeqIO.parse(f, 'fasta'):
        s = str(record.seq)
        overlap = sum(s[i:i+8] in kmers for i in range(len(s)-7))
        if overlap:
            hits.append((overlap, record))
    print('\n', f.name, f.parent.name, 'sequences with 8mer hits', len(hits))
    for overlap, record in sorted(hits, key=lambda x:x[0], reverse=True)[:8]:
        alignment = aligner.align(str(pdb.seq), str(record.seq))[0]
        print(overlap, round(alignment.score, 1), len(record.seq), record.id, record.description[:150])
