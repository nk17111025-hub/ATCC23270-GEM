from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
GFF = ROOT / '01_原始数据/01_基因组与官方注释/GCF_049532655.1_20260922/解压/ncbi_dataset/data/GCF_049532655.1/genomic.gff'
CROSS = ROOT / 'D3_全基因组候选功能发现/D3-1-R_全基因组重扫与修正/D3-1-R_全基因跨数据库交叉引用.tsv'

def attrs(value: str) -> dict[str, str]:
    return {key: unquote(raw) for key, raw in (item.split('=', 1) for item in value.split(';') if '=' in item)}

def read_gff():
    genes = []
    cds = {}
    with GFF.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cells = line.rstrip('\n').split('\t')
            if len(cells) != 9 or cells[0] != 'NZ_CP136162.1':
                continue
            a = attrs(cells[8])
            if cells[2] == 'gene' and a.get('gene_biotype') == 'protein_coding':
                genes.append((cells, a))
            elif cells[2] == 'CDS':
                cds.setdefault(a.get('locus_tag'), []).append((cells, a))
    assert len(genes) == 2929, len(genes)
    return genes, cds

def write_tsv(path, rows, fieldnames):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def main():
    genes, cds = read_gff()
    with CROSS.open(encoding='utf-8-sig', newline='') as f:
        cross_rows = {r['当前locus']: r for r in csv.DictReader(f, delimiter='\t')}
    group = list(enumerate(genes[879:1172], 880))
    assert len(group) == 293
    assert group[0][1][1]['locus_tag'] == 'RU820_RS04570'
    assert group[-1][1][1]['locus_tag'] == 'RU820_RS06225'
    index = []
    for number, (g, a) in group:
        locus = a['locus_tag']
        c = cds.get(locus, [])
        x = cross_rows.get(locus, {})
        wp = ';'.join(sorted({v['protein_id'] for _, v in c if v.get('protein_id')}))
        product = ';'.join(sorted({v['product'] for _, v in c if v.get('product')}))
        exact = bool(x.get('KEGG reaction') or x.get('BioCyc reaction ID') or x.get('Rhea精确交叉引用') or x.get('2016 GPR反应'))
        mapping_status = x.get('B1映射状态', '')
        mapped_afe = x.get('旧AFE locus', '') if mapping_status == '完全匹配' else '未确认'
        index.append({
            '总序号': number, '当前 RU820 locus': locus, '当前 WP protein_id': wp,
            '染色体坐标': f'{g[0]}:{g[3]}-{g[4]}({g[6]})',
            '当前 NCBI product': product, '当前 gene 名': a.get('gene', ''),
            '旧 RefSeq locus': a.get('old_locus_tag', ''),
            'CDS 条数': len(c),
            'gene/CDS 坐标链向一致': '是' if len(c) == 1 and c[0][0][0] == g[0] and c[0][0][3:5] == g[3:5] and c[0][0][6] == g[6] else '需复核',
            '假基因标记': a.get('pseudo', '') or a.get('pseudogene', ''),
            '旧 AFE locus': mapped_afe,
            'AFE 候选': x.get('旧AFE locus', '') if mapping_status != '完全匹配' else '',
            'AFE 映射依据': x.get('B1映射方法', ''), '映射状态': mapping_status,
            'KEGG KO': x.get('KEGG KO', ''), 'KEGG EC': x.get('KEGG EC', ''),
            'KEGG reaction': x.get('KEGG reaction', ''), 'BioCyc reaction ID': x.get('BioCyc reaction ID', ''),
            'Rhea 交叉引用': x.get('Rhea精确交叉引用', ''), '2016 GPR反应': x.get('2016 GPR反应', ''),
            '是否有反应ID线索': '是' if exact else '否', '是否找到精确 reaction': '待核实' if exact else '未找到',
            '动作': '待逐基因证据裁决', '关联行 ID': f'G54-{number:04d}',
            '审查状态': '身份与本地交叉引用已核；原始论文、蛋白域、化学及2016全字段对照未完成',
            '未决点': '核验 UniProt/InterPro/Pfam/KEGG/Rhea/BioCyc/BRENDA/TCDB、原始论文及反应定义',
            '身份来源': 'GCF_049532655.1 genomic.gff; CDS; 2026-09-22 本地下载',
            '交叉引用来源': 'D3-1-R_全基因跨数据库交叉引用.tsv; 2026-09-23',
        })
    assert len({r['当前 RU820 locus'] for r in index}) == 293
    assert all(re.fullmatch(r'WP_\d+\.\d+(;WP_\d+\.\d+)*', r['当前 WP protein_id']) for r in index)
    wp_counts = Counter(r['当前 WP protein_id'] for r in index)
    for r in index:
        r['本组同 WP 拷贝数'] = wp_counts[r['当前 WP protein_id']]
    write_tsv(OUT / 'G54_基因索引_初核.tsv', index, list(index[0]))
    alpha = ['Reaction ID', 'Reaction Name', 'Reaction Formula', 'Confidence Level', 'EC Number', 'PMID',
             'Subsystem', 'Gene-Reaction Association', 'Gene-Protein-Reaction Association',
             'Protein-Reaction-Association', 'Fe2+ lb', 'Fe2+ ub', 'tetrathionate lb', 'tetrathionate ub',
             'sulfur lb', 'sulfur ub', 'Record class', 'Candidate ID', 'Candidate group',
             '2016 linked Reaction ID', '2026 reaction object', 'Legacy AFE locus', 'Current locus',
             'Primary evidence DOI', 'Evidence type', 'Evidence boundary', '2026 action',
             '2026 reaction score', 'Source ID', 'Score scope']
    review = []
    for r in index:
        row = dict.fromkeys(alpha, '')
        row.update({'Record class': 'gene audit - incomplete', 'Candidate ID': r['关联行 ID'],
                    'Candidate group': 'G54', '2016 linked Reaction ID': r['2016 GPR反应'],
                    'Legacy AFE locus': r['旧 AFE locus'], 'Current locus': r['当前 RU820 locus'],
                    'Evidence type': 'current RefSeq CDS; inherited database cross-reference',
                    'Evidence boundary': '已核当前 GFF/CDS 身份；尚未完成独立蛋白功能、原始论文、精确反应和模型全字段复核。',
                    '2026 action': '待裁决', 'Source ID': 'G54-S001;G54-S002',
                    'Score scope': '未定义精确反应；不评分'})
        review.append(row)
    write_tsv(OUT / 'G54_Alpha结构审查表_初核.tsv', review, alpha)
    sources = [
        {'Source ID':'G54-S001','Full citation':'NCBI RefSeq GCF_049532655.1, AFEATCC23270_v5.2, genomic.gff and CDS','Year':'2026','DOI':'','Source type':'official genome annotation','Strain':'ATCC 23270','Experimental type':'genome annotation','Genes / proteins studied':'G54 293 protein-coding genes','Related Reaction IDs / candidate IDs':'G54-0880–G54-1172','What was directly measured':'No gene-specific enzyme reaction measured','Applicable confidence level':'not scored','Evidence limitation':'Current product is computational annotation; frozen local download 2026-09-22.'},
        {'Source ID':'G54-S002','Full citation':'D3-1-R_全基因跨数据库交叉引用.tsv, local project snapshot','Year':'2026','DOI':'','Source type':'derived cross-reference table','Strain':'ATCC 23270 mapping; mixed reaction databases','Experimental type':'none','Genes / proteins studied':'G54 293 protein-coding genes','Related Reaction IDs / candidate IDs':'G54-0880–G54-1172','What was directly measured':'No new experimental measurement','Applicable confidence level':'not scored','Evidence limitation':'Inherited database links are screening leads; exact chemistry and original records require independent verification.'}
    ]
    write_tsv(OUT / 'G54_来源表_初核.tsv', sources, list(sources[0]))
    print('G54 genes:', len(index), 'with reaction-ID leads:', sum(r['是否有反应ID线索']=='是' for r in index),
          'with 2016 GPR:', sum(bool(r['2016 GPR反应']) for r in index))
    print('mapping:', Counter(r['映射状态'] for r in index))

if __name__ == '__main__':
    main()
