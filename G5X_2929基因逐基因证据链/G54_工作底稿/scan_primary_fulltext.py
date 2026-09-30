import csv
import pathlib
import re
import urllib.request
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).parent
rows = list(csv.DictReader((ROOT / 'G54_基因索引_初核.tsv').open(encoding='utf-8-sig'), delimiter='\t'))
by_afe = {r['旧 AFE locus']: r for r in rows if re.fullmatch(r'AFE_\d{4}', r['旧 AFE locus'])}
papers = [
    ('Q2009', 'PMC2754497', '10.1186/1471-2164-10-394'),
    ('B2019', 'PMC6450195', '10.3389/fmicb.2019.00592'),
]
out = []
for source, pmcid, doi in papers:
    url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'
    tree = ET.fromstring(urllib.request.urlopen(url, timeout=40).read())
    body = tree.find('body')
    if body is None:
        raise ValueError(f'No body in {pmcid}')
    for node in body.iter():
        if node.tag not in {'p', 'caption', 'table-wrap'}:
            continue
        text = ' '.join(''.join(node.itertext()).split())
        for afe in sorted(set(re.findall(r'AFE_\d{4}', text)) & by_afe.keys()):
            r = by_afe[afe]
            out.append({
                '总序号': r['总序号'],
                'RU820': r['当前 RU820 locus'],
                'WP': r['当前 WP protein_id'],
                '旧 AFE': afe,
                '来源 ID': source,
                'DOI': doi,
                'PMCID': pmcid,
                '正文命中段落': text[:1800],
                '证据边界': '全文文本命中；需人工判读实验对象与反应，不自动算直接功能证据',
            })
with (ROOT / 'G54_两篇同株全文AFE命中.tsv').open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(out)
print(f'{len(out)} passages, {len(set((x["总序号"], x["来源 ID"]) for x in out))} gene-paper pairs, {len(set(x["总序号"] for x in out))} genes')
