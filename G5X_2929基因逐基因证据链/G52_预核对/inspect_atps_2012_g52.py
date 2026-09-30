import urllib.request
from pathlib import Path
from lxml import etree

DIR=Path(__file__).resolve().parent
path=DIR/'Jaramillo2012_PMC3449377.xml'
if not path.exists():
    with urllib.request.urlopen('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3449377/fullTextXML',timeout=40) as r:
        path.write_bytes(r.read())
root=etree.parse(str(path))
for p in root.xpath('//p'):
    t=' '.join(' '.join(p.itertext()).split())
    if any(s in t.lower() for s in ['primer','1674 bp','activity','2327613','atps gene']):
        print(t[:1800], '\n---')
