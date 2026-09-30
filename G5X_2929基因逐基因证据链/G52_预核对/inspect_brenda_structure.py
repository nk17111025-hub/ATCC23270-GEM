import urllib.request
from lxml import html

url='https://www.brenda-enzymes.org/enzyme.php?ecno=1.8.1.9&onlyTable=Reference&showtm=0'
with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'G52 evidence review/1.0'}),timeout=50) as r:
    body=r.read()
tree=html.fromstring(body)
print('body',len(body),'tr',len(tree.xpath('//tr')))
print('all div',len(tree.xpath('//div')),'root',tree.tag)
print('row-tab30',len(tree.xpath('//div[contains(@class, "row tab30")]')))
print('matching-row-tab30',len([x for x in tree.xpath('//div[contains(@class, "row tab30")]') if 'Acidithiobacillus ferrooxidans' in ' '.join(x.itertext())]))
for el in tree.iter():
    if el.text and 'Acidithiobacillus ferrooxidans' in el.text:
        print('tag',el.tag,'text',el.text[:150])
        for p in list(el.iterancestors())[:5]:
            print(' parent',repr(p.tag),repr(p.get('class')),p.get('id'),' '.join(' '.join(p.itertext()).split())[:250])
        break
