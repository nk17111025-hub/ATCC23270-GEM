import json, os, zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

root=Path(__file__).parent
data=json.loads((root/'evidence_data.json').read_text(encoding='utf-8'))
file=root/'outputs/01a0d100-c9f0-7122-903b-f8175e453668/ATCC23270-2026_最终模型证据表.xlsx'
tmp=file.with_suffix('.tmp.xlsx')
assert file.resolve().is_relative_to(root.resolve())

main='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
rel='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
pkg='http://schemas.openxmlformats.org/package/2006/relationships'
ET.register_namespace('',main)
ET.register_namespace('r',rel)

links1=[(f'X{618+i}','https://doi.org/'+x['doi']) for i,x in enumerate(data['candidates']) if x['doi']]
links2=[(f'D{2+i}','https://doi.org/'+x['doi']) for i,x in enumerate(data['sources'])]
replacements={}
with zipfile.ZipFile(file,'r') as zin:
    for sheetnum,links in [(1,links1),(2,links2)]:
        sheetpath=f'xl/worksheets/sheet{sheetnum}.xml'
        relpath=f'xl/worksheets/_rels/sheet{sheetnum}.xml.rels'
        ws=ET.fromstring(zin.read(sheetpath))
        relroot=ET.fromstring(zin.read(relpath)) if relpath in zin.namelist() else ET.Element(f'{{{pkg}}}Relationships')
        existing={int(e.attrib['Id'][3:]) for e in relroot if e.attrib.get('Id','').startswith('rId') and e.attrib['Id'][3:].isdigit()}
        start=max(existing,default=0)+1
        hyperlinks=ws.find(f'{{{main}}}hyperlinks')
        if hyperlinks is None:
            hyperlinks=ET.Element(f'{{{main}}}hyperlinks')
            # OOXML puts hyperlinks after conditional formatting and data validations, before print settings.
            before={'printOptions','pageMargins','pageSetup','headerFooter','rowBreaks','colBreaks','drawing','legacyDrawing','tableParts','extLst'}
            pos=next((i for i,e in enumerate(ws) if e.tag.rsplit('}',1)[-1] in before),len(ws))
            ws.insert(pos,hyperlinks)
        for i,(cell,url) in enumerate(links):
            rid=f'rId{start+i}'
            ET.SubElement(hyperlinks,f'{{{main}}}hyperlink',{'ref':cell,f'{{{rel}}}id':rid})
            ET.SubElement(relroot,f'{{{pkg}}}Relationship',{'Id':rid,'Type':'http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink','Target':url,'TargetMode':'External'})
        replacements[sheetpath]=ET.tostring(ws,encoding='utf-8',xml_declaration=True)
        replacements[relpath]=ET.tostring(relroot,encoding='utf-8',xml_declaration=True)
    with zipfile.ZipFile(tmp,'w') as zout:
        for info in zin.infolist():
            zout.writestr(info,replacements.pop(info.filename,zin.read(info.filename)))
        for name,blob in replacements.items():
            zout.writestr(name,blob,compress_type=zipfile.ZIP_DEFLATED)
os.replace(tmp,file)
print('native links',len(links1),len(links2),'file',file)
