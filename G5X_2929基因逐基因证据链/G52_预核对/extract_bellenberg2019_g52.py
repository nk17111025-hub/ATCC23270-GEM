import json
import re
from pathlib import Path
from pypdf import PdfReader

DIR=Path(__file__).resolve().parent
records=json.load((DIR/'g52_precheck.json').open(encoding='utf-8'))['records']
by_afe={r['旧AFE']:r for r in records if r['旧AFE']}
pdf=PdfReader(DIR/'Bellenberg2019_Data_Sheet_2.pdf')
out=[]
for page_index,page in enumerate(pdf.pages):
    text=page.extract_text()
    hits=list(re.finditer(r'AFE_\d{4}',text))
    for i,hit in enumerate(hits):
        afe=hit.group()
        if afe not in by_afe:
            continue
        stop=hits[i+1].start() if i+1<len(hits) else min(len(text),hit.end()+500)
        snippet=' '.join(text[hit.start():stop].split())
        record=by_afe[afe]
        out.append({'afe':afe,'locus':record['当前locus'],'wp':record['WP'],'page':page_index+1,'wp_exact':record['WP'] in snippet,'snippet':snippet[:500]})
with (DIR/'g52_bellenberg2019_hits.json').open('w',encoding='utf-8') as f:
    json.dump(out,f,ensure_ascii=False,indent=2)
print('rows',len(out),'unique_loci',len({x['locus'] for x in out}),'WP_exact',sum(x['wp_exact'] for x in out))
for x in out:
    print(x['page'],x['afe'],x['locus'],x['wp_exact'],x['snippet'][:150])
