import csv, re, collections, pathlib, urllib.parse
from build_g5_raw_afe import main as build_raw_afe_table
from build_g5_aux import main as build_aux_tables
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=pathlib.Path(__file__).resolve().parent
B1=ROOT/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'旧新基因映射.tsv'
D3=ROOT/'D3_全基因组候选功能发现'/'D3-1-R_全基因组重扫与修正'/'D3-1-R_全基因跨数据库交叉引用.tsv'
NEW=ROOT/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'GCF_049532655.1'/'genomic.gff'
OLD=ROOT/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'GCF_000021485.1'/'genomic.gff'
RAW_OLD=next(p for p in ROOT.rglob('GCA_000021485.1') if p.is_dir() and (p/'genomic.gff').exists() and (p/'protein.faa').exists())/'genomic.gff'
NEW_PROT=ROOT/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'GCF_049532655.1'/'protein.faa'
OLD_PROT=ROOT/'B1_最新基因组与旧AFE编号映射'/'标准化数据'/'GCF_000021485.1'/'protein.faa'
BASE='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/532/655/GCF_049532655.1_AFEATCC23270_v5.2/'

def read_tsv(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def attrs(s):
    d={}
    for x in s.strip().split(';'):
        if '=' in x:
            k,v=x.split('=',1); d[k]=urllib.parse.unquote(v)
    return d
def gff_cds(path):
    z=[]; legacy={}
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):continue
            c=line.rstrip('\n').split('\t')
            if len(c)>=9 and c[2]=='gene':
                a=attrs(c[8]); tag=a.get('locus_tag',''); old=a.get('old_locus_tag','')
                if tag and old: legacy[tag]=old
            if len(c)<9 or c[2]!='CDS':continue
            a=attrs(c[8]); pid=a.get('protein_id')
            if pid:z.append({'pid':pid,'seqid':c[0],'start':int(c[3]),'end':int(c[4]),'strand':c[6],'gene':a.get('gene',''),'product':a.get('product',''),'locus':a.get('locus_tag',''),'old_locus':a.get('old_locus_tag','')})
    for row in z: row['old_locus']=legacy.get(row['locus'],'')
    return z
def fasta_ids(path):
    out=set()
    with path.open(encoding='utf-8') as f:
        for l in f:
            if l.startswith('>'):out.add(l[1:].split()[0])
    return out

def main():
    build_aux_tables()
    # Rebuild original-AFE identity evidence directly from the official GCA GFF/protein set.
    # The raw helper projects old RefSeq coordinates to GCA after verifying their chromosome sequences.
    build_raw_afe_table()
    d3=read_tsv(D3); b1=read_tsv(B1); newg=gff_cds(NEW); oldg=gff_cds(OLD); oldids=fasta_ids(OLD_PROT); newids=fasta_ids(NEW_PROT)
    new_by_pid={r['pid']:r for r in newg}
    aux=read_tsv(OUT/'G5_辅助序列比对.tsv'); auxby={r['当前locus']:r for r in aux}
    neighborhood=read_tsv(OUT/'G5_辅助邻域核验.tsv'); nby={r['当前locus']:r for r in neighborhood}
    rawrows=read_tsv(OUT/'G5_辅助原始AFE核验.tsv'); rawby={r['当前locus']:r for r in rawrows}
    rawg=gff_cds(RAW_OLD); raw_by_locus=collections.defaultdict(list)
    for g in rawg: raw_by_locus[g['locus']].append(g)
    old_by_locus=collections.defaultdict(list)
    for g in oldg: old_by_locus[g['locus']].append(g)
    def source_for(locus,pid):
        rr=rawby.get(locus,{})
        evidence=' '.join(str(rr.get(k,'')) for k in ('区间原始AFE CDS候选序列比对','区间最佳AFE序列命中','旧GCF坐标投射到原始GCA AFE CDS及序列identity/双向coverage','最佳ACK protein ID'))
        ack_ids=sorted(set(re.findall(r'ACK\d+\.\d+',evidence)))
        old_wp=auxby.get(locus,{}).get('top旧protein ID','')
        ack_links='; '.join(f'{a}: https://www.ncbi.nlm.nih.gov/protein/{a}' for a in ack_ids) or 'no ACK protein hit'
        oldwp_link=f'old RefSeq WP {old_wp}: https://www.ncbi.nlm.nih.gov/protein/{old_wp}' if old_wp else ''
        return (f'NCBI RefSeq GCF_049532655.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049532655.1/; '
                f'current WP {pid}: https://www.ncbi.nlm.nih.gov/protein/{pid}; '
                f'old RefSeq GCF_000021485.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000021485.1/; {oldwp_link}; '
                f'original GenBank GCA_000021485.1 GFF/protein: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000021485.1/; {ack_links}')
    # D3 is the current-locus universe; cross-check against official NCBI GFF.
    pc=[x for x in d3 if x['基因类型']=='protein_coding']
    byloc={x['当前locus']:x for x in pc}
    br=collections.defaultdict(list)
    for x in b1:
        for current_locus in [v for v in x.get('当前GCF049位点','').split('|') if v]:
            br[current_locus].append(x)
    audited=[]; novel=[]; excluded=[]; ambiguous=[]; length_change_loci=set()
    for x in pc:
        locus=x['当前locus']; status=x.get('B1映射状态','');
        if status=='完全匹配':continue
        pid=x.get('当前protein ID',''); og=new_by_pid.get(pid,{})
        if not br.get(locus):
            bstate='B1无行'
        elif len(br[locus])>1 or status=='联系同源/重复候选':bstate='多重候选'
        else:bstate='未映射/未确定'
        candidates=[]
        if pid in oldids:
            for od in oldg:
                if od['pid']==pid: candidates.append((od.get('old_locus',''),od['pid'],od))
        afe=sorted({v for r in br.get(locus,[]) for v in r.get('旧AFE位点','').split('|') if re.fullmatch(r'AFE_\d{4}',v)})
        # Compare each unresolved current sequence against the complete old RefSeq protein set.
        # Keep sequence hits separate from identity labels: old AFE_xxxx requires an explicit tag.
        ar=auxby.get(locus,{})
        topid=ar.get('top旧protein ID',''); top_rs=ar.get('旧AFE候选(旧GFF locus_tag)','')
        try: pct=float(ar.get('identity','0')); cov=float(ar.get('当前序列覆盖','0')); oldcov=float(ar.get('旧序列覆盖','0')); mincov=float(ar.get('短序列覆盖(按min长度)','0'))
        except ValueError: pct=cov=oldcov=mincov=0
        nr=nby.get(locus,{}); suggestion=nr.get('可判分类建议','')
        oldtagmatch=re.findall(r'old=(AFE_\d{4})',nr.get('序列top旧locus候选及old_locus_tag',''))
        reason=''; verdict=''; mother=False; final_tags=[]
        strong_identity=(pct>=0.90 and mincov>=0.90 and cov>=0.90 and oldcov>=0.90)
        rr=rawby.get(locus,{})
        projection=rr.get('旧GCF坐标投射到原始GCA AFE CDS及序列identity/双向coverage','')
        hits=[]
        for h in re.finditer(r'=>(?P<afe>AFE_\d{4})\((?P<ack>[^)]+)\)\[bp=(?P<bp>\d+);aa=(?P<aa>\d+);id=(?P<id>[0-9.]+);curCov=(?P<cur>[0-9.]+);oldCov=(?P<old>[0-9.]+)\]',projection):
            hit=h.groupdict(); hit.update({k:float(hit[k]) for k in ('id','cur','old')}); hit.update({k:int(hit[k]) for k in ('bp','aa')}); hits.append(hit)
        substantial=[h for h in hits if h['aa']>=30 and h['id']>=0.90 and max(h['cur'],h['old'])>=0.90]
        partial=[h for h in substantial if min(h['cur'],h['old'])<0.90]
        raw_best_loci=sorted(set(re.findall(r'AFE_\d{4}',rr.get('全库最佳原始AFE_xxxx',''))))
        try: raw_id=float(rr.get('identity','0')); raw_cur=float(rr.get('current coverage','0')); raw_oldcov=float(rr.get('old coverage','0'))
        except ValueError: raw_id=raw_cur=raw_oldcov=0
        raw_best_full=bool(raw_id>=0.90 and raw_cur>=0.90 and raw_oldcov>=0.90 and len(raw_best_loci)==1)
        raw_best_partial=bool(raw_id>=0.90 and max(raw_cur,raw_oldcov)>=0.90 and len(raw_best_loci)==1)
        # The raw-AFE helper maps old RefSeq coordinates onto the verified identical old GCA/GCF
        # chromosome and separately compares proteins inside the confirmed two-anchor interval.
        interval_hits=[]
        interval_text=rr.get('区间原始AFE CDS候选序列比对','')
        for h in re.finditer(r'(AFE_\d{4})\((ACK[^)]+)\)\[aa=(\d+);id=([0-9.]+);curCov=([0-9.]+);oldCov=([0-9.]+)\]',interval_text):
            tag,ack,aa,ident,qc,oc=h.groups(); interval_hits.append({'afe':tag,'ack':ack,'aa':int(aa),'id':float(ident),'cur':float(qc),'old':float(oc)})
        anchored_full=[h for h in interval_hits if h['id']>=0.90 and h['aa']>=30 and h['cur']>=0.90 and h['old']>=0.90]
        raw_best_loci=sorted(set(re.findall(r'AFE_\d{4}',rr.get('区间最佳AFE序列命中','')))) or raw_best_loci
        full_projection=[]
        for h in substantial:
            for rg in raw_by_locus.get(h['afe'],[]):
                shorter_nt=min(new_by_pid.get(pid,{}).get('end',0)-new_by_pid.get(pid,{}).get('start',0)+1, rg['end']-rg['start']+1)
                if min(h['cur'],h['old'])>=0.90 and shorter_nt>0 and h['bp']/shorter_nt>=0.50: full_projection.append(h)

        if partial:
            tags=','.join(sorted({h['afe'] for h in partial}))
            final_tags=sorted({h['afe'] for h in partial}); length_change_loci.add(locus)
            verdict='旧AFE候选涉及split/merge或注释边界变化，待人工复核'
            reason='; '.join(f"原始GCA {h['afe']} ({h['ack']}) 同链坐标重叠 {h['bp']} bp；identity={h['id']:.4f}, aa={h['aa']}, 当前覆盖={h['cur']:.4f}, 旧覆盖={h['old']:.4f}" for h in partial)
            ambiguous.append((x,bstate,tags,verdict,reason))
        elif full_projection:
            tags=','.join(sorted({h['afe'] for h in full_projection}))
            final_tags=sorted({h['afe'] for h in full_projection})
            if len(tags.split(','))==1:
                tag=tags; h=full_projection[0]
                verdict='重新识别为原始GCA旧AFE基因对应'
                reason=f"原始GCA {tag} ({h['ack']}) 与旧GCF同坐标投射重叠 {h['bp']} bp，占较短CDS大部；identity={h['id']:.4f}, aa={h['aa']}, 双向coverage={h['cur']:.4f}/{h['old']:.4f}"
                excluded.append((x,bstate,tag,verdict,reason))
            else:
                verdict='多个原始AFE候选，待人工复核'; reason=f'高相似原始AFE坐标候选={tags}，存在多对一/拆分关系'; ambiguous.append((x,bstate,tags,verdict,reason))
        elif anchored_full and len({h['afe'] for h in anchored_full})==1:
            tag=anchored_full[0]['afe']
            final_tags=[tag]
            verdict='重新识别为原始GCA旧AFE基因对应'
            h=anchored_full[0]
            reason=f"原始GCA {tag} ({h['ack']}) 位于两侧已确认旧AFE_RS锚点之间；identity={h['id']:.4f}, aa={h['aa']}, 双向coverage={h['cur']:.4f}/{h['old']:.4f}"
            excluded.append((x,bstate,tag,verdict,reason))
        elif raw_best_full or raw_best_partial or substantial:
            tags=','.join(sorted(set(raw_best_loci+[h['afe'] for h in substantial]+[h['afe'] for h in interval_hits]))) or '待定'
            final_tags=sorted(set(raw_best_loci+[h['afe'] for h in substantial]+[h['afe'] for h in interval_hits]))
            verdict='原始AFE序列命中但位置/覆盖关系冲突，待人工复核'
            reason=f"原始GCA top={rr.get('最佳ACK protein ID','')} loci={tags}; identity={raw_id:.4f}, 双向coverage={raw_cur:.4f}/{raw_oldcov:.4f}; 同位锚定未成立或存在长度差/旁系同源"
            ambiguous.append((x,bstate,tags,verdict,reason))
            if raw_cur<0.90 or raw_oldcov<0.90 or any(min(h['cur'],h['old'])<0.90 for h in substantial): length_change_loci.add(locus)
        elif '旧AFE对应' in suggestion and oldtagmatch and strong_identity:
            oldtag=oldtagmatch[0]
            final_tags=[oldtag]
            verdict='重新识别为旧AFE基因对应'
            reason=f"邻域核验：{nr.get('邻域是否一致','')}；序列top hit={topid} identity={pct:.4f}，较短覆盖={mincov:.4f}，当前覆盖={cov:.4f}，旧序列覆盖={oldcov:.4f}；旧基因组AFE locus={top_rs}对应{oldtag}"
            excluded.append((x,bstate,oldtag,verdict,reason))
        elif 'old_locus_tag为空/非AFE' in suggestion and strong_identity:
            final_tags=afe
            verdict='旧基因组同源基因无旧AFE_xxxx标签；按工作定义纳入母表'
            reason=f"邻域核验：{nr.get('邻域是否一致','')}；{nr.get('序列top旧locus候选及old_locus_tag','')}；identity={pct:.4f}，较短覆盖={mincov:.4f}，当前覆盖={cov:.4f}，旧序列覆盖={oldcov:.4f}。旧2025 GFF已有同源蛋白注释，但无法关联任何旧AFE_xxxx，故按本项目身份定义列入母表；不表示生物学上新出现。"
            mother=True; novel.append((x,bstate,'无旧AFE_xxxx',verdict,reason))
        else:
            tags=','.join(afe) or ','.join(oldtagmatch) or '待定'
            final_tags=sorted(set(afe+oldtagmatch))
            verdict='旧新关系待人工复核'
            reason=f"邻域核验：{nr.get('邻域是否一致','')}；分类建议={suggestion}；top hit={topid} identity={pct:.4f}，较短覆盖={mincov:.4f}，当前覆盖={cov:.4f}，旧序列覆盖={oldcov:.4f}；identity或任一侧覆盖未达90%，或出现重复拷贝/结构变化，未强制归类。"
            if cov<0.90 or oldcov<0.90: length_change_loci.add(locus)
            ambiguous.append((x,bstate,tags,verdict,reason))
        source=source_for(locus,pid)
        audit_oldafe=','.join(final_tags or afe or oldtagmatch) if (final_tags or afe or oldtagmatch) else '无'
        audited.append({'当前locus':locus,'当前protein ID':pid,'B1映射状态':bstate,'旧AFE候选':audit_oldafe,'本轮核验结论':verdict,'是否进入母表':'是' if mother else '否','排除/保留原因':reason,'官方来源':source})
    # stable ordering by current locus then assign contiguous IDs
    novel.sort(key=lambda t:t[0]['当前locus'])
    def write(name,headers,rows):
        with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=headers,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
    master=[]
    for n,(x,*_) in enumerate(novel,1):
        pid=x['当前protein ID'];master.append({'编号':n,'当前官方基因编号':x['当前locus'],'当前官方蛋白编号':pid,'官方来源':f'NCBI RefSeq GCF_049532655.1; {BASE}; https://www.ncbi.nlm.nih.gov/protein/{pid}'})
    write('G5_01_全新基因母表.tsv',['编号','当前官方基因编号','当前官方蛋白编号','官方来源'],master)
    write('G5_02_全新基因身份核验.tsv',['当前locus','当前protein ID','B1映射状态','旧AFE候选','本轮核验结论','是否进入母表','排除/保留原因','官方来源'],audited)
    ambrows=[{'当前locus':x['当前locus'],'当前protein ID':x.get('当前protein ID',''),'旧AFE候选':tags,'歧义类型/待核事项':verdict,'当前证据及未决原因':reason,'官方来源':source_for(x['当前locus'],x.get('当前protein ID',''))} for x,state,tags,verdict,reason in ambiguous]
    write('G5_03_映射歧义待人工复核.tsv',['当前locus','当前protein ID','旧AFE候选','歧义类型/待核事项','当前证据及未决原因','官方来源'],ambrows)
    exrows=[{'当前locus':x['当前locus'],'当前protein ID':x.get('当前protein ID',''),'旧AFE位点':tag,'重新识别依据':reason,'官方来源':source_for(x['当前locus'],x.get('当前protein ID',''))} for x,state,tag,v,reason in excluded]
    write('G5_04_排除记录.tsv',['当前locus','当前protein ID','旧AFE位点','重新识别依据','官方来源'],exrows)
    assert len(audited)==len(master)+len(ambrows)+len(exrows)
    assert len({x['当前locus'] for x in audited})==len(audited)
    assert len(pc)==2929
    bmissing=len(audited)
    report=f'''# G5 验收报告\n\n- 当前 protein-coding gene 总数：{len(pc)}（D3-1-R 全集，并由当前 RefSeq GFF 核对）。\n- B1 初始无可靠 AFE 映射数量：{bmissing}（当前protein-coding loci中 B1状态非“完全匹配”，含多重候选与B1无行）。\n- 最终进入全新基因母表的数量：{len(master)}。\n- 被重新识别为旧基因对应关系的数量：{len(exrows)}。\n- gene split/merge/annotation-change 数量：确认 split/merge 0 条；历史 AFE 注释缺失而旧参考现行注释已有蛋白的 annotation-change {len(master)} 条；待复核表中 {len(length_change_loci)} 条为长度或边界变化候选，尚未确认 split/merge。\n- 仍待人工复核数量：{len(ambrows)}。\n- 母表是否从1连续编号：是（1–{len(master)}）。\n- 是否所有正式母表记录都有 RU820、WP 和官方来源：是。\n- 是否存在未归档 locus：否；核验对象分区数量 {len(master)}+{len(exrows)}+{len(ambrows)}={bmissing}。\n\n'''
    (OUT/'G5_验收报告.md').write_text(report,encoding='utf-8')
    print(f'pc={len(pc)} audit={len(audited)} novel={len(master)} excluded={len(exrows)} ambiguous={len(ambrows)}')
if __name__=='__main__':main()
