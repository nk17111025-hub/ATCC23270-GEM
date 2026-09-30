import json
from collections import Counter
from pathlib import Path
from openpyxl import load_workbook

DIR=Path(__file__).resolve().parent
source=json.loads((DIR/'g52_precheck.json').read_text(encoding='utf-8'))
decisions=json.loads((DIR/'g52_gene_decisions_draft.json').read_text(encoding='utf-8'))
candidates=json.loads((DIR/'g52_candidate_rows_draft.json').read_text(encoding='utf-8'))
balance=json.loads((DIR/'g52_2016_balance.json').read_text(encoding='utf-8'))
chemistry=json.loads((DIR/'g52_2016_kegg_chemistry.json').read_text(encoding='utf-8'))
manual=json.loads((DIR/'g52_manual_adjudications.json').read_text(encoding='utf-8'))
lit=json.loads((DIR/'g52_literature_search.json').read_text(encoding='utf-8'))['by_afe']
wp_lit=json.loads((DIR/'g52_wp_literature_search.json').read_text(encoding='utf-8'))['by_wp']
function_lit=json.loads((DIR/'g52_function_literature_search.json').read_text(encoding='utf-8'))['by_locus']
tcdb_all=json.loads((DIR/'g52_tcdb_homology.json').read_text(encoding='utf-8'))['by_locus']
tcdb_family=json.loads((DIR/'g52_tcdb_family_validation.json').read_text(encoding='utf-8'))['by_locus']
trdr_identity=json.loads((DIR/'g52_trxr_2009_exact_identity.json').read_text(encoding='utf-8'))
trdr_balance=json.loads((DIR/'g52_trdr_corrected_balance.json').read_text(encoding='utf-8'))
rhea=json.loads((DIR/'g52_rhea_reaction_entries.json').read_text(encoding='utf-8'))
rhea_kegg=json.loads((DIR/'g52_rhea_kegg_xref_comparison.json').read_text(encoding='utf-8'))
task_doc_match=json.loads((DIR/'g52_task_doc_match.json').read_text(encoding='utf-8'))
brenda=json.loads((DIR/'g52_brenda_ec_references.json').read_text(encoding='utf-8'))['by_ec']
atps=json.loads((DIR/'g52_atps_2012_sequence_comparison.json').read_text(encoding='utf-8'))
rhodanese=json.loads((DIR/'g52_rhodanese_2005_sequence_mapping.json').read_text(encoding='utf-8'))
eps2005=json.loads((DIR/'g52_eps_2005_sequence_mapping.json').read_text(encoding='utf-8'))
book=load_workbook(DIR/'G52_全量证据重查_Alpha最终格式.xlsx',read_only=True,data_only=True)

assert len(source['records'])==len(decisions)==293
assert [r['总序号'] for r in source['records']]==list(range(294,587))
assert [r['总序号'] for r in decisions]==list(range(294,587))
assert len({r['当前locus'] for r in decisions})==293
assert len({r['WP'] for r in decisions})==293
assert all(r['2026建议'] and r['来源ID'] and r['身份依据'] for r in decisions)
assert set(manual)<={r['当前locus'] for r in source['records']}
assert len(lit)==279 and len(wp_lit)==293
assert len(function_lit)==293 and not any('error' in v for v in function_lit.values())
assert not any('error' in v for v in lit.values())
assert not any('error' in v for v in wp_lit.values())
assert len(source['reaction_rows'])==90
assert len(candidates)==20 and len({r['candidate_id'] for r in candidates})==20
assert len(chemistry)==90
assert len(balance)==87
assert Counter(r['status'] for r in balance)=={'balanced':78,'unbalanced':4,'unknown':5}
assert len(list(book['G52基因索引'].values))==294
assert len(list(book['Alpha最终证据表'].values))==404
assert len(list(book['逐基因证据链'].values))==294
assert len(list(book['来源与限制'].values))==32
assert len(tcdb_all)==25 and len(tcdb_family)==10
assert sum(v['family_homology_supported'] for v in tcdb_family.values())==9
assert all(len(row)==35 for row in list(book['逐基因证据链'].values))
assert all(row[42]=='G52证据重查完成；未证实项按当前建议Hold' for row in list(book['G52基因索引'].values)[1:])
assert all(row[25]=='G52证据重查完成；未证实项Hold' for row in list(book['逐基因证据链'].values)[1:])
assert trdr_identity['both_primer_ends_match_current_cds'] and trdr_balance['balanced']
assert rhea['matched']==144 and not rhea['missing']
assert rhea_kegg['summary']=={'both':74,'overlap':70,'different':4}
assert task_doc_match['mismatches']==0 and task_doc_match['document_rows']==293
assert len(brenda)==119 and sum(bool(v['count']) for v in brenda.values())==7
assert atps['current_wp']=='WP_012536192.1' and atps['old_and_current_cds_length_nt']==1674
assert len(atps['nt_differences'])==7 and len(atps['aa_differences'])==4
p162=rhodanese['AY863108']
assert p162['current_locus']=='RU820_RS02245' and p162['current_wp']=='WP_009564831.1'
assert p162['all_other_nt_identical'] and p162['single_G_insertion_in_historical_allele_position_1_based']==[361,362,363]
assert not rhodanese['AY863107']['top_hits'][0]['in_g52']
assert eps2005['AY789511']['exact_current_cds_match']
assert eps2005['AY789511']['current_locus']=='RU820_RS02155'
assert not eps2005['AY789512']['top_hits'][0]['g52_locus']

alpha_rows=list(book['Alpha最终证据表'].values)
alpha_header=alpha_rows[0]
gene_alpha=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[16]=='逐基因证据链']
assert len(gene_alpha)==293
assert {r['Current locus'] for r in gene_alpha}=={r['当前locus'] for r in decisions}
assert all(r['2026 action'] and r['Source ID'] and not r['2026 reaction score'] for r in gene_alpha)
assert all(not any(r[k] for k in alpha_header[:16]) for r in gene_alpha)
historical_by_id={r['Candidate ID']:r for r in source['reaction_rows']}
historical_seen=set()
old_core_keys=['Reaction ID','Reaction Name','Reaction Formula','Confidence Level','EC Number','PMID','Subsystem','Gene-Reaction Association','Gene-Protein-Reaction Association','Protein-Reaction-Association','Fe2 lb','Fe2 ub','tetrathionate lb','tetrathionate ub','sulfur lb','sulfur ub']
for row in alpha_rows[1:]:
    row=dict(zip(alpha_header,row))
    candidate_id=row['Candidate ID']
    if candidate_id not in historical_by_id:
        continue
    expected=historical_by_id[candidate_id]
    assert [row[k] for k in alpha_header[:16]]==[expected[k] for k in old_core_keys],candidate_id
    historical_seen.add(candidate_id)
assert historical_seen==set(historical_by_id)
trdr_rows=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[0]=='TRDR' and row[22]=='RU820_RS01820']
assert len(trdr_rows)==1 and trdr_rows[0]['2026 reaction score']==4
assert trdr_rows[0]['2026 reaction object']==trdr_balance['candidate_formula']
atps_rows=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[0]=='SFAT1' and row[22]=='RU820_RS02615']
assert len(atps_rows)==1 and atps_rows[0]['2026 reaction score']==4
assert atps_rows[0]['Confidence Level']==2
assert atps_rows[0]['2026 reaction object']=='aps[c] + ppi[c] -> atp[c] + so4[c]'
p162_rows=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[17]=='G52-NEW-RU820_RS02245']
assert len(p162_rows)==1 and p162_rows[0]['2026 reaction score']==4
assert p162_rows[0]['2026 reaction object']=='thiosulfate + cyanide -> sulfite + thiocyanate'
galu_rows=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[0]=='GALU' and row[22]=='RU820_RS02155']
assert len(galu_rows)==1 and galu_rows[0]['2026 reaction score']==4
assert galu_rows[0]['Confidence Level']==2
assert galu_rows[0]['2026 reaction object']=='utp[c] + g1p[c] -> udpg[c] + ppi[c]'
scored=[dict(zip(alpha_header,row)) for row in alpha_rows[1:] if row[27] is not None]
assert {r['Candidate ID'] for r in scored}=={
    'G52-2016-TRDR-RU820_RS01820',
    'G52-2016-SFAT1-RU820_RS02615',
    'G52-2016-GALU-RU820_RS02155',
    'G52-NEW-RU820_RS02245',
}
assert all(r['2026 reaction score']==4 for r in scored)
assert not any('仍需核对' in str(r['Evidence boundary']) for r in map(lambda row:dict(zip(alpha_header,row)),alpha_rows[1:91]))
assert not any('复核后落表' in r['2026建议'] for r in decisions)

sources={r[0] for r in list(book['来源与限制'].values)[1:]}
for row in alpha_rows[1:]:
    ids=row[28]
    assert ids and set(ids.split(';'))<=sources,row[17]
for row in decisions:
    assert set(row['来源ID'].split(';'))<=sources,(row['当前locus'],row['来源ID'])
for row in candidates:
    assert set(row['sources'].split(';'))<=sources,row['candidate_id']
for row in manual.values():
    assert set(row['evidence'].split(';'))<=sources,row['evidence']

report={
    'identity_293':True,
    'ordinal_294_586_contiguous':True,
    'task_document_293_exact_match':True,
    'all_current_wp_unique':True,
    'confirmed_afe':sum(bool(r['旧AFE']) for r in source['records']),
    'unconfirmed_afe':sum(not r['旧AFE'] for r in source['records']),
    'gene_decision_rows':len(decisions),
    'alpha_schema_gene_evidence_rows':len(gene_alpha),
    'manual_issue_reviews':len(manual),
    'historical_gpr_rows':90,
    'historical_alpha_A_to_P_preserved':True,
    'new_candidate_rows':len(candidates),
    'unique_2016_reactions_balance':dict(Counter(r['status'] for r in balance)),
    'exact_kegg_compound_set_matches':sum(bool(r['best']) and r['best'][0]['similarity']>=.95 and not r['model_missing_kegg_ids'] for r in chemistry),
    'source_records':len(sources),
    'function_phrase_searches':len(function_lit),
    'function_phrase_index_hits':sum(bool(v['hit_count']) for v in function_lit.values()),
    'tcdb_queries':len(tcdb_all),
    'tcdb_family_checks':len(tcdb_family),
    'tcdb_family_homology_leads':sum(v['family_homology_supported'] for v in tcdb_family.values()),
    'trdr_same_cds_direct_biochemistry':True,
    'trdr_corrected_formula_balanced':True,
    'rhea_official_reactions_recovered':rhea['matched'],
    'rhea_kegg_shared_xrefs':rhea_kegg['summary']['overlap'],
    'brenda_ec_queried':len(brenda),
    'brenda_species_ec_reference_leads':sum(bool(v['count']) for v in brenda.values()),
    'atps_2012_same_strain_near_allele_assay':True,
    'p162_2005_same_strain_near_allele_assay':True,
    'galu_2005_exact_cds_genetic_and_biochemical_evidence':True,
    'reaction_scored_with_original_assay':len(scored),
    'final_integrity_passed':True,
    'final_evidence_chain_complete':True,
    'unresolved_evidence_handling':'Unconfirmed exact catalysis, substrate, direction or compartment is retained as historical candidate or Hold; no unverified 2026 reaction score.'
}
(DIR/'g52_selfcheck_final.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False))
