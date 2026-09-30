# E1 donor reaction pool

- Host: `D:\嗜酸氧化亚铁硫杆菌\before final\updatedv3.1_verified\updatedv3.1.xml`
- Observed host SHA256: `026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3` (expected hash matched)
- Host reactions/metabolites: 620/575
- Donor files and sources: iML1515: `D:\嗜酸氧化亚铁硫杆菌\ATCC23270-GEM\external\reference\donor_models\iML1515.json` (SHA256 `b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5`; Reused local file; reference page https://bigg.ucsd.edu/models/iML1515); iJN1463: `D:\嗜酸氧化亚铁硫杆菌\ATCC23270-GEM\external\reference\donor_models\iJN1463.xml` (SHA256 `d573833328ffae0dfa752a1fa3262ed939ed5862288beab287fca30d0fefb4a1`; Downloaded source https://bigg.ucsd.edu/static/models/iJN1463.xml); iCN1361: `D:\嗜酸氧化亚铁硫杆菌\ATCC23270-GEM\external\reference\donor_models\iCN1361.xml` (SHA256 `c4e418d831f428a93702f2d1b6c62d6fe3538545474f12397640028169d4ce0c`; Reused local file; author model repository https://github.com/SBRCNottingham/CnecatorGSM/tree/main/Model)

## Counts

| Donor | Read | Pseudo/biomass removed | Cross-compartment transport | Single non-cytosolic | Intracellular checked | Host equivalent rows | Direction difference source rows |
|---|---:|---:|---:|---:|---:|---:|---:|
| iML1515 | 2712 | 339 | 831 | 203 | 1339 | 322 | 9 |
| iJN1463 | 2927 | 383 | 821 | 149 | 1574 | 319 | 8 |
| iCN1361 | 1292 | 99 | 14 | 0 | 1179 | 263 | 60 |

- Unique donor reactions after donor-to-donor deduplication and host-equivalence removal: 2424
- Intracellular donor source rows before deduplication: 4092
- Unique donor reactions after donor-to-donor deduplication, before host filtering: 2768
- Unique donor reactions equivalent to host and removed: 344
- Host-equivalent intracellular source rows with donor direction capacity beyond host bounds: 77
- Shared BiGG reaction xrefs associated with multiple canonical signatures (diagnostic only): 20
- Transitive metabolite alias components rejected for same-model, same-compartment conflicts: 9
- Source rows collapsed by donor-to-donor chemical deduplication: 1324
- Deferred unique donor rows before host filtering: 1519
- Deferred rows equivalent to host and removed: 34
- Deferred unique rows retained (transport and non-cytosolic): 1485

## Mapping and boundaries

Metabolite identifiers use BiGG, MetaNetX, KEGG, and BioCyc annotations as a conflict-checked alias graph: co-annotations connect aliases only when each alias identifies one chemical per model. A transitive component is disabled if it maps to multiple distinct metabolites in the same model and compartment. The graph uses a deterministic BiGG/MNX/KEGG/BioCyc canonical label. Legacy host `_DASH_` IDs are normalized; generic `M_*` donor IDs are not used as fallback aliases. Names and formulas do not force matches. Reaction stoichiometry is reduced to primitive integer ratios and canonicalized over forward/reverse orientation. Host equivalence is therefore assessed without direction; differing donor/host bounds are counted above for later review. Shared BiGG reaction xrefs are counted as split-signature diagnostics and never force a merge.

Unmapped intracellular reaction-source rows: 82; after donor deduplication and host filtering, 82 final candidate rows still contain unresolved metabolite mappings. Unmapped deferred source rows: 4; these yield 4 unresolved deferred rows. The exact donor/reaction IDs and equations are listed in `unresolved_mapping.tsv`.

Exchanges, demands, sinks, boundary reactions, and biomass reactions are excluded. Cross-compartment reactions are categorized as transport; reactions wholly within one non-cytosolic compartment are categorized separately. Host-equivalent deferred signatures are removed. This output is a reaction candidate pool only; no strain design or model edit was run.

## Reproduction

Run with Python and COBRApy 0.32.1 installed (for example `pip install cobra==0.32.1`): `python "D:\嗜酸氧化亚铁硫杆菌\ATCC23270-GEM\scripts\build_E1_donor_reaction_pool.py"`.
