# -*- coding: utf-8 -*-
"""Curated universal central-metabolism reaction set (BiGG-style) for the
Phase 5B-6 large-universe search. Reactions are written in the model's own
BiGG metabolite namespace with compartment suffix [c].

Each entry: id -> dict(name, reversible, stoich, ec, source, tier)
"""

# (id, name, reversible, stoich, ec, source_db)
CURATED = [
    # --- Entner-Doudoroff / oxidative PPP ---
    ("G6PDH2r_fwd", "glucose-6-phosphate dehydrogenase (oxidative)", False,
     {"g6p-B[c]": -1, "nadp[c]": -1, "6pgl[c]": 1, "nadph[c]": 1, "h[c]": 1},
     "1.1.1.49", "BiGG"),
    ("G6PDH_NAD", "glucose-6-phosphate dehydrogenase (NAD)", False,
     {"g6p-B[c]": -1, "nad[c]": -1, "6pgl[c]": 1, "nadh[c]": 1, "h[c]": 1},
     "1.1.1.363", "BiGG"),
    ("GND", "6-phosphogluconate dehydrogenase", False,
     {"6pgc[c]": -1, "nadp[c]": -1, "ru5p-D[c]": 1, "co2[c]": 1, "nadph[c]": 1},
     "1.1.1.44", "BiGG"),
    ("EDD", "6-phosphogluconate dehydratase", False,
     {"6pgc[c]": -1, "2ddg6p[c]": 1, "h2o[c]": 1},
     "4.2.1.12", "BiGG"),
    ("EDA", "2-dehydro-3-deoxy-phosphogluconate aldolase", False,
     {"2ddg6p[c]": -1, "g3p[c]": 1, "pyr[c]": 1},
     "4.1.2.14", "BiGG"),
    ("RPIr", "ribose-5-phosphate isomerase (R5P->Ru5P)", True,
     {"ru5p-D[c]": -1, "r5p[c]": 1},
     "5.3.1.6", "BiGG"),
    # --- TCA completion ---
    ("ACONTa", "aconitase (citrate->cis-aconitate)", True,
     {"cit[c]": -1, "acon-C[c]": 1, "h2o[c]": 1},
     "4.2.1.3", "BiGG"),
    ("ACONTb", "aconitase (cis-aconitate->isocitrate)", True,
     {"acon-C[c]": -1, "h2o[c]": -1, "icit[c]": 1},
     "4.2.1.3", "BiGG"),
    ("ICDH", "isocitrate dehydrogenase (NADP)", False,
     {"icit[c]": -1, "nadp[c]": -1, "akg[c]": 1, "co2[c]": 1, "nadph[c]": 1},
     "1.1.1.42", "BiGG"),
    ("AKGD", "2-oxoglutarate dehydrogenase", False,
     {"akg[c]": -1, "coa[c]": -1, "nad[c]": -1, "succoa[c]": 1, "co2[c]": 1, "nadh[c]": 1},
     "1.2.1.105", "BiGG"),
    # --- anaplerotic ---
    ("PC", "pyruvate carboxylase", False,
     {"pyr[c]": -1, "co2[c]": -1, "atp[c]": -1, "oaa[c]": 1, "adp[c]": 1, "pi[c]": 1, "h[c]": 1},
     "6.4.1.1", "BiGG"),
    ("ME1", "malic enzyme (NAD)", False,
     {"mal-L[c]": -1, "nad[c]": -1, "pyr[c]": 1, "co2[c]": 1, "nadh[c]": 1},
     "1.1.1.38", "BiGG"),
    ("ME2", "malic enzyme (NADP)", False,
     {"mal-L[c]": -1, "nadp[c]": -1, "pyr[c]": 1, "co2[c]": 1, "nadph[c]": 1},
     "1.1.1.40", "BiGG"),
    ("PEPCK", "phosphoenolpyruvate carboxykinase (ATP)", True,
     {"oaa[c]": -1, "atp[c]": -1, "pep[c]": 1, "co2[c]": 1, "adp[c]": 1},
     "4.1.1.49", "BiGG"),
    ("PPS", "phosphoenolpyruvate synthase", False,
     {"pyr[c]": -1, "atp[c]": -1, "h2o[c]": -1, "pep[c]": 1, "amp[c]": 1, "pi[c]": 1, "h[c]": 2},
     "2.7.9.2", "BiGG"),
    ("PPDK", "pyruvate phosphate dikinase", False,
     {"pyr[c]": -1, "atp[c]": -1, "pi[c]": -1, "pep[c]": 1, "amp[c]": 1, "ppi[c]": 1},
     "2.7.9.1", "BiGG"),
    # --- glyoxylate shunt ---
    ("ICL", "isocitrate lyase", False,
     {"icit[c]": -1, "glx[c]": 1, "succ[c]": 1},
     "4.1.3.1", "BiGG"),
    ("MALS", "malate synthase", False,
     {"accoa[c]": -1, "glx[c]": -1, "h2o[c]": -1, "coa[c]": 1, "mal-L[c]": 1, "h[c]": 1},
     "2.3.3.9", "BiGG"),
    # --- nitrogen assimilation ---
    ("GDH_NADP", "glutamate dehydrogenase (NADP)", False,
     {"akg[c]": -1, "nh4[c]": -1, "nadph[c]": -1, "h[c]": -1, "glu-L[c]": 1, "nadp[c]": 1, "h2o[c]": 1},
     "1.4.1.4", "BiGG"),
    ("GDH_NAD", "glutamate dehydrogenase (NAD)", False,
     {"akg[c]": -1, "nh4[c]": -1, "nadh[c]": -1, "h[c]": -1, "glu-L[c]": 1, "nad[c]": 1, "h2o[c]": 1},
     "1.4.1.2", "BiGG"),
    # --- cofactor shuttles ---
    ("THD", "NAD(P) transhydrogenase", True,
     {"nadph[c]": -1, "nad[c]": -1, "nadp[c]": 1, "nadh[c]": 1},
     "1.6.1.2", "BiGG"),
    ("NADH16", "NADH dehydrogenase (forward ETC)", False,
     {"nadh[c]": -1, "q8[c]": -1, "h[c]": -1, "nad[c]": 1, "q8h2[c]": 1},
     "1.6.5.3", "BiGG"),
    ("NOX", "NADH oxidase (H2O-forming, soluble)", False,
     {"nadh[c]": -1, "o2[c]": -0.5, "h[c]": -1, "nad[c]": 1, "h2o[c]": 1},
     "1.6.3.4", "BiGG"),
    ("NOXP", "NADPH oxidase (H2O-forming, soluble)", False,
     {"nadph[c]": -1, "o2[c]": -0.5, "h[c]": -1, "nadp[c]": 1, "h2o[c]": 1},
     "1.6.3.1", "BiGG"),
    # --- pyruvate / formate ---
    ("PFL", "pyruvate formate lyase", False,
     {"pyr[c]": -1, "coa[c]": -1, "accoa[c]": 1, "for[c]": 1},
     "2.3.1.54", "BiGG"),
    ("FDHr", "formate dehydrogenase (reverse)", False,
     {"for[c]": -1, "nad[c]": -1, "co2[c]": 1, "nadh[c]": 1},
     "1.2.1.2", "BiGG"),
    # --- gluconeogenesis bypass (PEP synthesis from OAA) ---
    ("PEPCK_GTP", "PEP carboxykinase (GTP)", True,
     {"oaa[c]": -1, "gtp[c]": -1, "pep[c]": 1, "co2[c]": 1, "gdp[c]": 1},
     "4.1.1.32", "BiGG"),
]


def curated_dict():
    out = {}
    for rid, name, rev, stoich, ec, src in CURATED:
        out[rid] = {
            "name": name,
            "reversible": rev,
            "stoich": stoich,
            "ec": ec,
            "source_db": src,
            "tier": "U2",
        }
    return out
