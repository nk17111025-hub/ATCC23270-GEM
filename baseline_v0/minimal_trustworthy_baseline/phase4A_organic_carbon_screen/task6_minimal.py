# -*- coding: utf-8 -*-
"""Task 6: minimal functional-completion enumeration (hypothetical, model not modified)."""
import csv
from pathlib import Path

import phase4a_lib as lib

OUT = Path(r"D:\嗜酸氧化亚铁硫杆菌\baseline_v0\minimal_trustworthy_baseline\phase4A_organic_carbon_screen")

# (min_interventions, missing_functions, nearest_central_metabolite, predicted_length, intervention_type, notes)
MINIMAL = {
    "glucose": (1, "ADD_TRANSPORTER (glucose outer+inner-membrane uptake)",
                "g6p-B[c]", "1 (glc->g6p via BDGK); to pyruvate needs +2 reversals (GAPD1, PGK) or the PPP/RUBISCO route",
                "ADD_TRANSPORTER; native BDGK low-confidence",
                "BDGK confidence 1 / no GPR; lower EMP written gluconeogenic (GAPD1/PGK irreversible)"),
    "glycerol": (2, "ADD_TRANSPORTER + reverse G3PD2 (glyc3p->DHAP oxidation)",
                 "glyc3p[c]", "3 (glyc->glyc3p->DHAP->central) after G3PD2 reversal",
                 "ADD_TRANSPORTER + CHANGE_DIRECTIONALITY",
                 "GLYK reversible (conf1, no GPR); G3PD2 written DHAP->glyc3p (irreversible)"),
    "acetate": (1, "ADD_TRANSPORTER (acetate uptake)",
                "accoa[c]", "1 (ac->accoa via ACS)",
                "ADD_TRANSPORTER",
                "ACS present (conf2, AFE_1969); acid-stress toxicity adjudicated externally"),
    "pyruvate": (1, "ADD_TRANSPORTER (pyruvate uptake)",
                 "pyr[c]", "0 (already a central metabolite)",
                 "ADD_TRANSPORTER",
                 "pyruvate is already central; PDH->acetyl-CoA present"),
    "formate": (3, "ADD_TRANSPORTER + reverse FTHFL + build C1 assimilation module (major)",
                "10fthf[c]", "major (no serine/RuMP C1 assimilation in model)",
                "ADD_TRANSPORTER + CHANGE_DIRECTIONALITY + MULTI_STEP_HETEROLOGOUS_PATHWAY",
                "FTHFL written 10fthf->formate (irreversible); formate only -> CO2 (FDH) or purine formyl (GART)"),
    "ethanol": (3, "ADD_TRANSPORTER + ethanol dehydrogenase (ethanol->acetaldehyde) + aldehyde dehydrogenase (acetaldehyde->acetyl-CoA)",
                "accoa[c]", "3 (etoh->acald->accoa)",
                "ADD_TRANSPORTER + ADD_ONE_ENZYME x2",
                "no ethanol/acetaldehyde metabolite or enzyme in model; redox consequences (2 NADH)"),
    "fructose": (2, "ADD_TRANSPORTER + fructokinase (fructose->F6P)",
                 "f6p-B[c]", "2 (fru->f6p)",
                 "ADD_TRANSPORTER + ADD_ONE_ENZYME",
                 "no free-fructose metabolite; f6p-B present"),
    "sucrose": (2, "ADD_TRANSPORTER + invertase/sucrose phosphorylase (sucrose->glucose+fructose)",
                "glc-B[c]", "3 (sucr->glc->g6p via BDGK)",
                "ADD_TRANSPORTER + ADD_ONE_ENZYME",
                "no sucrose metabolite or cleavage enzyme"),
    "lactate": (2, "ADD_TRANSPORTER + lactate dehydrogenase (lactate->pyruvate)",
                "pyr[c]", "2 (lac->pyr)",
                "ADD_TRANSPORTER + ADD_ONE_ENZYME",
                "no lactate metabolite or LDH; GLYCH is D-lactate oxidase acting on glycolate not free lactate"),
    "methanol": (3, "ADD_TRANSPORTER + methanol dehydrogenase (methanol->formaldehyde/formate) + full C1 assimilation (major)",
                 "for[c]", "major (no C1 assimilation module)",
                 "ADD_TRANSPORTER + ADD_ONE_ENZYME + MULTI_STEP_HETEROLOGOUS_PATHWAY",
                 "no methanol metabolite or methanol dehydrogenase"),
    "4hba": (None, "native exchange+transport present but no 4hba->central-carbon catabolism (multi-enzyme, impractical)",
             "none", "n/a",
             "MULTI_STEP_HETEROLOGOUS_PATHWAY",
             "4hba is a thiamine-biosynthesis intermediate; not a realistic carbon substrate"),
}


def main():
    metas, rxns = lib.parse_sbml_meta()
    rows = []
    for spec in lib.SUBSTRATE_SPECS:
        sub = spec["substrate"]
        d = MINIMAL.get(sub, (None, "", "", "", "", ""))
        rows.append({
            "substrate": sub,
            "minimum_number_of_functional_interventions": d[0] if d[0] is not None else "n/a",
            "candidate_missing_functions": d[1],
            "nearest_existing_central_metabolite": d[2],
            "predicted_pathway_length_after_completion": d[3],
            "intervention_type": d[4],
            "notes": d[5],
        })

    fields = ["substrate", "minimum_number_of_functional_interventions", "candidate_missing_functions",
              "nearest_existing_central_metabolite", "predicted_pathway_length_after_completion",
              "intervention_type", "notes"]
    with (OUT / "minimal_completion_candidates.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote minimal completion rows", len(rows))
    for r in rows:
        print(r["substrate"], "min", r["minimum_number_of_functional_interventions"], "nearest", r["nearest_existing_central_metabolite"])


if __name__ == "__main__":
    main()
