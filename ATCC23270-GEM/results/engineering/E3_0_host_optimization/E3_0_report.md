# E3.0 whole-host optimization summary

Host SHA256: `026352372d0b04d2cc1518b0e37d92075f7af7f5112b94eb48d6302fc680f3a3`  
E1 pool SHA256: `23b45be4333e32f2b6535813d9bb017ef019afa417f3f494267145ad1b6df0d6`  
E2R3 fingerprint: `6a47044c69e9408352200ba87b45cf8a4e526822a2ecd4ab68c83f3f7940e52c`  
E3.0 script SHA256: `e124b2f647719022806312f205f5b1881a4c691072ba22687007341e14c0a8c1`  
Summary script SHA256: `a42af14082b6d1050dd2a375a67575e3897f48440040754121439ea46dee8d39`  
Solver: SCIP; Python 3.12.14; COBRApy 0.32.1; StrainDesign 1.19.1; SCIP seed 1.  
StrainDesign patch SHA256: `{"solver_interface.py": "e87fb3855d5269f0f2203a63c7aaee337125741f6f4d68c3868a7464024592c0", "strainDesignMILP.py": "1b76470e4a28874b3b1c5c0143aab1c60822c098082e3c92f1c8091616c088cc"}`

The branch jobs used separate fixed E2R3 reaction sets, references, candidate lists, checkpoints, and optimization fingerprints. Pareto metrics use reconstructed worst-case phenotype values and compare intervention count, Rubisco flux, absolute Fe2/O2/inorganic-carbon uptake, and biomass. The original optimizer TSVs remain unchanged. The unique frontier TSVs deduplicate identical KO plus regulatory constraint sets.

## S01

Fixed E1 reactions: `E1R01016; E1R01493; E1R01642`. Native candidate universe: 590 KO and 233 regulatory flux interventions.
Branch reference: biomass 0.0795062; glucose -0.5; Rubisco 0.0791345; Fe2 -127.901; O2 -31.2137; CO2/H2CO3 -0.0534467.

Budget results (statuses are shown as returned by the solver; ERROR and TIME_LIMIT remain distinct):

| Objective | Max cost | Solver status | Reconstruction | Solutions | Runtime (s) |
|---|---:|---|---|---:|---:|
| CBB_min | 1 | OPTIMAL | OPTIMAL | 1 | 9.07433 |
| CBB_min | 2 | OPTIMAL | OPTIMAL | 1 | 6.46946 |
| CBB_min | 4 | OPTIMAL | OPTIMAL | 1 | 4.27411 |
| CBB_min | 6 | OPTIMAL | OPTIMAL | 1 | 4.58402 |
| Fe2_efficiency | 1 | TIME_LIMIT |  | 0 | 313.487 |
| Fe2_efficiency | 2 | ERROR |  | 0 | 90.3631 |
| Fe2_efficiency | 4 | TIME_LIMIT |  | 0 | 305.66 |
| Fe2_efficiency | 6 | TIME_LIMIT |  | 0 | 305.989 |
| O2_efficiency | 1 | OPTIMAL | OPTIMAL | 5 | 19.6667 |
| O2_efficiency | 2 | TIME_LIMIT |  | 0 | 301.797 |
| O2_efficiency | 4 | TIME_LIMIT |  | 0 | 302.26 |
| O2_efficiency | 6 | TIME_LIMIT |  | 0 | 306.902 |
| CO2_efficiency | 1 | TIME_LIMIT_W_SOLS | OPTIMAL | 4 | 324.612 |
| CO2_efficiency | 2 | TIME_LIMIT_W_SOLS | OPTIMAL | 2 | 325.225 |
| CO2_efficiency | 4 | TIME_LIMIT_W_SOLS | OPTIMAL | 2 | 324.625 |
| CO2_efficiency | 6 | TIME_LIMIT_W_SOLS | OPTIMAL | 3 | 330.878 |

Unique valid designs after deduplicating KO/regulation sets: 6; nondominated designs: 6.
Repeated interventions across this branch's frontier: none observed.

Frontier phenotypes:

| ID | Native interventions | KO + OE | Biomass | Rubisco | Fe2 uptake | O2 uptake | CO2 uptake | Glucose off |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| S01_BASELINE_C0 | none | 0 + 0 | 0.0795062 | 0.0791345 | 157.785 | 38.4815 | 0.0534467 | optimal: 0.00206702 |
| S01_CO2_efficiency_C1_S1 | HCO3tex | 1 + 0 | 0.0781146 | 0.0791345 | 157.58 | 38.4815 | 0 | optimal: 6.13408e-14 |
| S01_CO2_efficiency_C1_S2 | HCO3tpp | 1 + 0 | 0.0781146 | 0.0791345 | 157.58 | 38.4815 | 0 | optimal: 6.13408e-14 |
| S01_CO2_efficiency_C1_S3 | TKT2 >= 0.097719444694 | 0 + 1 | 0.0782203 | 0.0791345 | 157.529 | 38.4815 | 0.0040619 | infeasible |
| S01_CO2_efficiency_C1_S4 | TKT1 <= -0.118536711739 | 0 + 1 | 0.0783929 | 0.0791345 | 157.563 | 38.4815 | 0.0106893 | infeasible |
| S01_CO2_efficiency_C6_S3 | HCO3E >= 0.144901475541;HYD3pp >= 0.429635297224 | 0 + 2 | 0.0781836 | 0.0791345 | 156.731 | 38.4815 | 0.00265013 | infeasible |

## S03

Fixed E1 reactions: `E1R00732; E1R00741`. Native candidate universe: 590 KO and 234 regulatory flux interventions.
Branch reference: biomass 0.0800899; glucose -0.5; Rubisco 0.0791345; Fe2 -130.795; O2 -31.9095; CO2/H2CO3 -0.0758606.

Budget results (statuses are shown as returned by the solver; ERROR and TIME_LIMIT remain distinct):

| Objective | Max cost | Solver status | Reconstruction | Solutions | Runtime (s) |
|---|---:|---|---|---:|---:|
| CBB_min | 1 | OPTIMAL | OPTIMAL | 1 | 4.37766 |
| CBB_min | 2 | OPTIMAL | OPTIMAL | 1 | 4.37271 |
| CBB_min | 4 | OPTIMAL | OPTIMAL | 1 | 4.62147 |
| CBB_min | 6 | OPTIMAL | OPTIMAL | 1 | 4.53896 |
| Fe2_efficiency | 1 | ERROR |  | 0 | 221.572 |
| Fe2_efficiency | 2 | TIME_LIMIT |  | 0 | 316.798 |
| Fe2_efficiency | 4 | TIME_LIMIT |  | 0 | 305.99 |
| Fe2_efficiency | 6 | TIME_LIMIT |  | 0 | 306.005 |
| O2_efficiency | 1 | OPTIMAL | OPTIMAL | 5 | 20.8061 |
| O2_efficiency | 2 | TIME_LIMIT |  | 0 | 302.668 |
| O2_efficiency | 4 | TIME_LIMIT |  | 0 | 302.193 |
| O2_efficiency | 6 | TIME_LIMIT |  | 0 | 305.966 |
| CO2_efficiency | 1 | TIME_LIMIT_W_SOLS | OPTIMAL | 2 | 320.617 |
| CO2_efficiency | 2 | TIME_LIMIT_W_SOLS | OPTIMAL | 2 | 326.112 |
| CO2_efficiency | 4 | TIME_LIMIT_W_SOLS | OPTIMAL | 2 | 324.539 |
| CO2_efficiency | 6 | TIME_LIMIT_W_SOLS | OPTIMAL | 4 | 335.622 |

Unique valid designs after deduplicating KO/regulation sets: 5; nondominated designs: 5.
Repeated interventions across this branch's frontier: none observed.

Frontier phenotypes:

| ID | Native interventions | KO + OE | Biomass | Rubisco | Fe2 uptake | O2 uptake | CO2 uptake | Glucose off |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| S03_BASELINE_C0 | none | 0 + 0 | 0.0800898 | 0.0791345 | 158.022 | 38.4815 | 0.0758606 | optimal: 0.00208257 |
| S03_CO2_efficiency_C1_S1 | HCO3tex | 1 + 0 | 0.0781146 | 0.0791345 | 157.727 | 38.4815 | 0 | optimal: 0 |
| S03_CO2_efficiency_C1_S2 | HCO3tpp | 1 + 0 | 0.0781146 | 0.0791345 | 157.727 | 38.4815 | 0 | optimal: 0 |
| S03_CO2_efficiency_C6_S3 | THRS >= 0.31893076874 | 0 + 1 | 0.0790496 | 0.0791345 | 157.867 | 38.4815 | 0.0359089 | infeasible |
| S03_CO2_efficiency_C6_S4 | TPI <= -0.318012781737 | 0 + 1 | 0.0785096 | 0.0791345 | 157.786 | 38.4815 | 0.0151721 | infeasible |

## S01 vs S03

| Branch | Unique valid intervention sets | Within-branch Pareto designs | Fewest positive host interventions |
|---|---:|---:|---:|
| S01 | 6 | 6 | 1 |
| S03 | 5 | 5 | 1 |

The unmodified S03 branch has the higher biomass optimum (0.0800899 versus 0.0795062).
Both branches have a reconstructed single native KO design that reaches zero external inorganic-carbon uptake at the biomass threshold. S01 uses 157.58 Fe2; S03 uses 157.727 under the common worst-case reconstruction. S03 therefore needs no fewer observed host interventions for this endpoint. The CO2-mode MILPs timed out with incumbents, so these designs are feasible findings rather than proven budget-wide optima.
No intervention was selected by the proven-optimal CBB runs, and the Fe2/O2 higher-budget runs did not resolve a new design. HCO3tex and HCO3tpp recur across CO2 budgets within each branch; no endogenous reaction is repeatedly selected across different objective modes. The current frontier leaves an overall backbone choice open: S03 has higher unmodified biomass, while S01 has lower Fe2 demand in the observed zero-carbon single-KO state.

Combined nondominated solutions retained: 9. Review `comparison/S01_vs_S03.tsv` for all branch Pareto designs and `comparison/nondominated_combined.tsv` for designs that remain nondominated across both branches.

Interpretation is limited to the model-level intervention hypotheses. Regulatory flux constraints do not imply a gene-expression fold change. A missing mode/budget row means the branch output did not contain that solve; inspect branch checkpoint status before interpreting it. No gene or experimental design was selected.
