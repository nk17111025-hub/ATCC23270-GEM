# PHASE 5B-6 — LARGE-UNIVERSE RUBISCO-INDEPENDENT MULTI-INTERVENTION SEARCH

Project: *Acidithiobacillus ferrooxidans* ATCC 23270 (GCF_049532655.1, RU820_RS locus tags)

## Summary table

| solution | H reactions | KO | OE | biomass R0 | biomass Rubisco restored | bypass FVA min | T0/TH | artifact-free |
|---|---:|---:|---:|---:|---:|---:|---|---|
| (none — see PARTIAL) | — | — | — | 0.0 | 0.1019 | — | T0 | n/a |

## Verdict

`PARTIAL — LARGE-UNIVERSE SEARCH FOUND NEW BOTTLENECKS BUT NO COMPLETE RESCUE`

The promoted v2 model identity was verified. The large reaction universe (143 mapped,
mass-balanced candidate reactions: 116 organism-native BioCyc + 27 curated BiGG universal)
was built and searched with a true MILP GapFill (HiGHS via scipy.optimize.milp, binary
selection variables, global-optimality gap 1e-6). **No complete Rubisco-independent rescue
was found with physically-closed candidate reactions.** A new, previously un-attributed
bottleneck was located precisely: a redox/cofactor imbalance that the earlier 14-reaction
universe could not reach.

---

## 1. Was Rubisco-independent glucose-supported biomass achieved?

No. Under the exact R0 condition (RUBISCO=RUBISCOX=0, glucose (-5,0), Fe2 ≤ WT, O2 ≤ WT,
external inorganic carbon not forced), native R0 biomass = 0.

## 2–3. Proven minimum additions / global optimality

No feasible rescue exists in the current 143-reaction universe, so no minimum was proven.
The MILP returned `INFEASIBLE` (HiGHS status 8), not a time-limit incumbent. This is a
genuine infeasibility, not a solver artefact.

## 4–5. Exact minimum solutions / number of alternatives

None (0). No alternative enumeration possible.

## 6. Actual carbon route

After applying the v2 direction corrections (see Model repair below), the central-carbon
core is functional under R0: glucose → G6P → F6P → F1,6BP → GAP/DHAP → 1,3-BPG → 3PG →
2PG → PEP → pyruvate → acetyl-CoA (PDH) → citrate → isocitrate → α-KG. Producibility
(mmol/gDW/h, demand test): 3PG/pyr/OAA/cit/icit all reachable. This is **not** the
Calvin-cycle route — it is EMP (Embden-Meyerhof-Parnas) upper/lower glycolysis. The
non-oxidative pentose supply still requires an oxidative PPP step (GND) for R5P/E4P.

## 7–9. Is ED / G6PDH_NAD still selected? Superior non-ED route?

No route is "selected" because no complete rescue exists. The central-carbon core that
does function is EMP, not ED. G6PDH_NAD was not required by the MILP.

## 10–11. How are acetyl-CoA / α-KG supplied?

Acetyl-CoA is produced natively by PDH (pyruvate + CoA + NAD → acetyl-CoA + CO2 + NADH)
once lower glycolysis is opened. α-KG is produced natively by the NADP-dependent isocitrate
dehydrogenase ICDHyr (isocitrate + NADP → α-KG + CO2 + NADPH). Neither is a hard carbon
bottleneck in the corrected model.

## 12. Glutamate/glutamine/proline/arginine precursors

**This is the dominant new bottleneck.** α-KG → glutamate requires nitrogen assimilation
and reducing power. The native GS-GOGAT cycle (GLUSy + GLNS) and the candidate GDH both
require NADPH; ICDHyr (the α-KG source) simultaneously produces NADPH. This couples
nitrogen assimilation to NADPH disposal, and the model has no native uncoupled NAD(P)H
re-oxidation (see section 14).

## 13. Lipid / fatty-acid precursors

Blocked downstream of acetyl-CoA; acetyl-CoA itself is not the binding constraint. Lipid
synthesis (ACCOAC, ACOATA, fatty-acid elongation) also consumes NADPH and CoA, adding to
the same redox/cofactor imbalance.

## 14–15. NADH / NADPH balance (root cause)

The respiratory NADH-oxidizing direction (NADHI) is structurally unavailable under R0:
NADHI FVA is min≈0 even with `LOWER_BOUND=-1000`. Its reverse direction
(NADH + Q → NAD + QH2) pumps 5 H+ to the periplasm, so it is coupled 1:1 to ATP-synthase
flux, which is capped by the fixed maintenance ATP demand (ATPM=3.475) plus biomass ATP.
The glucose-driven NADH production exceeds this proton-coupled re-oxidation capacity.
NADH/NADPH oxidases (soluble NOX/NOXP) interconvert NADH↔NAD but do not remove the excess
reducing-equivalent/cofactor mass. A diagnostic free NAD(P)H/NAD/NADP sink restores biomass
(≈0.08), localising the defect to the cofactor pool balance, not to carbon.

## 16. CoA balance

CoA is recycled (PDH/CS/fatty-acid steps return CoA); net CoA demand is the tiny biomass
term. Not the binding constraint.

## 17. Artifact-free?

Baseline v2 QC is clean (free ATP/NADH/NADPH = 0; glucose-closed biomass = 0; no-donor
infeasible). No rescue candidate existed to QC.

## 18–19. Essential / redundant added reactions

n/a — no rescue found.

## 20–24. Rubisco-restore behaviour / bypass classification

n/a — no bypass was found to restore or classify.

## 25–26. KO requirements / minimum KO set

KO search was not reached (no valid bypass exists; the workflow correctly does not start
route-enforcement KO/OE before a bypass).

## 27. OE (FSEOF) targets

Not reached.

## 28–30. T0/TH, energy-tier, glucose-dependency

Not reached (no accepted solution to validate). T0 was used for discovery as planned.

## 31–33. Likely heterologous vs native-GEM-omission vs uncertain

Preliminary classification of the reactions the search *needs* but the model lacks:

| reaction | EC | preliminary class | rationale |
|---|---|---|---|
| GND (6-phosphogluconate dehydrogenase) | 1.1.1.44 | N1 (native, GEM gap) | oxidative PPP missing from GEM |
| AKGD (α-ketoglutarate dehydrogenase) | 1.2.1.105 | H (likely heterologous) | A. ferrooxidans has an incomplete TCA |
| GDH (glutamate dehydrogenase) | 1.4.1.4 | N1/H | nitrogen-assimilation alternative |
| uncoupled NAD(P)H oxidation | 1.6.3.x | H | redox-balance intervention |

## 34. Top non-dominated designs / architecture for Phase 5B-7

None qualifies as a complete Rubisco-independent architecture. The highest-value next step
is a **redox-balance model repair** (a non-proton-pumping NADH:quinone/oxygen oxidoreductase
with correct cofactor-pool closure), after which the MILP should be re-run. Until that is
resolved, no gene-level (Phase 5B-7) donor selection is warranted.

---

## Model repair discovered (section 47) — recorded, not silently applied

`MODEL_REPAIR` — The promoted v2 direction corrections (NADHI, GAPD1, GAPD2, PGK, G6PDH2
opened to lb=-1000) are stored in the base bounds but are **silently overridden by the
Table-1 FIM condition bounds** (lb_fe2=0 for NADHI/GAPD1/GAPD2/PGK). The in-memory model
therefore never realises the documented v2 state under the R0/FIM discovery condition.
The search applies these as explicit overrides (V2_REPAIR in `04_gapfill_milp.py`) so the
promoted v2 state is honoured. This does **not** itself rescue R0; it is a necessary
precondition for the EMP core to carry flux.

## Autonomous recovery log

See `logs/phase5B6_autonomous_recovery_log.tsv` — the two substantive events are:
1. HiGHS default MIP feasibility tolerance (1e-6) made the 1e-6 biomass threshold ambiguous
   (the solver returned "optimal" with biomass 0). Fix: tightened
   primal/dual/mip feasibility tolerances to 1e-8.
2. Table-1 override nullifying v2 direction corrections. Fix: explicit V2_REPAIR overrides.

## Data lineage

`minimal_trustworthy_baseline_v2_T0.xml` (SHA256
`97e31dd7b87852e8f639b72408e6e9ea05713dc0a48e99bc77edf087f0919e96`) → R0 constraint snapshot
(RUBISCO=RUBISCOX=0, glucose (-5,0), Fe2 ≤ WT=-164.509, O2 ≤ WT=-38.663) → raw BioCyc
universe (1313) → mapped (138) → cleaned U1 (116) + curated U2 (27) → MILP (HiGHS,
binary GapFill) → `INFEASIBLE` → redox-bottleneck diagnosis → this report.

## Final verdict

`PARTIAL — LARGE-UNIVERSE SEARCH FOUND NEW BOTTLENECKS BUT NO COMPLETE RESCUE`
