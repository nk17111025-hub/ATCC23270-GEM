# PHASE 5B-1C — CO2 SUBSTITUTION MINIMAL CLOSEOUT

Question: with the forced `Ex_h2co3 = -2.0` constraint removed, does glucose reduce external
inorganic-carbon requirement at exactly matched biomass, Fe/S, and O2?

## Frozen model integrity

- Frozen v1 SHA256 = `bd9715e6d2419a5bd718ae728721626f9d2ecb0bf566b9bf7198fe5d3689af8c`
  (UNCHANGED).
- Tested models: existing Phase 4.1G T0 and TH candidates (not rebuilt).

## CO2 bound change

In test scenarios only: `Ex_h2co3[e]` set from forced `(-2.0, -2.0)` to `(-1000, 0)`
(allow uptake, forbid net export). Uptake is negative; the optimization maximizes the CO2
exchange flux toward 0, i.e. minimizes external CO2 uptake magnitude.

## Result — CONDITIONAL_CO2_SUBSTITUTION

| condition | donor | max CO2 saving (T0) | max CO2 saving (TH) | glucose uptake |
|---|---|---|---|---|
| FIM | Fe2 | **72.0%** | 54.2% | −0.995 / −0.730 |
| TTM | tetrathionate | 0% | 0% | 0 |
| TSM | thiosulfate | 0% | 0% | 0 |

In FIM, glucose reduces external CO2 from 2.0 to 0.559 mmol/gDW/h (72% saving, T0) at exactly
matched biomass (0.052076), Fe2 (−164.509) and O2 (−38.663). The saving first appears at
glucose cap 0.05 (15%) and saturates near cap 1.

In TTM and TSM, glucose uptake is 0 at every cap and CO2 stays at 2.0 (0% saving).

## Anti-artifact verification

For every FIM CO2-saving solution:
- Fe2 uptake = −164.509 (exact match, Δ=0);
- O2 uptake = −38.663 (exact match, Δ=0);
- biomass = 0.052076 (exact match, Δ=0);
- no new energy uptake: H2 is PRODUCED (+0.278, a redox byproduct), not consumed;
- free ATP = 0, free NADH = 0, free NADPH = 0.

The CO2 saving is therefore not purchased by additional Fe/S or O2 throughput, and not by a
hidden energy source.

## Sanity check (energy unfixed)

When Fe/S and O2 are not pinned, glucose lowers CO2 to 0.0 in ALL conditions while also
lowering Fe/S and O2 (~45%: FIM Fe2 −164.5→−89.5, O2 −38.7→−21.9). This shows the matched
(energy-fixed) case is the conservative one; when energy is allowed to fall, the CO2 saving
is even larger. (This is not the primary metric.)

## Mechanism

In FIM, glucose supplies carbon skeletons via the PPP→RuBP→RUBISCO route, so the external CO2
requirement falls (RUBISCO still active, ~2.07). The Fe2 oxidation's energy/redox yield makes
this substitution favorable. In TTM/TSM the sulfur-donor energy/redox balance does not allow
the same substitution, so glucose is not taken up.

## Conclusion

**CONDITIONAL_CO2_SUBSTITUTION.** Real CO2 substitution exists under Fe2 (FIM) growth
(up to 72% at matched resources), but not under tetrathionate (TTM) or thiosulfate (TSM).
The result is robust to T0/TH qualitatively (both save CO2 in FIM), though the magnitude
depends on the transport assumption (72% vs 54%).

## Files

co2_reference_states.tsv, co2_matched_resource_test.tsv, co2_substitution_summary.tsv,
co2_sanity_unfixed_energy.tsv, artifact_checks.tsv, validation_summary.json,
test_co2_substitution.py (reproducible).
