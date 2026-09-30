# PHASE 3.2-R2 FINAL REPORT — TRUE LOOPLESS FVA

审计对象：`minimal_trustworthy_baseline_v1`（SHA256 `bd9715e6…3689af8c`）。只计算、不修复。

## 最终状态

`PHASE3_2R2_TECHNICAL_INCOMPLETE — FASTSNP LOOPLESS BOUNDS NOT FULLY VERIFIED`

原因：`fastSNP`（MILP 型无环 FVA）在可用求解器 GLPK 下计算不可行（单反应无环 FVA 即需数分钟，6 个 context 在 40 分钟内未完成）。原始 2016 论文使用商业求解器 Gurobi 5.5。

`cycleFreeFlux`（Saa-Nielsen 凸无环公式）在 6 个 context 全部完成（约 5 秒），给出了无环 min/max 边界，可作为无环边界的实际证据。

## 完整性

- v1 SHA256 校验通过，与冻结一致。
- 参考生长（非 null）：FIM=0.052076、TTM=0.052076、TSM=0.052076。

## Q1 — ACCOAC / BIOC1 / BIOC2（确认 LOOP-DRIVEN）

cycleFreeFlux 无环 FVA（FIM 100%）：

| 反应 | 标准 FVA | 无环 FVA | 判定 |
|---|---|---|---|
| ACCOAC | [-999.93, 0.068] | [0.0, 0.068] | LOOP-DRIVEN FVA FREEDOM |
| BIOC1 | [0, 1000] | [0, 0.068] | LOOP-DRIVEN FVA FREEDOM |
| BIOC2 | [0, 1000] | [0, 0.068] | LOOP-DRIVEN FVA FREEDOM |

ACCOAC 的极端逆向范围（-999.93）与 BIOC1/BIOC2 的极端正向范围（1000）在无环 FVA 下消失，确认 `ACCOAC(reverse)+BIOC1+BIOC2` 是净零内部循环（前已证净化学计量为零）。

## Q2 — BDGK（正向 only，无逆向）

无环 FVA：BDGK = [0, 55.7]（FIM 100%）。三条件/两 fraction 均无逆向（loopless_reverse_possible=False）。

结论：BDGK 无环下仅正向（葡萄糖→G6P），逆向（G6P→葡萄糖+ATP）不可行。因此 BDGK 的“方向性隐患”经无环 FVA 后降级为“低置信度证据问题”，而非反向方向问题。仍为 WATCH（Confidence=1、无 GPR）。

## Q3 — NADTRHD

无环 FVA：NADTRHD = [0, 836.2]（正向 only）。其容量在无环后保留，非循环驱动。
修正措辞：NADTRHD KO 引起分布式全网络重路由（ETC、叶酸、核苷酸、甘氨酸/丝氨酸/谷氨酸、中央碳），无单一人工救援。

## Q4 — 乙醛酸循环

无环 FVA（FIM 100%）：MALS [0, 0.158]（正向 only）、GLXCL [-0.079, 7.14]（逆向仍可能）、GLYCK [-0.079, 7.14]（逆向仍可能）、GLXR [0, 64.27]（正向 only）。

GLXCL/GLYCK 的逆向（-0.079）在无环后仍存在 → DIRECTIONALITY WATCH。MALS/GLXR 为无环休眠容量（WATCH）。

## Q5 — MDH / FUM

无环 FVA：MDH [0.028, 835.5]、FUM [-835.5, -0.028]。其大 FVA 宽度在无环后完整保留，说明 MDH/FUM 的灵活性是真实化学计量灵活性（不完整 TCA 的可逆段），非循环驱动。FUM 的负向即 fumarate→malate 方向，为真实可逆。仍 WATCH（低置信度乙醛酸路线提供部分 KO 补偿）。

## 数值一致性

7 处“bound_consistency_failures”均为 TSM 100% 下 ~1e-7 级数值噪声（无环 max 略超标准 max 5e-8~1.5e-7），非实质违反（这些反应不在内部循环，无环=标准）。

## 结论

- fastSNP 未完成（GLPK MILP 限制）→ 按规则 TECHNICAL INCOMPLETE。
- cycleFreeFlux 无环边界全部完成，科学结论：ACCOAC/BIOC1/BIOC2 为 LOOP-DRIVEN；BDGK 无逆向；NADTRHD/MALS/GLXR 为无环休眠容量；GLXCL/GLYCK 逆向仍可能（DIRECTIONALITY WATCH）；MDH/FUM 为真实可逆灵活性。
- 无新 BLOCKER。
