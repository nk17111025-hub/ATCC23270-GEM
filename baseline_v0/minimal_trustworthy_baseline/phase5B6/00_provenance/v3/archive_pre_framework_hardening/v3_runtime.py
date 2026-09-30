# -*- coding: utf-8 -*-
"""v3 executable runtime.

The v3 execution architecture is:

    load v3 (base bounds authoritative)
        -> apply ENVIRONMENTAL condition bounds (exchange reactions only)
        -> apply experiment/search overrides

Historical Table-1 internal-reaction bounds are NEVER applied. Any attempt to apply a
legacy internal bound over v3 raises BOUND_PRECEDENCE_VIOLATION.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import p5b6_common as C
import phase4a_lib as lib  # noqa: E402  (historical parsing helpers only)
from build_glucose_model import build_memory  # noqa: E402


class BoundPrecedenceViolation(Exception):
    """Raised when a legacy internal bound is attempted over v3."""


# v3 base corrections applied on top of build_memory (the promoted v2 lineage).
V3_BASE_CORRECTIONS = {
    "NADHI": (-1000.0, 1000.0),
}

# Reactions whose approved directionality must survive every condition.
PROTECTED = ["BDGK", "G6PDH2", "GAPD1", "GAPD2", "PGK", "NADHI"]


def load_v3(transport: str = "T0"):
    """Load the in-memory v3 model (base bounds authoritative)."""
    metas, rxns = build_memory(transport)
    for rid, (lb, ub) in V3_BASE_CORRECTIONS.items():
        if rid in rxns:
            rxns[rid]["lb"] = lb
            rxns[rid]["ub"] = ub
    return metas, rxns


def condition_environmental_bounds(rxns, cond):
    """Return {reaction: (lb, ub)} of EXCHANGE-only condition bounds for `cond`."""
    lbc, ubc = lib.COND_COL[cond]
    out = {}
    for r, rx in rxns.items():
        lv = rx["table1_" + lbc]
        uv = rx["table1_" + ubc]
        if lv == "" and uv == "":
            continue
        if not r.startswith("Ex_"):
            continue  # internal reaction -> not environmental, never applied in v3
        out[r] = (float(lv), float(uv))
    return out


def legacy_internal_bounds(rxns, cond):
    """Return {reaction: (base_lb, base_ub, table_lb, table_ub)} of internal reactions whose
    legacy Table-1 bound would overwrite the v3 base bound. Used only for reporting/guard."""
    lbc, ubc = lib.COND_COL[cond]
    out = {}
    for r, rx in rxns.items():
        if r == lib.BIOMASS_RXN:
            continue
        lv = rx["table1_" + lbc]
        uv = rx["table1_" + ubc]
        if lv == "" and uv == "":
            continue
        if r.startswith("Ex_"):
            continue
        bl, bu = rx["lb"], rx["ub"]
        if float(lv) != bl or float(uv) != bu:
            out[r] = (bl, bu, float(lv), float(uv))
    return out


def guard_no_internal_override(rxns, cond):
    """Abort if a legacy internal bound would overwrite an approved corrected reaction."""
    lbc, ubc = lib.COND_COL[cond]
    for r in PROTECTED:
        rx = rxns[r]
        lv = rx["table1_" + lbc]
        uv = rx["table1_" + ubc]
        if lv == "" and uv == "":
            continue
        bl, bu = rx["lb"], rx["ub"]
        if float(lv) != bl or float(uv) != bu:
            raise BoundPrecedenceViolation(
                f"legacy Table-1 {cond} bound ({lv},{uv}) would overwrite corrected "
                f"reaction {r} base bound ({bl},{bu})"
            )


def setup_v3(metas, rxns, cond, override=None):
    """Build FBA matrices for v3.

    Base bounds are authoritative for internal reactions; only EXCHANGE reactions receive
    condition-specific bounds from Table-1.
    """
    S, met_list, rxn_list = lib.build_matrix(metas, rxns)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = rxns[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]

    # environmental condition bounds (exchanges only)
    env = condition_environmental_bounds(rxns, cond)
    for r, (el, eu) in env.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = el
            ub[rxn_idx[r]] = eu

    # project-internal patches (direction corrections already in base bounds; patches
    # are the accepted knockout/activation set from Phase 4A)
    for r, (pl, pu) in lib.PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl
            ub[rxn_idx[r]] = pu

    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    return S, met_list, rxn_list, rxn_idx, lb, ub


def max_biomass_v3(metas, rxns, cond, override=None):
    S, ml, rl, ri, lb, ub = setup_v3(metas, rxns, cond, override)
    c = np.zeros(len(rl))
    c[ri[lib.BIO_EX]] = -1.0
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    if not res.success:
        return None
    return {r: res.x[ri[r]] for r in rl}


def biomass_v3(metas, rxns, cond, override=None):
    v = max_biomass_v3(metas, rxns, cond, override)
    return v[lib.BIO_EX] if v else 0.0


def pfba_v3(metas, rxns, cond, override=None):
    """Parsimonious FBA: maximize biomass, then minimize L1 flux norm (cycle-free)."""
    S, ml, rl, ri, lb, ub = setup_v3(metas, rxns, cond, override)
    n = len(rl)
    c = np.zeros(n)
    c[ri[lib.BIO_EX]] = -1.0
    res1 = linprog(c, A_eq=S, b_eq=np.zeros(len(ml)), bounds=list(zip(lb, ub)), method="highs")
    if not res1.success:
        return None
    opt = -res1.fun
    # stage 2: minimize sum |v| subject to biomass >= opt - eps
    eps = 1e-6 * (1.0 + abs(opt))
    c2 = np.zeros(2 * n)
    c2[n:] = 1.0
    A_eq = np.zeros((len(ml), 2 * n))
    A_eq[:, :n] = S
    A_ub = np.zeros((2 * n + 1, 2 * n))
    b_ub = np.zeros(2 * n + 1)
    for i in range(n):
        A_ub[i, i] = 1.0
        A_ub[i, n + i] = -1.0
        A_ub[n + i, i] = -1.0
        A_ub[n + i, n + i] = -1.0
    A_ub[2 * n, :n] = c
    b_ub[2 * n] = -(opt - eps)
    bnd2 = [(lb[j], ub[j]) for j in range(n)] + [(0.0, None) for _ in range(n)]
    res2 = linprog(c2, A_eq=A_eq, b_eq=np.zeros(len(ml)), A_ub=A_ub, b_ub=b_ub, bounds=bnd2, method="highs")
    if not res2.success:
        return None
    return {r: res2.x[ri[r]] for r in rl}
