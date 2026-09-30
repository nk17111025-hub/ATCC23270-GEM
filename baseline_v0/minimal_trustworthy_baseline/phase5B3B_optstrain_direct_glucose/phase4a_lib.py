# -*- coding: utf-8 -*-
"""Phase 4A shared library: SBML/XLS parsing + FBA (scipy HiGHS) + graph helpers."""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import lxml.etree as ET
import numpy as np
import xlrd
from scipy.optimize import linprog

ROOT = Path(r"D:\嗜酸氧化亚铁硫杆菌")
V1_XML = ROOT / "baseline_v0" / "minimal_trustworthy_baseline" / "freeze" / "minimal_trustworthy_baseline_v1.xml"
MMC1_XLS = ROOT / "01_原始数据" / "02_2016_iMC507原始模型" / "mmc1.xls"

NS = "{http://www.sbml.org/sbml/level2}"

PATCHES = {
    "MACPD": (0.0, 0.0),
    "ACOATA": (0.0, 1000.0),
    "Htpp": (0.0, 0.0),
}

BIO_EX = "Ex_bio[e]"
BIOMASS_RXN = "Afe_biomass_mc507_WT_139p0M"

COND_COL = {"FIM": ("lb_fe2", "ub_fe2"), "TTM": ("lb_ttton", "ub_ttton"), "TSM": ("lb_tsul", "ub_tsul")}
DONOR_EX = {"FIM": "Ex_fe2[e]", "TTM": "Ex_ttton[e]", "TSM": "Ex_tsul[e]"}

CENTRAL_TARGETS = {
    "G6P_F6P": ["g6p-A[c]", "g6p-B[c]", "f6p-B[c]"],
    "triose_P": ["dhap[c]", "g3p[c]"],
    "pyruvate": ["pyr[c]"],
    "acetyl_CoA": ["accoa[c]"],
    "TCA_OAA": ["oaa[c]", "mal-L[c]", "fum[c]", "succ[c]", "cit[c]", "icit[c]", "akg[c]"],
    "pentose_P": ["r5p[c]", "ru5p-D[c]", "xu5p-D[c]", "e4p[c]", "s7p[c]"],
}

CURRENCY_BASES = {
    "h", "h2o", "o2", "h2o2", "co2", "hco3", "h2co3", "nh4",
    "pi", "ppi", "pppi", "coa", "nad", "nadh", "nadp", "nadph",
    "fad", "fadh2", "atp", "adp", "amp", "gtp", "gdp", "gmp",
    "ctp", "cdp", "cmp", "utp", "udp", "ump",
}


def is_currency(met: str) -> bool:
    base = met[: met.index("[")] if "[" in met else met
    return base in CURRENCY_BASES

SUBSTRATE_SPECS = [
    {
        "substrate": "glucose",
        "ex_met": None,
        "p_met": None,
        "c_met": "glc-B[c]",
        "alt_c_met": ["glc-A[c]"],
        "transport_rxns": [],
        "entry_rxns": ["BDGK"],
        "notes": "intracellular glucose present; BDGK (glucokinase, conf1, no GPR) links glc-B to g6p-B; no glucose transporter or exchange",
    },
    {
        "substrate": "glycerol",
        "ex_met": None,
        "p_met": None,
        "c_met": "glyc[c]",
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": ["GLYK", "G3PD2"],
        "notes": "glycerol->glyc3p via GLYK (reversible, conf1, no GPR); glyc3p->dhap via G3PD2 is irreversible in reduction direction so oxidation needs reversal",
    },
    {
        "substrate": "ethanol",
        "ex_met": None,
        "p_met": None,
        "c_met": None,
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "no ethanol/acetaldehyde metabolite in model; requires full pathway construction",
    },
    {
        "substrate": "acetate",
        "ex_met": None,
        "p_met": None,
        "c_met": "ac[c]",
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": ["ACS"],
        "notes": "acetate->acetyl-CoA via ACS (irreversible forward, conf2, AFE_1969); no acetate transporter/exchange",
    },
    {
        "substrate": "fructose",
        "ex_met": None,
        "p_met": None,
        "c_met": None,
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "no free fructose metabolite (only f6p-B and fdp); needs fructose transport + kinase or isomerization",
    },
    {
        "substrate": "sucrose",
        "ex_met": None,
        "p_met": None,
        "c_met": None,
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "no sucrose metabolite; needs transporter + invertase/sucrose phosphorylase",
    },
    {
        "substrate": "lactate",
        "ex_met": None,
        "p_met": None,
        "c_met": None,
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "no lactate metabolite; needs transporter + lactate dehydrogenase",
    },
    {
        "substrate": "pyruvate",
        "ex_met": None,
        "p_met": None,
        "c_met": "pyr[c]",
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "pyruvate is already a central metabolite; only missing transport/exchange",
    },
    {
        "substrate": "formate",
        "ex_met": None,
        "p_met": None,
        "c_met": "for[c]",
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": ["FDH", "FTHFL"],
        "notes": "formate->CO2 via FDH (irreversible); FTHFL is written in deformylation direction (irreversible) so formate->10fthf is blocked; no C1 assimilation module",
    },
    {
        "substrate": "methanol",
        "ex_met": None,
        "p_met": None,
        "c_met": None,
        "alt_c_met": [],
        "transport_rxns": [],
        "entry_rxns": [],
        "notes": "no methanol metabolite; requires methanol dehydrogenase + C1 assimilation pathway",
    },
    {
        "substrate": "4hba",
        "ex_met": "4hba[e]",
        "p_met": "4hba[p]",
        "c_met": "4hba[c]",
        "alt_c_met": [],
        "transport_rxns": ["4HBAtex", "4HBAtpp"],
        "entry_rxns": [],
        "notes": "only organic compound with complete native exchange+transport; a thiamine-biosynthesis intermediate, not a central-carbon substrate",
    },
]


def biomass_precursors(reactions):
    """Organic-carbon metabolites consumed by the biomass reaction (excludes currency/ions)."""
    if BIOMASS_RXN not in reactions:
        return []
    ions = {"ca2", "cu2", "fe2", "fe3", "k", "mg2", "mn2", "na1", "zn2", "so4", "so4aa"}
    out = []
    for m, c in reactions[BIOMASS_RXN]["stoich"].items():
        if c >= 0:
            continue
        base = m[: m.index("[")] if "[" in m else m
        if base in CURRENCY_BASES or base in ions:
            continue
        out.append(m)
    return sorted(out)


def local(tag: str) -> str:
    return tag.split("}", 1)[1] if "}" in tag else tag


def enc_met(base: str, comp: str) -> str:
    return "M_" + base.replace("-", "_DASH_") + "_" + comp


def enc_rxn(key: str) -> str:
    return "R_" + key.replace("[", "_LSQBKT_").replace("]", "_RSQBKT_").replace("-", "_DASH_")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_table2() -> Dict[str, dict]:
    book = xlrd.open_workbook(str(MMC1_XLS), on_demand=True)
    s = book.sheet_by_name("Table 2")
    out: Dict[str, dict] = {}
    for r in range(1, s.nrows):
        abbr = str(s.cell_value(r, 0)).strip()
        if not abbr:
            continue
        out[abbr] = {
            "name": str(s.cell_value(r, 1)).strip(),
            "formula": str(s.cell_value(r, 2)).strip(),
            "charge": s.cell_value(r, 3),
            "location": str(s.cell_value(r, 4)).strip(),
            "kegg": str(s.cell_value(r, 6)).strip(),
            "pubchem": s.cell_value(r, 7),
        }
    book.release_resources()
    return out


def _read_table1() -> Dict[str, dict]:
    book = xlrd.open_workbook(str(MMC1_XLS), on_demand=True)
    s = book.sheet_by_name("Table 1")
    out: Dict[str, dict] = {}
    for r in range(2, s.nrows):
        rid = str(s.cell_value(r, 0)).strip()
        if not rid:
            continue

        def cv(c):
            v = s.cell_value(r, c)
            if v == "":
                return ""
            if isinstance(v, float) and v.is_integer():
                return int(v)
            return v

        out[rid] = {
            "name": str(s.cell_value(r, 1)).strip(),
            "formula": str(s.cell_value(r, 2)).strip(),
            "confidence": cv(3),
            "ec": str(s.cell_value(r, 4)).strip(),
            "pmid": cv(5),
            "subsystem": str(s.cell_value(r, 6)).strip(),
            "gpr": str(s.cell_value(r, 7)).strip(),
            "gpr2": str(s.cell_value(r, 8)).strip(),
            "protein": str(s.cell_value(r, 9)).strip(),
            "lb_fe2": cv(10),
            "ub_fe2": cv(11),
            "lb_ttton": cv(12),
            "ub_ttton": cv(13),
            "lb_tsul": cv(14),
            "ub_tsul": cv(15),
        }
    book.release_resources()
    return out


def parse_sbml_meta() -> Tuple[Dict[str, dict], Dict[str, dict]]:
    t2 = _read_table2()
    t1 = _read_table1()
    sp2met: Dict[str, str] = {}
    for abbr in t2:
        base = abbr
        comp = ""
        if "[" in abbr and abbr.endswith("]"):
            base = abbr[: abbr.index("[")]
            comp = abbr[abbr.index("[") + 1 : -1]
        sp2met[enc_met(base, comp)] = abbr
    sbml2rxn: Dict[str, str] = {}
    for key in t1:
        sbml2rxn[enc_rxn(key)] = key
    tree = ET.parse(str(V1_XML))
    root = tree.getroot()
    model = [ch for ch in root if local(ch.tag) == "model"][0]
    metabolites: Dict[str, dict] = {}
    for abbr, meta in t2.items():
        base = abbr[: abbr.index("[")] if "[" in abbr else abbr
        comp = abbr[abbr.index("[") + 1 : -1] if "[" in abbr else ""
        metabolites[abbr] = dict(meta)
        metabolites[abbr]["base"] = base
        metabolites[abbr]["compartment"] = comp
    reactions: Dict[str, dict] = {}
    for r in model.iter():
        if local(r.tag) != "reaction":
            continue
        sid = r.get("id")
        key = sbml2rxn.get(sid)
        if key is None:
            raise ValueError("unmapped sbml reaction " + sid)
        lb = ub = 0.0
        for p in r.iter():
            if local(p.tag) == "parameter":
                if p.get("id") == "LOWER_BOUND":
                    lb = float(p.get("value"))
                if p.get("id") == "UPPER_BOUND":
                    ub = float(p.get("value"))
        stoich: Dict[str, float] = {}
        for sr in r.iter():
            if local(sr.tag) != "speciesReference":
                continue
            sp = sr.get("species")
            met = sp2met.get(sp)
            if met is None:
                raise ValueError("unmapped sbml species " + sp)
            parent = sr.getparent().tag
            coeff = float(sr.get("stoichiometry"))
            if parent == NS + "listOfReactants":
                stoich[met] = stoich.get(met, 0.0) - coeff
            else:
                stoich[met] = stoich.get(met, 0.0) + coeff
        meta = t1.get(key, {})
        reactions[key] = {
            "sbml_id": sid,
            "name": r.get("name") or meta.get("name", ""),
            "reversible": r.get("reversible") == "true",
            "lb": lb,
            "ub": ub,
            "stoich": stoich,
            "confidence": meta.get("confidence", ""),
            "ec": meta.get("ec", ""),
            "pmid": meta.get("pmid", ""),
            "subsystem": meta.get("subsystem", ""),
            "gpr": meta.get("gpr", ""),
            "gpr2": meta.get("gpr2", ""),
            "protein": meta.get("protein", ""),
            "table1_lb_fe2": meta.get("lb_fe2", ""),
            "table1_ub_fe2": meta.get("ub_fe2", ""),
            "table1_lb_ttton": meta.get("lb_ttton", ""),
            "table1_ub_ttton": meta.get("ub_ttton", ""),
            "table1_lb_tsul": meta.get("lb_tsul", ""),
            "table1_ub_tsul": meta.get("ub_tsul", ""),
        }
    return metabolites, reactions


def equation_str(stoich: Dict[str, float]) -> str:
    lhs = []
    rhs = []
    for met, c in stoich.items():
        if c < 0:
            lhs.append((met, -c))
        elif c > 0:
            rhs.append((met, c))

    def fmt(term):
        met, c = term
        cs = "" if abs(c - 1.0) < 1e-12 else (str(int(c)) if abs(c - round(c)) < 1e-12 else str(round(c, 6)))
        return (cs + " " + met).strip()

    return " + ".join(fmt(t) for t in sorted(lhs)) + " -> " + " + ".join(fmt(t) for t in sorted(rhs))


def build_matrix(metabolites, reactions):
    met_list = sorted(metabolites.keys())
    rxn_list = list(reactions.keys())
    met_idx = {m: i for i, m in enumerate(met_list)}
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    S = np.zeros((len(met_list), len(rxn_list)))
    for j, r in enumerate(rxn_list):
        for met, c in reactions[r]["stoich"].items():
            S[met_idx[met], j] = c
    return S, met_list, rxn_list


def run_fba(metabolites, reactions, condition, objective=BIO_EX, override=None, with_table1=True):
    S, met_list, rxn_list = build_matrix(metabolites, reactions)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = reactions[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]
    if with_table1:
        lbc, ubc = COND_COL[condition]
        for r in rxn_list:
            if r == BIOMASS_RXN:
                continue
            rx = reactions[r]
            lv = rx["table1_" + lbc]
            uv = rx["table1_" + ubc]
            if lv != "" and uv != "":
                lb[rxn_idx[r]] = float(lv)
                ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl
            ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    c = np.zeros(n)
    if objective in rxn_idx:
        c[rxn_idx[objective]] = -1.0
    bounds = [(lb[j], ub[j]) for j in range(n)]
    res = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=bounds, method="highs")
    if not res.success:
        return None, {}, res.status
    fluxes = {r: res.x[rxn_idx[r]] for r in rxn_list}
    return fluxes, dict(zip(rxn_list, bounds)), res.fun


def run_pfba(metabolites, reactions, condition, objective=BIO_EX, override=None, with_table1=True):
    """Two-stage parsimonious FBA: maximize objective, then minimize L1 flux norm."""
    S, met_list, rxn_list = build_matrix(metabolites, reactions)
    rxn_idx = {r: j for j, r in enumerate(rxn_list)}
    n = len(rxn_list)
    lb = np.zeros(n)
    ub = np.zeros(n)
    for j, r in enumerate(rxn_list):
        rx = reactions[r]
        lb[j], ub[j] = rx["lb"], rx["ub"]
    if with_table1:
        lbc, ubc = COND_COL[condition]
        for r in rxn_list:
            if r == BIOMASS_RXN:
                continue
            rx = reactions[r]
            lv = rx["table1_" + lbc]
            uv = rx["table1_" + ubc]
            if lv != "" and uv != "":
                lb[rxn_idx[r]] = float(lv)
                ub[rxn_idx[r]] = float(uv)
    for r, (pl, pu) in PATCHES.items():
        if r in rxn_idx:
            lb[rxn_idx[r]] = pl
            ub[rxn_idx[r]] = pu
    if override:
        for r, (ol, ou) in override.items():
            if r in rxn_idx:
                lb[rxn_idx[r]] = ol
                ub[rxn_idx[r]] = ou
    c = np.zeros(n)
    if objective in rxn_idx:
        c[rxn_idx[objective]] = -1.0
    bounds = [(lb[j], ub[j]) for j in range(n)]
    res1 = linprog(c, A_eq=S, b_eq=np.zeros(len(met_list)), bounds=bounds, method="highs")
    if not res1.success:
        return None, {}, res1.status
    opt = -res1.fun
    # stage 2: minimize sum |v| subject to objective >= opt - eps
    eps = 1e-6 * (1.0 + abs(opt))
    # variables: v (n), t (n)
    c2 = np.zeros(2 * n)
    c2[n:] = 1.0
    A_eq = np.zeros((len(met_list), 2 * n))
    A_eq[:, :n] = S
    # |v_i| <= t_i  ->  v_i - t_i <= 0 ; -v_i - t_i <= 0
    A_ub = np.zeros((2 * n + 1, 2 * n))
    b_ub = np.zeros(2 * n + 1)
    for i in range(n):
        A_ub[i, i] = 1.0
        A_ub[i, n + i] = -1.0
        A_ub[n + i, i] = -1.0
        A_ub[n + i, n + i] = -1.0
    # objective constraint: -c^T v <= -(opt-eps)
    A_ub[2 * n, :n] = c
    b_ub[2 * n] = -(opt - eps)
    bnd2 = [(lb[j], ub[j]) for j in range(n)] + [(0.0, None) for _ in range(n)]
    res2 = linprog(c2, A_eq=A_eq, b_eq=np.zeros(len(met_list)), A_ub=A_ub, b_ub=b_ub, bounds=bnd2, method="highs")
    if not res2.success:
        return None, dict(zip(rxn_list, bounds)), res2.status
    fluxes = {r: res2.x[rxn_idx[r]] for r in rxn_list}
    return fluxes, dict(zip(rxn_list, bounds)), opt


def graph_edges(metabolites, reactions, all_reversible=False):
    """Directed metabolite graph edges (src, dst, rxn, is_reverse).

    Forward direction always allowed; reverse allowed if reaction reversible (lb<0)
    or when all_reversible=True.
    """
    edges = {}
    for m in metabolites:
        edges[m] = []
    for r, rx in reactions.items():
        stoich = rx["stoich"]
        rev = all_reversible or rx["lb"] < 0
        reactants = [m for m, c in stoich.items() if c < 0]
        products = [m for m, c in stoich.items() if c > 0]
        for src in reactants:
            for dst in products:
                edges[src].append((dst, r, False))
                if rev:
                    edges[dst].append((src, r, True))
    return edges


def shortest_path_rxn(edges, start, targets, max_depth=40):
    """BFS shortest path over metabolite graph. Returns (target, met_path, rxn_path, n_reverse)."""
    from collections import deque
    targets = set(targets)
    if start in targets:
        return start, [start], [], 0
    q = deque([(start, [start], [], 0)])
    seen = {start}
    while q:
        cur, mpath, rpath, nrev = q.popleft()
        if len(mpath) - 1 > max_depth:
            continue
        for (nxt, rxn, is_rev) in edges.get(cur, []):
            if nxt in seen:
                continue
            if is_currency(nxt) and nxt not in targets:
                continue
            nm = mpath + [nxt]
            nr = rpath + [rxn]
            nn = nrev + (1 if is_rev else 0)
            if nxt in targets:
                return nxt, nm, nr, nn
            seen.add(nxt)
            q.append((nxt, nm, nr, nn))
    return None, [], [], 0


def adjacency(metabolites, reactions, directed=True):
    adj: Dict[str, Dict[str, List[str]]] = {m: {} for m in metabolites}
    for r, rx in reactions.items():
        stoich = rx["stoich"]
        for src, cs in stoich.items():
            if cs >= 0:
                continue
            for dst, cd in stoich.items():
                if cd <= 0 or src == dst:
                    continue
                adj.setdefault(src, {}).setdefault(dst, []).append(r)
                if not directed and rx["reversible"]:
                    adj.setdefault(dst, {}).setdefault(src, []).append(r)
    return adj


def shortest_met_path(metabolites, reactions, start, targets, directed=True, max_depth=30):
    adj = adjacency(metabolites, reactions, directed=directed)
    targets = set(targets)
    from collections import deque
    q = deque([(start, [start])])
    seen = {start}
    while q:
        cur, path = q.popleft()
        if len(path) - 1 > max_depth:
            continue
        if cur in targets and cur != start:
            return cur, path
        for nxt in adj.get(cur, {}):
            if nxt not in seen:
                seen.add(nxt)
                q.append((nxt, path + [nxt]))
    return None, []
