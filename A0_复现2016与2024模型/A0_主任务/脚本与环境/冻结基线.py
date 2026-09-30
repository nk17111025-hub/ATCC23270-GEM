from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import traceback
import xml.etree.ElementTree as ET
from collections import OrderedDict, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path


AGENT_THREAD = "01a0c8ec-839a-7b42-9dc9-9cb5a5a7ef32"
ACTUAL_MODEL = "未核实（当前工具没有返回实际模型字段）"
ACTUAL_REASONING = "未核实（当前工具没有返回实际 reasoning effort 字段）"
START_TIME = "2026-09-22T20:00:00+08:00"
EXPECTED_MODEL = "gpt-5.6-luna"
EXPECTED_REASONING = "high"

try:
    import xlrd
    import openpyxl
except Exception as exc:  # pragma: no cover - 给无依赖环境保留诊断路径
    xlrd = None
    openpyxl = None
    IMPORT_ERROR = repr(exc)
else:
    IMPORT_ERROR = ""


XML_NS = "http://www.sbml.org/sbml/level2"
NS = {"s": XML_NS}
ARROW_RE = re.compile(r"(<=>|<->|=>|->|→|↔|=)")
TERM_RE = re.compile(r"^\s*(?:(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+)?(.+?)\s*$")
SPECIES_RE = re.compile(r"^(.+?)\[([A-Za-z0-9_]+)\]$")
GENE_RE = re.compile(r"\b[A-Za-z][A-Za-z0-9_]*\b")


def now_iso() -> str:
    return dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def scalar(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer():
            return str(int(value))
        return format(value, ".15g")
    return str(value)


def clean(value) -> str:
    return scalar(value).replace("\r", " ").replace("\n", " ").strip()


def numeric(value: str) -> str:
    value = clean(value)
    if not value:
        return ""
    try:
        d = Decimal(value)
    except InvalidOperation:
        return value
    if d == d.to_integral():
        return str(d.quantize(Decimal(1)))
    text = format(d.normalize(), "f")
    text = text.rstrip("0").rstrip(".")
    return text or "0"


def num_float(value: str, default=None):
    value = clean(value)
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def write_tsv(path: Path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean(row.get(key, "")) for key in headers})


def write_json(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_xls(path: Path):
    if xlrd is None:
        raise RuntimeError("未安装 xlrd：" + IMPORT_ERROR)
    book = xlrd.open_workbook(str(path), on_demand=True)
    output = OrderedDict()
    for sheet in book.sheets():
        rows = []
        for r in range(sheet.nrows):
            values = {}
            for c in range(sheet.ncols):
                value = sheet.cell_value(r, c)
                if value != "":
                    values[c + 1] = scalar(value)
            values["__row__"] = r + 1
            rows.append(values)
        output[sheet.name] = rows
    return output


def read_xlsx(path: Path):
    if openpyxl is None:
        raise RuntimeError("未安装 openpyxl：" + IMPORT_ERROR)
    book = openpyxl.load_workbook(str(path), data_only=True, read_only=True)
    output = OrderedDict()
    for sheet in book.worksheets:
        rows = []
        for row_number, values in enumerate(sheet.iter_rows(values_only=True), 1):
            row = {c + 1: scalar(value) for c, value in enumerate(values) if value is not None and value != ""}
            row["__row__"] = row_number
            rows.append(row)
        output[sheet.title] = rows
    book.close()
    return output


def parse_formula(formula: str):
    """返回方向、带符号化学计量；空侧用于交换反应。"""
    text = clean(formula)
    match = ARROW_RE.search(text)
    if not match:
        return None, {}, "缺少反应箭头"
    arrow = match.group(1)
    direction = "可逆" if arrow in ("<=>", "<->", "↔") else "正向"
    left, right = text[: match.start()], text[match.end() :]

    def side(value, sign):
        result = defaultdict(float)
        if not value.strip():
            return result, ""
        for term in value.split("+"):
            term = term.strip()
            if not term:
                continue
            m = TERM_RE.match(term)
            if not m:
                return result, f"无法解析项：{term}"
            coefficient = float(m.group(1) or 1.0)
            species = m.group(2).strip()
            sm = SPECIES_RE.match(species)
            if not sm:
                return result, f"缺少区室标记：{species}"
            met = f"{sm.group(1)}[{sm.group(2)}]"
            result[met] += sign * coefficient
        return result, ""

    stoich, warning = side(left, -1)
    if warning:
        return direction, {}, warning
    right_stoich, warning = side(right, 1)
    if warning:
        return direction, {}, warning
    for key, value in right_stoich.items():
        stoich[key] += value
    return direction, {key: value for key, value in stoich.items() if abs(value) > 1e-15}, ""


def canonical_formula(formula: str) -> str:
    direction, stoich, warning = parse_formula(formula)
    if warning:
        return "解析失败：" + warning
    parts = []
    for key in sorted(stoich):
        value = stoich[key]
        number = numeric(str(value))
        parts.append(f"{key}:{number}")
    return f"{direction}|" + ";".join(parts)


def norm_text(value: str) -> str:
    return re.sub(r"\s+", " ", clean(value)).strip()


def norm_gpr(value: str) -> str:
    return re.sub(r"\s+", " ", clean(value)).strip().lower()


def xls_table(root: Path):
    book = read_xls(root / "01_原始资料/2016_Campodonico/mmc1.xls")
    sheet = book["Table 1"]
    records = OrderedDict()
    for row in sheet[2:]:  # 原始 Excel 第 3 行开始是反应数据
        reaction_id = clean(row.get(1, ""))
        if not reaction_id:
            continue
        records[reaction_id] = {
            "id": reaction_id,
            "name": clean(row.get(2, "")),
            "formula": clean(row.get(3, "")),
            "confidence": clean(row.get(4, "")),
            "ec": clean(row.get(5, "")),
            "pmid": clean(row.get(6, "")),
            "subsystem": clean(row.get(7, "")),
            "gpr": clean(row.get(8, "")),
            "lb": clean(row.get(11, "")),
            "ub": clean(row.get(12, "")),
            "cell": f"A{row['__row__']}:P{row['__row__']}",
            "gpr_cell": f"H{row['__row__']}",
            "source": "2016 mmc1.xls Table 1",
        }
    return records, {"sheet_rows": len(sheet), "data_rows": len(records), "blank_rows_after_header": sum(not any(v for k, v in row.items() if k != "__row__") for row in sheet[2:])}


def xlsx_tables(root: Path):
    book = read_xlsx(root / "01_原始资料/2024_Khaleque/from2024-mmc1.xlsx")
    s1 = book["Table S1"]
    reactions = OrderedDict()
    for row in s1[1:]:
        reaction_id = clean(row.get(1, ""))
        if not reaction_id:
            continue
        reactions[reaction_id] = {
            "id": reaction_id,
            "name": clean(row.get(2, "")),
            "formula": clean(row.get(3, "")),
            "gpr": clean(row.get(4, "")),
            "genes": clean(row.get(5, "")),
            "proteins": clean(row.get(6, "")),
            "subsystem": clean(row.get(7, "")),
            "reversible": clean(row.get(8, "")),
            "lb": clean(row.get(9, "")),
            "ub": clean(row.get(10, "")),
            "objective": clean(row.get(11, "")),
            "confidence": clean(row.get(12, "")),
            "ec": clean(row.get(13, "")),
            "notes": clean(row.get(14, "")),
            "references": clean(row.get(15, "")),
            "cell": f"A{row['__row__']}:O{row['__row__']}",
            "gpr_cell": f"D{row['__row__']}",
            "source": "2024 from2024-mmc1.xlsx Table S1",
        }
    s2 = book["Table S2"]
    metabolites = OrderedDict()
    for row in s2[1:]:
        metabolite_id = clean(row.get(1, ""))
        if not metabolite_id:
            continue
        metabolites[metabolite_id] = {
            "id": metabolite_id,
            "name": clean(row.get(2, "")),
            "formula_neutral": clean(row.get(3, "")),
            "formula_charged": clean(row.get(4, "")),
            "charge": clean(row.get(5, "")),
            "compartment": clean(row.get(6, "")),
            "kegg": clean(row.get(7, "")),
            "pubchem": clean(row.get(8, "")),
            "chebi": clean(row.get(9, "")),
            "cell": f"A{row['__row__']}:K{row['__row__']}",
            "source": "2024 from2024-mmc1.xlsx Table S2",
        }
    scenarios = {sheet: book[sheet] for sheet in ("Table S3", "Table S4")}
    return reactions, metabolites, scenarios, {
        "s1_sheet_rows": len(s1),
        "s1_data_rows": len(reactions),
        "s2_sheet_rows": len(s2),
        "s2_data_rows": len(metabolites),
        "s2_blank_rows": sum(not any(v for k, v in row.items() if k != "__row__") for row in s2[1:]),
    }


def sbml_stats(path: Path):
    root = ET.parse(path).getroot()
    model = root.find("s:model", NS)
    reactions = model.findall("s:listOfReactions/s:reaction", NS)
    species = model.findall("s:listOfSpecies/s:species", NS)
    compartments = model.findall("s:listOfCompartments/s:compartment", NS)
    objectives = []
    for reaction in reactions:
        params = {p.attrib.get("id"): p.attrib.get("value") for p in reaction.findall(".//s:parameter", NS)}
        if params.get("OBJECTIVE_COEFFICIENT") not in (None, "0", "0.0"):
            objectives.append({"reaction": reaction.attrib.get("id"), "coefficient": params.get("OBJECTIVE_COEFFICIENT")})
    gpr_tags = sum(1 for element in model.iter() if element.tag.lower().endswith("geneproductassociation"))
    boundary = [s.attrib.get("id") for s in species if s.attrib.get("boundaryCondition") == "true"]
    return {
        "sbml_level": root.attrib.get("level"),
        "sbml_version": root.attrib.get("version"),
        "model_id": model.attrib.get("id"),
        "compartment_count": len(compartments),
        "compartments": [{k: v for k, v in c.attrib.items()} for c in compartments],
        "species_count": len(species),
        "reaction_count": len(reactions),
        "gene_product_association_tag_count": gpr_tags,
        "boundary_species_count": len(boundary),
        "boundary_species": boundary,
        "fbc_objective_node_count": len(model.findall(".//{http://www.sbml.org/sbml/level3/version1/fbc/version2}objective")),
        "legacy_objective_coefficients": objectives,
    }


def sbml_warnings(path: Path):
    try:
        import libsbml
        # Windows下libSBML按文件名读取中文路径会误报File unreadable；
        # 按字节解码后从字符串读取，仍然检查同一份SBML内容。
        document = libsbml.readSBMLFromString(path.read_text(encoding="utf-8"))
        errors = []
        for i in range(document.getNumErrors()):
            error = document.getError(i)
            errors.append({"line": error.getLine(), "severity": error.getSeverityAsString(), "message": error.getMessage()})
        return {"available": True, "error_count": len(errors), "errors": errors[:100]}
    except Exception as exc:
        return {"available": False, "error_count": None, "errors": [repr(exc)]}


def sanitize(raw: str, prefix: str, used: set[str]) -> str:
    base = re.sub(r"[^A-Za-z0-9_]", "_", raw)
    if not base or base[0].isdigit():
        base = "x_" + base
    candidate = prefix + base
    index = 2
    while candidate in used:
        candidate = f"{prefix}{base}_{index}"
        index += 1
    used.add(candidate)
    return candidate


def build_cobra_model(reactions, metabolites, out_xml: Path, output_dir: Path):
    from cobra import Model, Reaction, Metabolite
    from cobra.io import write_sbml_model

    model = Model("iATF_ENG_2024_paper_reconstruction")
    met_map = {}
    reaction_map = {}
    used_met = set()
    compartment_map = {"cytosol": "c", "periplasm": "p", "extracellular": "e"}
    for raw, row in metabolites.items():
        sid = sanitize(raw, "M_", used_met)
        comp_name = row.get("compartment", "")
        comp = compartment_map.get(comp_name.lower(), comp_name[:1].lower() if comp_name else "c")
        met = Metabolite(sid, name=row.get("name", raw), compartment=comp)
        formula = row.get("formula_charged", "")
        if formula:
            met.formula = formula
        charge = num_float(row.get("charge", ""))
        if charge is not None and charge.is_integer():
            met.charge = int(charge)
        model.add_metabolites([met])
        met_map[raw] = met

    used_reactions = set()
    parse_warnings = []
    gpr_warnings = []
    for raw, row in reactions.items():
        rid = sanitize(raw, "R_", used_reactions)
        reaction_map[raw] = rid
        direction, stoich, warning = parse_formula(row.get("formula", ""))
        if warning:
            parse_warnings.append({"reaction": raw, "cell": row.get("cell"), "warning": warning, "formula": row.get("formula", "")})
        reaction = Reaction(rid, name=row.get("name", raw))
        lb = num_float(row.get("lb", ""), 0.0)
        ub = num_float(row.get("ub", ""), 0.0)
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        if stoich:
            for met_raw, coefficient in stoich.items():
                if met_raw not in met_map:
                    match = SPECIES_RE.match(met_raw)
                    comp = match.group(2) if match else "c"
                    met = Metabolite(sanitize(met_raw, "M_", used_met), name=met_raw, compartment=comp)
                    model.add_metabolites([met])
                    met_map[met_raw] = met
                reaction.add_metabolites({met_map[met_raw]: coefficient})
        gpr = row.get("gpr", "")
        if gpr:
            try:
                reaction.gene_reaction_rule = gpr
            except Exception as exc:
                gpr_warnings.append({"reaction": raw, "cell": row.get("gpr_cell"), "gpr": gpr, "warning": repr(exc)})
        model.add_reactions([reaction])

    objective_raw = "Ex_bio[e]" if "Ex_bio[e]" in reaction_map else "Afe_biomass_mc507_WT_139p0M"
    if objective_raw in reaction_map:
        model.objective = reaction_map[objective_raw]
        model.objective_direction = "max"
    # Windows 下 libSBML.writeSBMLToFile 对中文路径可能静默不落盘。
    # 先序列化为字符串，再由 Python 按 UTF-8 写入，并显式检查结果。
    buffer = io.StringIO()
    write_sbml_model(model, buffer)
    xml_text = buffer.getvalue()
    if not xml_text.strip():
        raise IOError(f"SBML输出为空：{out_xml}")
    out_xml.write_text(xml_text, encoding="utf-8", newline="\n")
    if not out_xml.is_file() or out_xml.stat().st_size == 0:
        raise IOError(f"SBML输出失败：{out_xml}")
    write_json(output_dir / "2024_ID映射.json", {"反应原始ID到SBML ID": reaction_map, "代谢物原始ID到SBML ID": {k: v.id for k, v in met_map.items()}})
    return model, reaction_map, met_map, {"reaction_formula_warnings": parse_warnings, "gpr_warnings": gpr_warnings}


def formula_diff(a: dict | None, b: dict | None) -> str:
    fields = []
    if a is None:
        return "2024新增"
    if b is None:
        return "2016存在而2024缺失"
    if canonical_formula(a.get("formula", "")) != canonical_formula(b.get("formula", "")):
        fields.append("化学计量/方向")
    if norm_text(a.get("name", "")) != norm_text(b.get("name", "")):
        fields.append("名称")
    if norm_gpr(a.get("gpr", "")) != norm_gpr(b.get("gpr", "")):
        fields.append("GPR")
    if numeric(a.get("lb", "")) != numeric(b.get("lb", "")) or numeric(a.get("ub", "")) != numeric(b.get("ub", "")):
        fields.append("上下界")
    if norm_text(a.get("subsystem", "")) != norm_text(b.get("subsystem", "")):
        fields.append("子系统")
    if a.get("ec", "") != b.get("ec", ""):
        fields.append("EC")
    if clean(a.get("reversible", "")) and clean(a.get("reversible", "")) != ("1" if parse_formula(a.get("formula", ""))[0] == "可逆" else "0"):
        fields.append("可逆性字段")
    return ";".join(fields) if fields else "无语义差异"


def make_diff(x16, x24):
    headers = ["反应ID", "状态", "差异字段", "2016反应式", "2024反应式", "2016规范化反应式", "2024规范化反应式", "2016上下界", "2024上下界", "2016GPR", "2024GPR", "2016单元格", "2024单元格", "层级提示"]
    rows = []
    for reaction_id in list(x16) + [key for key in x24 if key not in x16]:
        a, b = x16.get(reaction_id), x24.get(reaction_id)
        diff = formula_diff(a, b)
        if a is None:
            status = "2024新增"
            hint = classify_candidate(reaction_id)
        elif b is None:
            status = "2016独有"
            hint = "保留待核"
        elif diff == "无语义差异":
            status = "共同且无语义差异"
            hint = "共同基础"
        else:
            status = "共同且有语义差异"
            hint = classify_candidate(reaction_id)
        rows.append({
            "反应ID": reaction_id,
            "状态": status,
            "差异字段": diff,
            "2016反应式": a.get("formula", "") if a else "",
            "2024反应式": b.get("formula", "") if b else "",
            "2016规范化反应式": canonical_formula(a.get("formula", "")) if a else "",
            "2024规范化反应式": canonical_formula(b.get("formula", "")) if b else "",
            "2016上下界": f"{a.get('lb', '')}..{a.get('ub', '')}" if a else "",
            "2024上下界": f"{b.get('lb', '')}..{b.get('ub', '')}" if b else "",
            "2016GPR": a.get("gpr", "") if a else "",
            "2024GPR": b.get("gpr", "") if b else "",
            "2016单元格": a.get("cell", "") if a else "",
            "2024单元格": b.get("cell", "") if b else "",
            "层级提示": hint,
        })
    return headers, rows


def classify_candidate(reaction_id: str) -> str:
    if reaction_id.startswith("ZZ_Ecto"):
        return "异源工程：ectoine；不能自动视为ATCC23270野生型固有能力"
    if reaction_id.startswith("ZZ_Tre") or reaction_id.startswith("ZZ_glc") or reaction_id.startswith("ZZ_g1p") or reaction_id == "ZZ_glycogen_EX":
        return "工程/假设层：兼容溶质、转运或模拟交换；不能自动视为野生型固有能力"
    return "宿主本底更新候选：需保留2016证据并经后续人工审计"


def gpr_rows(x16, x24):
    headers = ["反应ID", "2016 XML GPR", "2016 Excel GPR", "2024 S1 GPR", "恢复GPR", "恢复来源", "2016 GPR单元格", "2024 GPR单元格", "原始逻辑保留说明"]
    rows = []
    for reaction_id in list(x16) + [key for key in x24 if key not in x16]:
        a, b = x16.get(reaction_id), x24.get(reaction_id)
        gpr16 = a.get("gpr", "") if a else ""
        gpr24 = b.get("gpr", "") if b else ""
        if gpr24:
            restored, source = gpr24, "2024 S1 D列"
        elif gpr16:
            restored, source = gpr16, "2016 mmc1.xls H列"
        else:
            restored, source = "", "无Excel证据"
        rows.append({
            "反应ID": reaction_id,
            "2016 XML GPR": "缺失（逐反应检查）",
            "2016 Excel GPR": gpr16,
            "2024 S1 GPR": gpr24,
            "恢复GPR": restored,
            "恢复来源": source,
            "2016 GPR单元格": a.get("gpr_cell", "") if a else "",
            "2024 GPR单元格": b.get("gpr_cell", "") if b else "",
            "原始逻辑保留说明": "仅在工作副本使用；不改2016原始XML；空白不补写" if restored else "保持空白；不能从无证据推断",
        })
    return headers, rows


def scenario_data(sheets, reactions):
    groups = OrderedDict([
        ("TreYZ", {"start": 3, "result_cols": (5, 6, 7), "product": "ZZ_Tre_EX", "inactivate": ["ZZ_TreT1_ADP", "ZZ_TreT2_UDP"]}),
        ("TreT-ADP", {"start": 8, "result_cols": (10, 11, 12), "product": "ZZ_Tre_EX", "inactivate": ["ZZ_TreT2_UDP", "ZZ_TreYZ"]}),
        ("TreT-UDP", {"start": 13, "result_cols": (15, 16, 17), "product": "ZZ_Tre_EX", "inactivate": ["ZZ_TreT1_ADP", "ZZ_TreYZ"]}),
        ("Ectoine", {"start": 18, "result_cols": (20, 21, 22), "product": "ZZ_Ecto4_EX", "inactivate": []}),
    ])
    output = []
    table_info = {}
    for sheet_name, rows in sheets.items():
        data_rows = {clean(row.get(1, "")): row for row in rows[2:] if clean(row.get(1, ""))}
        n2 = data_rows.get("Ex_n2[e]", {})
        nh4 = data_rows.get("Ex_nh4[e]", {})
        n2_open = num_float(n2.get(3, ""), 0) < 0 and num_float(n2.get(4, ""), 0) > 0
        nh4_open = num_float(nh4.get(3, ""), 0) < 0 and num_float(nh4.get(4, ""), 0) > 0
        nitrogen = "氮气（N2）" if n2_open and not nh4_open else "氨（NH4+）" if nh4_open and not n2_open else "未能从S3/S4边界唯一判断"
        table_info[sheet_name] = {
            "氮源判定": nitrogen,
            "证据": {"Ex_n2[e]": [n2.get(3, ""), n2.get(4, "")], "Ex_nh4[e]": [nh4.get(3, ""), nh4.get(4, "")]},
            "判定性质": "根据对应表的边界列判定；S3/S4本身没有显式氮源列",
        }
        for pathway, meta in groups.items():
            base = meta["start"]
            lb_col, ub_col = base, base + 1
            biomass_col, product_col, constrained_col = meta["result_cols"]
            for objective_name, result_col, constraint in (("最大生长", biomass_col, "无"), ("产物最大化", product_col, "无"), ("产物最大化_30%生长", constrained_col, "Ex_bio[e]下界=最大生长值×0.30")):
                expected = {}
                for raw_id, row in data_rows.items():
                    expected[raw_id] = {
                        "lb": clean(row.get(lb_col, "")),
                        "ub": clean(row.get(ub_col, "")),
                        "flux": clean(row.get(result_col, "")),
                        "cell": f"{chr(64 + result_col)}{row['__row__']}",
                    }
                output.append({
                    "场景ID": f"{sheet_name}_{pathway}_{objective_name}",
                    "来源表": sheet_name,
                    "氮源": nitrogen,
                    "通路": pathway,
                    "优化目标": "Ex_bio[e]" if objective_name == "最大生长" else meta["product"],
                    "生长约束": constraint,
                    "通路禁用反应": ";".join(meta["inactivate"]) or "无（按论文提供信息）",
                    "边界列": f"{lb_col}/{ub_col}",
                    "结果列": str(result_col),
                    "预期通量": expected,
                })
    return output, table_info


def json_scenario_config(scenarios, table_info, reactions, generated_model_stats):
    return {
        "文件性质": "2024论文模型复原工作副本的场景定义；工程和假设层不等同于野生型固有能力",
        "来源": {
            "论文": "Simulating compatible solute biosynthesis using a metabolic flux model of the biomining acidophile, Acidithiobacillus ferrooxidans ATCC 23270",
            "原始文件": "01_原始资料/2024_Khaleque/from2024-mmc1.xlsx",
            "论文证据页": [3, 5, 6, 7],
            "检索日期": "2026-09-22 Asia/Shanghai（项目内原始资料）",
        },
        "模型统计": generated_model_stats,
        "通量单位": {"反应": "mmol/gDCW/h", "细胞生物量交换": "h^-1", "论文产率表述": "mol/mol incoming carbon（论文正文）"},
        "S1_目标解释": {
            "Objective列实际": "631条S1反应均为空白",
            "执行规则": "不把空白解释为0；根据S3/S4结果列和论文正文重建目标",
            "最大生长": "目标 Ex_bio[e]",
            "产物最大化": "Trehalose 使用 ZZ_Tre_EX；Ectoine 使用 ZZ_Ecto4_EX",
            "30%生长约束": "先求最大 Ex_bio[e]，再将 Ex_bio[e] 下界设为最大生长值×0.30，随后最大化产物交换",
        },
        "氮源判定": table_info,
        "通路模式": {
            "TreYZ": "禁用 ZZ_TreT1_ADP 与 ZZ_TreT2_UDP；其余边界采用 S3/S4 对应列",
            "TreT-ADP": "禁用 ZZ_TreT2_UDP 与 ZZ_TreYZ；其余边界采用 S3/S4 对应列",
            "TreT-UDP": "禁用 ZZ_TreT1_ADP 与 ZZ_TreYZ；其余边界采用 S3/S4 对应列",
            "Ectoine": "使用 ZZ_Ecto1/2/3 与 ZZ_Ecto4_EX；论文提供的S3/S4未给出额外关闭Trehalose通路的文字证据，因此不擅自添加关闭约束",
        },
        "场景": scenarios,
        "裁决边界": [
            "ZZ_Ecto1/2/3/4_EX属于异源工程层，不能自动写入ATCC23270野生型本底",
            "ZZ_Tre*、ZZ_glc*、ZZ_g1p*及ZZ_glycogen_EX属于论文工程/假设/模拟交换层，不能单独作为菌株固有能力证据",
            "S3判定为氮气、S4判定为氨是对边界列的证据驱动判读，不是对S1 Objective空白的填充",
        ],
    }


def load_generated_cobra(path: Path):
    from cobra.io import read_sbml_model
    return read_sbml_model(str(path))


def set_bounds_safely(reaction, lower: float, upper: float):
    """避免固定负通量边界在逐字段更新中短暂变成非法区间。"""
    if lower > upper:
        raise ValueError(f"边界本身非法：{lower} > {upper}")
    if lower > reaction.upper_bound:
        reaction.upper_bound = upper
        reaction.lower_bound = lower
    else:
        reaction.lower_bound = lower
        reaction.upper_bound = upper


def run_scenarios(model, scenario_specs, x24, output_path: Path):
    rows = []
    solver_note = ""
    try:
        solver_note = f"solver={model.solver.interface.__name__ if hasattr(model.solver, 'interface') else model.solver}"
    except Exception:
        solver_note = "solver=未能读取"
    for spec in scenario_specs:
        scenario = model.copy()
        expected = spec["预期通量"]
        errors = []
        try:
            for raw_id, values in expected.items():
                if raw_id not in x24:
                    continue
                rid = next((r.id for r in scenario.reactions if r.name == x24[raw_id].get("name", "") and r.id.endswith(re.sub(r"[^A-Za-z0-9_]", "_", raw_id))), None)
                # 直接使用保存的映射更稳妥；备用匹配仅用于旧版脚本生成的模型。
                if rid is None:
                    rid = "R_" + re.sub(r"[^A-Za-z0-9_]", "_", raw_id)
                if rid not in scenario.reactions:
                    continue
                rxn = scenario.reactions.get_by_id(rid)
                set_bounds_safely(rxn, num_float(values["lb"], 0.0), num_float(values["ub"], 0.0))
            if spec["优化目标"] not in x24:
                raise RuntimeError("目标反应不在2024 S1：" + spec["优化目标"])
            target_id = "R_" + re.sub(r"[^A-Za-z0-9_]", "_", spec["优化目标"])
            if target_id not in scenario.reactions:
                raise RuntimeError("目标反应SBML ID不存在：" + target_id)
            if spec["优化目标"] != "Ex_bio[e]":
                biomass_id = "R_" + re.sub(r"[^A-Za-z0-9_]", "_", "Ex_bio[e]")
                if biomass_id not in scenario.reactions:
                    raise RuntimeError("缺少 Ex_bio[e]，无法应用30%生长约束")
                biomass = scenario.reactions.get_by_id(biomass_id)
                if spec["生长约束"].startswith("Ex_bio"):
                    scenario.objective = biomass
                    max_growth = scenario.optimize().objective_value
                    set_bounds_safely(biomass, 0.30 * max_growth, biomass.upper_bound)
            scenario.objective = scenario.reactions.get_by_id(target_id)
            solution = scenario.optimize()
            status = str(solution.status)
            flux = solution.fluxes
            max_residual = 0.0
            for metabolite in scenario.metabolites:
                residual = 0.0
                for reaction in metabolite.reactions:
                    residual += reaction.metabolites[metabolite] * float(flux[reaction.id])
                max_residual = max(max_residual, abs(residual))
            expected_target = num_float(expected.get(spec["优化目标"], {}).get("flux", ""))
            observed_target = float(solution.objective_value) if solution.objective_value is not None else None
            abs_error = abs(observed_target - expected_target) if expected_target is not None and observed_target is not None else None
            key_ids = ["Ex_bio[e]", spec["优化目标"], "Ex_fe2[e]", "Ex_n2[e]", "Ex_nh4[e]", "Ex_o2[e]"]
            key_parts = []
            for raw_id in key_ids:
                if raw_id not in expected:
                    continue
                rid = "R_" + re.sub(r"[^A-Za-z0-9_]", "_", raw_id)
                if rid not in flux.index:
                    continue
                observed = float(flux[rid])
                published = num_float(expected[raw_id].get("flux", ""))
                delta = observed - published if published is not None else None
                key_parts.append(f"{raw_id}:文献{expected[raw_id].get('flux','')}|复现{observed:.12g}|差{delta:.6g}")
            passed = status == "optimal" and abs_error is not None and abs_error <= 1e-7 and max_residual <= 1e-8
            rows.append({
                "场景ID": spec["场景ID"], "来源表": spec["来源表"], "氮源": spec["氮源"], "通路": spec["通路"],
                "优化目标": spec["优化目标"], "生长约束": spec["生长约束"], "求解状态": status,
                "文献目标值": expected_target, "复现目标值": observed_target, "绝对误差": abs_error,
                "关键交换通量对照": "; ".join(key_parts), "Sv最大残差": max_residual,
                "验收状态": "通过" if passed else "未通过：目标值或残差超出容差",
                "求解器诊断": solver_note,
            })
        except Exception as exc:
            rows.append({
                "场景ID": spec["场景ID"], "来源表": spec["来源表"], "氮源": spec["氮源"], "通路": spec["通路"],
                "优化目标": spec["优化目标"], "生长约束": spec["生长约束"], "求解状态": "未执行",
                "文献目标值": "", "复现目标值": "", "绝对误差": "", "关键交换通量对照": "",
                "Sv最大残差": "", "验收状态": "未通过：" + repr(exc), "求解器诊断": solver_note,
            })
    headers = ["场景ID", "来源表", "氮源", "通路", "优化目标", "生长约束", "求解状态", "文献目标值", "复现目标值", "绝对误差", "关键交换通量对照", "Sv最大残差", "验收状态", "求解器诊断"]
    write_tsv(output_path, headers, rows)
    return headers, rows


def raw_manifest(root: Path):
    manifest = root / "00_项目总控/原始资料SHA256清单.tsv"
    rows = []
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = root / row["相对路径"]
            observed = sha256(path) if path.is_file() else "缺失"
            rows.append({"相对路径": row["相对路径"], "清单字节数": row["字节数"], "实际字节数": path.stat().st_size if path.is_file() else "", "清单SHA256": row["SHA256"], "实际SHA256": observed, "状态": "通过" if observed.lower() == row["SHA256"].lower() else "未通过"})
    return rows


def main():
    parser = argparse.ArgumentParser(description="A0冻结基线、真实差异、2024论文模型复原和FBA诊断")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / "08_模型基线/A0"
    report = root / "12_验收报告/A0"
    env = root / "11_脚本与环境/A0"
    out.mkdir(parents=True, exist_ok=True)
    report.mkdir(parents=True, exist_ok=True)
    env.mkdir(parents=True, exist_ok=True)
    intermediate = env / "中间"
    intermediate.mkdir(parents=True, exist_ok=True)

    hashes = raw_manifest(root)
    x16, x16_info = xls_table(root)
    x24, metabolites, scenario_sheets, x24_info = xlsx_tables(root)
    original_xml = root / "01_原始资料/2016_Campodonico/mmc3.xml"
    original_stats = sbml_stats(original_xml)
    original_import = sbml_warnings(original_xml)
    generated_xml = out / "2024_论文复原.xml"
    model, reaction_map, met_map, build_warnings = build_cobra_model(x24, metabolites, generated_xml, out)
    generated_stats = sbml_warnings(generated_xml)
    generated_cobra_stats = {"reaction_count": len(model.reactions), "metabolite_count": len(model.metabolites), "gene_count": len(model.genes), "objective": str(model.objective.expression), "default_objective": "Ex_bio[e]（S1 Objective为空，按S3/S4和论文恢复）"}
    diff_headers, diff_rows = make_diff(x16, x24)
    write_tsv(out / "2016_到2024_真实差异.tsv", diff_headers, diff_rows)
    gpr_headers, gpr_data = gpr_rows(x16, x24)
    write_tsv(out / "恢复GPR及原始证据定位.tsv", gpr_headers, gpr_data)
    candidate_headers = ["反应ID", "候选层级", "2016单元格", "2024单元格", "差异字段", "2024反应式", "2024GPR", "处理结论"]
    candidate_rows = []
    for row in diff_rows:
        if row["状态"] in ("共同且无语义差异",):
            continue
        candidate_rows.append({"反应ID": row["反应ID"], "候选层级": row["层级提示"], "2016单元格": row["2016单元格"], "2024单元格": row["2024单元格"], "差异字段": row["差异字段"], "2024反应式": row["2024反应式"], "2024GPR": row["2024GPR"], "处理结论": "保留为候选；本轮不自动纳入ATCC23270野生型本底"})
    write_tsv(out / "宿主本底更新候选清单.tsv", candidate_headers, candidate_rows)

    scenario_specs, table_info = scenario_data(scenario_sheets, x24)
    config = json_scenario_config(scenario_specs, table_info, x24, generated_cobra_stats)
    write_json(out / "2024_情景配置.json", config)
    try:
        _, comparison_rows = run_scenarios(model, scenario_specs, x24, out / "复现比较.tsv")
    except Exception as exc:
        comparison_rows = [{"验收状态": "未通过：FBA总执行异常：" + repr(exc), "求解器诊断": traceback.format_exc()}]
        write_tsv(out / "复现比较.tsv", ["验收状态", "求解器诊断"], comparison_rows)

    snapshot = out / "2016原版只读快照"
    snapshot.mkdir(parents=True, exist_ok=True)
    snapshot_files = [
        root / "01_原始资料/2016_Campodonico/mmc3.xml",
        root / "01_原始资料/2016_Campodonico/mmc1.xls",
    ]
    snapshot_rows = []
    for source in snapshot_files:
        destination = snapshot / source.name
        if destination.exists():
            # 仅处理本脚本先前生成的A0快照，允许幂等重跑。
            try:
                os.chmod(destination, 0o666)
            except OSError:
                pass
        shutil.copy2(source, destination)
        try:
            os.chmod(destination, 0o444)
        except OSError:
            pass
        snapshot_rows.append({"快照文件": destination.name, "来源": str(source.relative_to(root)), "SHA256": sha256(destination), "状态": "只读意图；原始文件未改动"})
    write_tsv(snapshot / "快照SHA256.tsv", ["快照文件", "来源", "SHA256", "状态"], snapshot_rows)

    passed_hash = all(row["状态"] == "通过" for row in hashes)
    passed_parse = generated_stats.get("available") and generated_stats.get("error_count") == 0 and not build_warnings["reaction_formula_warnings"]
    passed_fba = bool(comparison_rows) and all(row.get("验收状态") == "通过" for row in comparison_rows)
    baseline = {
        "任务": "A0_冻结基线与2024论文模型复原",
        "agent_thread标识": AGENT_THREAD,
        "执行类型": "Codex",
        "目标模型": EXPECTED_MODEL,
        "目标推理等级": EXPECTED_REASONING,
        "实际模型": ACTUAL_MODEL,
        "实际reasoning等级": ACTUAL_REASONING,
        "开始时间": START_TIME,
        "完成记录时间": now_iso(),
        "原始资料校验": hashes,
        "原始资料冻结状态": "通过" if passed_hash else "未通过",
        "2016_XML_SBML与导入": {"统计": original_stats, "导入检查": original_import, "Excel统计": x16_info},
        "2024补充表统计": x24_info,
        "2024模型工作副本": {"SBML导入检查": generated_stats, "COBRA统计": generated_cobra_stats, "构建警告": build_warnings},
        "差异统计": {"总行数": len(diff_rows), "共同无语义差异": sum(r["状态"] == "共同且无语义差异" for r in diff_rows), "共同有语义差异": sum(r["状态"] == "共同且有语义差异" for r in diff_rows), "2024新增": sum(r["状态"] == "2024新增" for r in diff_rows), "2016独有": sum(r["状态"] == "2016独有" for r in diff_rows)},
        "GPR恢复": {"2016 XML逐反应检查": "缺失", "2016 Excel非空": sum(bool(r.get("gpr")) for r in x16.values()), "2024 S1非空": sum(bool(r.get("gpr")) for r in x24.values()), "唯一基因数（2024 GPR粗计）": len({g for row in x24.values() for g in GENE_RE.findall(row.get("gpr", "")) if g.lower() not in {"and", "or"}})},
        "场景复现": {"场景数": len(comparison_rows), "通过数": sum(r.get("验收状态") == "通过" for r in comparison_rows), "未通过数": sum(r.get("验收状态") != "通过" for r in comparison_rows)},
        "A0状态": {"原始资料已冻结": "通过" if passed_hash else "未通过", "2024基线已复现": "通过" if passed_fba else "未通过", "模型可解析": "通过" if passed_parse else "未通过"},
    }
    write_json(out / "2016_基线登记.json", baseline)

    unresolved = [
        "本轮没有修改01_原始资料和00_项目总控共享表。",
        "2016 mmc3.xml为SBML Level 2 Version 1；逐反应检查未发现GPR标签，反应级GPR仅能从2016 mmc1.xls H列及2024 S1 D列恢复。",
        "2016 XML没有FBC Objective节点，但有一个反应参数 OBJECTIVE_COEFFICIENT=1，对应 biomass exchange；工作副本按论文S3/S4将 Ex_bio[e] 作为默认最大生长目标。",
        "2024 S1 Objective列631条均为空；本轮根据论文第3、5—7页和S3/S4结果列恢复最大生长、产物最大化、30%生长约束。",
        "S3通过 Ex_n2[e] 开放、Ex_nh4[e] 关闭判定为氮气；S4通过 Ex_nh4[e] 开放、Ex_n2[e] 关闭判定为氨。此判定是边界证据推断，已写入2024_情景配置.json。",
        "异源 ectoine 反应、Trehalose路径、葡萄糖/糖原相关反应及模拟专用交换均保留为工程/假设层，未自动命名为已验证野生型。",
        "模型重建使用2024 S1/S2直接构建，SBML ID经过合法化映射；原始ID与SBML ID见2024_ID映射.json。",
        "复现比较.tsv中的每个场景均需同时检查目标绝对误差≤1e-7和Sv最大残差≤1e-8；若有求解器或公式解析警告，状态以未通过为准。",
    ]
    write_text = out / "未解决问题.md"
    write_text.write_text("# A0未解决问题与边界\n\n" + "\n".join(f"- {item}" for item in unresolved) + "\n", encoding="utf-8")

    (env / "requirements.txt").write_text("cobra==0.32.1\npython-libsbml==5.21.2\nxlrd==2.0.2\nopenpyxl==3.1.5\npypdf==6.19.0\n", encoding="utf-8")
    (env / "运行说明.md").write_text(
        "# A0运行说明\n\n"
        "项目隔离环境位于 `11_脚本与环境/A0/.venv`。\n\n"
        "```powershell\n"
        ".\\11_脚本与环境\\A0\\.venv\\Scripts\\python.exe .\\11_脚本与环境\\A0\\冻结基线.py --root .\n"
        "```\n\n"
        "脚本只写 `08_模型基线/A0`、`11_脚本与环境/A0` 和 `12_验收报告/A0`；原始资料和00共享表仅读。\n",
        encoding="utf-8",
    )

    report_lines = [
        "# A0验收报告",
        "",
        f"完成记录时间：{now_iso()}",
        f"agent/thread：{AGENT_THREAD}",
        f"实际模型：{ACTUAL_MODEL}",
        f"实际reasoning等级：{ACTUAL_REASONING}",
        f"目标配置：{EXPECTED_MODEL} / {EXPECTED_REASONING}",
        "",
        "## 状态",
        "",
        f"- 原始资料已冻结：{'通过' if passed_hash else '未通过'}（7项逐项SHA256）",
        f"- 2016原始SBML解析：{'通过' if original_import.get('error_count') == 0 else '见登记中的导入警告'}",
        f"- 2024工作副本SBML解析：{'通过' if passed_parse else '未通过'}",
        f"- 2024场景复现：{'通过' if passed_fba else '未通过'}；通过 {sum(r.get('验收状态') == '通过' for r in comparison_rows)}/{len(comparison_rows)}",
        "",
        "## 核心事实",
        "",
        f"- 2016 XML：SBML Level {original_stats.get('sbml_level')} Version {original_stats.get('sbml_version')}，{original_stats.get('compartment_count')}区室、{original_stats.get('species_count')}代谢物、{original_stats.get('reaction_count')}反应、GPR标签 {original_stats.get('gene_product_association_tag_count')}，非零旧式目标系数 {len(original_stats.get('legacy_objective_coefficients', []))}。",
        f"- 2016 MMC1：{x16_info['data_rows']}条反应；2024 S1：{x24_info['s1_data_rows']}条反应；2024 S2：{x24_info['s2_data_rows']}个非空代谢物ID。",
        f"- 真实差异行：{len(diff_rows)}；共同无语义差异 {sum(r['状态'] == '共同且无语义差异' for r in diff_rows)}；共同有语义差异 {sum(r['状态'] == '共同且有语义差异' for r in diff_rows)}；2024新增 {sum(r['状态'] == '2024新增' for r in diff_rows)}。",
        "- S3按边界恢复为氮气场景，S4按边界恢复为氨场景；两组均分别包含最大生长、产物最大化和30%生长约束。",
        "",
        "## 输出",
        "",
        "详见 `08_模型基线/A0/2016_基线登记.json`、`2016_到2024_真实差异.tsv`、`2024_论文复原.xml`、`2024_情景配置.json`、`复现比较.tsv`、`恢复GPR及原始证据定位.tsv`、`宿主本底更新候选清单.tsv`、`未解决问题.md`。",
        "",
        "## 阻塞点",
        "",
        "若复现比较存在未通过项，以该表中的求解状态、目标误差、Sv最大残差和求解器诊断为准；没有通过的场景不得被称为已验证野生型。",
    ]
    (report / "A0验收报告.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    write_json(report / "A0_完成记录.json", {
        "任务": "A0_冻结基线与2024论文模型复原",
        "agent/thread标识": AGENT_THREAD,
        "实际模型": ACTUAL_MODEL,
        "实际reasoning等级": ACTUAL_REASONING,
        "目标模型": EXPECTED_MODEL,
        "目标reasoning等级": EXPECTED_REASONING,
        "开始时间": START_TIME,
        "完成时间": now_iso(),
        "写入范围": ["08_模型基线/A0", "11_脚本与环境/A0", "12_验收报告/A0"],
        "未写入": ["01_原始资料", "00_项目总控共享表", "B1", "C2"],
        "验收状态": baseline["A0状态"],
        "实际模型和推理等级说明": "当前工具未返回实际模型和reasoning字段，已按要求标记未核实；未以目标配置冒充实际配置。",
    })
    print(json.dumps({"原始资料已冻结": passed_hash, "模型可解析": passed_parse, "2024基线已复现": passed_fba, "场景数": len(comparison_rows), "通过场景数": sum(r.get('验收状态') == '通过' for r in comparison_rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
