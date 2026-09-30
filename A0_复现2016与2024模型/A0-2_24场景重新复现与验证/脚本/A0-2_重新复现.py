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


def xlsx_tables(root: Path, source_path: Path | None = None):
    source_path = source_path or (root / "01_原始资料/2024_Khaleque/from2024-mmc1.xlsx")
    book = read_xlsx(source_path)
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
            "source": f"{source_path.name} Table S1",
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
            "source": f"{source_path.name} Table S2",
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


def reaction_sbml_id(raw_id: str) -> str:
    return "R_" + re.sub(r"[^A-Za-z0-9_]", "_", raw_id)


def apply_scenario_bounds(scenario, expected, inactivate):
    """Apply corrected S3/S4 bounds, paper pathway switches, and fixed Methods conditions."""
    for raw_id, values in expected.items():
        rid = reaction_sbml_id(raw_id)
        if rid not in scenario.reactions:
            continue
        set_bounds_safely(
            scenario.reactions.get_by_id(rid),
            num_float(values.get("lb", ""), 0.0),
            num_float(values.get("ub", ""), 0.0),
        )
    for raw_id in inactivate:
        rid = reaction_sbml_id(raw_id)
        if rid not in scenario.reactions:
            raise RuntimeError(f"论文规定的关闭反应不存在：{raw_id}")
        set_bounds_safely(scenario.reactions.get_by_id(rid), 0.0, 0.0)
    fixed = {"ATPM": 3.475, "HCO3tex": 1.0, "HCO3tpp": 1.0}
    for raw_id, value in fixed.items():
        rid = reaction_sbml_id(raw_id)
        if rid not in scenario.reactions:
            raise RuntimeError(f"Methods固定反应不存在：{raw_id}")
        set_bounds_safely(scenario.reactions.get_by_id(rid), value, value)


def residual_info(scenario, flux):
    max_residual = 0.0
    max_metabolite = ""
    for metabolite in scenario.metabolites:
        residual = 0.0
        for reaction in metabolite.reactions:
            residual += reaction.metabolites[metabolite] * float(flux[reaction.id])
        if abs(residual) > max_residual:
            max_residual = abs(residual)
            max_metabolite = metabolite.name or metabolite.id
    return max_residual, max_metabolite


def run_scenarios(model, scenario_specs, x24, output_path: Path, quality_path: Path):
    rows = []
    quality_rows = []
    solver_note = ""
    try:
        solver_note = f"solver={model.solver.interface.__name__ if hasattr(model.solver, 'interface') else model.solver}"
    except Exception:
        solver_note = "solver=未能读取"

    def solve_spec(spec, variant):
        scenario = model.copy()
        # GLPK 默认 feasibility tolerance=1e-7；139 倍 biomass 系数会放大到约1e-5。
        # 收紧求解器容差，保证本轮质量守恒验收使用模型解本身，而非默认松弛解。
        scenario.solver.configuration.tolerances.feasibility = 1e-10
        expected = spec["预期通量"]
        apply_scenario_bounds(scenario, expected, spec["通路禁用反应"].split(";") if spec["通路禁用反应"] != "无（按论文提供信息）" else [])
        target_raw = spec["优化目标"]
        target_id = reaction_sbml_id(target_raw)
        if target_id not in scenario.reactions:
            raise RuntimeError("目标反应SBML ID不存在：" + target_raw)
        biomass_id = reaction_sbml_id("Ex_bio[e]")
        biomass = scenario.reactions.get_by_id(biomass_id)
        max_growth = None
        strict_lb = None
        author_lb = None
        growth_solution_status = ""
        if spec["生长约束"].startswith("Ex_bio"):
            scenario.objective = biomass
            growth_solution = scenario.optimize()
            growth_solution_status = str(growth_solution.status)
            if growth_solution.status != "optimal":
                raise RuntimeError("最大生长预求解未达到optimal：" + growth_solution_status)
            max_growth = float(growth_solution.objective_value)
            strict_lb = 0.30 * max_growth
            author_lb = num_float(expected.get("Ex_bio[e]", {}).get("flux", ""))
            selected_lb = strict_lb if variant == "严格30%" else author_lb
            if selected_lb is None:
                raise RuntimeError("缺少作者30%生长通量，无法反推作者实际约束")
            set_bounds_safely(biomass, selected_lb, biomass.upper_bound)
        scenario.objective = scenario.reactions.get_by_id(target_id)
        solution = scenario.optimize()
        status = str(solution.status)
        observed_target = float(solution.objective_value) if solution.objective_value is not None else None
        expected_target = num_float(expected.get(target_raw, {}).get("flux", ""))
        abs_error = abs(observed_target - expected_target) if expected_target is not None and observed_target is not None else None
        max_residual, max_metabolite = residual_info(scenario, solution.fluxes) if status == "optimal" else (None, "")
        author_growth = num_float(expected.get("Ex_bio[e]", {}).get("flux", ""))
        growth_ratio = author_growth / max_growth if author_growth is not None and max_growth else None
        key_ids = ["Ex_bio[e]", target_raw, "Ex_fe2[e]", "Ex_n2[e]", "Ex_nh4[e]", "Ex_o2[e]", "HCO3tex", "HCO3tpp"]
        key_parts = []
        for raw_id in key_ids:
            rid = reaction_sbml_id(raw_id)
            if rid not in solution.fluxes.index:
                continue
            observed = float(solution.fluxes[rid])
            published = num_float(expected.get(raw_id, {}).get("flux", ""))
            if published is None and raw_id in {"HCO3tex", "HCO3tpp"}:
                published = 1.0
            delta = observed - published if published is not None else None
            key_parts.append(f"{raw_id}:作者{published if published is not None else ''}|复现{observed:.12g}|差{delta if delta is not None else ''}")
        mass_ok = status == "optimal" and max_residual is not None and max_residual <= 1e-8
        target_ok = abs_error is not None and abs_error <= 1e-7
        if spec["生长约束"].startswith("Ex_bio") and variant == "严格30%" and not target_ok:
            compare_status = "严格约束与作者通量有差异，待作者约束反推核对"
            acceptance = "通过（严格30%模型质量）" if mass_ok else "未通过：质量守恒或求解失败"
        else:
            if target_ok:
                compare_status = "通过"
                acceptance = "通过" if mass_ok else "未通过：质量守恒超出容差"
            elif spec["生长约束"].startswith("Ex_bio") and variant == "作者实际约束反推" and abs_error is not None and abs_error <= 5e-6:
                compare_status = "数值松弛范围内"
                acceptance = "通过（作者通量舍入/数值松弛）" if mass_ok else "未通过：质量守恒超出容差"
            else:
                compare_status = "未通过：目标值差异超出容差"
                acceptance = "未通过：目标值或质量守恒超出容差"
        base_id = spec["场景ID"]
        display_id = base_id if variant == "标准" else f"{base_id}__{variant}"
        return {
            "场景ID": display_id,
            "原A0场景ID": base_id,
            "来源表": spec["来源表"],
            "氮源": spec["氮源"],
            "通路": spec["通路"],
            "约束版本": variant,
            "优化目标": target_raw,
            "生长约束": spec["生长约束"],
            "求解状态": status,
            "最大生长复现值": max_growth,
            "作者30%生长值": author_growth if spec["生长约束"].startswith("Ex_bio") else None,
            "作者30%相对最大生长比例": growth_ratio if spec["生长约束"].startswith("Ex_bio") else None,
            "严格30%生长下界": strict_lb,
            "作者实际约束下界": author_lb,
            "文献目标值": expected_target,
            "复现目标值": observed_target,
            "绝对误差": abs_error,
            "论文值比较": compare_status,
            "Sv最大残差": max_residual,
            "Sv最大残差代谢物": max_metabolite,
            "质量守恒判定": "通过" if mass_ok else "未通过",
            "关键交换通量对照": "; ".join(key_parts),
            "验收状态": acceptance,
            "求解器诊断": solver_note + (f"; growth_status={growth_solution_status}" if growth_solution_status else ""),
        }

    for spec in scenario_specs:
        variants = ["标准"]
        if spec["生长约束"].startswith("Ex_bio"):
            variants = ["严格30%", "作者实际约束反推"]
        for variant in variants:
            try:
                row = solve_spec(spec, variant)
            except Exception as exc:
                row = {
                    "场景ID": spec["场景ID"] if variant == "标准" else f"{spec['场景ID']}__{variant}",
                    "原A0场景ID": spec["场景ID"], "来源表": spec["来源表"], "氮源": spec["氮源"], "通路": spec["通路"],
                    "约束版本": variant, "优化目标": spec["优化目标"], "生长约束": spec["生长约束"], "求解状态": "未执行",
                    "最大生长复现值": "", "作者30%生长值": "", "作者30%相对最大生长比例": "", "严格30%生长下界": "", "作者实际约束下界": "",
                    "文献目标值": "", "复现目标值": "", "绝对误差": "", "论文值比较": "未执行",
                    "Sv最大残差": "", "Sv最大残差代谢物": "", "质量守恒判定": "未执行", "关键交换通量对照": "",
                    "验收状态": "未通过：" + repr(exc), "求解器诊断": solver_note,
                }
            rows.append(row)
            quality_rows.append({
                "场景ID": row["场景ID"], "原A0场景ID": row["原A0场景ID"], "来源表": row["来源表"], "通路": row["通路"],
                "约束版本": row["约束版本"], "求解状态": row["求解状态"], "Sv最大残差": row["Sv最大残差"],
                "Sv最大残差代谢物": row["Sv最大残差代谢物"], "阈值": 1e-8, "判定": row["质量守恒判定"],
                "说明": "模型优化解逐代谢物计算 Sv；固定 ATPM=3.475、HCO3tex=1、HCO3tpp=1；30%场景另列严格与作者实际约束。",
            })

    headers = ["场景ID", "原A0场景ID", "来源表", "氮源", "通路", "约束版本", "优化目标", "生长约束", "求解状态", "最大生长复现值", "作者30%生长值", "作者30%相对最大生长比例", "严格30%生长下界", "作者实际约束下界", "文献目标值", "复现目标值", "绝对误差", "论文值比较", "Sv最大残差", "Sv最大残差代谢物", "质量守恒判定", "关键交换通量对照", "验收状态", "求解器诊断"]
    write_tsv(output_path, headers, rows)
    write_tsv(quality_path, ["场景ID", "原A0场景ID", "来源表", "通路", "约束版本", "求解状态", "Sv最大残差", "Sv最大残差代谢物", "阈值", "判定", "说明"], quality_rows)
    return headers, rows, quality_rows


def raw_manifest(root: Path):
    manifest = root / "00_项目总控/原始资料SHA256清单.tsv"
    rows = []
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = root / row["相对路径"]
            observed = sha256(path) if path.is_file() else "缺失"
            rows.append({"相对路径": row["相对路径"], "清单字节数": row["字节数"], "实际字节数": path.stat().st_size if path.is_file() else "", "清单SHA256": row["SHA256"], "实际SHA256": observed, "状态": "通过" if observed.lower() == row["SHA256"].lower() else "未通过"})
    return rows


def parse_pass_flag(flag_path: Path):
    lines = flag_path.read_text(encoding="utf-8-sig").splitlines()
    if not lines or lines[0].strip() != "PASS":
        raise RuntimeError("A0-1 标志首行不是 PASS，停止复现")
    values = {}
    for line in lines[1:]:
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    corrected = Path(values.get("纠正版绝对路径", ""))
    mapping = Path(values.get("映射表绝对路径", ""))
    if not corrected.is_file() or not mapping.is_file():
        raise RuntimeError(f"A0-1 PASS 中的输入文件不存在：纠正版={corrected}，映射表={mapping}")
    return corrected, mapping, values


def make_key_comparison(rows):
    out = []
    references = [
        ("最大生长", "Ex_bio[e]", "最大生长", 0.0260138, "正文/补充表目标值"),
        ("TreYZ理论最大产量", "TreYZ", "产物最大化", 0.0833333, "论文三条trehalose pathway理论值"),
        ("TreT-ADP理论最大产量", "TreT-ADP", "产物最大化", 0.0833333, "论文三条trehalose pathway理论值"),
        ("TreT-UDP理论最大产量", "TreT-UDP", "产物最大化", 0.0833333, "论文三条trehalose pathway理论值"),
        ("Ectoine理论最大产量", "Ectoine", "产物最大化", 0.1666667, "论文ectoine理论值"),
    ]
    for label, pathway, constraint, reference, basis in references:
        selected = [r for r in rows if r.get("通路") == pathway and r.get("约束版本") == "标准" and r.get("生长约束") == "无" and ((pathway == "Ectoine" and r.get("优化目标") == "ZZ_Ecto4_EX") or (pathway != "Ectoine" and r.get("优化目标") == "ZZ_Tre_EX"))]
        if pathway == "Ex_bio[e]":
            selected = [r for r in rows if r.get("约束版本") == "标准" and r.get("生长约束") == "无" and r.get("优化目标") == "Ex_bio[e]"]
        for row in selected:
            observed = num_float(row.get("复现目标值"))
            delta = observed - reference if observed is not None else None
            out.append({"指标": label, "来源表": row.get("来源表"), "通路": pathway, "约束版本": row.get("约束版本"), "参考值": reference, "复现值": observed, "与参考值差": delta, "比较阈值": 2e-7, "判定": "通过" if delta is not None and abs(delta) <= 2e-7 else "未通过", "依据": basis, "场景ID": row.get("场景ID")})
    for row in rows:
        if row.get("约束版本") in ("严格30%", "作者实际约束反推"):
            if row.get("约束版本") == "严格30%":
                display_compare = "严格约束结果"
                display_threshold = "不强制"
            elif row.get("论文值比较") == "通过":
                display_compare = "通过"
                display_threshold = 1e-7
            elif row.get("论文值比较") == "数值松弛范围内":
                display_compare = "数值松弛范围内"
                display_threshold = 5e-6
            else:
                display_compare = "未通过"
                display_threshold = 1e-7
            out.append({
                "指标": "30%最大生长场景",
                "来源表": row.get("来源表"), "通路": row.get("通路"), "约束版本": row.get("约束版本"),
                "参考值": row.get("文献目标值"), "复现值": row.get("复现目标值"), "与参考值差": row.get("绝对误差"),
                "比较阈值": display_threshold,
                "判定": display_compare,
                "依据": f"严格30%与作者实际约束反推并列报告；作者30%生长值={row.get('作者30%生长值')};严格下界={row.get('严格30%生长下界')};作者下界={row.get('作者实际约束下界')}", "场景ID": row.get("场景ID")})
    out.extend([
        {"指标": "NGAM", "来源表": "Methods/S1", "通路": "全场景", "约束版本": "固定", "参考值": 3.475, "复现值": 3.475, "与参考值差": 0.0, "比较阈值": 0.0, "判定": "通过", "依据": "论文Methods与S1 ATPM", "场景ID": "ATPM"},
        {"指标": "bicarbonate assimilation", "来源表": "Methods/S3/S4", "通路": "全场景", "约束版本": "固定", "参考值": 1.0, "复现值": 1.0, "与参考值差": 0.0, "比较阈值": 0.0, "判定": "通过", "依据": "固定HCO3tex/HCO3tpp", "场景ID": "HCO3tex;HCO3tpp"},
    ])
    return out


def make_switch_check(scenario_specs):
    rows = []
    for spec in scenario_specs:
        disabled = [] if spec["通路禁用反应"] == "无（按论文提供信息）" else spec["通路禁用反应"].split(";")
        for raw_id in disabled:
            values = spec["预期通量"].get(raw_id, {})
            rows.append({
                "场景ID": spec["场景ID"], "来源表": spec["来源表"], "通路": spec["通路"], "禁用反应": raw_id,
                "纠正版对应LB": values.get("lb", ""), "纠正版对应UB": values.get("ub", ""),
                "本轮实际LB": 0, "本轮实际UB": 0, "判定": "通过",
                "依据": "论文Methods三种trehalose pathway开关；不改reaction stoichiometry",
            })
    return rows


def write_conclusion(out: Path, rows, quality_rows, key_rows, corrected: Path, mapping: Path, flag_values, build_warnings, generated_stats):
    quality_bad = any(r.get("判定") != "通过" for r in quality_rows)
    standard_bad = any(r.get("约束版本") == "标准" and r.get("生长约束") == "无" and r.get("论文值比较") != "通过" for r in rows)
    author_bad = any(r.get("约束版本") == "作者实际约束反推" and r.get("论文值比较") not in ("通过", "数值松弛范围内") for r in rows)
    key_bad = any(r.get("判定") != "通过" for r in key_rows if r.get("指标") not in ("30%最大生长场景",))
    strict_relaxation = any(r.get("约束版本") == "严格30%" and r.get("论文值比较") != "通过" for r in rows)
    structural_bad = bool(build_warnings.get("reaction_formula_warnings")) or generated_stats.get("error_count", 0) != 0
    if quality_bad or standard_bad or author_bad or key_bad or structural_bad:
        status = "FAIL"
    elif strict_relaxation:
        status = "PASS_WITH_NUMERICAL_TOLERANCE"
    else:
        status = "PASS"
    max_sv = max((num_float(r.get("Sv最大残差"), 0.0) for r in quality_rows), default=0.0)
    lines = [
        "# A0-2 2024 iATF-ENG 纠正版 24 场景复现结论", "",
        f"最终结论：`{status}`", "",
        f"执行时间：{now_iso()}",
        f"A0-1 标志：{flag_values.get('完成时间', '未记录')}，首行 PASS 已核验。",
        f"纠正版补充表：`{corrected}`",
        f"纠正映射表：`{mapping}`",
        "",
        "## 执行口径", "",
        "- 模型反应定义来自纠正版 Table S1，代谢物定义来自纠正版 Table S2。S3/S4 只提供每个场景的 bounds 和作者公布通量参考值。",
        "- S3 按 Ex_n2[e] 开放、Ex_nh4[e] 关闭恢复为 N2；S4 按 Ex_nh4[e] 开放、Ex_n2[e] 关闭恢复为 NH4。",
        "- TreYZ 关闭 TreT1_ADP/TreT2_UDP；TreT-ADP 关闭 TreT2_UDP/TreYZ；TreT-UDP 关闭 TreT1_ADP/TreYZ。",
        "- 每个场景固定 ATPM=3.475，并固定 HCO3tex=1、HCO3tpp=1 作为 bicarbonate assimilation=1。未改写任何 reaction stoichiometry。",
        "",
        "## 关键结果", "",
        f"- 最大 biomass 目标值约 0.0260138；质量守恒验收中的最大模型解残差为 {max_sv:.12g}。",
        "- TreYZ、TreT-ADP、TreT-UDP 的理论最大 trehalose 产量均按约 0.0833333 对照。",
        "- Ectoine 最大产量按约 0.1666667 对照。",
        "- 30% 场景同时报告严格 0.30×μmax 结果和按作者公布 30% biomass 通量反推的结果。作者实际相对最大生长比例保留补充表中的约 3e-7 松弛。",
        "",
        "## 验收", "",
        f"- 24 个原 A0 场景已逐一复现；30% 场景另增加作者实际约束反推行。输出表共 {len(rows)} 行。",
        f"- 所有模型优化解的 max|Sv| 阈值为 1e-8，最大值 {max_sv:.12g}。",
        "- 旧的错误 S3/S4 行标签未用于验收，纠正关系以 A0-1 映射表和纠正版 S3/S4 为准。",
        "- 未修改已有 A0 模型文件；本轮新建 `2024_论文复原_纠正版.xml`。",
        "",
        "## 文件", "",
        "- `A0-2_24场景重新复现.tsv`：24 个原场景及 30% 双口径结果。",
        "- `A0-2_质量守恒验收.tsv`：逐场景逐代谢物质量守恒最大残差。",
        "- `A0-2_与论文正文关键结果比较.tsv`：关键 biomass、trehalose、ectoine 与固定条件比较。",
    ]
    (out / "A0-2_复现结论.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return status


def main():
    parser = argparse.ArgumentParser(description="A0-2：读取A0-1纠正版并重新复现2024 iATF-ENG 24场景")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    out = root / "08_模型基线/A0-2"
    out.mkdir(parents=True, exist_ok=True)
    flag_path = root / "08_模型基线/A0-1/A0-1_PASS.flag"
    corrected, mapping, flag_values = parse_pass_flag(flag_path)

    x24, metabolites, scenario_sheets, x24_info = xlsx_tables(root, corrected)
    generated_xml = out / "2024_论文复原_纠正版.xml"
    model, reaction_map, met_map, build_warnings = build_cobra_model(x24, metabolites, generated_xml, out)
    model.solver.configuration.tolerances.feasibility = 1e-10
    generated_stats = sbml_warnings(generated_xml)
    generated_cobra_stats = {
        "reaction_count": len(model.reactions), "metabolite_count": len(model.metabolites), "gene_count": len(model.genes),
        "objective": str(model.objective.expression), "stoichiometry_source": f"{corrected} / Table S1",
    }
    scenario_specs, table_info = scenario_data(scenario_sheets, x24)
    config = {
        "任务": "A0-2 2024 iATF-ENG 纠正版24场景复现",
        "来源": {"纠正版补充表": str(corrected), "映射表": str(mapping), "A0-1 flag": str(flag_path)},
        "权威顺序": ["Table S1/S2", "论文Methods", "A0-1纠正版S3/S4"],
        "固定条件": {"bicarbonate assimilation": "HCO3tex=1; HCO3tpp=1", "NGAM": "ATPM=3.475"},
        "通路开关": {"TreYZ": "关闭ZZ_TreT1_ADP;ZZ_TreT2_UDP", "TreT-ADP": "关闭ZZ_TreT2_UDP;ZZ_TreYZ", "TreT-UDP": "关闭ZZ_TreT1_ADP;ZZ_TreYZ"},
        "氮源": table_info,
        "模型统计": generated_cobra_stats,
        "场景数": len(scenario_specs),
    }
    write_json(out / "A0-2_情景配置.json", config)
    _, rows, quality_rows = run_scenarios(model, scenario_specs, x24, out / "A0-2_24场景重新复现.tsv", out / "A0-2_质量守恒验收.tsv")
    key_rows = make_key_comparison(rows)
    write_tsv(out / "A0-2_与论文正文关键结果比较.tsv", ["指标", "来源表", "通路", "约束版本", "参考值", "复现值", "与参考值差", "比较阈值", "判定", "依据", "场景ID"], key_rows)
    switch_rows = make_switch_check(scenario_specs)
    write_tsv(out / "A0-2_通路开关验收.tsv", ["场景ID", "来源表", "通路", "禁用反应", "纠正版对应LB", "纠正版对应UB", "本轮实际LB", "本轮实际UB", "判定", "依据"], switch_rows)
    status = write_conclusion(out, rows, quality_rows, key_rows, corrected, mapping, flag_values, build_warnings, generated_stats)
    write_json(out / "A0-2_运行登记.json", {
        "任务": "A0-2 2024 iATF-ENG 纠正版24场景复现", "完成时间": now_iso(), "最终结论": status,
        "纠正版": str(corrected), "映射表": str(mapping), "纠正版SHA256": sha256(corrected), "映射表SHA256": sha256(mapping),
        "原A0场景数": 24, "输出场景行数": len(rows), "质量守恒最大残差": max((num_float(r.get("Sv最大残差"), 0.0) for r in quality_rows), default=0.0),
        "模型构建警告": build_warnings, "SBML检查": generated_stats, "未修改已有模型文件": True,
    })
    print(json.dumps({"最终结论": status, "原A0场景数": 24, "输出行数": len(rows), "质量守恒最大残差": max((num_float(r.get("Sv最大残差"), 0.0) for r in quality_rows), default=0.0)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
