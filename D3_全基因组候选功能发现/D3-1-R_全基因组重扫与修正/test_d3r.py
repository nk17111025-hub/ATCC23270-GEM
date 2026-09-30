import csv
import importlib.util
import json
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("d3r", HERE / "build_d3r.py")
d3r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d3r)


def rows(name):
    with (HERE / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class D3RRegression(unittest.TestCase):
    def test_exact_id_never_accepts_substring(self):
        reaction = {"对象ID": "PHOSNACMURPENTATRANS-RXN", "外部ID": "MetaCyc:ABC;RPE2|RPE-X"}
        self.assertFalse(d3r.exact_biocyc_match("RPE", reaction))
        self.assertTrue(d3r.exact_biocyc_match("RPE2", reaction))
        self.assertTrue(d3r.exact_biocyc_match("PHOSNACMURPENTATRANS-RXN", reaction))
        self.assertFalse(d3r.exact_biocyc_match("TRANS-RXN", reaction))
        for reaction in d3r.load_biocyc()[3].values():
            self.assertFalse(d3r.exact_biocyc_match("RPE", reaction))

    def test_current_genome_is_starting_population(self):
        genome, _ = d3r.parse_gff()
        all_rows = rows("D3-1-R_全基因跨数据库交叉引用.tsv")
        self.assertEqual(len(genome), 3021)
        self.assertEqual({r["当前locus"] for r in all_rows}, set(genome))
        self.assertEqual(len(rows("D3-1-R_模型反应精确匹配.tsv")), 631)

    def test_afe_2841_chain_and_evidence_limits(self):
        row = next(r for r in rows("D3-1-R_全基因跨数据库交叉引用.tsv") if r["当前locus"] == "RU820_RS13140")
        self.assertEqual(row["当前protein ID"], "WP_012537424.1")
        self.assertEqual(row["旧AFE locus"], "AFE_2841")
        self.assertEqual(row["KEGG KO"], "K25026")
        self.assertEqual(row["KEGG EC"], "2.7.1.2")
        self.assertEqual(row["2024 GPR反应"], "BDGK")
        self.assertEqual(row["2016 GPR反应"], "")
        self.assertEqual(row["BioCyc reaction ID"], "")
        self.assertEqual(row["BRENDA EC记录"], "")
        self.assertIn("RHEA:17825", row["Rhea精确交叉引用"])

    def test_focus_tables_do_not_include_keyword_false_positives(self):
        organic = rows("D3-1-R_有机碳转运证据.tsv")
        self.assertTrue(organic)
        for row in organic:
            name = row["NCBI当前功能"].lower()
            self.assertRegex(name, r"transporter|permease|porter|porin")
            self.assertNotIn("transferase", name)
            self.assertNotIn("synthase", name)
        central = rows("D3-1-R_中央碳代谢证据.tsv")
        self.assertTrue(central)
        self.assertFalse(any("c-terminal-like domain-containing" in r["NCBI当前功能"].lower() for r in central))
        fes = rows("D3-1-R_FeS_ETC证据.tsv")
        self.assertTrue(fes)
        self.assertFalse(any("alcohol oxidase" in r["NCBI当前功能"].lower() for r in fes))
        all_rows = {r["当前locus"]: r for r in rows("D3-1-R_全基因跨数据库交叉引用.tsv")}
        self.assertEqual(all_rows["RU820_RS11605"]["有机碳转运表入选"], "否")
        self.assertIn("无机底物", all_rows["RU820_RS11605"]["重点分类冲突"])
        for locus in ("RU820_RS07440", "RU820_RS12715"):
            self.assertEqual(all_rows[locus]["转运候选"], "否")
            self.assertEqual(all_rows[locus]["旧模型未覆盖功能候选"], "否")
        self.assertEqual(all_rows["RU820_RS00635"]["FeS_ETC表入选"], "否")
        for locus in ("RU820_RS12075", "RU820_RS12090", "RU820_RS12120", "RU820_RS12125", "RU820_RS12135", "RU820_RS09965", "RU820_RS09970"):
            self.assertEqual(all_rows[locus]["FeS_ETC表入选"], "是", locus)

    def test_multimapping_is_not_confirmed_gpr(self):
        all_rows = {r["当前locus"]: r for r in rows("D3-1-R_全基因跨数据库交叉引用.tsv")}
        for locus, rid in (("RU820_RS03065", "PHFT"), ("RU820_RS03835", "PGDH")):
            row = all_rows[locus]
            self.assertEqual(row["映射冲突"], "是")
            self.assertEqual(row["旧模型GPR覆盖"], "歧义候选")
            self.assertEqual(row["旧模型未覆盖功能候选"], "待定")
            self.assertIn(rid, row["2024 GPR歧义候选反应"])

    def test_brenda_details_resolve_to_raw_enzyme_reaction(self):
        path = d3r.ROOT / "05_BRENDA/D3-1/原始响应/BRENDA_FeS_reaction_compounds.json"
        raw = json.loads(path.read_text(encoding="utf-8"))["results"]["bindings"]
        valid = {(r["enzyme"]["value"].rsplit("/", 1)[-1], r["rxn"]["value"], r["orgName"]["value"]) for r in raw}
        evidence = rows("D3-1-R_注释证据表.tsv")
        for row in evidence:
            if row["数据库"] != "BRENDA":
                continue
            if row["证据字段"] == "enzyme/EC":
                self.assertIn("不据EC聚合赋值", row["记录内容"])
            if row["证据字段"] == "reaction/compound":
                record = row["记录ID"].removeprefix("enzyme:")
                enzyme, reaction = record.split(";reaction:", 1)
                self.assertIn((enzyme, reaction, row["菌株"]), valid)

    def test_checksums_and_missing_data(self):
        manifest = rows("D3-1-R_原始数据校验清单.tsv")
        self.assertTrue(manifest)
        self.assertFalse(any(r["校验"] == "FAIL" for r in manifest))
        self.assertTrue(any(r["校验"] == "PASS" for r in manifest))
        missing = rows("D3-1-R_访问缺失表.tsv")
        self.assertTrue(any(r["来源"] == "MetaCyc" and r["状态"] == "未独立取得" for r in missing))
        self.assertTrue(any(r["来源"] == "UniProt" and r["状态"] == "未检查" for r in missing))
        for old_gene in ("AFE_2291", "AFE_3147"):
            self.assertTrue(any(r["当前locus"] == "未映射:" + old_gene and r["状态"] == "模型GPR旧locus未在B1映射表" for r in missing))
        model = {r["模型反应ID"]: r for r in rows("D3-1-R_模型反应精确匹配.tsv")}
        self.assertIn("AFE_3147", model["CYTAA31"]["GPR旧locus未在B1映射表"])
        self.assertIn("AFE_2291", model["PSSA120"]["GPR旧locus未在B1映射表"])

    def test_rhea_full_ec_query_scope(self):
        manifest = rows("Rhea_新增查询清单.tsv")
        self.assertEqual(len(manifest), 46)
        self.assertTrue(all((d3r.ROOT / r["文件"]).exists() for r in manifest))
        self.assertEqual(sum(len(d3r.tokens(r["查询EC"])) for r in manifest), 552)
        metrics = json.loads((HERE / "D3-1-R_metrics.json").read_text(encoding="utf-8"))
        self.assertEqual(metrics["rhea_full_ec"], metrics["rhea_ec_result"] + metrics["rhea_ec_no_return"])


if __name__ == "__main__":
    unittest.main()
