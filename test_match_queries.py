import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import Stressbiomarkers_insects as extractor


class MatchQueryTests(unittest.TestCase):
    def setUp(self):
        extractor.POSITIVE_MATCH_CACHE.clear()
        extractor.NEGATIVE_MATCH_CACHE.clear()
        self.cfg = {
            "token_min_len": 4,
            "max_n": 1,
            "fuzzy_min_len": 15,
            "fuzzy_max_len": 30,
            "fuzzy_k": 0,
            "fuzzy_m": 30,
            "fuzzy_threshold": 80,
            "fuzzy_len_diff_ratio": 0.2,
        }

    def extract(self, text, db):
        extractor.METABOLITE_KEYS = list(db)
        with patch.object(extractor, "extract_text", return_value=text), \
             patch.object(extractor, "extract_insect_species", return_value=["Apis mellifera"]), \
             patch.object(extractor, "get_logp", return_value=None):
            return extractor.process_pdf("paper.pdf", db, set(), {}, self.cfg)

    def test_exact_queries_include_all_paper_spellings_and_survive_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            hmdb_path = Path(directory) / "hmdb.xml"
            hmdb_path.write_text(
                '<hmdb xmlns="http://www.hmdb.ca"><metabolite>'
                '<accession>HMDB1</accession><name>Canonical</name>'
                '<chemical_formula>C3H7NO2</chemical_formula>'
                '<synonyms><synonym>betaalanine</synonym>'
                '<synonym>Alternate</synonym></synonyms>'
                '</metabolite></hmdb>', encoding="utf-8"
            )
            lipid_path = Path(directory) / "lipids.sdf"
            lipid_path.write_text(
                '> <LM_ID>\nLMP1\n\n> <NAME>\nLipidname\n\n'
                '> <FORMULA>\nC4H8O2\n\n> <SYSTEMATIC_NAME>\nSystematic\n\n'
                '> <ABBREVIATION>\nLipabbr\n\n> <SYNONYMS>\nLipsyn; Lipalias\n$$$$\n',
                encoding="utf-8"
            )
            db = extractor.build_combined_db(
                extractor.load_hmdb(hmdb_path), extractor.load_lipidmaps(lipid_path)
            )

        text = "βalanine Betaalanine Alternate Canonical Lipsyn Systematic Lipabbr"
        for _ in range(2):
            results = {row["HMDB ID/LMP"]: row for row in self.extract(text, db)}
            self.assertEqual(set(results), {"HMDB1", "LMP1"})
            self.assertEqual(
                set(results["HMDB1"]["Fuzzy query"].splitlines()),
                {"βalanine", "Betaalanine", "Alternate", "Canonical"},
            )
            self.assertEqual(
                set(results["LMP1"]["Fuzzy query"].splitlines()),
                {"Lipsyn", "Systematic", "Lipabbr"},
            )
            for row in results.values():
                self.assertEqual(row["Match type"], "EXACT")
                self.assertFalse(row["Fuzzy match"])

        result = self.extract("Alternate", db)[0]
        self.assertEqual(result["Fuzzy query"], "Alternate")

    def test_cached_fuzzy_match_remains_fuzzy(self):
        db = {
            "abcdefghijklmnop": {
                "NAME": "abcdefghijklmnop",
                "HMDB/LMP_ID": "HMDB2",
                "CHEMICAL_FORMULA": "C3H7NO2",
                "DATABASE_SOURCE": "HMDB",
            }
        }
        self.cfg["fuzzy_k"] = 10
        for _ in range(2):
            result = self.extract("abcdefghijklmnoz", db)[0]
            self.assertEqual(result["Match type"], "FUZZY")
            self.assertEqual(result["Fuzzy query"], "abcdefghijklmnoz")
            self.assertEqual(result["Fuzzy match"], "abcdefghijklmnop")

    def test_exact_queries_are_written_to_workbook(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "paper.pdf").touch()
            hmdb_path = root / "hmdb.xml"
            hmdb_path.write_text(
                '<hmdb xmlns="http://www.hmdb.ca"><metabolite>'
                '<accession>HMDB1</accession><name>Canonical</name>'
                '<chemical_formula>C3H7NO2</chemical_formula>'
                '<synonyms><synonym>Alternate</synonym></synonyms>'
                '</metabolite></hmdb>', encoding="utf-8"
            )
            taxon_path = root / "taxon.txt"
            taxon_path.write_text(
                "class\tscientificName\nInsecta\tApis mellifera\n", encoding="utf-8"
            )
            output_path = root / "results.xlsx"
            args = Namespace(
                pdf_dir=str(root), hmdb_path=str(hmdb_path), lipidmaps_path="",
                taxon_path=str(taxon_path), log_path=str(root / "test.log"),
                output_xlsx=str(output_path), **self.cfg
            )
            with patch.object(extractor, "parse_args", return_value=args), \
                 patch.object(extractor, "extract_text", return_value="Apis mellifera Canonical Alternate"), \
                 patch.object(extractor, "get_logp", return_value=None):
                extractor.main()
            rows = extractor.pd.read_excel(output_path, sheet_name="metabolite_hits")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows.iloc[0]["Match type"], "EXACT")
            self.assertEqual(
                set(rows.iloc[0]["Fuzzy query"].splitlines()), {"Canonical", "Alternate"}
            )


if __name__ == "__main__":
    unittest.main()
