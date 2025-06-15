import os
import json
import sys
import unittest
from pathlib import Path
import pandas as pd
from unittest.mock import patch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from services.remote.blast_service.find import (
    validate_variant,
    validate_disease_name,
    validate_significance,
    extract_clinvar_diseases,
    load_variants_from_json,
    process_variants,
    shrink_and_count_duplicates,
    file_controller,
)


class TestFind(unittest.TestCase):
    def setUp(self):
        self.test_files_dir = Path(__file__).parent / "test_files"
        self.variations_json = (
            self.test_files_dir
            / "unknown_transcript1_1732746177_leukema_vs_GRCh38_direct_variations.json"
        )
        self.disease_associations_csv = (
            self.test_files_dir / "disease_associations_20250615-144526.csv"
        )

        self.sample_variant = {
            "chromosome": "1",
            "position": 123456,
            "reference_allele": "A",
            "query_allele": "G",
        }

        self.sample_clinvar_response = [
            {
                "_id": "chr1:g.123456A>G",
                "clinvar": {
                    "rcv": [
                        {
                            "clinical_significance": "Pathogenic",
                            "conditions": {
                                "name": "Test Disease",
                                "synonyms": ["Test Disease Synonym"],
                            },
                        }
                    ]
                },
            }
        ]

    def test_validate_variant(self):
        self.assertTrue(validate_variant("1", 123456, "A", "G"))
        self.assertTrue(validate_variant("X", 123456, "AT", "GC"))
        self.assertTrue(validate_variant("MT", 123456, "A", "-"))

        self.assertFalse(validate_variant("25", 123456, "A", "G"))
        self.assertFalse(validate_variant("1", -1, "A", "G"))
        self.assertFalse(validate_variant("1", 123456, "X", "G"))
        self.assertFalse(validate_variant("1", 123456, "A", "X"))

    def test_validate_disease_name(self):
        self.assertTrue(validate_disease_name("Cancer"))
        self.assertTrue(validate_disease_name("Diabetes Type 2"))

        self.assertFalse(validate_disease_name("unknown"))
        self.assertFalse(validate_disease_name("not specified"))
        self.assertFalse(validate_disease_name(""))
        self.assertFalse(validate_disease_name("   "))

    def test_validate_significance(self):
        self.assertTrue(validate_significance("Pathogenic"))
        self.assertTrue(validate_significance("Likely pathogenic"))

        self.assertFalse(validate_significance("unknown"))
        self.assertFalse(validate_significance("not specified"))

    @patch("services.remote.blast_service.find.query_into_myvariant")
    def test_extract_clinvar_diseases(self, mock_query):
        mock_query.return_value = self.sample_clinvar_response

        diseases = extract_clinvar_diseases(self.sample_clinvar_response)

        self.assertEqual(len(diseases), 1)
        disease = diseases[0]
        self.assertEqual(disease["disease_name"], "Test Disease")
        self.assertEqual(disease["clinical_significance"], "Pathogenic")
        self.assertEqual(disease["synonyms"], ["Test Disease Synonym"])

    def test_load_variants_from_json(self):
        test_data = [self.sample_variant]
        test_json_path = self.test_files_dir / "test_variants.json"

        try:
            with open(test_json_path, "w") as f:
                json.dump(test_data, f)

            variants_df = load_variants_from_json(str(test_json_path))

            self.assertFalse(variants_df.empty)
            self.assertEqual(len(variants_df), 1)
            self.assertEqual(variants_df.iloc[0]["chromosome"], "1")
            self.assertEqual(variants_df.iloc[0]["position"], 123456)
            self.assertEqual(variants_df.iloc[0]["reference"], "A")
            self.assertEqual(variants_df.iloc[0]["alternate"], "G")

        finally:
            if test_json_path.exists():
                test_json_path.unlink()

    def test_shrink_and_count_duplicates(self):
        data = {
            "chromosome": ["1", "1", "2"],
            "position": [123, 123, 456],
            "reference": ["A", "A", "C"],
            "alternate": ["G", "G", "T"],
        }
        df = pd.DataFrame(data)

        result = shrink_and_count_duplicates(df)

        self.assertEqual(len(result), 2)  # Should have 2 unique rows
        self.assertTrue("count" in result.columns)
        self.assertEqual(result[result["chromosome"] == "1"]["count"].iloc[0], 2)

    @patch("services.remote.blast_service.find.query_into_myvariant")
    @patch("services.remote.blast_service.find.save_to_csv")
    @patch("services.remote.blast_service.find.file_controller")
    def test_process_variants(self, mock_file_controller, mock_save_csv, mock_query):
        mock_query.return_value = self.sample_clinvar_response

        test_data = [self.sample_variant]
        test_json_path = self.test_files_dir / "test_variants.json"

        try:
            with open(test_json_path, "w") as f:
                json.dump(test_data, f)

            disease_data = [
                {
                    "clinical_significance": "Pathogenic",
                    "disease_name": "Test Disease",
                    "synonyms": ["Test Disease Synonym"],
                    "chromosome": "1",
                    "position": 123456,
                    "reference": "A",
                    "alternate": "G",
                    "hgvs_id": "chr1:g.123456A>G",
                }
            ]

            with patch(
                "services.remote.blast_service.find.extract_clinvar_diseases",
                return_value=disease_data,
            ):
                result = process_variants(str(test_json_path))

                mock_query.assert_called_once()
                mock_save_csv.assert_called_once()
                mock_file_controller.add_mutation_entry.assert_called_once_with(
                    file_name=test_json_path.name,
                    chromosome="1",
                    position=123456,
                    reference="A",
                    alternate="G",
                    clinical_significance="Pathogenic",
                    disease_name="Test Disease",
                    synonyms=["Test Disease Synonym"],
                    hgvs_id="chr1:g.123456A>G",
                )

                self.assertIsNotNone(result)

        finally:
            if test_json_path.exists():
                test_json_path.unlink()


if __name__ == "__main__":
    unittest.main()
