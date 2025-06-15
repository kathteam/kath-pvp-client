import unittest
import sys
import os

# Add the src directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from services.helpers.hgvs_parser import (
    parse_hgvs,
    format_to_hgvs,
    vcf_to_c_hgvs,
    determine_variant_type,
)


class TestHGVS(unittest.TestCase):
    def test_parse_snv(self):
        """Test parsing of single nucleotide variants"""
        hgvs = "chr1:g.123A>G"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["reference"], "A")
        self.assertEqual(result["alternate"], "G")
        self.assertEqual(result["type"], "SNV")

    def test_parse_mnv(self):
        """Test parsing of multi-nucleotide variants"""
        hgvs = "chr1:g.123AT>GC"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["reference"], "AT")
        self.assertEqual(result["alternate"], "GC")
        self.assertEqual(result["type"], "MNV")
        self.assertEqual(result["end_position"], 124)

    def test_parse_deletion(self):
        """Test parsing of deletion variants"""
        hgvs = "chr1:g.123_125del"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["end_position"], 125)
        self.assertEqual(result["alternate"], "-")
        self.assertEqual(result["type"], "DEL")

    def test_parse_insertion(self):
        """Test parsing of insertion variants"""
        hgvs = "chr1:g.123_124insAT"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["end_position"], 124)
        self.assertEqual(result["reference"], "-")
        self.assertEqual(result["alternate"], "AT")
        self.assertEqual(result["type"], "INS")

    def test_parse_delins(self):
        """Test parsing of deletion-insertion variants"""
        hgvs = "chr1:g.123_125delinsAT"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["end_position"], 125)
        self.assertEqual(result["alternate"], "AT")
        self.assertEqual(result["type"], "DELINS")

    def test_parse_duplication(self):
        """Test parsing of duplication variants"""
        hgvs = "chr1:g.123_125dup"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["end_position"], 125)
        self.assertEqual(result["type"], "DUP")

    def test_format_to_hgvs(self):
        """Test formatting variants to HGVS notation"""
        # Test SNV
        self.assertEqual(format_to_hgvs("1", 123, "A", "G"), "chr1:g.123A>G")

        # Test deletion
        self.assertEqual(format_to_hgvs("1", 123, "AT", "-"), "chr1:g.123_124del")

        # Test insertion
        self.assertEqual(format_to_hgvs("1", 123, "-", "AT"), "chr1:g.123_124insAT")

        # Test MNV
        self.assertEqual(format_to_hgvs("1", 123, "AT", "GC"), "chr1:g.123AT>GC")

    def test_vcf_to_c_hgvs(self):
        """Test conversion of VCF format to coding HGVS"""
        # Test SNV
        self.assertEqual(
            vcf_to_c_hgvs("NM_000314.8", 395, "G", "T"), "NM_000314.8:c.395G>T"
        )

        # Test deletion
        self.assertEqual(
            vcf_to_c_hgvs("NM_000314.8", 395, "GT", ""), "NM_000314.8:c.395_396del"
        )

        # Test insertion
        self.assertEqual(
            vcf_to_c_hgvs("NM_000314.8", 395, "", "GT"), "NM_000314.8:c.395_396insGT"
        )

    def test_determine_variant_type(self):
        """Test variant type determination"""
        self.assertEqual(determine_variant_type("A", "G"), "SNV")
        self.assertEqual(determine_variant_type("AT", "GC"), "MNV")
        self.assertEqual(determine_variant_type("AT", "-"), "DEL")
        self.assertEqual(determine_variant_type("-", "AT"), "INS")
        self.assertEqual(determine_variant_type("AT", "GC"), "MNV")
        self.assertEqual(determine_variant_type("", ""), "UNKNOWN")

    def test_invalid_hgvs(self):
        """Test handling of invalid HGVS strings"""
        invalid_cases = [
            "",  # Empty string
            "invalid",  # Invalid format
            "chr1:g.",  # Incomplete
            "chr1:g.123",  # Missing variant information
        ]

        for hgvs in invalid_cases:
            result = parse_hgvs(hgvs)
            self.assertIsNone(result)

    def test_non_standard_format(self):
        """Test parsing of non-standard genomic coordinate formats"""
        hgvs = "1-123-125-AT-GC"
        result = parse_hgvs(hgvs)
        self.assertIsNotNone(result)
        self.assertEqual(result["chromosome"], "1")
        self.assertEqual(result["position"], 123)
        self.assertEqual(result["end_position"], 125)
        self.assertEqual(result["reference"], "AT")
        self.assertEqual(result["alternate"], "GC")


if __name__ == "__main__":
    unittest.main()
