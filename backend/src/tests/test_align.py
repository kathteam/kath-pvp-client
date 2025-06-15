import unittest
import sys
import os
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA,
)
from services.remote.blast_service.align import process_single_fasta


class TestAlign(unittest.TestCase):
    def setUp(self):
        self.test_files_dir = Path(__file__).parent / "test_files"
        self.query_fasta = (
            self.test_files_dir / "unknown_transcript1_1732746177_leukema.fasta"
        )
        self.reference_fasta_dir = PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA / "reference"
        self.blast_output_dir = self.test_files_dir / "blast_output"

    def test_process_single_fasta_with_valid_input(self):
        self.blast_output_dir.mkdir(parents=True, exist_ok=True)

        print(self.query_fasta)
        print(self.reference_fasta_dir)
        print(self.blast_output_dir)

        result_path = process_single_fasta(
            query_fasta_path=str(self.query_fasta),
            reference_fasta_dir=self.reference_fasta_dir,
            blast_output_dir=self.blast_output_dir,
        )

        self.assertIsNotNone(result_path, "Result path should not be None")
        self.assertTrue(Path(result_path).exists(), "Result file should exist")
        self.assertTrue(result_path.endswith(".json"), "Result should be a JSON file")

        if Path(result_path).exists():
            Path(result_path).unlink()

    def test_process_single_fasta_with_invalid_query(self):
        invalid_path = self.test_files_dir / "nonexistent.fasta"

        result = process_single_fasta(
            query_fasta_path=str(invalid_path),
            reference_fasta_dir=self.reference_fasta_dir,
            blast_output_dir=self.blast_output_dir,
        )

        self.assertIsNone(result, "Result should be None for invalid query file")

    def test_process_single_fasta_with_invalid_reference_dir(self):
        invalid_ref_dir = Path("/nonexistent/directory")

        result = process_single_fasta(
            query_fasta_path=str(self.query_fasta),
            reference_fasta_dir=invalid_ref_dir,
            blast_output_dir=self.blast_output_dir,
        )

        self.assertIsNone(
            result, "Result should be None for invalid reference directory"
        )

    def tearDown(self):
        if self.blast_output_dir.exists():
            for file in self.blast_output_dir.glob("*"):
                if file.is_file():
                    file.unlink()
            self.blast_output_dir.rmdir()


if __name__ == "__main__":
    unittest.main()
