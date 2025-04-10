import os
from logging import Logger

from utils.logger import get_logger
from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS,
)

from .align import perform_blast_aligning
from .analyze import analyze_blast_xml


class BlastService:
    def __init__(self):
        self.logger: Logger = get_logger(__name__)
        self.data_dir: str = PROGRAM_STORAGE_DIR_SHARED_BLAST
        self.uploads_dir: str = PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS
        os.makedirs(self.data_dir, exist_ok=True)
        os.makedirs(self.uploads_dir, exist_ok=True)

    def align_mutations(self):
        try:
            result_file = perform_blast_aligning()

            return {"status": "success", "result_file": result_file}
        except Exception as e:
            self.logger.error(f"Error performing blast analysis: {str(e)}")
            return {"status": "error", "result_file": "Failed to perform blast analysis"}

    def perform_blast_analysis(self):
        try:
            result_file = analyze_blast_xml()

            return {"status": "success", "result_file": result_file}
        except Exception as e:
            self.logger.error(f"Error analyzing blast XML: {str(e)}")
            return {"status": "error", "result_file": "Failed to analyze blast XML"}
