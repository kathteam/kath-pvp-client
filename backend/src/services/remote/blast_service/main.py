from logging import Logger

from utils.logger import get_logger
from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS,
)

from .align import process_single_fasta
from .find import process_variants


class BlastService:
    def __init__(self):
        self.logger: Logger = get_logger(__name__)
        self.data_dir: str = PROGRAM_STORAGE_DIR_SHARED_BLAST
        self.uploads_dir: str = PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS

    def align_mutations(self, fasta_file: str):
        try:
            aligned_file = process_single_fasta(fasta_file)

            return {"status": "success", "result_file": aligned_file}
        except Exception as e:
            self.logger.error(f"Error performing blast analysis: {str(e)}")
            return {"status": "error", "result_file": "Failed to perform blast analysis"}

    def disease_extraction(self, fasta_file: str):

        try:
            print(f"Processing file: {fasta_file}")
            # Perform blast aligning
            result_file = process_single_fasta(fasta_file)

            if not result_file:
                raise Exception("Failed to perform blast aligning")


            disease_file = process_variants(result_file)

            if not disease_file:
                raise Exception("No diseases found or failed to process variants")

            return {"status": "success", "result_file": disease_file}

        except Exception as e:
            self.logger.error(f"Error performing thorough analysis: {str(e)}")
            return {"status": "error", "result_file": "Failed to perform thorough analysis"}
