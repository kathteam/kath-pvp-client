import os
from logging import Logger

from utils import get_logger
from shared import REF_GENOME, PROGRAM_STORAGE_DIR_SHARED_DATA
from .downloader import (
    download_disease_related_genes_clinvar,
    download_reference_genome_direct,
    reference_genome_exists,
)


class FastaService:
    def __init__(self):
        self.logger: Logger = get_logger(__name__)
        self.data_dir: str = PROGRAM_STORAGE_DIR_SHARED_DATA
        os.makedirs(self.data_dir, exist_ok=True)

    def create_disease_download(self, disease_term: str, max_results: int = 20):
        try:
            downloaded_files = download_disease_related_genes_clinvar(
                disease_term=disease_term,
                max_results=max_results,
                # output_dir=data_dir # Currently we dont need this
            )

            return {
                "status": "success",
                "disease_term": disease_term,
                "max_results": max_results,
                "downloaded_files": downloaded_files,
                "count": len(downloaded_files),
            }
        except Exception as e:
            self.logger.error(f"Error downloading disease genes for {disease_term}: {str(e)}")
            return {"error": f"Download failed: {str(e)}"}, 500

    def download_reference_genome_grch38(self, ref_genome: str = REF_GENOME):
        if reference_genome_exists(ref_genome):
            return {"status": "success", "reference_genome": ref_genome}

        ref_file = download_reference_genome_direct(version=ref_genome)

        if ref_file is None:
            return {"error": "Failed to download reference genome"}, 500

        return {"status": "success", "reference_genome": ref_file}
