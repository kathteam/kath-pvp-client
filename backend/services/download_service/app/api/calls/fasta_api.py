from typing import List, Dict, Any, Optional
from pathlib import Path
import os
import logging
import json
import sys

# Add the correct paths to import modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from fasta_downloader import download_disease_related_genes_clinvar, download_reference_genome_direct
from constant_variables import PROGRAM_STORAGE_FOLDER

logger = logging.getLogger("download-service.fasta")

data_dir = os.path.join(PROGRAM_STORAGE_FOLDER, "shared", "data")
os.makedirs(data_dir, exist_ok=True)

class FastaAPI:
    def list_downloads(self):
        """
        List all available FASTA downloads.
        """
        try:
            files = []
            if data_dir.exists():
                files = [f.name for f in data_dir.glob("*.fasta") if f.is_file()]
            
            return {"downloads": files}
        except Exception as e:
            logger.error(f"Error listing downloads: {str(e)}")
            return {"error": f"Failed to list downloads: {str(e)}"}, 500

    def get_download(self, file_name: str):
        """
        Get details about a specific downloaded FASTA file.
        """
        file_path = data_dir / file_name
        
        if not file_path.exists() or not file_path.is_file():
            return {"error": f"File {file_name} not found"}, 404
        
        try:
            stats = file_path.stat()
            return {
                "file_name": file_name,
                "size_bytes": stats.st_size,
                "last_modified": stats.st_mtime,
                "path": str(file_path)
            }
        except Exception as e:
            logger.error(f"Error getting file details for {file_name}: {str(e)}")
            return {"error": f"Failed to get file details: {str(e)}"}, 500

    def create_disease_download(self, disease_term: str, max_results: int = 20):
        """
        Download FASTA files for genes related to a specific disease from ClinVar.
        """
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
                "count": len(downloaded_files)
            }
        except Exception as e:
            logger.error(f"Error downloading disease genes for {disease_term}: {str(e)}")
            return {"error": f"Download failed: {str(e)}"}, 500

    def download_file(self, file_name: str):
        """
        Return the path to a specific FASTA file to download.
        Note: In pywebview, you'll typically return the file path and handle
        the actual download in JavaScript or open the file directly.
        """
        file_path = data_dir / file_name
        
        if not file_path.exists() or not file_path.is_file():
            return {"error": f"File {file_name} not found"}, 404
        
        try:
            return {
                "file_path": str(file_path),
                "file_name": file_name
            }
        except Exception as e:
            logger.error(f"Error downloading file {file_name}: {str(e)}")
            return {"error": f"Failed to download file: {str(e)}"}, 500
    
    def download_reference_genome_grch38(self, ref_gene: str = "GRCh38"):
        """
        Downloads the reference genome for GRCh38/GRCh37.
        """
        ref_file = download_reference_genome_direct(version=ref_gene)
        
        if ref_file is None:
            return {"error": "Failed to download reference genome"}, 500
        
        return {
            "status": "success",
            "reference_genome": ref_file
        }