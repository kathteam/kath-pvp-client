from dotenv import load_dotenv
import os
from pathlib import Path
from Bio import Entrez
from src.shared.constants import *

from utils.logger import get_logger

logger = get_logger(__name__)


class EnvSetup:
    def __init__(self):

        env_file = Path(PROGRAM_STORAGE_DIR_ENVIRONMENT, ".env")
        load_dotenv(env_file)

        Entrez.email = os.getenv("NCBI_API_EMAIL", "")
        if not Entrez.email:
            logger.warning(
                "Entrez email is not set in the environment variables. Please set it in the .env file."
            )
            return

        Entrez.api_key = os.getenv("NCBI_API_KEY", "")
        if not Entrez.api_key:
            logger.warning(
                "Entrez API key is not set in the environment variables. Please set it in the .env file."
            )
            return
        
        print("Entrez email and API key loaded successfully.")


class FolderSetup:
    def __init__(self):
        self.create_directory(PROGRAM_STORAGE_DIR)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_DATA)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_BLAST)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_SAMPLES)
        self.create_directory(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS)
        self.create_directory(PROGRAM_STORAGE_DIR_ENVIRONMENT)

    def create_directory(self, path: str):
        try:
            os.makedirs(path, exist_ok=True)
        except Exception as e:
            logger.warning(f"Error creating directory {path}: {e}")
