import os
import sqlite3
from logging import Logger

from utils import get_logger

WORKDIR = os.path.join(os.path.expanduser("~"), "Documents")


class FileController:
    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def list_files(self, path=".") -> list[dict]:
        if not os.path.exists(os.path.join(WORKDIR, path)):
            raise ValueError(f"Invalid path: {path}")
        files = os.listdir(os.path.join(WORKDIR, path))
        file_info = []
        for file in files:
            file_path = os.path.join(WORKDIR, path, file)
            try:
                if os.path.isfile(file_path):
                    size_kb = os.path.getsize(file_path) / 1024
                    file_type = self._determine_file_type(file)
                    file_info.append(
                        {
                            "filename": file,
                            "type": file_type,
                            "size_kb": size_kb,
                            "item_count": None,  # Not applicable for files
                        }
                    )
                elif os.path.isdir(file_path):
                    item_count = len(os.listdir(file_path))
                    file_info.append(
                        {
                            "filename": file,
                            "type": "folder",
                            "size_kb": None,  # Not applicable for folders
                            "item_count": item_count,
                        }
                    )
            except PermissionError:
                self.logger.warning(f"Permission denied for file or folder: {file_path}")
                file_info.append(
                    {"filename": file, "type": "unknown", "size_kb": None, "item_count": None}
                )
        files_string = ", ".join([info["filename"] for info in file_info])
        self.logger.info(f"Listing files in {path}")
        self.logger.info(f"Files: {files_string}")
        return file_info

    def upload_file(self, path: str, file_name: str, file_content: bytes) -> None:
        target_path = os.path.join(WORKDIR, path, file_name)
        with open(target_path, "wb") as f:
            f.write(bytes(file_content))  # Convert list to bytes
        self.logger.info(f"Uploaded file: {file_name} to {path}")

    def delete_file(self, path: str, file_name: str) -> None:
        target_path = os.path.join(WORKDIR, path, file_name)
        os.remove(target_path)
        self.logger.info(f"Deleted file: {file_name} from {path}")

    def rename_file(self, path: str, old_name: str, new_name: str) -> None:
        old_path = os.path.join(WORKDIR, path, old_name)
        new_path = os.path.join(WORKDIR, path, new_name)
        os.rename(old_path, new_path)
        self.logger.info(f"Renamed file: {old_name} to {new_name} in {path}")

    def _determine_file_type(self, filename: str) -> str:
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
        file_types = {
            "text file": [".txt"],
            "CSV": [".csv", ".tsv", ".xls", ".xlsx", ".ods"],
            "fasta": [".fasta", ".fa"],
            "VCF": [".vcf"],
            "database": [".db", ".sqlite"],
            "PDF": [".pdf", ".doc", ".docx", ".odt", ".rtf"],
            "executable": [".exe", ".sh", ".bat", ".py", ".pl", ".rb", ".jar"],
            "image": [".jpg", ".jpeg", ".png", ".gif", ".bmp", ".svg"],
            "video": [".mp4", ".mkv", ".avi", ".mov", ".wmv", ".flv"],
            "audio": [".mp3", ".wav", ".flac", ".ogg", ".m4a"],
            "archive": [".zip", ".tar", ".gz", ".rar", ".7z", ".iso"],
        }
        for file_type, extensions in file_types.items():
            if ext in extensions:
                return file_type
        return "file"

    def get_kath_directory(self) -> str:
        """
        Determines the ~/.kath directory, creating it if it doesn't exist.

        Returns:
            str: The path to the ~/.kath directory.
        """
        kath_dir = os.path.join(os.path.expanduser("~"), ".kath")
        if not os.path.exists(kath_dir):
            os.makedirs(kath_dir)
        return kath_dir

    def create_vcf_database(self, db_path: str) -> None:
        """
        Creates a SQLite database with a table to store mutation and disease data.
        Only creates the database if it doesn't already exist.

        Args:
            db_path (str): Path to the SQLite database file.
        """
        if os.path.exists(db_path):
            self.logger.info(f"Database already exists at {db_path}")
            return

        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS mutation_data (
                mutation_id INTEGER PRIMARY KEY AUTOINCREMENT,
                file_name TEXT NOT NULL,
                clinical_significance TEXT,
                disease_name TEXT,
                synonyms TEXT,
                chromosome TEXT NOT NULL,
                position INTEGER NOT NULL,
                reference TEXT NOT NULL,
                alternate TEXT NOT NULL,
                hgvs_id TEXT,
                UNIQUE (file_name, clinical_significance, disease_name)
            )
        """
        )

        connection.commit()
        connection.close()
        self.logger.info(f"Created new mutation database at {db_path}")

    def add_mutation_entry(
        self,
        file_name: str,
        chromosome: str,
        position: int,
        reference: str,
        alternate: str,
        clinical_significance: str = None,
        disease_name: str = None,
        synonyms: list = None,
        hgvs_id: str = None,
        db_path: str = None,
    ) -> int:
        """
        Adds a mutation entry to the database.

        Args:
            file_name: Source file name for the mutation
            chromosome: Chromosome identifier (e.g., '1', 'X', 'MT')
            position: Genomic position of the mutation
            reference: Reference nucleotide(s)
            alternate: Alternate nucleotide(s)
            clinical_significance: Clinical significance (e.g., 'Pathogenic')
            disease_name: Associated disease name
            synonyms: List of disease synonyms
            hgvs_id: HGVS identifier for the mutation
            db_path: Optional custom path to database file. If None, uses default in ~/.kath/

        Returns:
            int: ID of the newly inserted mutation entry
        """
        if db_path is None:
            kath_dir = self.get_kath_directory()
            db_path = os.path.join(kath_dir, "mutations.db")

        self.create_vcf_database(db_path)

        if synonyms is not None:
            import json

            synonyms_str = json.dumps(synonyms)
        else:
            synonyms_str = None

        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()

        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO mutation_data
                (file_name, clinical_significance, disease_name, synonyms,
                chromosome, position, reference, alternate, hgvs_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    file_name,
                    clinical_significance,
                    disease_name,
                    synonyms_str,
                    chromosome,
                    position,
                    reference,
                    alternate,
                    hgvs_id,
                ),
            )

            mutation_id = cursor.lastrowid

            connection.commit()
            self.logger.info(f"Added mutation entry (ID: {mutation_id}) to database")
            return mutation_id

        except sqlite3.Error as e:
            self.logger.error(f"Error adding mutation to database: {e}")
            connection.rollback()
            raise
        finally:
            connection.close()

    def get_mutation_entries(
        self,
    ) -> list[dict]:

        kath_dir = self.get_kath_directory()
        db_path = os.path.join(kath_dir, "mutations.db")

        if not os.path.exists(db_path):
            return []

        connection = sqlite3.connect(db_path)
        cursor = connection.cursor()

        query = "SELECT * FROM mutation_data ORDER BY mutation_id DESC"
        params = []

        try:
            cursor.execute(query, params)
            rows = cursor.fetchall()

            # Map rows to dictionaries
            columns = [column[0] for column in cursor.description]
            entries = [dict(zip(columns, row)) for row in rows]

            self.logger.info(f"Retrieved {len(entries)} mutation entries from database")
            return entries

        except sqlite3.Error as e:
            self.logger.error(f"Error retrieving mutation entries: {e}")
            raise
        finally:
            connection.close()
