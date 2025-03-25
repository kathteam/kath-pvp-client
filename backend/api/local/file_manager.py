import os
from logging import Logger
from utils.logger import get_logger

WORKDIR = os.path.join(os.path.expanduser("~"), "Documents")

class FileManager:
    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def list_files(self, path=".") -> list[dict]:
        files = os.listdir(os.path.join(WORKDIR, path))
        file_info = []
        for file in files:
            file_path = os.path.join(WORKDIR, path, file)
            if os.path.isfile(file_path):
                size_kb = os.path.getsize(file_path) / 1024
                file_type = self._determine_file_type(file)
                file_info.append({
                    "filename": file,
                    "type": file_type,
                    "size_kb": size_kb,
                    "item_count": None  # Not applicable for files
                })
            elif os.path.isdir(file_path):
                item_count = len(os.listdir(file_path))
                file_info.append({
                    "filename": file,
                    "type": "folder",
                    "size_kb": None,  # Not applicable for folders
                    "item_count": item_count
                })
        files_string = ", ".join([
            f"{info['filename']} ({info['type']}, {info['size_kb']:.2f} KB)" if info['size_kb'] is not None
            else f"{info['filename']} ({info['type']}, {info['item_count']} items)"
            for info in file_info
        ])
        self.logger.info(f"Listing files in {path}")
        self.logger.info(f"Files: {files_string}")
        return file_info

    def _determine_file_type(self, filename: str) -> str:
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
        if ext in [".txt"]:
            return "text file"
        elif ext in [".csv"]:
            return "CSV"
        elif ext in [".fasta", ".fa"]:
            return "fasta"
        elif ext in [".vcf"]:
            return "VCF"
        elif ext in [".db", ".sqlite"]:
            return "database"
        else:
            return "file"
