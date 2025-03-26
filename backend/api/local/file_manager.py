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
            try:
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
            except PermissionError:
                self.logger.warning(f"Permission denied for file or folder: {file_path}")
                file_info.append({
                    "filename": file,
                    "type": "unknown",
                    "size_kb": None,
                    "item_count": None
                })
        files_string = ", ".join([info["filename"] for info in file_info])
        self.logger.info(f"Listing files in {path}")
        self.logger.info(f"Files: {files_string}")
        return file_info

    def _determine_file_type(self, filename: str) -> str:
        _, ext = os.path.splitext(filename)
        ext = ext.lower()
        if ext in [".txt"]:
            return "text file"
        if ext in [".csv", ".tsv", ".xls", ".xlsx", ".ods"]:
            return "CSV"
        if ext in [".fasta", ".fa"]:
            return "fasta"
        if ext in [".vcf"]:
            return "VCF"
        if ext in [".db", ".sqlite"]:
            return "database"
        if ext in [".pdf", ".doc", ".docx", ".odt", ".rtf"]:
            return "PDF"
        if ext in [".exe", ".sh", ".bat", ".py", ".pl", ".rb", ".jar"]:
            return "executable"
        if ext in [".jpg", ".jpeg", ".png", ".gif", ".bmp", ".svg"]:
            return "image"
        if ext in [".mp4", ".mkv", ".avi", ".mov", ".wmv", ".flv"]:
            return "video"
        if ext in [".mp3", ".wav", ".flac", ".ogg", ".m4a"]:
            return "audio"
        if ext in [".zip", ".tar", ".gz", ".rar", ".7z", ".iso"]:
            return "archive"
        return "file"
