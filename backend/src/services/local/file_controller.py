import os
from logging import Logger

from utils.logger import get_logger

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
