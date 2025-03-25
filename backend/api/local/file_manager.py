import os
from logging import Logger
from utils.logger import get_logger

WORKDIR = os.path.join(os.path.expanduser("~"), "Documents")

class FileManager:
    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def list_files(self, path=".") -> list[tuple[str, float]]:
        files = os.listdir(os.path.join(WORKDIR, path))
        file_info = []
        for file in files:
            file_path = os.path.join(WORKDIR, path, file)
            if os.path.isfile(file_path):
                size_kb = os.path.getsize(file_path) / 1024
                file_info.append((file, size_kb))
            else:
                file_info.append((file, None))  # None for directories or non-files
        files_string = ", ".join([f"{name} ({size:.2f} KB)" if size is not None
                                  else f"{name} (dir)" for name, size in file_info])
        self.logger.info(f"Listing files in {path}")
        self.logger.info(f"Files: {files_string}")
        return file_info
