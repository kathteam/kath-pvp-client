from logging import Logger

from utils.logger import get_logger

import os

WORKDIR = os.path.join(os.path.expanduser("~"), "Documents")


class FileManager:
    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def list_files(self, path=".") -> list[str]:
        files = os.listdir(os.path.join(WORKDIR, path))
        self.logger.info(f"Listing files in {path}")
        self.logger.info(f"Files: {', '.join(files)}")
        return files
