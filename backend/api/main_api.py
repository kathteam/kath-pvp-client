import logging
from .local import UiController
from .local import FileManager
from .remote import HttpClient, FastaService

class Api:
    def __init__(self) -> None:
        # Local functionality
        self.ui_controller: UiController = UiController()
        self.file_manager: FileManager = FileManager()


        # Remote functionality
        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()
        # Initialize fasta service
        self.fasta_service = FastaService()

    #
    # Local operations
    #
    def fullscreen(self) -> None:
        self.ui_controller.fullscreen()

    #
    # Remote operations
    #
    def list_files(self, path=".") -> list[dict]:
        try:
            return self.file_manager.list_files(path)
        except Exception as e:
            logging.error(f"Failed to list files in path '{path}': {str(e)}")
            return [{"error": str(e)}]

    def upload_file(self, path: str, file_name: str, file_content: bytes) -> None:
        try:
            self.file_manager.upload_file(path, file_name, file_content)
        except Exception as e:
            logging.error(f"Failed to upload file '{file_name}' to path '{path}': {str(e)}")
            raise ValueError(f"Failed to upload file: {str(e)}") from e
        
    def delete_file(self, path: str, file_name: str) -> None:
        try:
            self.file_manager.delete_file(path, file_name)
        except Exception as e:
            logging.error(f"Failed to delete file '{file_name}' from path '{path}': {str(e)}")
            raise ValueError(f"Failed to delete file: {str(e)}") from e
        
    def rename_file(self, path: str, old_name: str, new_name: str) -> None:
        try:
            self.file_manager.rename_file(path, old_name, new_name)
        except Exception as e:
            logging.error(f"Failed to rename file '{old_name}' to '{new_name}' in path '{path}': {str(e)}")
            raise ValueError(f"Failed to rename file: {str(e)}") from e
