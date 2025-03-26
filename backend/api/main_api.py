from .local import UiController
from .local import FileManager
from .remote import HttpClient


class Api:
    def __init__(self) -> None:
        # Local functionality
        self.ui_controller: UiController = UiController()

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()

        self.file_manager: FileManager = FileManager()

        # Remove functionality

    #
    # Local operations
    #

    # UI operations
    def fullscreen(self) -> None:
        self.ui_controller.fullscreen()

    #
    # Remote operations
    #
    def list_files(self, path=".") -> list[dict]:
        try:
            return self.file_manager.list_files(path)
        except Exception as e:
            # Log the error
            return [{"error": str(e)}]

    def upload_file(self, path: str, file_name: str, file_content: bytes) -> None:
        try:
            self.file_manager.upload_file(path, file_name, file_content)
        except Exception as e:
            # Log the error and re-raise or handle appropriately
            raise Exception(f"Failed to upload file: {str(e)}")
