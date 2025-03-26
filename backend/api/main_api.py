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
        return  self.file_manager.list_files(path)
