from .local import UiController
from .remote import HttpClient, FastaService


class Api:
    def __init__(self) -> None:
        #
        # Local operations
        #

        # Initialize ui controller
        self.ui_controller: UiController = UiController()

        #
        # Remote operations
        #

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()

        # Initialize fasta service
        self.fasta_service = FastaService()
