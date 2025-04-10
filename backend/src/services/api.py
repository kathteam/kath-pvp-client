from .local import UiController, FileController
from .remote import HttpClient, FastaService, BlastService


class Api:
    def __init__(self) -> None:

        # -------------------
        # Local functionality
        # -------------------

        self.ui_controller: UiController = UiController()
        self.file_controller: FileController = FileController()

        # --------------------
        # Remote functionality
        # --------------------

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()
        # Initialize fasta service
        self.fasta_service = FastaService()
        # Initialize blast service
        self.blast_service = BlastService()
