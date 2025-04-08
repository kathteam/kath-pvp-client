from .local import UiController, FileManager
from .remote import HttpClient, FastaService, BlastService, DiseaseInformationService


class Api:
    def __init__(self) -> None:

        # -------------------
        # Local functionality
        # -------------------

        self.ui_controller: UiController = UiController()
        self.file_manager: FileManager = FileManager()

        # --------------------
        # Remote functionality
        # --------------------

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()
        # Initialize fasta service
        self.fasta_service = FastaService()
        # Initialize blast service
        self.blast_service = BlastService()

        self.data_service = DiseaseInformationService()
