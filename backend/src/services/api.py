from .local import UiController, FileController
from .remote import HttpClient, FastaService, BlastService, DiseaseService


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

        self.http_client: HttpClient = HttpClient()
        self.fasta_service = FastaService()
        self.blast_service = BlastService()
        self.disease_service = DiseaseService()
