from .local import UiController
from .remote import HttpClient

from services.download_service.app.api.calls.fasta_api import FastaAPI

class Api:
    def __init__(self) -> None:
        # Local functionality
        self.ui_controller: UiController = UiController()

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()

        # Initialize FastaAPI
        self.fasta_api = FastaAPI()

    #
    # Local operations
    #

    # UI operations
    def fullscreen(self) -> None:
        self.ui_controller.fullscreen()

    #
    # Remote operations
    #
    
    # Fasta operations
    
    def list_downloads(self):
        return self.fasta_api.list_downloads()
    
    def get_download(self, file_name: str):
        return self.fasta_api.get_download(file_name)
    
    def create_disease_download(self, disease_term: str, max_results: int = 20):
        return self.fasta_api.create_disease_download(disease_term, max_results)
    
    def download_file(self, file_name: str):
        return self.fasta_api.download_file(file_name)
    
    def download_reference_genome_grch38(self, ref_gene: str = "GRCh38"):
        return self.fasta_api.download_reference_genome_grch38(ref_gene)