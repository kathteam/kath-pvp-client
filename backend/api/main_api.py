from .local import UiController
from .remote import HttpClient

from services.download_service.app.api.calls.fasta_api import FastaAPI

FastaAPI = FastaAPI()

class Api:
    def __init__(self) -> None:
        # Local functionality
        self.ui_controller: UiController = UiController()

        # Initialize HTTP client that will be shared across all services
        self.http_client: HttpClient = HttpClient()

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
    
    # Fasta operations
    
    def list_downloads(self):
        return self.FastaAPI.list_downloads()
    
    def get_download(self, file_name: str):
        return self.FastaAPI.get_download(file_name)
    
    def create_disease_download(self, disease_term: str, max_results: int = 20):
        return self.FastaAPI.create_disease_download(disease_term, max_results)
    
    def download_file(self, file_name: str):
        return self.FastaAPI.download_file(file_name)
    
    def download_reference_genome_grch38(self, ref_gene: str = "GRCh38"):
        return self.FastaAPI.download_reference_genome_grch38(ref_gene)