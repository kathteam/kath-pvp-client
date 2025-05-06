from logging import Logger

from utils.logger import get_logger
from services.utils.prepare_disease import generate_json_diseases


class DiseaseService:
    def __init__(self):
        self.logger: Logger = get_logger(__name__)

    def get_disease_data(self, file_path: str):
        """
        Get disease data from a CSV file and prepare it for display.

        :param file_path: Path to the CSV file containing disease data.
        :return: List of dictionaries with prepared disease data.
        """
        try:
            # Read and prepare disease data
            prepared_data = generate_json_diseases(file_path)
            json_format = []
            for _, row in enumerate(prepared_data):
                json_format.append(
                    {
                        "clinical_significance": row["clinical_significance"],
                        "disease_name": row["disease_name"],
                    }
                )
            return json_format

        except Exception as e:
            self.logger.error(f"Error preparing disease data: {str(e)}")
            return {"error": f"Preparation failed: {str(e)}"}, 500
