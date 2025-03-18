from logging import Logger
from requests import Session, Response, exceptions

from utils.logger import get_logger


class HttpClient:
    def __init__(self, base_url: str = None) -> None:
        self.base_url: str = base_url or "http://localhost:5000/"
        self.session: Session = Session()
        self.logger: Logger = get_logger(__name__)

    def set_auth_token(self, token: str) -> None:
        """Set authentication token for all requests"""
        self.session.headers.update({"Authorization": f"Bearer {token}"})

    def get(self, endpoint: str, params: dict = None) -> dict:
        """Make a GET request to the remote API"""
        url: str = f"{self.base_url}/{endpoint}"
        self.logger.info(f"Making GET request to {url}")
        try:
            response: Response = self.session.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except exceptions.RequestException as e:
            self.logger.error(f"HTTP request failed: {e}")
            return {"error": str(e)}

    def post(self, endpoint: str, data: dict = None, json: dict = None) -> dict:
        """Make a POST request to the remote API"""
        url: str = f"{self.base_url}/{endpoint}"
        self.logger.info(f"Making POST request to {url}")
        try:
            response: Response = self.session.post(url, data=data, json=json)
            response.raise_for_status()
            return response.json()
        except exceptions.RequestException as e:
            self.logger.error(f"HTTP request failed: {e}")
            return {"error": str(e)}
