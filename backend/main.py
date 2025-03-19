from logging import Logger
from webview import Window, create_window, start

from api import Api

from utils.logger import init_logging, get_logger
from utils.helpers import get_entrypoint

# Initialize logging
init_logging()
logger: Logger = get_logger(__name__)


# Entry point for the application
if __name__ == "__main__":
    try:
        entrypoint: str = get_entrypoint()
        logger.info("Starting the application")
        window: Window = create_window(title="Kath", url=entrypoint, js_api=Api())
        start(icon="assets/logo.png", debug=True)
        logger.info("Application closed successfully")
    except Exception as e:
        logger.exception(f"Failed to start the application", exc_info=e)
