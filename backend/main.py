import webview
from utils.logger import get_logger
from utils.helpers import get_entrypoint

logger = get_logger(__name__)


# Move elsewhere
class Api:
    def fullscreen(self) -> None:
        window = webview.windows[0] if webview.windows else None
        if not window:
            logger.error("No window found")
            return

        window.toggle_fullscreen()
        logger.info("Toggled fullscreen mode")


if __name__ == "__main__":
    try:
        entrypoint = get_entrypoint()
        logger.info("Starting the application")
        window = webview.create_window(title="Kath", url=entrypoint, js_api=Api())
        webview.start(icon="assets/logo.png")
        logger.info("Application started successfully")
    except Exception as e:
        logger.exception(f"Failed to start the application: {e}")
        raise e
