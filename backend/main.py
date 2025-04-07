import os
import platform
from logging import Logger
from webview import Window, create_window, start

from api import Api

from utils.logger import init_logging, get_logger
from utils.helpers import get_entrypoint, initial_route

# Initialize logging
init_logging()
logger: Logger = get_logger(__name__)


# Entry point for the application
if __name__ == "__main__":
    try:
        entrypoint: str = get_entrypoint()
        logger.info("Starting the application")
        window: Window = create_window(
            title="Kath",
            url=entrypoint,
            js_api=Api(),
            min_size=(640, 480),
            maximized=True,
            draggable=True,
        )

        is_wsl = "microsoft-standard" in platform.uname().release.lower() or os.path.exists(
            "/proc/sys/fs/binfmt_misc/WSLInterop"
        )

        # If in WSL, don't pass icon to prevent creating a file
        if is_wsl:
            start(func=initial_route, args=window, debug=True)
        else:
            start(func=initial_route, args=window, icon="assets/logo.png", debug=True)

        logger.info("Application closed successfully")
    except FileNotFoundError as e:
        logger.exception("Failed to find required files to start the application", exc_info=e)
    except ImportError as e:
        logger.exception("Failed to import required modules", exc_info=e)
    except Exception as e:
        logger.exception("Unexpected error when starting the application", exc_info=e)
