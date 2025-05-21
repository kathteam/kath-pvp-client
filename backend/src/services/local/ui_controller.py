from logging import Logger
from webview import Window, windows

from utils import get_logger


class UiController:
    def __init__(self) -> None:
        self.logger: Logger = get_logger(__name__)

    def fullscreen(self) -> None:
        window: Window = windows[0] if windows else None
        if not window:
            self.logger.error("No window found")
            return

        window.toggle_fullscreen()
        self.logger.info("Toggled fullscreen mode")
