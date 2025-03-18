import os
from .logger import get_logger

logger = get_logger(__name__)


# TODO Might need to change this to be more generic, paths are hardcoded
def get_entrypoint() -> str:
    base_dir = os.path.dirname(__file__)

    possible_paths = {
        os.path.join(base_dir, "../../gui/index.html"): "../gui/index.html",
        os.path.join(base_dir, "../../Resources/gui/index.html"): "../Resources/gui/index.html",
        os.path.join(base_dir, "./../gui/index.html"): "./gui/index.html",
    }

    for key, value in possible_paths.items():
        logger.debug(f"Checking path: {key}")
        if os.path.exists(key):
            return value

    logger.error(f"No index.html found. Checked paths: {possible_paths}")
    raise FileNotFoundError(f"No index.html found. Checked paths: {possible_paths}")
