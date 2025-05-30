import os
from logging import Logger
from webview import Window

from .logger import get_logger

logger: Logger = get_logger(__name__)


# TODO Might need to change this to be more generic, paths are hardcoded
def get_entrypoint() -> str:
    base_dir: str = os.path.dirname(__file__)

    # TODO Fix at some point (first three work for build exe, other three for development)
    possible_paths: dict[str, str] = {
        os.path.join(base_dir, "../../gui/index.html"): "../gui/index.html",
        os.path.join(base_dir, "../../Resources/gui/index.html"): "../Resources/gui/index.html",
        os.path.join(base_dir, "./../gui/index.html"): "./gui/index.html",
        os.path.join(base_dir, "../../../gui/index.html"): "../../gui/index.html",
        os.path.join(
            base_dir, "../../../Resources/gui/index.html"
        ): "../../Resources/gui/index.html",
        os.path.join(base_dir, "./../../gui/index.html"): "./../gui/index.html",
    }

    for key, value in possible_paths.items():
        logger.debug(f"Checking path: {key}")
        if os.path.exists(key):
            return value

    logger.error(f"No index.html found. Checked paths: {possible_paths}")
    raise FileNotFoundError("No index.html found.")


def initial_func(window: Window) -> None:
    _initial_route(window)
    _initial_setup(window)


def _initial_setup(window: Window) -> None:
    window.evaluate_js(
        """
        window.dispatchEvent(new CustomEvent('setup-status', {
            detail: { status: 'started' }
        }))
        """
    )
    from services.remote.fasta_service import FastaService

    fastaService = FastaService()
    tasks: list[tuple[function, str]] = [
        (fastaService.download_reference_genome_grch38, "Downloading reference genome"),
    ]

    for i, (task, message) in enumerate(tasks):
        progress = i / len(tasks)
        window.evaluate_js(
            f"""
            window.dispatchEvent(new CustomEvent('setup-status', {{
                detail: {{
                    status: 'progress', 
                    progress: {progress}, 
                    message: "{message}"
                }}
            }}))
            """
        )
        task()

    window.evaluate_js(
        """
        window.dispatchEvent(new CustomEvent('setup-status', {
            detail: { status: 'completed' }
        }))
        """
    )


def _initial_route(window: Window) -> None:
    window.evaluate_js("window.location.href = '#/dashboard'")
