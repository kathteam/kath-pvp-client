import logging
import os


def init_logging() -> None:
    # Create logs directory in the backend directory if it doesn't exist
    log_dir = os.path.join(os.path.dirname(__file__), "..", "logs")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(log_dir, "backend.log")),
            logging.StreamHandler(),
        ],
    )


# Initialize logging when this module is imported
init_logging()


def get_logger(name) -> logging.Logger:
    return logging.getLogger(name)
