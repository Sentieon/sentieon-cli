import logging

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))


def get_logger(name: str) -> logging.Logger:
    """Return a logger with a StreamHandler."""
    logger = logging.getLogger(name)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def set_level(level: int | str) -> None:
    """Set the level of the package logger."""
    logging.getLogger(__package__).setLevel(level)
