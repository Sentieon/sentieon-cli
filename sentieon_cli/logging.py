import logging

handler = logging.StreamHandler()
handler.setFormatter(
    logging.Formatter("%(levelname)s:%(name)s:%(message)s")
)


def get_logger(name: str):
    """Return a logger with a StreamHandler."""
    logger = logging.getLogger(name)
    logger.addHandler(handler)
    logger.propagate = False
    return logger
