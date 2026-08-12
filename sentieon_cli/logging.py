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
    """Set the level of the package logger.

    The package's module loggers (``sentieon_cli.dag``,
    ``sentieon_cli.executor``, ...) are all at NOTSET and inherit their
    effective level from the ``sentieon_cli`` package logger, so setting it
    here scopes verbosity to this package without touching the root logger
    (and thus third-party libraries).
    """
    logging.getLogger(__package__).setLevel(level)
