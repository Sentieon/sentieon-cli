"""
Logging configuration

Handlers live only on the package logger; the loggers returned by
``get_logger`` are plain children that propagate to it. The package logger is
always at DEBUG so verbosity can be enforced per handler: the console handler
follows ``-v``/``-q``/``-d`` while the ``run.log`` file handler always records
DEBUG.
"""

import logging
import os
from typing import Optional

LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"

_package_logger = logging.getLogger(__package__)
_package_logger.setLevel(logging.DEBUG)
_package_logger.propagate = False  # avoid duplicates through the root logger

_console_handler = logging.StreamHandler()
_console_handler.setFormatter(logging.Formatter(LOG_FORMAT))
_console_handler.setLevel(logging.INFO)
_package_logger.addHandler(_console_handler)

# The file handler installed by `add_file_handler`, kept so repeated setup in
# a single process replaces it rather than accumulating handlers.
_file_handler: Optional[logging.FileHandler] = None


def get_logger(name: str) -> logging.Logger:
    """Return a module logger propagating to the package logger."""
    return logging.getLogger(name)


def set_console_level(level: int | str) -> None:
    """Set the verbosity of the console handler."""
    _console_handler.setLevel(level)


def add_file_handler(path: "os.PathLike[str] | str") -> logging.FileHandler:
    """Attach a DEBUG file handler, replacing any previously added one."""
    global _file_handler
    if _file_handler is not None:
        remove_file_handler(_file_handler)
    handler = logging.FileHandler(path, mode="w")
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    handler.setLevel(logging.DEBUG)
    _package_logger.addHandler(handler)
    _file_handler = handler
    return handler


def remove_file_handler(handler: logging.FileHandler) -> None:
    """Detach and close a handler added by ``add_file_handler``."""
    global _file_handler
    _package_logger.removeHandler(handler)
    handler.close()
    if _file_handler is handler:
        _file_handler = None
