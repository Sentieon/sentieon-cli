import colorlog

handler = colorlog.StreamHandler()
handler.setFormatter(
    colorlog.ColoredFormatter(
        "%(log_color)s%(levelname)s:%(name)s:%(message)s"
    )
)


def get_logger(name):
    """Return a logger with a colorlog handler."""
    logger = colorlog.getLogger(name)
    logger.addHandler(handler)
    return logger
