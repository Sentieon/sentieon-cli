import argh

from . import dnascope_longread


def main():
    """main entry point for this project"""
    parser = argh.ArghParser()
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging",
        action="store_const",
        dest="loglevel",
        const="INFO",
        default="WARNING",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Print debugging info",
        action="store_const",
        dest="loglevel",
        const="DEBUG",
    )

    parser.add_commands([dnascope_longread.dnascope_longread])
    parser.dispatch()


if __name__ == "__main__":
    main()
