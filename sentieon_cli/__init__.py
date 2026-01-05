from . import argh_parser
from .dnascope import DNAscopePipeline
from .dnascope_hybrid import DNAscopeHybridPipeline
from .dnascope_longread import DNAscopeLRPipeline
from .pangenome import PangenomePipeline
from .sentieon_pangenome import SentieonPangenome


def main():
    """main entry point for this project"""
    parser = argh_parser.CustomArgparseParser()
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
    subparsers = parser.add_subparsers(required=True)

    # DNAscope parser
    pipeline = DNAscopePipeline()
    dnascope_subparser = subparsers.add_parser("dnascope")
    pipeline.add_arguments(dnascope_subparser)
    dnascope_subparser.set_defaults(pipeline=pipeline.main)

    # DNAscope LongRead parser
    pipeline = DNAscopeLRPipeline()
    dnascopelr_subparser = subparsers.add_parser("dnascope-longread")
    pipeline.add_arguments(dnascopelr_subparser)
    dnascopelr_subparser.set_defaults(pipeline=pipeline.main)

    pipeline = DNAscopeHybridPipeline()
    dnascope_hybrid_subparser = subparsers.add_parser("dnascope-hybrid")
    pipeline.add_arguments(dnascope_hybrid_subparser)
    dnascope_hybrid_subparser.set_defaults(pipeline=pipeline.main)

    # Pangenome parser
    pipeline = PangenomePipeline()
    pangenome_subparser = subparsers.add_parser("pangenome")
    pipeline.add_arguments(pangenome_subparser)
    pangenome_subparser.set_defaults(pipeline=pipeline.main)

    # Sentieon pangenome
    pipeline = SentieonPangenome()
    sentieon_pangenome_subparser = subparsers.add_parser("sentieon-pangenome")
    pipeline.add_arguments(sentieon_pangenome_subparser)
    sentieon_pangenome_subparser.set_defaults(pipeline=pipeline.main)

    args = parser.parse_args()
    args.pipeline(args)


if __name__ == "__main__":
    main()
