import argparse

from . import argh_parser
from .dnascope import DNAscopePipeline
from .dnascope_hybrid import DNAscopeHybridPipeline
from .dnascope_longread import DNAscopeLRPipeline
from .hybrid_pangenome import HybridPangenome
from .job import Job
from .sentieon_pangenome import SentieonPangenome
from .util import __version__


def add_logging_args(parser: argparse.ArgumentParser) -> None:
    """Add the console verbosity flags"""
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging (the default)",
        action="store_true",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        help="Only log warnings and errors",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Print debugging info",
        action="store_true",
    )


def resolve_loglevel(args: argparse.Namespace) -> str:
    """The console log level implied by the verbosity flags"""
    if args.debug:
        return "DEBUG"
    if args.quiet and not args.verbose:
        return "WARNING"
    return "INFO"


def main():
    """main entry point for this project"""
    parser = argh_parser.CustomArgparseParser()
    add_logging_args(parser)
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
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

    # DNAscope pangenome
    pipeline = SentieonPangenome()
    dnascope_pangenome_subparser = subparsers.add_parser(
        "dnascope-pangenome", aliases=["sentieon-pangenome"]
    )
    pipeline.add_arguments(dnascope_pangenome_subparser)
    dnascope_pangenome_subparser.set_defaults(pipeline=pipeline.main)

    # Hybrid pangenome
    pipeline = HybridPangenome()
    hybrid_pangenome_subparser = subparsers.add_parser("hybrid-pangenome")
    pipeline.add_arguments(hybrid_pangenome_subparser)
    hybrid_pangenome_subparser.set_defaults(pipeline=pipeline.main)

    args = parser.parse_args()
    args.loglevel = resolve_loglevel(args)
    # Job ids must be unique for the whole run, which may execute more than
    # one DAG, so numbering restarts here rather than per DAG.
    Job.reset_ids()
    args.pipeline(args)


if __name__ == "__main__":
    main()
