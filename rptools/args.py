from argparse import ArgumentParser
from collections.abc import Callable

from brs_utils import add_logger_args

from rptools._version import __version__


def build_args_parser(
    prog: str,
    description: str = "",
    epilog: str = "",
    m_add_args: Callable | None = None,
) -> ArgumentParser:

    parser = ArgumentParser(prog=prog, description=description, epilog=epilog)

    # Build Parser with rptools common arguments
    parser = _add_arguments(parser)

    # Add module specific arguments
    if m_add_args is not None:
        parser = m_add_args(parser)

    return parser


def _add_arguments(parser: ArgumentParser) -> ArgumentParser:
    # Add arguments related to the logger
    parser = add_logger_args(parser)

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="show the version number and exit",
    )

    return parser
