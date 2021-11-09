from glob import glob
from os import path as os_path
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger
from typing import List
from rptools.Args import build_args_parser
from colored import fg, bg, attr


def init(
    parser: ArgumentParser,
    args: Namespace
) -> Logger:
    from brs_utils import create_logger
    from rptools._version import __version__

    if args.log.lower() in ['silent', 'quiet'] or args.silent:
        args.log = 'CRITICAL'

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{color}{typo}rptools {version}{rst}{color} ({prog}){rst}\n'.format(
            prog = logger.name,
            version = __version__,
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    logger.debug(args)

    return logger


def get_modules(path: str) -> List[str]:
    paths = [
        os_path.abspath(
            os_path.join(f + os_path.pardir)
        ) for f
        in glob(
            os_path.join(
                path,
                '*',
                '__init__.py'
            )
        )
    ]
    return [
        os_path.basename(
            os_path.dirname(f)
        ) for f in paths
    ]

def entry_point():
  
    modules = get_modules(os_path.dirname(os_path.abspath(__file__)))

    description = '\nWelcome to rpTools!\n'
    description += '\n\'rptools\' is a package to process rpSBML files but cannot be directly run. Runnable tools are:\n'
    for module in modules:
        description += '   - '+module+'\n'
    description += '\nTo find help for a specific tool, please type:\n'
    description += '   python -m rptools.<tool_name> --help\n\n'

    print(description)

    parser = build_args_parser(
        prog = 'rptools',
        description = 'Package to process rpSBML files'
    )
    args = parser.parse_args()

    logger = init(parser, args)


def _cli():

    entry_point()

    return 0


if __name__ == '__main__':
    _cli()
