from os import path as os_path
from argparse import (
    ArgumentParser,
    Namespace
)
from logging import Logger


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
        'rptools {version} ({prog})\n'.format(
            prog = logger.name,
            version = __version__
        )
    )
    logger.debug(args)

    return logger


def _cli():

    with open(os_path.join(os_path.dirname(os_path.abspath(__file__)), '.env'), 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('MODULES='):
                modules = line.splitlines()[0].split('=')[1].lower().split(',')

    print()
    print('Welcome to rpTools!')
    print()
    print('\'rptools\' is a package and cannot be directly executed. Executable tools are:')
    for module in modules:
        print('   - '+module)
    print()
    print('To find help for a specific tool, please type:')
    print('   python -m rptools.<tool_name> --help')
    print()

    return 0


if __name__ == '__main__':
    _cli()
