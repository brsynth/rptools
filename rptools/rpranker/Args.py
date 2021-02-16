from argparse  import ArgumentParser
from rptools._version import __version__


def add_arguments(parser):
    parser.add_argument(
        'pathways',
        nargs = '+',
        type = str,
        default = [],
        help = 'Paths to pathways with global score'
    )
    parser.add_argument(
        '--outfile',
        type = str,
        default = '',
        help = 'Path to write out file with the ranking'
    )

    return parser
