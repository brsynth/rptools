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
        '--top',
        type = int,
        default = '10',
        help = 'Numbers of pathways to be selected'
    )
    parser.add_argument(
        '--outfile',
        type = str,
        default = '',
        help = 'Path to write out file with the ranking'
    )
    parser.add_argument(
        '--outdir',
        type = str,
        default = '',
        help = 'Path to write out selected pathways'
    )

    return parser
