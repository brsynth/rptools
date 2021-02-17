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
        help = 'Path to store selected pathways as a .tar.gz archive, unless --light is passed.'
    )
    parser.add_argument(
        '--light',
        action = 'store_true',
        default = False,
        help = 'If set, selected pathway filenames are stored into a file (default: False). Needs --outfile argument.'
    )
    # parser.add_argument(
    #     '--outdir',
    #     type = str,
    #     default = '',
    #     help = 'Path to write out selected pathways'
    # )

    return parser
