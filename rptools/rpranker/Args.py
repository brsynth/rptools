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
        '--rank_outfile',
        type = str,
        default = '',
        help = 'Path to store the ranking result.'
    )
    parser.add_argument(
        '--data_outfile',
        type = str,
        default = '',
        help = 'Path to store selected pathways as a .tar.gz archive.'
    )
    parser.add_argument(
        '--data_outdir',
        type = str,
        default = '',
        help = 'Path to store selected pathways within a folder.'
    )
    # parser.add_argument(
    #     '--outdir',
    #     type = str,
    #     default = '',
    #     help = 'Path to write out selected pathways'
    # )

    return parser
