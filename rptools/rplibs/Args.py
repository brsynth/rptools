from os import path as os_path
from argparse  import ArgumentParser
from rptools._version import __version__

def add_arguments(parser):
    parser.add_argument(
        '--pathways',
        required=True,
        type=str,
        nargs='+',
        help='Pathways (rpSBML) to make statistics on'
    )
    return parser
