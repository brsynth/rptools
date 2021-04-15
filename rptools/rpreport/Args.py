from argparse import ArgumentParser
from rptools._version import __version__

def add_arguments(parser):
    parser.add_argument(
        '-d', '--input_dir',
        action='store_true',
        help='source_path is a folder containing standalone rpSBML file(s).'
    )
    parser.add_argument(
        'source_path',
        type=str,
        help='Path to a tar archive (default) or folder (using \'-d\' option) containing rpSBML file(s).')
    parser.add_argument(
        'output_folder',
        type=str,
        help='Output folder where report file(s) will be generated.')
    parser.add_argument(
        '--dev',
        action='store_true',
        help='For dev purpose only : create supplementary files into a dev folder'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Turn on console verbose mode.'
    )
    parser.add_argument(
        '--standalone',
        action='store_true',
        help='if set will output an autonomous HTML containing all css and js files dependencies.'
    )
    return parser
