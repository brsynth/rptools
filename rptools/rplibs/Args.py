from argparse  import ArgumentParser


def build_args_parser():
    parser = ArgumentParser(prog='rpcache', description='Python script to pre-compute data')
    parser = _add_arguments(parser)
    return parser


def _add_arguments(parser):
    # parser.add_argument('-sm', '--store_mode', type=str, default='file',
    #                     help='data storage mode: file or db')
    parser.add_argument('--gen_cache', default=None, type=str, dest='cache_dir',
                        help='generate the cache and exits')
    parser.add_argument('--log', metavar='ARG',
                        type=str, choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='error',
                        help='Adds a console logger for the specified level (default: error)')
    return parser

