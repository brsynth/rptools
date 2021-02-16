from argparse  import ArgumentParser


def add_arguments(parser):
    # parser.add_argument('-sm', '--store_mode', type=str, default='file',
    #                     help='data storage mode: file or db')
    parser.add_argument('--gen_cache', default=None, type=str, dest='cache_dir',
                        help='generate the cache and exits')
    return parser

