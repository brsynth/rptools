#!/usr/bin/env python

from rr_cache import rrCache
from rptools.rplibs.Args import add_arguments
from rptools import build_args_parser


def gen_cache(outdir, logger):
    rrCache.generate_cache(outdir, logger)


def _cli():
    parser = build_args_parser(
        prog = 'rpcache',
        description = 'Pre-compute data',
        m_add_args = add_arguments
    )
    args  = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    if args.cache_dir:
        print("rrCache is going to be generated into " + args.cache_dir)
        gen_cache(args.cache_dir, logger)


if __name__ == '__main__':
    _cli()
