#!/usr/bin/env python

from rptools.rplibs      import rpCache
from rptools.rplibs.Args import build_args_parser


def gen_cache(outdir, logger):
    rpCache.generate_cache(outdir, logger)


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    if args.cache_dir:
        print("rpCache is going to be generated into " + args.cache_dir)
        gen_cache(args.cache_dir, logger)


if __name__ == '__main__':
    _cli()
