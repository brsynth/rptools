#!/usr/bin/env python

from rptools.rpextractsink  import genSink, build_args_parser
from rptools.rplibs         import rpCache


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    rpcache = rpCache('file', ['cid_strc'], logger=logger)
    genSink(rpcache,
            args.input_sbml,
            args.output_sink,
            args.remove_dead_end,
            args.compartment_id,
            logger=logger)


if __name__ == '__main__':
    _cli()
