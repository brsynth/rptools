#!/usr/bin/env python

from rr_cache import rrCache
from rptools.rpextractsink  import genSink
from rptools.rpextractsink.Args import add_arguments
from rptools import build_args_parser


def _cli():
    parser = build_args_parser(
        prog = 'rpextractsink',
        description = 'Generate the sink from a model SBML by specifying the compartment',
        m_add_args = add_arguments
    )
    args  = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    cache = rrCache('file', ['cid_strc'], logger=logger)
    genSink(cache,
            args.input_sbml,
            args.output_sink,
            args.remove_dead_end,
            args.compartment_id,
            logger=logger)


if __name__ == '__main__':
    _cli()
