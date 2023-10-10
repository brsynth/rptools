#!/usr/bin/env python

from rr_cache import rrCache
from .rpextractsink import genSink
from .Args import add_arguments
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

    cache = rrCache(
        attrs=['cid_strc'],
        cache_dir=args.cache_dir,
        logger=logger
    )
    sink = genSink(
        cache,
        args.input_sbml,
        args.remove_dead_end,
        args.compartment_id,
        args.standalone,
        logger=logger
    )
    
    logger.debug(f'Writing sink to {args.output_sink}...')
    # Write the sink file
    with open(args.output_sink, 'w', encoding='utf-8') as outS:
        write(outS, ['Name', 'InChI'])
        for _mnx_id, _inchi in sink.items():
            write(outS, [_mnx_id, _inchi])


def write(outFile, elts, delimiter=',', quotechar='"'):
    """
    Write elements of elts list into file as 'csv' would do

    :param file: file to write into
    :param elts: list of elements to write
    :param delimiter: character to insert between each element
    :param quotechar: character to put around each element
    """
    if elts:
        outFile.write(quotechar+elts[0]+quotechar)
        for elt in elts[1:]:
            outFile.write(delimiter+quotechar+elt+quotechar)
    outFile.write('\n')


if __name__ == '__main__':
    _cli()
