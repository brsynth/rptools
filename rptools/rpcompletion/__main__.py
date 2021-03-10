#!/usr/bin/env python

from rptools.rplibs       import rpCache
from rptools.rpcompletion import rp_completion
from rptools.rpcompletion.Args import add_arguments
from rptools import build_args_parser


def _cli():
    parser = build_args_parser(
        prog = 'rpcompletion',
        description = 'Parse RP2 pathways to generate rpSBML collection of unique and complete (cofactors) pathways',
        m_add_args = add_arguments
    )
    args  = parser.parse_args()

    args.pubchem_search = args.pubchem_search.lower() in ['true', 't']

    from rptools.__main__ import init
    logger = init(parser, args)

    cache = rpCache(db='file', logger=logger)

    try:
        result = rp_completion(
            cache,
            args.rp2_pathways,
            args.rp2paths_compounds,
            args.rp2paths_pathways,
            args.outdir,
            int(args.upper_flux_bound),
            int(args.lower_flux_bound),
            int(args.max_subpaths_filter),
            args.pathway_id,
            args.compartment_id,
            args.species_group_id,
            args.sink_species_group_id,
            args.pubchem_search,
            logger=logger
        )
        return result
    except ValueError as e:
        logger.error(str(e))
        return 2



if __name__ == '__main__':
    _cli()
