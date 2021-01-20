#!/usr/bin/env python

import logging
from rptools.rplibs       import rpCache
from rptools.rpcompletion import rp_completion, build_args_parser


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    args.pubchem_search = args.pubchem_search.lower() in ['true', 't']

    # Create logger
    logger = logging.getLogger('rptools - rpCompletion')
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, args.log.upper()))

    cache = rpCache(db='file', logger=logger)

    try:
        result = rp_completion(cache,
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
                               logger=logger)
        return result
    except ValueError as e:
        logging.error(str(e))
        return 2



if __name__ == '__main__':
    _cli()
