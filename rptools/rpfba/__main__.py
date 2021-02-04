#!/usr/bin/env python

from rptools.rpfba import (
    runFBA,
    build_args_parser
)
from brs_utils import create_logger

def entry_point():
    parser = build_args_parser()
    args   = parser.parse_args()

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '\n{prog} {version}\n'.format(
            prog = logger.name,
            version = __version__
        )
    )
    logger.debug(args)

    result = runFBA(
        args.input_sbml,
        args.gem_sbml,
        args.outfile,
        args.sim_type,
        args.source_reaction,
        args.target_reaction,
        args.source_coefficient,
        args.target_coefficient,
        args.is_max,
        args.fraction_of,
        args.dont_merge,
        args.pathway_id,
        args.objective_id,
        args.compartment_id,
        args.species_group_id,
        args.sink_species_group_id,
        logger=logger
    )

    return result


if __name__ == '__main__':
    entry_point()
