from os import (
  path as os_path,
  makedirs as os_makedirs
)
from errno import (
  EEXIST as errno_EEXIST
)
from logging import (
  Logger,
  getLogger
)
from typing import(
  List,
  Dict,
  Tuple
)
from copy import deepcopy
from rptools import build_args_parser
from rptools.rpfba.Args import add_arguments
from rptools.rpfba import runFBA
from rptools.rplibs import (
  rpSBML,
  rpPathway
)

def entry_point():
    parser = build_args_parser(
        prog = 'rpfba',
        description='Process to Flux Balance Analysis',
        m_add_args=add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    pathway = rpSBML(
      inFile=args.pathway,
      logger=logger
    ).to_Pathway()

    runFBA(
      pathway=pathway,
      gem_sbml_path=args.model,
      sim_type=args.sim,
      src_rxn_id=args.source_reaction,
      tgt_rxn_id=args.target_reaction,
      src_coeff=args.source_coefficient,
      tgt_coeff=args.target_coefficient,
      is_max=args.is_max,
      frac_of_src=args.fraction_of,
      merge=args.merge,
      pathway_id=args.pathway_id,
      objective_id=args.objective_id,
      compartment_id=args.compartment_id,
      ignore_orphan_species=not args.dont_ignore_orphan_species,
      species_group_id=args.species_group_id,
      sink_species_group_id=args.sink_species_group_id,
      upper_flux_bound=float(args.upper_flux_bound),
      lower_flux_bound=float(args.lower_flux_bound),
      logger=logger
    )

    if pathway is None:
      logger.info('No results written. Exiting...')
    else:
      logger.info('Writing into file...')
      if not os_path.exists(os_path.dirname(args.outfile)):
          try:
              os_makedirs(os_path.dirname(args.outfile))
          except OSError as exc: # Guard against race condition
              if exc.errno != errno_EEXIST:
                  raise
      rpSBML.from_Pathway(pathway).write_to_file(args.outfile)
      logger.info('   |--> written in ' + args.outfile)


if __name__ == '__main__':
    entry_point()
