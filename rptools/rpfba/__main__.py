import pandas as pd

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
from rptools.rpfba.medium import (
    __MEDIUM_PATH,
    load_compounds,
)
from rptools.rpfba import (
    runFBA
)
from rptools.rplibs.rpPathway import rpPathway

def entry_point():
    # Args.
    parser = build_args_parser(
        prog = 'rpfba',
        description='Process to Flux Balance Analysis',
        m_add_args=add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    pathway = rpPathway.from_rpSBML(
      infile=args.pathway_file,
      logger=logger
    )

    # Load compounds.
    df_medium_base, df_medium_user = load_compounds(
        medium_id=args.medium_id,
        file_base=__MEDIUM_PATH,
        file_user=args.medium_file,
        logger=logger
    )

    results = runFBA(
      pathway=pathway,
      gem_sbml_path=args.model_file,
      compartment_id=args.compartment_id,
      biomass_rxn_id=args.biomass_rxn_id,
      objective_rxn_id=args.objective_rxn_id,
      sim_type=args.sim,
      fraction_coeff=args.fraction_of,
      merge=args.merge,
      ignore_orphan_species=args.ignore_orphan_species,
      medium_compartment_id=args.medium_compartment_id,
      df_medium_base=df_medium_base,
      df_medium_user=df_medium_user,
      logger=logger
    )

    if results is None:
      logger.info('No results written. Exiting...')
    else:
      logger.info('Writing into file...')
      dirname = os_path.dirname(args.outfile)
      if dirname != '' and not os_path.exists(dirname):
          try:
              os_makedirs(dirname)
          except OSError as exc: # Guard against race condition
              if exc.errno != errno_EEXIST:
                  raise
      pathway.to_rpSBML().write_to_file(args.outfile)
      logger.info('   |--> written in ' + args.outfile)


if __name__ == '__main__':
    entry_point()
