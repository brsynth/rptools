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
        description = 'Process to Flux Balance Analysis',
        m_add_args = add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    pathway = rpSBML(
      inFile=args.pathway,
      logger=logger
    ).to_Pathway()

    results = runFBA(
                      pathway = pathway,
                gem_sbml_path = args.model,
                     sim_type = args.sim,
                   src_rxn_id = args.source_reaction,
                   tgt_rxn_id = args.target_reaction,
                    src_coeff = args.source_coefficient,
                    tgt_coeff = args.target_coefficient,
                       is_max = args.is_max,
                  frac_of_src = args.fraction_of,
                        merge = args.merge,
                   pathway_id = args.pathway_id,
                 objective_id = args.objective_id,
               compartment_id = args.compartment_id,
        ignore_orphan_species = not args.dont_ignore_orphan_species,
             species_group_id = args.species_group_id,
        sink_species_group_id = args.sink_species_group_id,
                       logger = logger
    )

    # results = {'species': {'MNXM15': {'biomass_shadow_price': {'value': -0.0}, 'fraction_shadow_price': {'value': -0.0}}, 'MNXM4': {'biomass_shadow_price': {'value': -0.0}, 'fraction_shadow_price': {'value': -0.0}}, 'MNXM5': {'biomass_shadow_price': {'value': -0.3745738207611868}, 'fraction_shadow_price': {'value': -2.6225806451612894}}, 'TARGET_0000000001': {'biomass_shadow_price': {'value': 0.0}, 'fraction_shadow_price': {'value': -1.0}}, 'MNXM6': {'biomass_shadow_price': {'value': -0.3838016323334834}, 'fraction_shadow_price': {'value': -2.6870967741935474}}, 'MNXM13': {'biomass_shadow_price': {'value': -0.0}, 'fraction_shadow_price': {'value': -0.0}}, 'MNXM188': {'biomass_shadow_price': {'value': -0.13370910645572764}, 'fraction_shadow_price': {'value': -0.9419354838709677}}, 'MNXM1': {'biomass_shadow_price': {'value': 0.0009416134257445724}, 'fraction_shadow_price': {'value': 0.0064516129032255835}}, 'CMPD_0000000003': {'biomass_shadow_price': {'value': 0.0}, 'fraction_shadow_price': {'value': -0.9870967741935488}}}, 'reactions': {'rxn_1': {'biomass': {'value': 0.7638744755010182, 'units': 'milimole / gDW / hour'}, 'fraction': {'value': 1.3296695186776557, 'units': 'milimole / gDW / hour'}}, 'rxn_2': {'biomass': {'value': 0.7638744755010182, 'units': 'milimole / gDW / hour'}, 'fraction': {'value': 1.3296695186776557, 'units': 'milimole / gDW / hour'}}}, 'pathway': {'biomass': {'value': 0.7638744755010182, 'units': 'milimole / gDW / hour'}, 'fraction': {'value': 1.3296695186776557, 'units': 'milimole / gDW / hour'}}}
    if pathway is None:
      logger.info('No results written. Exiting...')
    else:
      logger.info('Writing into file...')
      # Write results into the pathway
      write_results(
        pathway,
        results,
        args.sim,
        logger
      )
      rpSBML.from_Pathway(pathway).write_to_file(args.outfile)
      logger.info('   |--> written in ' + args.outfile)


def write_results(
  pathway: rpPathway,
  results: Dict,
  sim: str,
  logger: Logger = getLogger(__name__)
) -> None:

  # Write species results
  for spe_id, score in results['species'].items():
    for k, v in score.items():
      pathway.get_specie(spe_id).set_fba_info(
        key=k,
        value=v
      )
  # Write reactions results
  for rxn_id, score in results['reactions'].items():
    for k, v in score.items():
      pathway.get_reaction(rxn_id).set_fba_info(
        key=k,
        value=v
      )
  # Write pathway result
  for k, v in results['pathway'].items():
    pathway.set_fba_info(
      key=k,
      value=v
    )

if __name__ == '__main__':
    entry_point()
