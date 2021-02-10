from rptools.rpfba import (
    runFBA,
    build_args_parser
)
from brs_utils import create_logger
from rptools._version import __version__

def entry_point():
    parser = build_args_parser()
    args   = parser.parse_args()

    # Create logger
    logger = create_logger(parser.prog, args.log)

    logger.info(
        '{prog} {version}\n'.format(
            prog = logger.name,
            version = __version__
        )
    )
    logger.debug(args)

    rpsbml = runFBA(
                  rpsbml_path = args.pathway,
                gem_sbml_path = args.model,
                     sim_type = args.sim,
                   src_rxn_id = args.source_reaction,
                   tgt_rxn_id = args.target_reaction,
                    src_coeff = args.source_coefficient,
                    tgt_coeff = args.target_coefficient,
                       is_max = args.is_max,
                  frac_of_src = args.fraction_of,
                   dont_merge = args.dont_merge,
                   pathway_id = args.pathway_id,
                 objective_id = args.objective_id,
               compartment_id = args.compartment_id,
             species_group_id = args.species_group_id,
        sink_species_group_id = args.sink_species_group_id,
                       logger = logger
    )

    if not rpsbml is None:
      rpsbml.writeSBML(args.outfile)
      logger.info('  |--> written in ' + args.outfile)



if __name__ == '__main__':
    entry_point()
