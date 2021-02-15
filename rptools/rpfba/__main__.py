from rptools.rpfba import (
    runFBA,
    build_args_parser
)


def entry_point():
    parser = build_args_parser()
    args   = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

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
      logger.info('Writing into file...')
      rpsbml.writeToFile(args.outfile)
      logger.info('  |--> written in ' + args.outfile)


if __name__ == '__main__':
    entry_point()
