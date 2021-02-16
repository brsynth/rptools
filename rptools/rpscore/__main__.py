from rptools.rpscore import (
    compute_globalscore,
)
from rptools.rpscore.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import rpSBML


def entry_point():
  
    parser = build_args_parser(
        prog = 'rpscore',
        description = 'Calculate global score by combining all scores (rules, FBA, Thermo)',
        m_add_args = add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    rpsbml = rpSBML(
        inFile = args.pathway_file,
        logger = logger
    )

    rpsbml_dict = compute_globalscore(
                   rpsbml = rpsbml,
          weight_rp_steps = args.weight_rp_steps,
        weight_rule_score = args.weight_rule_score,
               weight_fba = args.weight_fba,
            weight_thermo = args.weight_thermo,
             max_rp_steps = args.max_rp_steps,
              thermo_ceil = args.thermo_ceil,
             thermo_floor = args.thermo_floor,
                 fba_ceil = args.fba_ceil,
                fba_floor = args.fba_floor,
               pathway_id = args.pathway_id,
             objective_id = args.objective_id,
                thermo_id = args.thermo_id,
                   logger = logger
    )

    score = rpsbml_dict['pathway']['brsynth']['global_score']

    if args.outfile != '':
        rpsbml.updateBRSynthPathway(rpsbml_dict, args.pathway_id)
        rpsbml.writeToFile(args.outfile)

    if not args.silent:
        if args.log.lower() in ['critical', 'error', 'warning']:
            print(score)
        else:
            logger.info('\nGlobal Score = ' + str(score))


if __name__ == '__main__':
    entry_point()
