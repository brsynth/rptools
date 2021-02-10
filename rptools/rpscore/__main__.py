from rptools.rpscore import (
    compute_globalscore,
    build_args_parser
)
from brs_utils import create_logger
from rptools.rplibs import rpSBML
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
                thermo_id = args.thermo_id
    )

    score = rpsbml_dict['pathway']['brsynth']['global_score']

    # rpsbml.setGlobalScore(score)

    if args.outfile != '':
        rpsbml.updateBRSynthPathway(rpsbml_dict, args.pathway_id)
        rpsbml.writeSBML(args.outfile)

    print(score)


if __name__ == '__main__':
    entry_point()
