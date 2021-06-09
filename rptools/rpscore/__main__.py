from json import (
  load as json_load,
  dump as json_dump
)
from typing import (
  Dict,
  List
)
from logging import (
  Logger,
  getLogger
)
from rptools.rpscore import (
    compute_globalscore,
)
from rptools.rpscore.Args import add_arguments
from rptools.rplibs import rpSBML
from rptools import build_args_parser


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

    with open(args.pathway, 'r') as fp:
      pathway = json_load(fp)

    global_score = compute_globalscore(
                  pathway = pathway,
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

    if pathway is None:
      logger.info('No results written. Exiting...')
    else:
      logger.info('Writing into file...')
      # Write results into the pathway
      write_results(
        pathway,
        global_score,
        logger
      )
      # Write pathway into file
      if args.outfile is not None or args.outfile != '':
        with open(args.outfile, 'w') as fp:
            pathway = json_dump(pathway, fp, indent=4)
        logger.info('   |--> written in ' + args.outfile)

    if not args.silent:
        if args.log.lower() in ['critical', 'error', 'warning']:
            print(global_score)
        else:
            logger.info('\nGlobal Score = ' + str(global_score))


def write_results(
  pathway: Dict,
  score: float,
  logger: Logger = getLogger(__name__)
) -> None:
  # Write pathway result
  if 'scores' not in pathway:
    pathway['scores'] = {}
  pathway['scores']['global'] = score


if __name__ == '__main__':
    entry_point()
