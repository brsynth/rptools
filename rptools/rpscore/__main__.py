from typing import (
  Dict,
  List
)
from logging import (
  Logger,
  getLogger
)
from rptools.rpscore import (
    load_training_data,
    predict_score
)
from rptools.rpscore.Args import (
  add_arguments,
  __MODELS_PATH as models_path
)
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

    # rpsbml = rpSBML(
    #     inFile = args.pathway_file,
    #     logger = logger
    # )

    features_dset_train = load_training_data(args.data_train_file)
    predict_score(
      args.test_data_file,
      args.test_score_file,
      args.data_predict_file,
      models_path,
      features_dset_train,
      no_of_rxns_thres=10
    )


    # if pathway is None:
    #   logger.info('No results written. Exiting...')
    # else:
    #   logger.info('Writing into file...')
    #   # Write results into the pathway
    #   write_results(
    #     pathway,
    #     global_score,
    #     logger
    #   )
    #   # Write pathway into file
    #   if args.outfile is not None or args.outfile != '':
    #     with open(args.outfile, 'w') as fp:
    #         pathway = json_dump(pathway, fp, indent=4)
    #     logger.info('   |--> written in ' + args.outfile)

    # if not args.silent:
    #     if args.log.lower() in ['critical', 'error', 'warning']:
    #         print(global_score)
    #     else:
    #         logger.info('\nGlobal Score = ' + str(global_score))


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
