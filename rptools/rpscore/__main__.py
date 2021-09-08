from typing import (
    Dict,
    List
)
from logging import (
    Logger,
    getLogger
)
from os import (
    path as os_path,
    makedirs
)
from tempfile import NamedTemporaryFile
from rptools.rpscore import predict_score
from rptools.rpscore.Args import add_arguments
from rptools.rplibs import rpPathway
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

    # if len(args.pathways) == 1:
    #   if args.outfile is None or args.outfile == '':
    #     logger.error('Option --outfile has to be set in case of single input pathway, exiting...')
    #     exit(1)

    # pathways = []
    # for pathway in args.pathways:
    #     pathways.append(
    #         rpPathway.from_rpSBML(
    #             infile=pathway,
    #             logger=logger
    #         )
    #     )

    pathway = rpPathway.from_rpSBML(
        infile=args.infile,
        logger=logger
    )

    score = predict_score(
        pathway=pathway,
        # data_train_file=args.data_train_file,
        # models_path=models_path,
        no_of_rxns_thres=args.no_of_rxns_thres
    )

    # if len(pathways) > 1:
    #     if not os_path.exists(args.outdir):
    #         makedirs(args.outdir)
    #     for i in range(len(pathways)):
    #         # Write results into the pathway
    #         pathways[i].set_global_score(
    #             scores[i]
    #         )
    #         # Write pathway into file
    #         pathways[i].to_rpSBML().write_to_file(
    #             os_path.join(
    #                 args.outdir,
    #                 os_path.basename(args.pathways[i])
    #             )
    #         )
    #     else:
    # Write results into the pathway
    pathway.set_global_score(score)
    # Write pathway into file
    pathway.to_rpSBML().write_to_file(
        args.outfile
    )


if __name__ == '__main__':
    entry_point()
