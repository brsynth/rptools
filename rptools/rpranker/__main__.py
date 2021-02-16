from rptools.rpranker import (
    rank,
)
from rptools.rpranker.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import rpSBML


def entry_point():
  
    parser = build_args_parser(
        prog = 'rpranker',
        description = 'Rank pathways accoring to their global score',
        m_add_args = add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    ranked_pathways = rank(
        pathways = args.pathways,
        logger = logger
    )

    if args.outfile != '':
        with open(args.outfile, 'w') as fp:
            fp.write(str(ranked_pathways))

    if not args.silent:
        if args.log.lower() in ['critical', 'error', 'warning']:
            print(ranked_pathways)
        else:
            logger.info('\nRanked Pathways')
            logger.info('   |-' + '\n   |-'.join('{}: {}'.format(*k) for k in enumerate(ranked_pathways)))


if __name__ == '__main__':
    entry_point()
