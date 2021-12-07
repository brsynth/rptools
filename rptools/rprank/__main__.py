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
from .rpRank import rank
from .Args import add_arguments
from rptools.rplibs import rpPathway
from rptools import build_args_parser


def entry_point():
  
    parser = build_args_parser(
        prog='rprank',
        description='Rank pathways',
        m_add_args=add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    # Build the list of pathways to rank
    pathways = [
        rpPathway.from_rpSBML(
            infile=pathway_filename,
            logger=logger
        ) for pathway_filename in args.pathways
    ]

    # Rank pathways
    sorted_pathways = rank(pathways)

    sorted_pathways_str = '\n'.join(
        args.delimiter.join(item) for item
        in sorted_pathways.items()
    )
    print(f'#Name{args.delimiter}Score')
    print(sorted_pathways_str)


if __name__ == '__main__':
    entry_point()
