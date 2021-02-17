from rptools.rpranker import (
    rank,
)
from rptools.rpranker.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import rpSBML
from os import (
    path as os_path,
    mkdir
)
from shutil import copy
from typing import (
    List,
    Dict,
    Tuple
)
from tempfile import TemporaryDirectory
from brs_utils import compress_tar_gz
from tarfile import open as tf_open


def entry_point():
  
    parser = build_args_parser(
        prog = 'rpranker',
        description = 'Rank pathways accoring to their global score',
        m_add_args = add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    # Process to the ranking
    ranked_pathways = rank(
        pathways = args.pathways,
        logger = logger
    )

    # Write into a file the list of top ranked pathways
    if args.rank_outfile != '':
        store_paths_into_file(
            ranked_pathways[:args.top],
            args.rank_outfile
        )
        logger.info(
            '\nThe {top} pathways with highest scores are available in file {file}.'.format(
                top = args.top,
                file = args.rank_outfile
            )
        )

    # Write into an archive the top ranked pathways
    if args.data_outfile != '':
        store_into_tar_gz_file(
            ranked_pathways[:args.top],
            args.data_outfile
        )
        logger.info(
            '\nThe {top} pathway paths with highest scores are available in file {file}.'.format(
                top = args.top,
                file = args.data_outfile
            )
        )

    # # Copy top ranked rpsbml files into a folder
    # if args.outdir != '':
    #     copy_into_folder(
    #         ranked_pathways[:args.top],
    #         args.outdir
    #     )
    #     logger.info(
    #         '\nTop {top} pathways are available in folder {folder}.'.format(
    #             top = args.top,
    #             folder = args.outdir
    #         )
    #     )

    if not args.silent:
        if args.log.lower() in ['critical', 'error', 'warning']:
            print(ranked_pathways)
        else:
            logger.info('\nRanked Pathways')
            logger.info('   |-' + '\n   |-'.join('{}: {}'.format(*k) for k in enumerate(ranked_pathways)))


def store_paths_into_file(
    pathways: List[ Tuple[float, str] ],
    outfile: str
) -> None:
    with open(outfile, 'w') as fp:
        for item in pathways:
            fp.write(
                '{score} {filename}\n'.format(
                    score  = item[0],
                    filename = item[1]
                )
            )


def store_into_tar_gz_file(
    pathways: List[ Tuple[float, str] ],
    outfile: str
) -> None:
    with TemporaryDirectory() as temp_d:
        with tf_open(outfile, 'w:gz') as tar:
            for item in pathways:
                path = item[1]
                tar.add(
                    path,
                    arcname = os_path.basename(path)
                )


# def copy_into_folder(
#     pathways: List[ Tuple[float, str] ],
#     outdir: str
# ) -> None:
#     if not os_path.exists(outdir):
#         mkdir(outdir)
#     for item in pathways:
#         copy(item[1], outdir)


if __name__ == '__main__':
    entry_point()
