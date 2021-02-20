from rptools.rpranker import (
    rank,
)
from rptools.rpranker.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import rpSBML
from os import (
    path as os_path,
    mkdir,
    rename as os_rename
)
from shutil import (
    copy,
    copyfile
)
from typing import (
    List,
    Dict,
    Tuple
)
from tempfile import TemporaryDirectory
from brs_utils import compress_tar_gz
from tarfile import open as tf_open
from logging import (
    Logger,
    getLogger
)


def entry_point():
  
    parser = build_args_parser(
        prog = 'rpranker',
        description = 'Rank pathways accoring to their global score',
        m_add_args = add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    # Check arguments
    if args.rename:
        if not args.data_outdir and not args.data_outfile:
            logger.error('\'rename\' arg needs --data_outdir or --data_outfile arg. Exiting...')
            return
            
    # Process to the ranking
    ranked_pathways = rank(
        pathways = args.pathways,
        logger = logger
    )

    if not args.silent:
        logger.info('\nWriting results...')

    # Write into a file the list of top ranked pathways
    if args.rank_outfile != '':
        store_paths_into_file(
            ranked_pathways[:args.top],
            args.rank_outfile
        )
        if not args.silent:
            logger.info(
                '\nThe {top} pathways with highest scores are available in file {file}.'.format(
                    top = args.top,
                    file = args.rank_outfile
                )
            )

    # Write into a folder the top ranked pathways
    if args.data_outdir != '':
        store_paths_into_folder(
            pathways = ranked_pathways[:args.top],
            outdir = args.data_outdir,
            rename = args.rename
        )
        if not args.silent:
            logger.info(
                '\nThe {top} pathways with highest scores are available in folder {folder}.'.format(
                    top = args.top,
                    folder = args.data_outdir
                )
            )

    # Write into an archive the top ranked pathways
    if args.data_outfile != '':
        store_into_tar_gz_file(
            ranked_pathways[:args.top],
            args.data_outfile,
            args.rename
        )
        if not args.silent:
            logger.info(
                '\nThe {top} pathway paths with highest scores are available in file {file}.'.format(
                    top = args.top,
                    file = args.data_outfile
                )
            )

    if not args.silent:
        if args.log.lower() in ['critical', 'error', 'warning']:
            print(ranked_pathways)
        else:
            logger.info('\nRanked Pathways')
            logger.info('   |-' + '\n   |-'.join(
                '{}: {}'.format(*k) for k in enumerate(ranked_pathways))
            )


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


def store_paths_into_folder(
    pathways: List[ Tuple[float, str, str] ],
    outdir: str,
    rename: bool = False,
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Store pathways into a folder.

    Parameters
    ----------
    pathways: List[ Tuple[float, str, str] ]
        Pathway (score, filenames, rpSBML name).
    outdir: str
        Folder to store files into.
    rename: bool
        Have files to be renamed with rpsbml names
    logger : Logger
        The logger object.

    Returns
    -------
    None
    """
    if not os_path.exists(outdir):
        mkdir(outdir)

    for i in range(len(pathways)):
        orig_filename = pathways[i][1]
        sbml_filename = pathways[i][2]
        prefix = str(i) + ' - '

        if rename:
            outfile = os_path.join(
                outdir,
                prefix + sbml_filename
            )
        else:
            outfile = os_path.join(
                outdir,
                prefix + os_path.basename(orig_filename)
            )

        copy(orig_filename, outfile)

        # # Prepend ranking index to filename    
        # os_rename(
        #     outfile,
        #     os_path.join(
        #         os_path.dirname(outfile),
        #         str(i)+'-'+os_path.basename(outfile)
        #     )
        # )


def store_into_tar_gz_file(
    pathways: List[ Tuple[float, str, str] ],
    outfile: str,
    rename: bool = False,
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Store pathways into a folder.

    Parameters
    ----------
    pathways: List[ Tuple[float, str, str] ]
        Pathway (score, filenames, rpSBML name).
    rename: bool
        Have files to be renamed with rpsbml names
    logger : Logger
        The logger object.

    Returns
    -------
    None
    """
    with TemporaryDirectory() as temp_d:
        with tf_open(outfile, 'w:gz') as tar:
            for i in range(len(pathways)):
                path = pathways[i][1]
                prefix = str(i) + ' - '
                if rename:
                    arcname = os_path.basename(pathways[i][2])
                else:
                    arcname = os_path.basename(path)
                tar.add(
                    path,
                    arcname = prefix + arcname
                )


if __name__ == '__main__':
    entry_point()
