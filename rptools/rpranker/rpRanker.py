from logging import (
    Logger,
    getLogger
)
from rptools.rplibs import rpSBML
from typing import (
    List,
    Dict,
    Tuple
)


def rank(
    pathways: List[str],
    logger: Logger = getLogger(__name__)
) -> List[ Tuple[float, str] ]:
    """
    From a list of pathway (rpsbml) filenames, rank them according to global score.

    Parameters
    ----------
    pathways: List[str]
        Pathway (rpSBML) filenames
    logger : Logger
        The logger object.

    Returns
    -------
    pathways: List[ Tuple[float, str] ]
        List of tuple (global score, pathway filenames) sorted according to global score
    """

    sorted_pathways = []

    for pathway in pathways:
        rpsbml = rpSBML(
            inFile = pathway,
            logger = logger
        )
        logger.info('Pathway {rp_name}'.format(rp_name = rpsbml.getName()))
        logger.info('   |- loaded from ' + pathway)
        score = rpsbml.read_global_score()
        sorted_pathways += [ (score, pathway) ]
        logger.info('   |- ranked with score: {score}'.format(score = score))

    sorted_pathways.sort(
        reverse = True,
        key = lambda r: r[0]
    )

    return sorted_pathways