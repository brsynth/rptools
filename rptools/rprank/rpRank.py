from typing import List
from rptools.rplibs import rpPathway


def rank(pathways: List[rpPathway]) -> List[str]:
    _pathways = {}
    for pathway in pathways:
        _pathways[pathway.get_id().replace(' ', '_')] = str(pathway.get_global_score())
    sorted_pathways = dict(
        sorted(
            _pathways.items(),
            key=lambda item: float(item[1]),
            reverse=True
        )
    )
    return sorted_pathways
