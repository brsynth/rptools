def rank(pathways: dict) -> list[str]:
    _pathways = {}
    for pathway_name, pathway in pathways.items():
        _pathways[pathway_name] = str(pathway["pathway"].get_global_score())
    sorted_pathways = dict(
        sorted(_pathways.items(), key=lambda item: float(item[1]), reverse=True)
    )
    return sorted_pathways
