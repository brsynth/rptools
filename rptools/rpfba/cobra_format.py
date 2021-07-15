from logging import (
    Logger,
    getLogger
)
from typing import (
    Dict
)
from copy import deepcopy

at_pattern = '__64__'
bigg_prefix = 'M_'

def to_cobra(string: str) -> str:
    if string.startswith(bigg_prefix):
        string = string[len(bigg_prefix):]
    return string.replace(at_pattern, '@')

def cobra_suffix(comp_id: str) -> str:
    return at_pattern+comp_id

def cobraize(
    string: str,
    comp_id: str
) -> str:
    suffix = cobra_suffix(comp_id)
    if string.endswith(suffix):
        return string
    else:
        return string+cobra_suffix(comp_id)

def uncobraize(
    string: str
) -> str:
    return string.split(at_pattern)[0]

def uncobraize_results(
    results: Dict,
    cobra_suffix: str,
    logger: Logger = getLogger(__name__)
) -> None:

    res = {
        'species': {},
        'reactions': {},
        'pathway': {},
        'ignored_species': []
    }
    logger.debug(cobra_suffix)
    logger.debug(results)
    # Uncobraize species results
    for spe_id, score in results['species'].items():
        res['species'][spe_id.replace(cobra_suffix, '')] = score

    # Uncobraize rpfba_ignored_species results
    for spe_id in results['ignored_species']:
        res['ignored_species'] += [spe_id.replace(cobra_suffix, '')]

    # Copy other results
    for key in ['reactions', 'pathway']:
        res[key] = deepcopy(results[key])
    return res

