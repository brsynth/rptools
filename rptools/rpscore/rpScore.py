from logging import (
    Logger,
    getLogger
)
from rptools.rplibs      import rpSBML
from typing import (
    List,
    Dict,
    Tuple
)
from copy import deepcopy
from numpy import (
    average as np_avg,
    std as np_std
)


## Extract the reaction SMILES from an SBML, query rule_score and write the results back to the SBML
#
# Higher is better
#
# NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
# Higher is better
# TODO: try to standardize the values instead of normalisation.... Advantage: not bounded
def compute_globalscore(
    rpsbml: rpSBML,
    weight_rp_steps: float = 0.10002239003499142,
    weight_rule_score: float = 0.13346271414277305,
    weight_fba: float = 0.6348436269211155,
    weight_thermo: float = 0.13167126890112002,
    max_rp_steps: int = 15, # TODO: add this as a limit in RP2
    thermo_ceil: float = 5000.0,
    thermo_floor: float = -5000.0,
    fba_ceil: float = 5.0,
    fba_floor: float = 0.0,
    pathway_id: str = 'rp_pathway',
    objective_id: str = 'obj_fraction',
    thermo_id: str = 'dfG_prime_m',
    logger: Logger = getLogger(__name__)
) -> Dict:
    """From a rpsbml object, retreive the different characteristics of a pathway and combine them to calculate a global score.

    Note that the results are written to the rpsbml directly

    :param rpsbml: rpSBML object
    :param weight_rp_steps: The weight associated with the number of steps (Default: 0.10002239003499142)
    :param weight_rule_score: The weight associated with the mean of reaction rule scores (Default: 0.13346271414277305)
    :param weight_fba: The weight associated with the flux of the target (Default: 0.6348436269211155)
    :param weight_thermo: The weight associated with the sum of reaction Gibbs free energy (Default: 0.13167126890112002)
    :param max_rp_steps: The maximal number of steps are run in RP2 (Default: 15)
    :param thermo_ceil: The upper limit of Gibbs free energy for each reaction (Default: 5000.0)
    :param thermo_floor: The lower limit of Gibbs free energy for each reaction (Default: -5000.0)
    :param fba_ceil: The upper flux limit of the heterologous pathway (Default: 5.0)
    :param fba_floor: The lower flux limit of the heterologous pathway (Default: 5.0)
    :param pathway_id: The ID of the heterologous pathway (Default: rp_pathway)
    :param objective_id: The ID of the FBA objective (Default: obj_fraction)
    :param thermo_id: The ID of the Gibbs free energy that may be used. May be either dfG_prime_m or dfG_prime_o (Default: dfG_prime_m)

    :rtype: float
    :return: The global score
    """

    # WARNING: we do this because the list gets updated
    logger.debug('thermo_ceil:  '+str(thermo_ceil))
    logger.debug('thermo_floor: '+str(thermo_floor))
    logger.debug('fba_ceil:     '+str(fba_ceil))
    logger.debug('fba_floor:    '+str(fba_floor))

    bounds = {
        'thermo': {
            'floor': thermo_floor,
            'ceil' : thermo_ceil
        },
        'fba': {
            'floor': fba_floor,
            'ceil' : fba_ceil
        },
        'max_rp_steps': max_rp_steps
    }
    # Dict to store list of scores over reactions
    scores = {}
    # score_names = ['dfG_prime_m', 'dfG_uncert', 'dfG_prime_o', 'rule_score', 'fba_obj_biomass', 'fba_obj_fraction']
    rpsbml_dict = rpsbml.toDict(pathway_id)

    rpsbml_dict, scores = score_from_reactions(
        rpsbml_dict,
        scores,
        bounds,
        pathway_id,
        logger
    )

    rpsbml_dict = score_from_pathway(
        rpsbml_dict,
        scores,
        bounds,
        pathway_id,
        logger
    )

    #################################################
    ################# GLOBAL ########################
    #################################################

    ##### global score #########
    weights = {
        'norm_rule_score': weight_rule_score,
        'norm_'+str(thermo_id): weight_thermo,
        'norm_steps': weight_rp_steps,
        'norm_fba_'+str(objective_id): weight_fba
    }
    scores_l = weights_l = []
    logger.info('Checking scores...')
    for key in weights.keys():
        try:
            score = rpsbml_dict['pathway']['brsynth'][key]
            scores_l += [score]
            weights_l += [weights[key]]
            logger.info('   |- \'' + key + '\' added...')
        except KeyError as e:
            key = str(e)
            logger.debug('KeyError: ' + key)
            logger.warning('   |- ' + key + ' not found...')
            # msg = 'It seems that '
            # if key.startswith('\'norm_dfG_'):
            #     msg += 'Thermo '
            # elif key.startswith('\'norm_fba_'):
            #     msg += 'FBA '
            # msg += 'step has not been completed '
            # logger.warning('   |- ' + msg)

    globalScore = np_avg(
        scores_l,
        weights = weights_l
    )

    rpsbml_dict['pathway']['brsynth']['global_score'] = globalScore

    return rpsbml_dict


def score_from_reactions(
    rpsbml_dict: Dict,
    scores: Dict,
    bounds: Dict,
    pathway_id: str,
    logger: Logger = getLogger(__name__)
) -> Tuple[Dict, Dict]:

    rpsbml_dict_copy = deepcopy(rpsbml_dict)
    scores_copy      = deepcopy(scores)

    for reac_id in list(rpsbml_dict_copy['reactions'].keys()):

        for bd_id in list(
            rpsbml_dict_copy['reactions'][reac_id]['brsynth'].keys()
        ):

            thermo     = bd_id.startswith('dfG_')
            fba        = bd_id.startswith('fba_')
            rule_score = bd_id=='rule_score'

            if thermo or fba:

                try:
                    ####### Thermo ############
                    # lower is better -> -1.0 to have highest better
                    # WARNING: we will only take the dfG_prime_m value
                    ####### FBA ##############
                    # higher is better
                    # return all the FBA values
                    # ------- reactions ----------
                    if bd_id not in scores_copy:
                        scores_copy[bd_id] = []

                    value = rpsbml_dict_copy['reactions'][reac_id]['brsynth'][bd_id]['value']
                    floor = bounds['thermo']['floor'] if thermo else bounds['fba']['floor']
                    ceil  = bounds['thermo']['ceil']  if thermo else bounds['fba']['ceil']
                    norm_score = minmax_score(value, floor, ceil)

                except (KeyError, TypeError) as e:
                    logger.warning('Cannot find: '+str(bd_id)+' for the reaction: '+str(reac_id))
                    norm_score = 0.0

                if thermo:
                    norm_score = 1 - norm_score
                rpsbml_dict_copy['reactions'][reac_id]['brsynth']['norm_'+bd_id] = norm_score
                scores_copy[bd_id].append(norm_score)

            elif rule_score:
                if bd_id not in scores_copy:
                    scores_copy[bd_id] = []
                # rule score higher is better
                scores_copy[bd_id].append(
                    rpsbml_dict_copy['reactions'][reac_id]['brsynth'][bd_id]
                )
            
            else:
                logger.debug('Not normalising: '+str(bd_id))

    return rpsbml_dict_copy, scores_copy


def score_from_pathway(
    rpsbml_dict: Dict,
    scores: Dict,
    bounds: Dict,
    pathway_id: str,
    logger: Logger = getLogger(__name__)
) -> Dict:

    rpsbml_dict_copy = deepcopy(rpsbml_dict)

    for bd_id in scores:

        ############### FBA ################
        # higher is better
        if bd_id.startswith('fba_'):
            rpsbml_dict_copy['pathway']['brsynth']['norm_'+bd_id] = \
                minmax_score(
                    rpsbml_dict_copy['pathway']['brsynth'][bd_id]['value'],
                    bounds['fba']['floor'], bounds['fba']['ceil']
                )

        ############# Thermo ################
        elif bd_id.startswith('dfG_'):
            # here add weights based on std
            rpsbml_dict_copy['pathway']['brsynth']['norm_'+bd_id] = \
                np_avg(
                    [np_avg(scores[bd_id]),
                    1.0-np_std(scores[bd_id])],
                    weights = [0.5, 0.5]
                )
            # the score is higher is better - (-1 since we want lower variability)
            # rpsbml_dict['pathway']['brsynth']['var_'+bd_id] = 1.0-np.var(path_norm[bd_id])

    bd_id = 'rule_score'
    ############# rule score ############
    # higher is better
    if not bd_id in scores:
        logger.warning('Cannot detect rule_score: '+str(scores))
        rpsbml_dict_copy['pathway']['brsynth']['norm_'+bd_id] = 0.0
    else:
        rpsbml_dict_copy['pathway']['brsynth']['norm_'+bd_id] = np_avg(scores[bd_id])

    ##### length of pathway ####
    # lower is better -> -1.0 to reverse it
    norm_steps = 0.0
    if len(rpsbml_dict_copy['reactions']) > bounds['max_rp_steps']:
        logger.warning('There are more steps than specified')
        norm_steps = 0.0
    else:
        try:
            norm_steps = (
                float(len(rpsbml_dict_copy['reactions']))-1.0
            ) / (
                float(bounds['max_rp_steps'])-1.0
            )
            norm_steps = 1.0 - norm_steps
        except ZeroDivisionError:
            norm_steps = 0.0

    rpsbml_dict_copy['pathway']['brsynth']['norm_steps'] = norm_steps

    return rpsbml_dict_copy


def minmax_score(
    value: float,
    floor: float,
    ceil: float,
    logger: Logger = getLogger(__name__)
) -> float:
    """Compute and returns score
    """
    if ceil >= value >= floor:
        # min-max feature scaling
        norm = (value-floor) / (ceil-floor)
    elif value < floor:
        norm = 0.0
    # then value > ceil
    else:
        norm = 1.0
    # rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value'] = norm

    return norm


