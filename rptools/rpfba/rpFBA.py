import csv
import pickle
import sys

import pandas as pd
from copy import deepcopy
from logging import (
    Logger,
    getLogger
)
from os import remove
from pandas.core.series  import Series    as np_series
from typing import (
    List,
    Dict,
    Tuple
)
from tempfile import (
    NamedTemporaryFile
)
from json import dumps as json_dumps
from cobra.flux_analysis import pfba
from cobra               import io        as cobra_io
from cobra.io.sbml       import (
    validate_sbml_model,
    CobraSBMLError
)
from cobra.core.model    import Model     as cobra_model
from cobra.core.solution import Solution  as cobra_solution
 
from rptools.rpfba.medium import (
    add_missing_specie,
    build_minimal_medium,
    crossref_medium_id,
    df_to_medium,
    is_df_medium_defined,
    merge_medium,
    merge_medium_exchange
)
from rptools.rplibs.rpCompound import rpCompound
from rptools.rplibs.rpSBML import rpSBML
from rptools.rplibs.rpPathway import rpPathway
from rptools.rplibs.rpReaction import rpReaction
from .cobra_format import (
    cobraize,
    uncobraize_results,
    cobra_suffix,
    to_cobra
)

# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes

def runFBA(
    pathway: rpPathway,
    gem_sbml_path: str,
    compartment_id: str,
    objective_rxn_id: str = 'rxn_target',
    biomass_rxn_id: str = 'biomass',
    sim_type: str = 'fraction',
    fraction_coeff: float = 0.75,
    merge: bool = False,
    ignore_orphan_species: bool = True,
    medium_compartment_id: str = 'MNXC2',
    df_medium_base: pd.DataFrame = None,
    df_medium_user: pd.DataFrame = None,
    logger: Logger = getLogger(__name__)
) -> Dict:
    """Single rpSBML simulation

    :param file_name: The name of the model
    :param pathway_fn: Path to the pathway file (JSON)
    :param gem_sbml: Path to the GEM file
    :param sim_type: The type of simulation to use. Available simulation types include: fraction, fba, rpfba
    :param src_rxn_id: The reaction id of the source reaction.
    :param target_reaction: The reaction id of the target reaction. Note that if fba or rpfba options are used, then these are ignored
    :param source_coefficient: The source coefficient
    :param target_coefficient: The target coefficient
    :param is_max: Maximise or minimise the objective'last', i
    :param fraction_of: The fraction of the optimum. Note that this value is ignored is fba is used
    :param tmpOutputFolder: The path to the output document
    :param merge: Output the merged model (Default: False)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the auto-generated id of the results (Default: None)
    :param compartment_id: The SBML compartment id (Default: MNXC3)
    :param fill_orphan_species: Add pseudo reactions that consume/produce single parent species. Note in development
    :param species_group_id: The id of the central species (Default: central_species)
    :param sink_species_group_id: The id of the sink species (Default: rp_sink_species)
 not 
    :type inputTar: str
    :type gem_sbml_path: str
    :type sim_type: str
    :type src_rxn_id: str
    :type target_reaction: str
    :type source_coefficient: float
    :type target_coefficient: float
    :type is_max: bool
    :type fraction_of: float
    :type tmpOutputFolder: str
    :type merge: bool
    :type num_workers: int
    :type pathway_id: str
    :type objective_id: str
    :type compartment_id: str
    :type fill_orphan_species: bool
    :type species_group_id: str
    :type sink_species_group_id: str

    :return: Succcess or failure of the function
    :rtype: bool
    """

    logger.debug('           pathway_fn: ' + str(pathway))
    logger.debug('        gem_sbml_path: ' + str(gem_sbml_path))
    logger.debug('             sim_type: ' + str(sim_type))
    logger.debug('     objective_rxn_id: ' + objective_rxn_id)
    logger.debug('       biomass_rxn_id: ' + str(biomass_rxn_id))
    logger.debug('       fraction_coeff: ' + str(fraction_coeff))
    logger.debug('                merge: ' + str(merge))
    logger.debug('       compartment_id: ' + str(compartment_id))
    logger.debug('ignore_orphan_species: ' + str(ignore_orphan_species))


    ## MODEL
    # Create rpSBML object
    rpsbml_gem = rpSBML(inFile=gem_sbml_path, logger=logger)
    # Check compartment ID
    compartment_id = check_SBML_compartment(
        rpsbml=rpsbml_gem,
        compartment_id=compartment_id,
        logger=logger
    )
    if compartment_id is None:
        return None
    # Check biomass reaction ID
    biomass_rxn_id = check_SBML_rxnid(
        rpsbml=rpsbml_gem,
        rxn_id=biomass_rxn_id,
        logger=logger
    )
    if biomass_rxn_id is None:
        return None

    # PATHWAY
    rpsbml = build_rpsbml(
        pathway=pathway,
        logger=logger
    )
    # Check objective reaction ID
    objective_rxn_id = check_SBML_rxnid(
        rpsbml=rpsbml,
        rxn_id=objective_rxn_id,
        logger=logger
    )
    if objective_rxn_id is None:
        return None

    ## MERGE
    # Merge predicted pathway with the full model
    # missing_species are species that are not detected in the model
    logger.info('Merging rpSBML models: ' + rpsbml.getName() + ' and ' + rpsbml_gem.getName() + '...')
    (
        rpsbml_merged,
        reactions_in_both,
        missing_species,
        compartment_id
    ) = rpSBML.merge(
        pathway=rpsbml,
        model=rpsbml_gem,
        compartment_id=compartment_id,
        logger=logger
    )
    if rpsbml_merged is None:
        return None
    logger.debug('rpsbml_merged: ' + str(rpsbml_merged))
    logger.debug('reactions_in_both: ' + str(reactions_in_both))
    cobra_rpsbml_merged = rpSBML.cobraize(rpsbml_merged)

    ## MEDIUM
    df_medium = pd.DataFrame()
    if is_df_medium_defined(df_medium_base) or is_df_medium_defined(df_medium_user): 
        # Check medium compartment id.
        medium_compartment_id = check_SBML_compartment(
          rpsbml=cobra_rpsbml_merged,
          compartment_id=medium_compartment_id,
          logger=logger
        )
        if medium_compartment_id is None:
            logger.warning('Medium id not find in the model -> ignore modifications')
        else:
            # CrossRef ids
            df_medium_base = crossref_medium_id(
                df=df_medium_base,
                model=cobra_rpsbml_merged,
                compartment_id=medium_compartment_id,
                logger=logger
            )            
            df_medium_user = crossref_medium_id(
                df=df_medium_user,
                model=cobra_rpsbml_merged,
                compartment_id=medium_compartment_id,
                logger=logger
            )
            # Merge coumponds.
            df_medium = merge_medium(
                first=df_medium_base, 
                second=df_medium_user
            )

            # Select exchange reaction
            df_exchange_reaction = cobra_rpsbml_merged.build_exchange_reaction(
                compartment_id=medium_compartment_id
            )

            # Merge df medium with exchange reactions
            df_medium = merge_medium_exchange(
                medium=df_medium,
                exchange_reaction=df_exchange_reaction
            )

            # Add specie missing in the model
            cobra_rpsbml_merged = add_missing_specie(
                model=cobra_rpsbml_merged,
                df=df_medium,
                compartment_id=medium_compartment_id,
                logger=logger
            )

    # Detect orphan species among missing ones in the model,
    # i.e. that are only consumed or produced
    if ignore_orphan_species:
        hidden_species = build_hidden_species(
            rpsbml=rpsbml,
            missing_species=missing_species,
            compartment_id=compartment_id,
            logger=logger
        )
    else: hidden_species = []
    
    # NOTE: reactions is organised with key being the rpsbml reaction and value being the rpsbml_gem value`
    # BUG: when merging the rxn_sink (very rare cases) can be recognised if another reaction contains the same species as a reactant
    ## under such as scenario the algorithm will consider that they are the same -- TODO: overwrite it
    if objective_rxn_id in reactions_in_both:
        logger.warning(
            'The target_reaction ('+str(objective_rxn_id)+') ' \
          + 'has been detected in model ' + str(gem_sbml_path.getName()) + ', ' \
          + 'ignoring this model...'
        )
        return None

    ######## FBA ########
    results = {}
    if sim_type.lower() in ['fba', 'pfba']:
        objective_id = cobra_rpsbml_merged.find_or_create_objective(
            rxn_id=objective_rxn_id,
            obj_id=f'brs_obj_{objective_rxn_id}',
        )
        cobra_results, minimal_medium = runCobra(
            sim_type=sim_type,
            rpsbml=cobra_rpsbml_merged,
            hidden_species=hidden_species,
            objective_id=objective_id,
            fraction_coeff=fraction_coeff,
            medium=df_medium,
            logger=logger
        )
    else:
        (
            cobra_results,
            results_biomass,
            objective_id,
            minimal_medium
        ) = rp_fraction(
            rpsbml=cobra_rpsbml_merged,
            objective_rxn_id=objective_rxn_id,
            biomass_rxn_id=biomass_rxn_id,
            hidden_species=hidden_species,
            fraction_coeff=fraction_coeff,
            medium=df_medium,
            logger=logger
        )
        results['biomass'] = results_biomass

    # print(cobra_results.objective_value)
    # print(cobra_results.status)
    # print(cobra_results.fluxes)
    # print(cobra_results.shadow_prices)

    results[sim_type] = cobra_results

    # Write results for merged model
    write_results_to_rpsbml(
        rpsbml=cobra_rpsbml_merged,
        objective_id=objective_id,
        cobra_results=cobra_results,
        sim_type=sim_type,
        logger=logger
    )

    _results = build_results(
        results=results,
        pathway=pathway,
        compartment_id=compartment_id
    )
    _results['ignored_species'] = deepcopy(hidden_species)

    # Remove the Cobra standard ('compound@compartment') from all compounds
    pathway.uncobraize()
    _results = uncobraize_results(
        _results,
        cobra_suffix(compartment_id)
    )

    # Write results into the pathway
    write_results_to_pathway(
        pathway,
        _results,
        logger
    )

    _results['minimal_medium'] = minimal_medium
    return _results

    # if cobra_results is None:
    #     return None

    # '''
    # ###### multi objective #####
    # elif sim_type=='multi_fba':
    #     rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    # '''
    # if not merge:
    #     complete_heterologous_pathway(
    #         rpsbml = rpsbml,
    #         missing_species = missing_species,
    #         rpsbml_merged = rpsbml_merged,
    #         species_group_id = species_group_id,
    #         sink_species_group_id = sink_species_group_id,
    #         pathway_id = pathway_id,
    #         reactions_in_both = reactions_in_both,
    #         logger = logger
    #     )
    #     logger.info('Returning model with heterologous pathway only')
    #     pathway_fba = rpsbml.to_dict()
    # else:
    #     logger.info('Returning the full model')
    #     pathway_fba = rpsbml_merged.to_dict()

    # pathway_fba = rename_species(
    #     pathway_fba,
    #     logger
    # )

    # from json import dumps
    # print(dumps(pathway_fba, indent=4))
    # exit()


def build_rpsbml(
    pathway: rpPathway,
    logger: Logger = getLogger(__name__)
) -> rpSBML:
    # Create consumption of the target
    rxn_target = create_target_consumption_reaction(
        pathway.get_target_id(),
        logger
    )
    # Set Flux Bounds
    for rxn in pathway.get_list_of_reactions()+[rxn_target]:
        rxn.set_fbc(
            l_bound=0,
            u_bound=rpReaction.get_default_fbc_upper()
        )
        rxn.set_reversible(False)
    # Create rpSBML object
    rpsbml = pathway.to_rpSBML()
    # Create the target consumption reaction in the rpSBML
    rpsbml.createReaction(
        id=rxn_target.get_id(),
        reactants=rxn_target.get_reactants(),
        products=rxn_target.get_products(),
        smiles=rxn_target.get_smiles(),
        fbc_upper=rxn_target.get_fbc_upper(),
        fbc_lower=rxn_target.get_fbc_lower(),
        fbc_units=rxn_target.get_fbc_units(),
        reversible=rxn_target.reversible(),
        reacXref={'ec': rxn_target.get_ec_numbers()},
        infos=rxn_target._to_dict(full=False)
    )
    return rpsbml


def build_hidden_species(
    rpsbml: rpSBML,
    missing_species: List[str],
    compartment_id: bool,
    logger: Logger = getLogger(__name__)
) -> List[str]:
    rpsbml.search_isolated_species(missing_species)
    return [
        cobraize(spe_id, compartment_id) for
        spe_id in rpsbml.get_isolated_species()
    ]

def check_SBML_compartment(
    rpsbml: rpSBML,
    compartment_id: str,
    logger: Logger = getLogger(__name__)
) -> Tuple[str, str]:

    # Check model compartment ID
    # Set new compartment ID in case it exists under another ID
    _compartment_id = rpsbml.search_compartment_id(compartment_id)
    if _compartment_id is None:
        logger.debug(f'Compartment \'{compartment_id}\' not found in the model \'{rpsbml.getName()}\'')
    else:
        logger.debug(f'Compartment \'{compartment_id}\' found in the model \'{rpsbml.getName()}\'')
        return _compartment_id

    logger.error(f'Compartment \'{compartment_id}\' not found in the model \'{rpsbml.getName()}\'')
    logger.error(
        'Available compartments: {comp_ids}'.format(
            comp_ids=', '.join([
                comp.getId()
                for comp in list(rpsbml.getModel().getListOfCompartments())
            ])
        )
    )
    return None

def check_SBML_rxnid(
    rpsbml: rpSBML,
    rxn_id: str,
    logger: Logger = getLogger(__name__)
) -> Tuple[str, str]:

    # Check model reaction ID
    # Get reaction from the rpSBML
    _rxn = rpsbml.getModel().getReaction(rxn_id)
    if _rxn is None:
        logger.error(f'Reaction ID \'{rxn_id}\' not found in the model \'{rpsbml.getName()}\'')
        possible_rxn_ids = [
            rxn.getId()
            for rxn in list(rpsbml.getModel().getListOfReactions())
            if rxn_id in rxn.getId().lower()
        ]
        if possible_rxn_ids != []:
            logger.error(
                'Possible reactions: {rxn_ids}'.format(
                    rxn_ids=', '.join(possible_rxn_ids)
                )
            )
        return None

    return _rxn.getId()


def build_results(
    results: Dict,
    pathway: rpPathway,
    compartment_id: str,
) -> Dict:
    _results = {
        'species': {},
        'reactions': {},
        'pathway': {}
    }

    # SPECIES
    for spe_id in pathway.get_species_ids():
        _results['species'][spe_id] = {}
        for sim_type, cobra_r in results.items():
            value = cobra_r.shadow_prices.get(
                to_cobra(cobraize(spe_id, compartment_id))
            )
            _results['species'][spe_id][sim_type+'_shadow_price'] = {
                'value': value,
                # 'units': 'milimole / gDW / hour',
            }
    # REACTIONS
    for rxn_id in pathway.get_reactions_ids():
        _results['reactions'][rxn_id] = {}
        for sim_type, cobra_r in results.items():
            _results['reactions'][rxn_id][sim_type] = {
                'value': cobra_r.fluxes[rxn_id],
                'units': 'milimole / gDW / hour'
            }
            if sim_type == 'biomass':
                _results['reactions'][rxn_id][sim_type]['units'] = 'gDW / gDW / hour'
    # PATHWAY
    _results['pathway'] = {}
    for sim_type, cobra_r in results.items():
        _results['pathway'][sim_type] = {
            'value': cobra_r.objective_value,
            'units': 'milimole / gDW / hour'
        }
        if sim_type == 'biomass':
            _results['pathway'][sim_type]['units'] = 'gDW / gDW / hour'
    return _results

def create_target_consumption_reaction(
    target_id: str,
    logger: Logger = getLogger(__name__)
) -> 'rpReaction':
    rxn = rpReaction(
        id='rxn_target',
        logger=logger
    )
    rxn.add_reactant(
        compound_id=target_id,
        stoichio=1
    )
    return rxn


def write_results_to_pathway(
  pathway: rpPathway,
  results: Dict,
  logger: Logger = getLogger(__name__)
) -> None:
    # Write species results
    for spe_id, score in results['species'].items():
        for k, v in score.items():
            pathway.get_specie(spe_id).add_fba_info(
                key=k,
                value=v
            )
    # Write reactions results
    for rxn_id, score in results['reactions'].items():
        for k, v in score.items():
            pathway.get_reaction(rxn_id).add_fba_info(
                key=k,
                value=v
            )
    # Write pathway result
    for k, v in results['pathway'].items():
        pathway.add_fba_info(
            key=k,
            value=v
        )
    # Write ignored species
    pathway.add_species_group(
        'fba_ignored',
        results['ignored_species']
    )


def complete_heterologous_pathway(
    rpsbml: rpSBML,
    hidden_species: List[str],
    rpsbml_merged: rpSBML,
    species_group_id: str,
    sink_species_group_id: str,
    pathway_id: str,
    reactions_in_both: Dict,
    logger: Logger = getLogger(__name__)
) -> None:

    # Save the central species
    # groups = rpsbml.getPlugin('groups')
    central = rpsbml.getGroup(species_group_id)
    sink_group = rpsbml.getGroup(sink_species_group_id)
    rp_group = rpsbml.getGroup(pathway_id)
    cent_spe = [str(i.getIdRef()) for i in central.getListOfMembers()]
    sink_spe = [str(i.getIdRef()) for i in sink_group.getListOfMembers()]
    rp_reac  = [str(i.getIdRef()) for i in rp_group.getListOfMembers()]
    logger.debug('old central species: ' + str(cent_spe))
    logger.debug('old sink species:    ' + str(sink_spe))
    logger.debug('old rp reactions:    ' + str(rp_reac))

    rev_reactions = {v: k for k, v in reactions_in_both.items()}
    logger.debug('reactions_in_both: ' + str(reactions_in_both))
    logger.debug('rev_reactions:     ' + str(rev_reactions))
    logger.info('Building model with heterologous pathway only')
    groups = rpsbml_merged.getPlugin('groups')
    rp_pathway = rpsbml_merged.getGroup(pathway_id)
    logger.debug('---- Reactions ----')
    for member in rp_pathway.getListOfMembers():
        #### reaction annotation
        logger.debug(member.getIdRef())
        reacFBA = rpsbml_merged.getModel().getReaction(member.getIdRef())
        logger.debug(reacFBA)
        try:
            #reacIN = rpsbml.model.getReaction(reactions_convert[member.getIdRef()])
            reacIN = rpsbml.getModel().getReaction(rev_reactions[member.getIdRef()])
        except KeyError:
            reacIN = rpsbml.getModel().getReaction(member.getIdRef())
        logger.debug(reacIN)
        logger.debug(reacFBA.getAnnotation())
        reacIN.setAnnotation(reacFBA.getAnnotation())
        #### species TODO: only for shadow price
    #### add groups ####
    source_groups = rpsbml_merged.getPlugin('groups')
    target_groups = rpsbml.getPlugin('groups')
    target_groupsID = [i.getId() for i in rpsbml.getListOfGroups()]
    for source_group in rpsbml_merged.getListOfGroups():
        logger.debug('Replacing group id: '+str(source_group.getId()))
        if source_group.getId() == species_group_id:
            target_group = rpsbml.getGroup(source_group.getId())
            # TODO: #### replace the new potentially incorect central species with the normal ones #####
            # delete all the previous members
            logger.debug('Removing rp_core_species')
            for i in range(target_group.getNumMembers()):
                logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                target_group.removeMember(0)
            # add the new ones
            for cs in cent_spe:
                logger.debug('Creating new member: '+str(cs))
                newM = target_group.createMember()
                newM.setIdRef(cs)
        elif source_group.getId()==sink_species_group_id:
            target_group = rpsbml.getGroup(source_group.getId())
            logger.debug('Removing sink species')
            for i in range(target_group.getNumMembers()):
                logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                target_group.removeMember(0)
            #add the new ones
            for cs in sink_spe:
                logger.debug('Creating new member: '+str(cs))
                newM = target_group.createMember()
                newM.setIdRef(cs)
        elif source_group.getId() in target_groupsID:
            target_group = rpsbml.getGroup(source_group.getId())
            target_group.setAnnotation(source_group.getAnnotation())
        '''
        elif source_group.getId()==pathway_id:
            target_group = target_groups.getGroup(source_group.getId())
            logger.debug('Removing rp ractions')
            for i in range(target_group.getNumMembers()):
                logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                target_group.removeMember(0)
            #add the new ones
            for cs in rp_reac:
                logger.debug('Creating new member: '+str(cs))
                newM = target_group.createMember()
                newM.setIdRef(cs)
        '''
    # add ignored_species group
    rpsbml.set_isolated_species(rpsbml_merged.get_isolated_species())
    create_ignored_species_group(
        rpsbml,
        hidden_species,
        logger
    )
    #### add objectives ####
    source_fbc = rpsbml_merged.getPlugin('fbc')
    target_fbc = rpsbml.getPlugin('fbc')
    target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
    for source_obj in source_fbc.getListOfObjectives():
        source_obj_id = source_obj.getId()
        if source_obj.getId() in target_objID:
            target_obj = target_fbc.getObjective(source_obj.getId())
            target_obj.setAnnotation(source_obj.getAnnotation())
            for target_fluxObj in target_obj.getListOfFluxObjectives():
                for source_fluxObj in source_obj.getListOfFluxObjectives():
                    if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                        target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
        else:
            target_fbc.addObjective(source_obj)
    # rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
    target_fbc.setActiveObjectiveId(source_obj_id) # tmp random assignement of objective


# def rp_fba(
#           rpsbml: rpSBML,
#     objective_id: str,
#   hidden_species: List[str],
#           logger: Logger = getLogger(__name__)
# ) -> cobra_solution:
#     """Run FBA using a single objective

#     :param reaction_id: The id of the reactions involved in the objective
#     :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
#     :param is_max: Maximise or minimise the objective (Default: True)
#     :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
#     :param objective_id: Overwrite the default id (Default: None)

#     :type reaction_id: str
#     :type coefficient: float
#     :type is_max: bool
#     :type pathway_id: str
#     :type objective_id: str

#     :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
#     :rtype: tuple
#     """
#     logger.info('Running FBA...')
#     logger.debug('rpsbml:       ' + str(rpsbml))
#     logger.debug('hidden_species:  ' + str(hidden_species))

#     cobra_results = runCobra(
#         rpsbml=rpsbml,
#         objective_id=objective_id,
#         hidden_species=hidden_species,
#         logger=logger
#     )

#     return cobra_results


# def rp_pfba(
#           rpsbml: rpSBML,
#     objective_id: str,
#   hidden_species: List[str],
#      fraction_coeff: float = 0.95,
#           logger: Logger = getLogger(__name__)
# ) -> cobra_solution:
#     """Run parsimonious FBA using a single objective

#     :param reaction_id: The id of the reactions involved in the objective
#     :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
#     :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.95)
#     :param is_max: Maximise or minimise the objective (Default: True)
#     :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
#     :param objective_id: Overwrite the default id (Default: None)
#     :type reaction_id: str
#     :type coefficient: float
#     :type fraction_of_optimum: float
#     :type is_max: bool
#     :type pathway_id: str
#     :type objective_id: str

#     :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
#     :rtype: tuple
#     """
#     logger.info('Running FBA (parsimonious)...')
#     logger.debug('rpsbml:       ' + str(rpsbml))
#     logger.debug('fraction_coeff:  ' + str(fraction_coeff))

#     sim_type = 'pfba'
#     cobra_results = runCobra(
#         method=sim_type,
#         rpsbml=rpsbml,
#         hidden_species=hidden_species,
#         objective_id=objective_id,
#         fraction_coeff=fraction_coeff,
#         logger=logger
#     )
#     if cobra_results is None:
#         return None
#     return cobra_results


def rp_fraction(
    rpsbml: rpSBML,
    objective_rxn_id: str,
    biomass_rxn_id: str,
    hidden_species: List[str],
    fraction_coeff:  float = 0.75,
    medium: pd.DataFrame=None,
    logger: Logger = getLogger(__name__)
) -> cobra_solution:
    """Optimise for a target reaction while fixing a source reaction to the fraction of its optimum

    :param source_reaction: The id of the source reaction
    :param source_coefficient: The source coefficient associated with the source reaction id
    :param target_reaction: The id of the target reaction
    :param target_coefficient: The source coefficient associated with the target reaction id
    :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.75)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type source_reaction: str
    :type source_coefficient: float
    :type target_reaction: str
    :type target_coefficient: float
    :type fraction_of_optimum: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """

    logger.debug('rpsbml:       ' + str(rpsbml))
    logger.debug('hidden_species:   ' + str(hidden_species))
    logger.debug('objective_rxn_id:   ' + objective_rxn_id)
    logger.debug('biomass_rxn_id: ' + str(biomass_rxn_id))
    logger.debug('fraction_coeff:  ' + str(fraction_coeff))

    def get_annot_objective(
        rpsbml: rpSBML,
        objective_id: str
    ) -> str:

        fbc_plugin = rpsbml.getPlugin('fbc')
        objective = fbc_plugin.getObjective(objective_id)
        rpsbml.checklibSBML(
            objective,
            'Getting objective '+str(objective_id)
        )

        fbc_obj_annot = objective.getAnnotation()
        # TODO: if this is None need to set it up
        if fbc_obj_annot is None:
            return None

        logger.debug('Already calculated flux for '+str(objective_id))

        return fbc_obj_annot

    # retreive the biomass objective and flux results and set as maxima
    biomass_objective_id = rpsbml.find_or_create_objective(
        rxn_id=biomass_rxn_id,
        obj_id=f'brs_obj_{biomass_rxn_id}'
    )

    # fbc_plugin = rpsbml.getPlugin('fbc')
    # TODO: use the rpSBML BRSynth annotation parser
    # try:
    fbc_obj_annot = get_annot_objective(rpsbml, biomass_objective_id)

    # except (AttributeError, ValueError) as e:
    if fbc_obj_annot is None:
        # logger.debug(e)
        logger.debug('Performing FBA to calculate the source reaction')

        ### FBA ###
        # logger.info('Running the FBA (fraction of reaction)...')
        # rpsbml.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
        logger.info('Processing FBA (biomass)...')
        cobra_results, minimal_medium = runCobra(
            sim_type='biomass',
            rpsbml=rpsbml,
            hidden_species=hidden_species,
            objective_id=biomass_objective_id,
            medium=medium,
            logger=logger
        )

        # rxn>scores>fba>biomass = cobra_results.fluxes('rxn_X')
        # scores>fba>biomass = cobra_results.objective_value
        if cobra_results is None:
            return None, None, biomass_objective_id, minimal_medium

        results_biomass = cobra_results

        write_results_to_rpsbml(
            rpsbml=rpsbml,
            objective_id=biomass_objective_id,
            sim_type='biomass',
            cobra_results=cobra_results,
            logger=logger
        )

        fbc_obj_annot = get_annot_objective(rpsbml, biomass_objective_id)
        if fbc_obj_annot is None:
            logger.error('No annotation available for: '+str(biomass_objective_id))

    # print(cobra_results.objective_value)
    # rpsbml.write_to_file('joan.xml')
    # exit()

    flux = float(
        fbc_obj_annot.getChild(
            'RDF'
        ).getChild(
            'BRSynth'
        ).getChild(
            'brsynth'
        ).getChild(
            0
        ).getAttrValue(
            'value'
        )
    )

    # exit()
    # # TODO: add another to check if the objective id exists
    # logger.debug('FBA source flux ('+str(src_rxn_id)+') is: '+str(source_flux))
    # if not objective_id:
    #     objective_id = 'obj_'+str(tgt_rxn_id)+'__restricted_'+str(src_rxn_id)

    objective_id = rpsbml.find_or_create_objective(
        rxn_id=objective_rxn_id,
        obj_id=f'brs_obj_{objective_rxn_id}',
    )
    # objective_id = rpsbml.find_or_create_objective(
    #     reactions = [tgt_rxn_id],
    #     coefficients = [tgt_coeff],
    #     is_max = is_max,
    #     objective_id = objective_id
    # )

    logger.debug('Optimising the objective: '+str(biomass_rxn_id))
    logger.debug('     Setting upper bound: '+str(flux*fraction_coeff))
    logger.debug('     Setting lower bound: '+str(flux*fraction_coeff))

    old_upper_bound, old_lower_bound = rpsbml.getReactionConstraints(biomass_rxn_id)
    rpsbml.setReactionConstraints(
        biomass_rxn_id,
        flux*fraction_coeff,
        flux*fraction_coeff
    )

    logger.info('Processing FBA (fraction)...')
    sim_type = 'fraction'
    cobra_results, minimal_medium = runCobra(
        sim_type=sim_type,
        rpsbml=rpsbml,
        hidden_species=hidden_species,
        objective_id=objective_id,
        fraction_coeff=fraction_coeff,
        medium=medium,
        logger=logger
    )
    if cobra_results is None:
        return results_biomass, None, objective_id, minimal_medium

    # # rxn>scores>fba>fraction = cobra_results.fluxes('rxn_X')
    # # scores>fba>fraction = cobra_results.objective_value
    # results[sim_type] = cobra_results

    # # Write results for merged model
    # write_results_to_rpsbml(
    #     rpsbml=rpsbml,
    #     objective_id=objective_id,
    #     sim_type=sim_type,
    #     cobra_results=cobra_results,
    #     pathway_id=pathway_id,
    #     logger=logger
    # )

    # rpsbml.write_to_file('joan.xml')
    ##### print the biomass results ######
    logger.debug('Biomass: '+str(cobra_results.fluxes.get(biomass_objective_id)))
    logger.debug(' Target: '+str(cobra_results.fluxes.get(objective_id)))

    # reset the bounds to the original values for the target
    rpsbml.setReactionConstraints(
        biomass_rxn_id,
        old_upper_bound,
        old_lower_bound
    )

    logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))

    return cobra_results, results_biomass, objective_id, minimal_medium


def runCobra(
    sim_type: str,
    rpsbml: rpSBML,
    objective_id: str,
    hidden_species: List[str] = [],
    fraction_coeff: float = 0.95,
    medium: pd.DataFrame=None,
    logger: Logger = getLogger(__name__)
) -> Tuple[cobra_solution, pd.DataFrame]:
    """Run Cobra to optimize model.

    :param sim_type: The type of simulation to use. Available simulation types include: fraction, fba, rpfba
    :param rpsbml: The model to analyse.
    :param objective_id: Overwrite the auto-generated id of the results (Default: None)
    :param hidden_species: List of species to mask (Optional).
    :param fraction_coeff: The fraction of the optimum. Used in pfba simulation (Default: 0.95).
    :param medium: A DataFrame describing medium composition (Optional).
    :param logger: A logger (Optional).

    :type sim_type: str
    :type rpsbml: rpSBML
    :type objective_id: str
    :type hidden_species: List[str]
    :type fraction_coeff: float
    :type medium: pd.DataFrame
    :type logger: Logger

    :return: Results of the simulation.
    :rtype: cobra.Solution
    """

    cobraModel = build_cobra_model(
        rpsbml=rpsbml,
        objective_id=objective_id,
        hidden_species=hidden_species,
        logger=logger
    )
    if not cobraModel:
        return None

    if is_df_medium_defined(medium):
        cobraModel.medium = df_to_medium(medium)
        logger.debug('Medium modify')
        logger.debug(cobraModel.medium)

    cobra_results = None
    if sim_type.lower() == 'pfba':
        cobra_results = pfba(cobraModel, fraction_coeff)
    else:
        cobra_results = cobraModel.optimize(
            objective_sense='maximize',
            raise_error=True
        )

    logger.debug(cobra_results)

    # Minimal medium.
    logger.debug('Minimal medium')
    minimal_medium = build_minimal_medium(
        model=cobraModel,
        solution=cobra_results
    )
    logger.debug(minimal_medium)

    return cobra_results, minimal_medium

def build_cobra_model(
    rpsbml: rpSBML,
    objective_id: str,
    hidden_species: List[str] = [],
    logger: Logger = getLogger(__name__)
) -> cobra_model:
    """Convert the rpSBML object to cobra object

    :return: Success or failure of the function
    :rtype: bool
    """

    rpsbml.logger.info('Creating Cobra object from rpSBML...')

    rpsbml.activateObjective(
        objective_id = objective_id,
        plugin = 'fbc'
    )

    # To handle file removing (Windows)
    error = False
    with NamedTemporaryFile(delete=False) as temp_f:
        rpsbml.write_to_file(temp_f.name)
        temp_f.close()
        try:
            cobraModel = cobra_io.read_sbml_model(temp_f.name, use_fbc_package=True)
        except CobraSBMLError:
            logger.error('Something went wrong reading the SBML model')
            (model, errors) = validate_sbml_model(temp_f.name)
            logger.error(str(json_dumps(errors, indent=4)))
            error = True

    # To handle file removing (Windows)
    remove(temp_f.name)
    if error:
        return None

    # Hide to Cobra species that are isolated
    cobraModel.remove_metabolites(
        [
            cobraModel.metabolites.get_by_id(
                to_cobra(met)
            ) for met in hidden_species
        ]
    )

    logger.debug(cobraModel)

    return cobraModel


def write_results_to_rpsbml(
    rpsbml: rpSBML,
    objective_id: str,
    sim_type: str,
    cobra_results: cobra_solution,
    pathway_id: str = 'rp_pathway',
    logger: Logger = getLogger(__name__)
) -> None:
    """Method to hardcode into BRSynth annotations the results of a COBRA analysis

    :param objective_id: The id of the objective to optimise
    :param cobra_results: The cobrapy results object
    :param pathway_id: The id of the heterologous pathway group (Default: rp_pathway)

    :type cobra_results: cobra.ModelSummary
    :type objective_id: str
    :type pathway_id: str

    :return: None
    :rtype: None
    """
    rpsbml.logger.debug('----- Setting the results for '+str(objective_id)+ ' -----')
    rpsbml.logger.debug('rpsbml: ' + str(rpsbml))
    rpsbml.logger.debug('sim_type: ' + str(sim_type))
    rpsbml.logger.debug('objective_id: ' + str(objective_id))
    rpsbml.logger.debug('cobra_results: ' + str(cobra_results))
    rpsbml.logger.debug('cobra_results.objective_value: ' + str(cobra_results.objective_value))
    rpsbml.logger.debug('pathway_id: ' + str(pathway_id))

    write_objective_to_pathway(
        cobra_results.objective_value,
        rpsbml,
        pathway_id,
        sim_type,
        logger
    )

    # create_ignored_species_group(
    #     rpsbml,
    #     missing_species,
    #     logger
    # )

    write_fluxes_to_objectives(
        cobra_results.fluxes,
        cobra_results.objective_value,
        rpsbml,
        objective_id,
        logger
    )

    write_fluxes_to_reactions(
        cobra_results.fluxes,
        rpsbml,
        pathway_id,
        sim_type,
        logger
    )


def write_fluxes_to_objectives(
    fluxes: np_series,
    objective_value: float,
    rpsbml: rpSBML,
    objective_id: str,
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Write Cobra results to the reactions of pathway with id pathway_id in rpsbml object.

    Parameters
    ----------
    fluxes: np_series
        Fluxes.
    objective_value: float
        Value of the objective.
    rpsbml: rpSBML
        rpSBML object of which reactions will be updated with results
    objective_id: str
        The id of the objective to optimise
    """

    rpsbml.logger.debug(
        'Set the objective ' + str(objective_id) \
      + ' a flux_value of ' + str(objective_value)
    )

    # get the objective
    obj = rpsbml.getObjective(
        objective_id,
        objective_value
    )
    rpsbml.updateBRSynth(
        obj,
        'flux_value',
        str(objective_value)
    )

    Bigg_Reaction_Prefix = 'R_'

    for flux_obj in obj.getListOfFluxObjectives():

        rxn_id = flux_obj.getReaction()
        if rxn_id.startswith(Bigg_Reaction_Prefix):
            rxn_id = rxn_id[2:]
        flux = fluxes.get(rxn_id)

        # sometimes flux cannot be returned
        if flux is None:
            rpsbml.logger.warning(
                f'Cannot retrieve reaction ID {str(rxn_id)} flux from cobrapy... setting to 0.0'
            )
            flux = 0.0

        rpsbml.logger.debug(
            'Set the reaction ' + str(rxn_id) \
            + ' a flux_value of ' + str(flux)
        )
        rpsbml.updateBRSynth(
            flux_obj,
            'flux_value',
            str(flux)
        )


def write_fluxes_to_reactions(
    fluxes: np_series,
    rpsbml: rpSBML,
    pathway_id: str,
    sim_type: str,
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Write Cobra results to the reactions of pathway with id pathway_id in rpsbml object.

    Parameters
    ----------
    fluxes: np_series
        Fluxes.
    rpsbml: rpSBML
        rpSBML object of which reactions will be updated with results
    pathway_id: str
        The id of the pathway within reactions will be updated
    objective_id: str
        The id of the objective to optimise
    """

    rp_pathway = rpsbml.getGroup(pathway_id)

    for member in rp_pathway.getListOfMembers():

        rxn = rpsbml.getModel().getReaction(member.getIdRef())

        if rxn is None:
            logger.error(
                'Cannot retreive the following reaction: ' \
              + str(member.getIdRef())
            )
            #return False
            continue
        
        flux = fluxes.get(rxn.getId())

        logger.debug(
            'Set the reaction ' + str(member.getIdRef()) \
          + ' a ' + str('fba_' + str(sim_type)) \
          + ' of ' + str(flux))

        rpsbml.updateBRSynth(
            rxn,
            'fba_'+str(sim_type),
            str(flux)
        )


def write_objective_to_pathway(
    objective_value: float,
    rpsbml: rpSBML,
    pathway_id: str,
    sim_type: str,
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Write Cobra results to the pathway with id pathway_id in rpsbml object.

    Parameters
    ----------
    objective_value: float
        Value of the objective.
    rpsbml: rpSBML
        rpSBML object of which reactions will be updated with results
    pathway_id: str
        The id of the pathway within reactions will be updated
    objective_id: str
        The id of the objective to optimise
    """

    logger.debug(
        'Set ' + str(pathway_id) + ' with ' \
      + str('fba_'+str(sim_type)) + ' to ' \
      + str(objective_value)
    )

    rpsbml.updateBRSynth(
        rpsbml.getGroup(pathway_id),
        'fba_'+str(sim_type),
        str(objective_value)
    )


def create_ignored_species_group(
    rpsbml: rpSBML,
    hidden_species: List[str],
    logger: Logger = getLogger(__name__)
) -> None:
    """
    Write ignored species during FBA to the pathway with id pathway_id in rpsbml object.

    Parameters
    ----------
    rpsbml: rpSBML
        rpSBML object of which reactions will be updated with results
    pathway_id: str
        The id of the pathway within reactions will be updated
    """

    group_id = 'rpfba_ignored_species'

    rpsbml.createGroup(
        id=group_id,
        brs_annot=False
    )
    for spe in hidden_species:
        rpsbml.addMember(
            group_id=group_id,
            idRef=spe
        )
    
    logger.debug(
        'Create ' + str(group_id) + ' group with ' \
      + str(hidden_species)
    )

    # rpsbml.updateBRSynth(
    #     sbase_obj=rpsbml.getGroup(group_id),
    #     annot_header='ignored_species',
    #     value=str(missing_species),
    #     isList=True
    # )

########################################################################
############################### FBA pathway ranking ####################
########################################################################

# 1) Number of interventions
# need to calculate the number of steps that are not native to know the number of interventions

# 2) Maximal growth rate

# 3) Minimum product yeild at maximal growth rate

# 4) Minimum product yeild

# 5) Anaerobic condition

# 6) Number of potentially disruptive products

    # Toxicity?

# 7) Number of accessible metabolites (avoid intermediate accumulation)

# 8) Thermodynamics (MDF)

# 9) The overlap of the same changes --> might not be applicable in our case

# 10) Reduced model

# 11) ECM


# def runMultiObjective(rpsbml,
#                       reactions,
#                       coefficients,
#                       is_max=True,
#                       pathway_id='rp_pathway',
#                       objective_id=None):
#     """Run FBA using multiple objectives
#
#     :param reactions: The ids of the reactions involved in the objective
#     :param coefficients: The coefficients associated with the reactions id
#     :param is_max: Maximise or minimise the objective (Default: True)
#     :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
#     :param objective_id: Overwrite the default id (Default: None)
#
#     :type reactions: list
#     :type coefficients: list
#     :type is_max: bool
#     :type pathway_id: str
#     :type objective_id: str
#
#     :return: Success or failure of the function
#     :rtype: bool
#     """
#     fbc_plugin = rpsbml.getModel().getPlugin('fbc')
#     rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
#     objective_id = rpsbml.findCreateObjective(reactions, coefficients, is_max)
#     rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
#                         'Setting active objective '+str(objective_id))
#     cobraModel = rpsbml.convertToCobra()
#     if not cobraModel:
#         return False
#     cobra_results = cobraModel.optimize()
#     rpsbml.writeFBAResults(objective_id, cobra_results, pathway_id)
#     return rpsbml

