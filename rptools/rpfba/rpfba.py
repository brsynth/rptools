import pandas as pd
from logging import Logger, getLogger
from os import remove
from argparse import Namespace as arg_nspace
from pandas.core.series import Series as np_series
from typing import List, Dict, Tuple
from tempfile import NamedTemporaryFile
from json import dumps as json_dumps
from cobra.flux_analysis import pfba
from cobra import io as cobra_io
from cobra.io.sbml import validate_sbml_model, CobraSBMLError
from cobra.core.model import Model as cobra_model
from cobra.core.solution import Solution as cobra_solution

from rptools.rplibs import (
    rpSBML,
    rpPathway
)
from .cobra_format import cobraize, to_cobra
from .Args import DEFAULT_ARGS as DEFAULT_RPFBA_ARGS

# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes

class ModelError(Exception):
    pass


def preprocess(
    args: arg_nspace,
    logger: Logger = getLogger(__name__),
):
    pathway = rpPathway(args.pathway_file, logger=logger)
    pathway.setup_pathway_fba()
    model = rpSBML(inFile=args.model_file, logger=logger)

    try:
        ids = check_ids(
            pathway=pathway,
            model=model,
            objective_rxn_id=args.objective_rxn_id,
            biomass_rxn_id=args.biomass_rxn_id,
            compartment_id=args.compartment_id,
            logger=logger,
        )
    except ModelError as e:
        logger.error(e)
        return 1

    # MERGE
    (
        merged_model,
        reactions_in_both,
        missing_species,
        compartment_id
    ) = rpSBML.merge(
        pathway=pathway.get_rpsbml(),
        model=model,
        compartment_id=ids['comp_id'],
        logger=logger
    )
    logger.debug(f"model: {model}")
    logger.debug(f"reactions_in_both: {reactions_in_both}")
    logger.debug(f"missing_species: {missing_species}")

    # CHECKING
    # Detect orphan species among missing ones in the model,
    # i.e. that are only consumed or produced
    if not args.with_orphan_species:
        merged_model.search_isolated_species(missing_species)

    if args.merge != "":
        logger.info(f"Write merged rpSBML file to {args.merge}")
        merged_model.write_to_file(args.merge)

    return merged_model, pathway, ids


def check_ids(
    pathway: rpPathway,
    model: rpSBML,
    objective_rxn_id: str,
    biomass_rxn_id: str,
    compartment_id: str,
    logger: Logger = getLogger(__name__)
) -> Dict:
    '''Check the IDs of the pathway and the model.
    If the IDs are not found, then raise ModelError exception.
    
    :param pathway: The pathway rpSBML object
    :param model: The model rpSBML object
    :param objective_rxn_id: The objective reaction ID
    :param biomass_rxn_id: The biomass reaction ID
    :param compartment_id: The SBML compartment ID
    :param logger: The logger object
    
    :type pathway: rpPathway
    :type model: rpSBML
    :type objective_rxn_id: str
    :type biomass_rxn_id: str
    :type compartment_id: str
    :type logger: Logger
    
    :return: The objective reaction ID, the model compartment ID, and the biomass reaction ID
    :rtype: Dict
    '''

    # PATHWAY
    # Check objective reaction ID
    objective_rxn_id = pathway.get_rpsbml().check_SBML_rxnid(objective_rxn_id)
    if objective_rxn_id is None:
        raise ModelError(f"No objective reaction ID found in the pathway {pathway.get_id()}.")

    # MODEL
    # Check compartment ID
    compartment_id = model.check_SBML_compartment(compartment_id)
    if compartment_id is None:
        raise ModelError(f"No compartment ID found in the model {model.get_id()}.")
    # Check biomass reaction ID
    biomass_rxn_id = model.check_SBML_rxnid(biomass_rxn_id)
    if biomass_rxn_id is None:
        raise ModelError(f"No biomass reaction ID found in the model {model.get_id()}.")

    return {
    'obj_rxn_id': objective_rxn_id,
    'comp_id': compartment_id,
    'biomass_rxn_id': biomass_rxn_id
    }


def runFBA_fromFile(
    model_file: str,
    compartment_id: str,
    objective_rxn_id: str = DEFAULT_RPFBA_ARGS["objective_rxn_id"],
    biomass_rxn_id: str = DEFAULT_RPFBA_ARGS["biomass_rxn_id"],
    sim_type: str = DEFAULT_RPFBA_ARGS["sim"],
    fraction_coeff: float = DEFAULT_RPFBA_ARGS["fraction_coeff"],
    hidden_species: List[str] = [],
    logger: Logger = getLogger(__name__),
) -> Dict:
    """Single rpSBML simulation

    :param model_file: Path to the model file (SBML)
    :param compartment_id: The compartment ID
    :param objective_rxn_id: The objective reaction ID (Default: rxn_target)
    :param biomass_rxn_id: The biomass reaction ID (Default: biomass)
    :param sim_type: The simulation type (Default: fraction)
    :param fraction_coeff: The fraction coefficient (Default: 0.75)
    :param hidden_species: List of hidden species (Default: [])
    :param logger: The logger object

    :type model_file: str
    :type compartment_id: str
    :type objective_rxn_id: str
    :type biomass_rxn_id: str
    :type sim_type: str
    :type fraction_coeff: float
    :type hidden_species: List[str]
    :type logger: Logger

    :return: The results of the simulation
    :rtype: Dict
    """

    logger.debug("           model_file: " + str(model_file))
    logger.debug("             sim_type: " + str(sim_type))
    logger.debug("     objective_rxn_id: " + str(objective_rxn_id))
    logger.debug("       biomass_rxn_id: " + str(biomass_rxn_id))
    logger.debug("       fraction_coeff: " + str(fraction_coeff))
    logger.debug("       compartment_id: " + str(compartment_id))
    logger.debug("       hidden_species: " + str(hidden_species))

    return runFBA(
        model=rpSBML(inFile=model_file, logger=logger),
        compartment_id=compartment_id,
        objective_rxn_id=objective_rxn_id,
        biomass_rxn_id=biomass_rxn_id,
        sim_type=sim_type,
        fraction_coeff=fraction_coeff,
        hidden_species=hidden_species,
        logger=logger
    )


def runFBA(
    model: rpSBML,
    compartment_id: str,
    objective_rxn_id: str = DEFAULT_RPFBA_ARGS["objective_rxn_id"],
    biomass_rxn_id: str = DEFAULT_RPFBA_ARGS["biomass_rxn_id"],
    sim_type: str = DEFAULT_RPFBA_ARGS["sim"],
    fraction_coeff: float = DEFAULT_RPFBA_ARGS["fraction_coeff"],
    logger: Logger = getLogger(__name__),
) -> Dict:
    """Single rpSBML simulation

    :param model_file: Path to the model file (SBML)
    :param compartment_id: The compartment ID
    :param objective_rxn_id: The objective reaction ID (Default: rxn_target)
    :param biomass_rxn_id: The biomass reaction ID (Default: biomass)
    :param sim_type: The simulation type (Default: fraction)
    :param fraction_coeff: The fraction coefficient (Default: 0.75)
    :param hidden_species: List of hidden species (Default: [])
    :param logger: The logger object

    :type model_file: str
    :type compartment_id: str
    :type objective_rxn_id: str
    :type biomass_rxn_id: str
    :type sim_type: str
    :type fraction_coeff: float
    :type hidden_species: List[str]
    :type logger: Logger

    :return: The results of the simulation
    :rtype: Dict
    """

    logger.debug("                model: " + str(model))
    logger.debug("             sim_type: " + str(sim_type))
    logger.debug("     objective_rxn_id: " + str(objective_rxn_id))
    logger.debug("       biomass_rxn_id: " + str(biomass_rxn_id))
    logger.debug("       fraction_coeff: " + str(fraction_coeff))
    logger.debug("       compartment_id: " + str(compartment_id))

    # # Load the model
    # model_c = rpSBML.cobraize(model)

    # # NOTE: reactions is organised with key being the rpsbml reaction and value being the rpsbml_gem value`
    # # BUG: when merging the rxn_sink (very rare cases) can be recognised if another reaction contains the same species as a reactant
    # ## under such as scenario the algorithm will consider that they are the same -- TODO: overwrite it
    # if reactions_in_both is not None and pathway_obj_rxn_id in reactions_in_both:
    #     logger.warning(
    #         "The target_reaction ("
    #         + str(pathway_obj_rxn_id)
    #         + ") "
    #         + "has been detected in model "
    #         + str(model.getName())
    #         + ", "
    #         + "ignoring this model..."
    #     )
    #     return 2

    ######## FBA ########
    results = {}
    if sim_type.lower() in ["fba", "pfba"]:
        objective_id = model.find_or_create_objective(
            rxn_id=objective_rxn_id,
            obj_id=f"brs_obj_{objective_rxn_id}",
        )
        cobra_results = runCobra(
            sim_type=sim_type,
            rpsbml=model,
            objective_id=objective_id,
            fraction_coeff=fraction_coeff,
            logger=logger,
        )
    else:
        (cobra_results, results_biomass, objective_id) = rp_fraction(
            rpsbml=model,
            objective_rxn_id=objective_rxn_id,
            biomass_rxn_id=biomass_rxn_id,
            fraction_coeff=fraction_coeff,
            logger=logger,
        )

        results["biomass"] = results_biomass

    # print(cobra_results.objective_value)
    # print(cobra_results.status)
    # print(cobra_results.fluxes)
    # print(cobra_results.shadow_prices)

    results[sim_type] = cobra_results

    # Write results for merged model
    write_results_to_rpsbml(
        rpsbml=model,
        objective_id=objective_id,
        cobra_results=cobra_results,
        sim_type=sim_type,
        logger=logger,
    )

    return results

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


# def build_hidden_species(
#     rpsbml: rpSBML,
#     missing_species: List[str],
#     compartment_id: bool,
#     logger: Logger = getLogger(__name__),
# ) -> List[str]:
#     rpsbml.search_isolated_species(missing_species)
#     return rpsbml.get_isolated_species()
#     return [
#         cobraize(spe_id, compartment_id) for spe_id in rpsbml.get_isolated_species()
#     ]


def build_results(
    results: Dict,
    pathway: rpPathway,
    compartment_id: str,
    hidden_species: List[str],
    logger: Logger = getLogger(__name__),
) -> Dict:
    _results = {
        "species": {},
        "reactions": {},
        "pathway": {},
        "ignored_species": hidden_species,
    }

    # SPECIES
    for spe_id in pathway.get_species_ids():
        _results["species"][spe_id] = {}
        for sim_type, cobra_r in results.items():
            value = cobra_r.shadow_prices.get(
                to_cobra(cobraize(spe_id, compartment_id))
            )
            _results["species"][spe_id][sim_type + "_shadow_price"] = {
                "value": value,
                # 'units': 'milimole / gDW / hour',
            }
    # REACTIONS
    for rxn_id in pathway.get_reactions_ids():
        _results["reactions"][rxn_id] = {}
        for sim_type, cobra_r in results.items():
            _results["reactions"][rxn_id][sim_type] = {
                "value": cobra_r.fluxes[rxn_id],
                "units": "milimole / gDW / hour",
            }
            if sim_type == "biomass":
                _results["reactions"][rxn_id][sim_type]["units"] = "gDW / gDW / hour"
    # PATHWAY
    _results["pathway"] = {}
    for sim_type, cobra_r in results.items():
        _results["pathway"][sim_type] = {
            "value": cobra_r.objective_value,
            "units": "milimole / gDW / hour",
        }
        if sim_type == "biomass":
            _results["pathway"][sim_type]["units"] = "gDW / gDW / hour"

    return _results


def write_results_to_pathway(
    pathway: rpPathway, results: Dict, logger: Logger = getLogger(__name__)
) -> None:
    # Write species results
    for spe_id, score in results["species"].items():
        for k, v in score.items():
            pathway.get_specie(spe_id).add_fba_info(key=k, value=v)
    # Write reactions results
    for rxn_id, score in results["reactions"].items():
        for k, v in score.items():
            pathway.get_reaction(rxn_id).add_fba_info(key=k, value=v)
    # Write pathway result
    for k, v in results["pathway"].items():
        pathway.add_fba_info(key=k, value=v)
    # Write ignored species
    pathway.add_species_group("fba_ignored", results["ignored_species"])


# def complete_heterologous_pathway(
#     rpsbml: rpSBML,
#     hidden_species: List[str],
#     rpsbml_merged: rpSBML,
#     species_group_id: str,
#     sink_species_group_id: str,
#     pathway_id: str,
#     reactions_in_both: Dict,
#     logger: Logger = getLogger(__name__),
# ) -> None:
#     # Save the central species
#     # groups = rpsbml.getPlugin('groups')
#     central = rpsbml.getGroup(species_group_id)
#     sink_group = rpsbml.getGroup(sink_species_group_id)
#     rp_group = rpsbml.getGroup(pathway_id)
#     cent_spe = [str(i.getIdRef()) for i in central.getListOfMembers()]
#     sink_spe = [str(i.getIdRef()) for i in sink_group.getListOfMembers()]
#     rp_reac = [str(i.getIdRef()) for i in rp_group.getListOfMembers()]
#     logger.debug("old central species: " + str(cent_spe))
#     logger.debug("old sink species:    " + str(sink_spe))
#     logger.debug("old rp reactions:    " + str(rp_reac))

#     rev_reactions = {v: k for k, v in reactions_in_both.items()}
#     logger.debug("reactions_in_both: " + str(reactions_in_both))
#     logger.debug("rev_reactions:     " + str(rev_reactions))
#     logger.info("Building model with heterologous pathway only")
#     groups = rpsbml_merged.getPlugin("groups")
#     rp_pathway = rpsbml_merged.getGroup(pathway_id)
#     logger.debug("---- Reactions ----")
#     for member in rp_pathway.getListOfMembers():
#         #### reaction annotation
#         logger.debug(member.getIdRef())
#         reacFBA = rpsbml_merged.getModel().getReaction(member.getIdRef())
#         logger.debug(reacFBA)
#         try:
#             # reacIN = rpsbml.model.getReaction(reactions_convert[member.getIdRef()])
#             reacIN = rpsbml.getModel().getReaction(rev_reactions[member.getIdRef()])
#         except KeyError:
#             reacIN = rpsbml.getModel().getReaction(member.getIdRef())
#         logger.debug(reacIN)
#         logger.debug(reacFBA.getAnnotation())
#         reacIN.setAnnotation(reacFBA.getAnnotation())
#         #### species TODO: only for shadow price
#     #### add groups ####
#     source_groups = rpsbml_merged.getPlugin("groups")
#     target_groups = rpsbml.getPlugin("groups")
#     target_groupsID = [i.getId() for i in rpsbml.getListOfGroups()]
#     for source_group in rpsbml_merged.getListOfGroups():
#         logger.debug("Replacing group id: " + str(source_group.getId()))
#         if source_group.getId() == species_group_id:
#             target_group = rpsbml.getGroup(source_group.getId())
#             # TODO: #### replace the new potentially incorect central species with the normal ones #####
#             # delete all the previous members
#             logger.debug("Removing rp_core_species")
#             for i in range(target_group.getNumMembers()):
#                 logger.debug(
#                     "Deleting group member: "
#                     + str(target_group.getMember(0).getIdRef())
#                 )
#                 target_group.removeMember(0)
#             # add the new ones
#             for cs in cent_spe:
#                 logger.debug("Creating new member: " + str(cs))
#                 newM = target_group.createMember()
#                 newM.setIdRef(cs)
#         elif source_group.getId() == sink_species_group_id:
#             target_group = rpsbml.getGroup(source_group.getId())
#             logger.debug("Removing sink species")
#             for i in range(target_group.getNumMembers()):
#                 logger.debug(
#                     "Deleting group member: "
#                     + str(target_group.getMember(0).getIdRef())
#                 )
#                 target_group.removeMember(0)
#             # add the new ones
#             for cs in sink_spe:
#                 logger.debug("Creating new member: " + str(cs))
#                 newM = target_group.createMember()
#                 newM.setIdRef(cs)
#         elif source_group.getId() in target_groupsID:
#             target_group = rpsbml.getGroup(source_group.getId())
#             target_group.setAnnotation(source_group.getAnnotation())
#         """
#         elif source_group.getId()==pathway_id:
#             target_group = target_groups.getGroup(source_group.getId())
#             logger.debug('Removing rp ractions')
#             for i in range(target_group.getNumMembers()):
#                 logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
#                 target_group.removeMember(0)
#             #add the new ones
#             for cs in rp_reac:
#                 logger.debug('Creating new member: '+str(cs))
#                 newM = target_group.createMember()
#                 newM.setIdRef(cs)
#         """
#     # add ignored_species group
#     rpsbml.set_isolated_species(rpsbml_merged.get_isolated_species())
#     create_ignored_species_group(rpsbml, hidden_species, logger)
#     #### add objectives ####
#     source_fbc = rpsbml_merged.getPlugin("fbc")
#     target_fbc = rpsbml.getPlugin("fbc")
#     target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
#     for source_obj in source_fbc.getListOfObjectives():
#         source_obj_id = source_obj.getId()
#         if source_obj.getId() in target_objID:
#             target_obj = target_fbc.getObjective(source_obj.getId())
#             target_obj.setAnnotation(source_obj.getAnnotation())
#             for target_fluxObj in target_obj.getListOfFluxObjectives():
#                 for source_fluxObj in source_obj.getListOfFluxObjectives():
#                     if target_fluxObj.getReaction() == source_fluxObj.getReaction():
#                         target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
#         else:
#             target_fbc.addObjective(source_obj)
#     # rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
#     target_fbc.setActiveObjectiveId(
#         source_obj_id
#     )  # tmp random assignement of objective


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
    fraction_coeff: float = DEFAULT_RPFBA_ARGS["fraction_coeff"],
    logger: Logger = getLogger(__name__),
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

    logger.debug("rpsbml:       " + str(rpsbml))
    logger.debug("objective_rxn_id:   " + objective_rxn_id)
    logger.debug("biomass_rxn_id: " + str(biomass_rxn_id))
    logger.debug("fraction_coeff:  " + str(fraction_coeff))

    def get_annot_objective(rpsbml: rpSBML, objective_id: str) -> str:
        fbc_plugin = rpsbml.getPlugin("fbc")
        objective = fbc_plugin.getObjective(objective_id)
        rpsbml.checklibSBML(objective, "Getting objective " + str(objective_id))

        fbc_obj_annot = objective.getAnnotation()
        # TODO: if this is None need to set it up
        if fbc_obj_annot is None:
            return None

        logger.debug("Already calculated flux for " + str(objective_id))

        return fbc_obj_annot

    # retreive the biomass objective and flux results and set as maxima
    biomass_objective_id = rpsbml.find_or_create_objective(
        rxn_id=biomass_rxn_id, obj_id=f"brs_obj_{biomass_rxn_id}"
    )

    # fbc_plugin = rpsbml.getPlugin('fbc')
    # TODO: use the rpSBML BRSynth annotation parser
    # try:
    fbc_obj_annot = get_annot_objective(rpsbml, biomass_objective_id)

    # except (AttributeError, ValueError) as e:
    if fbc_obj_annot is None:
        # logger.debug(e)
        logger.debug("Performing FBA to calculate the source reaction")

        ### FBA ###
        # logger.info('Running the FBA (fraction of reaction)...')
        # rpsbml.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
        logger.info("Processing FBA (biomass)...")
        cobra_results = runCobra(
            sim_type="biomass",
            rpsbml=rpsbml,
            objective_id=biomass_objective_id,
            logger=logger,
        )

        # rxn>scores>fba>biomass = cobra_results.fluxes('rxn_X')
        # scores>fba>biomass = cobra_results.objective_value
        if cobra_results is None:
            return None, None, biomass_objective_id

        results_biomass = cobra_results

        write_results_to_rpsbml(
            rpsbml=rpsbml,
            objective_id=biomass_objective_id,
            sim_type="biomass",
            cobra_results=cobra_results,
            logger=logger,
        )

        fbc_obj_annot = get_annot_objective(rpsbml, biomass_objective_id)
        if fbc_obj_annot is None:
            logger.error("No annotation available for: " + str(biomass_objective_id))

    flux = float(
        fbc_obj_annot.getChild("RDF")
        .getChild("BRSynth")
        .getChild("brsynth")
        .getChild(0)
        .getAttrValue("value")
    )

    objective_id = rpsbml.find_or_create_objective(
        rxn_id=objective_rxn_id,
        obj_id=f"brs_obj_{objective_rxn_id}",
    )
    logger.debug(f"objective_id: {objective_id}")

    logger.debug(f"Optimising the objective: {biomass_rxn_id}")
    logger.debug(f"     Setting upper bound: {flux*fraction_coeff}")
    logger.debug(f"     Setting lower bound: {flux*fraction_coeff}")

    old_upper_bound, old_lower_bound = rpsbml.getReactionConstraints(biomass_rxn_id)
    rpsbml.setReactionConstraints(
        biomass_rxn_id, flux * fraction_coeff, flux * fraction_coeff
    )

    logger.info("Processing FBA (fraction)...")
    sim_type = "fraction"
    cobra_results = runCobra(
        sim_type=sim_type,
        rpsbml=rpsbml,
        objective_id=objective_id,
        fraction_coeff=fraction_coeff,
        logger=logger,
    )
    if cobra_results is None:
        return results_biomass, None, objective_id

    ##### print the biomass results ######
    logger.debug("Biomass: " + str(cobra_results.fluxes.get(biomass_objective_id)))
    logger.debug(" Target: " + str(cobra_results.fluxes.get(objective_id)))

    # reset the bounds to the original values for the target
    rpsbml.setReactionConstraints(biomass_rxn_id, old_upper_bound, old_lower_bound)

    logger.debug(
        "The objective "
        + str(objective_id)
        + " results "
        + str(cobra_results.objective_value)
    )

    return cobra_results, results_biomass, objective_id


def runCobra(
    sim_type: str,
    rpsbml: rpSBML,
    objective_id: str,
    fraction_coeff: float = 0.95,
    logger: Logger = getLogger(__name__),
) -> Tuple[cobra_solution, pd.DataFrame]:
    """Run Cobra to optimize model.

    :param sim_type: The type of simulation to use. Available simulation types include: fraction, fba, rpfba
    :param rpsbml: The model to analyse.
    :param objective_id: Overwrite the auto-generated id of the results (Default: None)
    :param hidden_species: List of species to mask (Optional).
    :param fraction_coeff: The fraction of the optimum. Used in pfba simulation (Default: 0.95).
    :param logger: A logger (Optional).

    :type sim_type: str
    :type rpsbml: rpSBML
    :type objective_id: str
    :type hidden_species: List[str]
    :type fraction_coeff: float
    :type logger: Logger

    :return: Results of the simulation.
    :rtype: cobra.Solution
    """

    cobraModel = build_cobra_model(
        rpsbml=rpsbml,
        objective_id=objective_id,
        logger=logger,
    )
    if not cobraModel:
        return None

    cobra_results = None
    # cobraModel.objective = {
    #     cobraModel.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M'): 1,
    #     cobraModel.reactions.get_by_id('PYRt2'): 2
    # }
    if sim_type.lower() == "pfba":
        cobra_results = pfba(cobraModel, fraction_coeff)
    else:
        cobra_results = cobraModel.optimize(
            objective_sense="maximize", raise_error=True
        )

    logger.debug(cobra_results)

    return cobra_results


def build_cobra_model(
    rpsbml: rpSBML,
    objective_id: str,
    logger: Logger = getLogger(__name__),
) -> cobra_model:
    """Convert the rpSBML object to cobra object

    :return: Success or failure of the function
    :rtype: bool
    """

    rpsbml.logger.info("Creating Cobra object from rpSBML...")
    logger.debug(f"objective_id: {objective_id}")

    rpsbml.activateObjective(objective_id=objective_id, plugin="fbc")

    # To handle file removing (Windows)
    error = False
    with NamedTemporaryFile(delete=False) as temp_f:
        rpsbml.write_to_file(temp_f.name)
        temp_f.close()
        try:
            cobraModel = cobra_io.read_sbml_model(temp_f.name, use_fbc_package=True)
        except CobraSBMLError:
            logger.error("Something went wrong reading the SBML model")
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
            cobraModel.metabolites.get_by_id(to_cobra(met))
            for met in rpsbml.get_isolated_species()
        ]
    )

    logger.debug(cobraModel)

    return cobraModel


def write_results_to_rpsbml(
    rpsbml: rpSBML,
    objective_id: str,
    sim_type: str,
    cobra_results: cobra_solution,
    pathway_id: str = "rp_pathway",
    logger: Logger = getLogger(__name__),
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
    rpsbml.logger.debug("----- Setting the results for " + str(objective_id) + " -----")
    rpsbml.logger.debug("rpsbml: " + str(rpsbml))
    rpsbml.logger.debug("sim_type: " + str(sim_type))
    rpsbml.logger.debug("objective_id: " + str(objective_id))
    rpsbml.logger.debug("cobra_results: " + str(cobra_results))
    rpsbml.logger.debug(
        "cobra_results.objective_value: " + str(cobra_results.objective_value)
    )
    rpsbml.logger.debug("pathway_id: " + str(pathway_id))

    write_objective_to_pathway(
        cobra_results.objective_value, rpsbml, pathway_id, sim_type, logger
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
        logger,
    )

    write_fluxes_to_reactions(
        cobra_results.fluxes, rpsbml, pathway_id, sim_type, logger
    )


def write_fluxes_to_objectives(
    fluxes: np_series,
    objective_value: float,
    rpsbml: rpSBML,
    objective_id: str,
    logger: Logger = getLogger(__name__),
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
        "Set the objective "
        + str(objective_id)
        + " a flux_value of "
        + str(objective_value)
    )

    # get the objective
    obj = rpsbml.getObjective(objective_id, objective_value)
    rpsbml.updateBRSynth(obj, "flux_value", str(objective_value))

    Bigg_Reaction_Prefix = "R_"

    for flux_obj in obj.getListOfFluxObjectives():
        rxn_id = flux_obj.getReaction()
        if rxn_id.startswith(Bigg_Reaction_Prefix):
            rxn_id = rxn_id[2:]
        flux = fluxes.get(rxn_id)

        # sometimes flux cannot be returned
        if flux is None:
            rpsbml.logger.warning(
                f"Cannot retrieve reaction ID {str(rxn_id)} flux from cobrapy... setting to 0.0"
            )
            flux = 0.0

        rpsbml.logger.debug(
            "Set the reaction " + str(rxn_id) + " a flux_value of " + str(flux)
        )
        rpsbml.updateBRSynth(flux_obj, "flux_value", str(flux))


def write_fluxes_to_reactions(
    fluxes: np_series,
    rpsbml: rpSBML,
    pathway_id: str,
    sim_type: str,
    logger: Logger = getLogger(__name__),
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
                "Cannot retreive the following reaction: " + str(member.getIdRef())
            )
            # return False
            continue

        flux = fluxes.get(rxn.getId())

        logger.debug(
            "Set the reaction "
            + str(member.getIdRef())
            + " a "
            + str("fba_" + str(sim_type))
            + " of "
            + str(flux)
        )

        rpsbml.updateBRSynth(rxn, "fba_" + str(sim_type), str(flux))


def write_objective_to_pathway(
    objective_value: float,
    rpsbml: rpSBML,
    pathway_id: str,
    sim_type: str,
    logger: Logger = getLogger(__name__),
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
        "Set "
        + str(pathway_id)
        + " with "
        + str("fba_" + str(sim_type))
        + " to "
        + str(objective_value)
    )

    rpsbml.updateBRSynth(
        rpsbml.getGroup(pathway_id), "fba_" + str(sim_type), str(objective_value)
    )


# def create_ignored_species_group(
#     rpsbml: rpSBML, hidden_species: List[str], logger: Logger = getLogger(__name__)
# ) -> None:
#     """
#     Write ignored species during FBA to the pathway with id pathway_id in rpsbml object.

#     Parameters
#     ----------
#     rpsbml: rpSBML
#         rpSBML object of which reactions will be updated with results
#     pathway_id: str
#         The id of the pathway within reactions will be updated
#     """

#     group_id = "rpfba_ignored_species"

#     rpsbml.createGroup(id=group_id, brs_annot=False)
#     for spe in hidden_species:
#         rpsbml.addMember(group_id=group_id, idRef=spe)

#     logger.debug("Create " + str(group_id) + " group with " + str(hidden_species))

#     # rpsbml.updateBRSynth(
#     #     sbase_obj=rpsbml.getGroup(group_id),
#     #     annot_header='ignored_species',
#     #     value=str(missing_species),
#     #     isList=True
#     # )


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
