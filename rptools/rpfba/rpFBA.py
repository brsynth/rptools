from logging import (
    Logger,
    getLogger
)
from cobra.flux_analysis import pfba
from cobra               import io        as cobra_io
from cobra.io.sbml       import (
    validate_sbml_model,
    CobraSBMLError
)
from cobra.core.model    import Model     as cobra_model
from cobra.core.solution import Solution  as cobra_solution
from pandas.core.series  import Series    as np_series
from libsbml             import Objective as sbml_objective
from rptools.rplibs      import rpSBML
from typing import (
    List,
    Dict,
    Tuple
)
from tempfile import (
    NamedTemporaryFile
)
from json import dumps as json_dumps



# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes


# TODO: do not use the species_group_id and the sink_species_group_id. Loop through all the groups (and if the same) and overwrite the annotation instead
def runFBA(
              rpsbml_path: str,
            gem_sbml_path: str,
                 sim_type: str,
               src_rxn_id: str,
                src_coeff: float,
               tgt_rxn_id: str,
                tgt_coeff: float,
                   is_max: bool = True,
              frac_of_src: float = 0.75,
                    merge: bool = False,
               pathway_id: str = 'rp_pathway',
             objective_id: str = None,
           compartment_id: str = 'MNXC3',
    ignore_orphan_species: bool = True,
         species_group_id: str = 'central_species',
    sink_species_group_id: str = 'rp_sink_species',
                   logger: Logger = getLogger(__name__)
) -> rpSBML:
    """Single rpSBML simulation

    :param file_name: The name of the model
    :param rpsbml_path: Path to the rpSBML file
    :param gem_sbml: Path to the GEM file
    :param sim_type: The type of simulation to use. Available simulation types include: fraction, fba, rpfba
    :param src_rxn_id: The reaction id of the source reaction.
    :param target_reaction: The reaction id of the target reaction. Note that if fba or rpfba options are used, then these are ignored
    :param source_coefficient: The source coefficient
    :param target_coefficient: The target coefficient
    :param is_max: Maximise or minimise the objective
    :param fraction_of: The fraction of the optimum. Note that this value is ignored is fba is used
    :param tmpOutputFolder: The path to the output document
    :param merge: Output the merged model (Default: False)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the auto-generated id of the results (Default: None)
    :param compartment_id: The SBML compartment id (Default: MNXC3)
    :param fill_orphan_species: Add pseudo reactions that consume/produce single parent species. Note in development
    :param species_group_id: The id of the central species (Default: central_species)
    :param sink_species_group_id: The id of the sink species (Default: rp_sink_species)

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

    logger.debug('          rpsbml_path: ' + str(rpsbml_path))
    logger.debug('        gem_sbml_path: ' + str(gem_sbml_path))
    logger.debug('             sim_type: ' + str(sim_type))
    logger.debug('           src_rxn_id: ' + src_rxn_id)
    logger.debug('            src_coeff: ' + str(src_coeff))
    logger.debug('           tgt_rxn_id: ' + tgt_rxn_id)
    logger.debug('            tgt_coeff: ' + str(tgt_coeff))
    logger.debug('          frac_of_src: ' + str(frac_of_src))
    logger.debug('               is_max: ' + str(is_max))
    logger.debug('                merge: ' + str(merge))
    logger.debug('           pathway_id: ' + pathway_id)
    logger.debug('         objective_id: ' + str(objective_id))
    logger.debug('       compartment_id: ' + str(compartment_id))
    logger.debug('ignore_orphan_species: ' + str(ignore_orphan_species))
    logger.debug('     species_group_id: ' + str(species_group_id))
    logger.debug('sink_species_group_id: ' + str(sink_species_group_id))

    rpsbml = rpSBML(rpsbml_path, logger=logger)
    logger.debug('input_sbml: ' + str(rpsbml))

    rpsbml_gem = rpSBML(gem_sbml_path, logger=logger)
    logger.debug('rpsbml_gem: ' + str(rpsbml_gem))

    logger.info('Merging rpSBML models: ' + rpsbml.getName() + ' and ' + rpsbml_gem.getName() + '...')

    rpsbml_merged, reactions_in_both = rpSBML.mergeModels(
        source_rpsbml = rpsbml,
        target_rpsbml = rpsbml_gem,
        logger = logger
    )

    if rpsbml_merged is None:
        return None
    
    if ignore_orphan_species:
        rpsbml_merged.search_isolated_species()

    logger.debug('rpsbml_merged: ' + str(rpsbml_merged))
    logger.debug('reactions_in_both: ' + str(reactions_in_both))

    # NOTE: reactions is organised with key being the rpsbml reaction and value being the rpsbml_gem value`
    # BUG: when merging the rxn_sink (very rare cases) can be recognised if another reaction contains the same species as a reactant
    ## under such as scenario the algorithm will consider that they are the same -- TODO: overwrite it
    if tgt_rxn_id in reactions_in_both:
        logger.warning(
            'The target_reaction ('+str(tgt_rxn_id)+') ' \
          + 'has been detected in model ' + str(gem_sbml_path.getName()) + ', ' \
          + 'ignoring this model...'
        )
        return None

    ## Check COBRA

    ######## FBA ########
    if sim_type == 'fraction':
        cobra_results, rpsbml_merged = rp_fraction(
                  rpsbml = rpsbml_merged,
              ignore_met = ignore_orphan_species,
              src_rxn_id = src_rxn_id,
               src_coeff = src_coeff,
              tgt_rxn_id = tgt_rxn_id,
               tgt_coeff = tgt_coeff,
             frac_of_src = frac_of_src,
                  is_max = is_max,
              pathway_id = pathway_id,
            objective_id = objective_id,
                  logger = logger
        )
    else:
        objective_id = rpsbml_merged.find_or_create_objective(
            reactions = [tgt_rxn_id],
            coefficients = [tgt_coeff],
            is_max = is_max,
            objective_id = objective_id
        )
        if sim_type == 'fba':
            cobra_results = rp_fba(
                      rpsbml = rpsbml_merged,
                objective_id = objective_id,
                  ignore_met = ignore_orphan_species,
                      logger = logger
            )
        ####### pFBA #######
        elif sim_type == 'pfba':
            cobra_results = rp_pfba(
                      rpsbml = rpsbml_merged,
                objective_id = objective_id,
                  ignore_met = ignore_orphan_species,
                 frac_of_opt = frac_of_src,
                      logger = logger
            )
        else:
            logger.error('Cannot recognise sim_type: ' + str(sim_type))
            return None

        if cobra_results is None:
            return None

        write_results(
            rpsbml = rpsbml_merged,
            objective_id = objective_id,
            cobra_results = cobra_results,
            pathway_id = pathway_id,
            logger = logger
        )
        # for group in rpsbml_merged.getModel().getPlugin('groups').getListOfGroups():
        #     print(group)

    if cobra_results is None:
        return None

    '''
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    '''
    if not merge:
        complete_heterologous_pathway(
            rpsbml = rpsbml,
            rpsbml_merged = rpsbml_merged,
            species_group_id = species_group_id,
            sink_species_group_id = sink_species_group_id,
            pathway_id = pathway_id,
            reactions_in_both = reactions_in_both,
            logger = logger
        )
        logger.info('Returning model with heterologous pathway only')
        return rpsbml
    else:
        logger.info('Returning the full model')
        return rpsbml_merged


def complete_heterologous_pathway(
    rpsbml: rpSBML,
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
            logger.debug('Removing central_species')
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
    create_ignored_species_group(rpsbml, logger)
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


def rp_fba(
          rpsbml: rpSBML,
    objective_id: str,
      ignore_met: bool = True,
          logger: Logger = getLogger(__name__)
) -> cobra_solution:
    """Run FBA using a single objective

    :param reaction_id: The id of the reactions involved in the objective
    :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type reaction_id: str
    :type coefficient: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """
    logger.info('Running FBA...')
    logger.debug('rpsbml:       ' + str(rpsbml))
    logger.debug('ignore_met:  ' + str(ignore_met))

    cobra_results = runCobra(
        rpsbml = rpsbml,
        objective_id = objective_id,
        ignore_met = ignore_met,
        logger = logger
    )

    return cobra_results


def rp_pfba(
          rpsbml: rpSBML,
    objective_id: str,
      ignore_met: bool = True,
     frac_of_opt: float = 0.95,
          logger: Logger = getLogger(__name__)
) -> cobra_solution:
    """Run parsimonious FBA using a single objective

    :param reaction_id: The id of the reactions involved in the objective
    :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
    :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.95)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type reaction_id: str
    :type coefficient: float
    :type fraction_of_optimum: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """
    logger.info('Running FBA (parsimonious)...')
    logger.debug('rpsbml:       ' + str(rpsbml))
    logger.debug('ignore_met:  ' + str(ignore_met))
    logger.debug('frac_of_opt:  ' + str(frac_of_opt))

    rpsbml.activateObjective(
        objective_id = objective_id,
        plugin = 'fbc'
    )

    cobraModel = cobra_model(
        rpsbml = rpsbml,
        logger = logger
    )
    if not cobraModel:
        return None

    if ignore_met:
        isolated_species = [cmp.replace('__64__', '@') for cmp in rpsbml.get_isolated_species()]
        # print([cobraModel.metabolites.get_by_id(met) for met in isolated_species])
        cobraModel.remove_metabolites([cobraModel.metabolites.get_by_id(met) for met in isolated_species])

    cobra_results = pfba(cobraModel, frac_of_opt)

    logger.debug(cobra_results)

    return cobra_results


def rp_fraction(
          rpsbml: rpSBML,
      src_rxn_id:    str,
       src_coeff:  float,
      tgt_rxn_id:    str,
       tgt_coeff:  float,
     frac_of_src:  float = 0.75,
          is_max:   bool = True,
      pathway_id:    str = 'rp_pathway',
    objective_id:    str = None,
      ignore_met:   bool = True,
          logger: Logger = getLogger(__name__)
) -> Tuple[cobra_solution, rpSBML]:
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
    
    logger.info('Processing FBA (fraction of reaction)...')
    logger.debug('rpsbml:       ' + str(rpsbml))
    logger.debug('ignore_met:   ' + str(ignore_met))
    logger.debug('src_rxn_id:   ' + src_rxn_id)
    logger.debug('src_coeff:    ' + str(src_coeff))
    logger.debug('tgt_rxn_id:   ' + tgt_rxn_id)
    logger.debug('tgt_coeff:    ' + str(tgt_coeff))
    logger.debug('frac_of_src:  ' + str(frac_of_src))
    logger.debug('is_max:       ' + str(is_max))
    logger.debug('pathway_id:   ' + pathway_id)
    logger.debug('objective_id: ' + str(objective_id))

    def get_annot_objective(
        rpsbml: rpSBML,
        objective_id: str
    ) -> str:

        fbc_plugin = rpsbml.getPlugin('fbc')
        fbc_obj = fbc_plugin.getObjective(objective_id)

        fbc_obj_annot = fbc_obj.getAnnotation()
        # TODO: if this is None need to set it up
        if fbc_obj_annot is None:
            return None

        logger.debug('Already calculated flux for '+str(objective_id))

        return fbc_obj_annot

    # retreive the biomass objective and flux results and set as maxima
    src_obj_id = rpsbml.find_or_create_objective(
        reactions = [src_rxn_id],
        coefficients = [src_coeff],
        is_max = is_max
    )

    # objective = rpsbml.getPlugin(
    #     'fbc'
    # ).getObjective(objective_id)

    # rpsbml.checklibSBML(
    #     objective,
    #     'Getting objective '+str(objective_id)
    # )



    # fbc_plugin = rpsbml.getPlugin('fbc')
    # TODO: use the rpSBML BRSynth annotation parser
    # try:
    fbc_obj_annot = get_annot_objective(rpsbml, src_obj_id)

    # except (AttributeError, ValueError) as e:
    if fbc_obj_annot is None:
        # logger.debug(e)
        logger.debug('Performing FBA to calculate the source reaction')

        ### FBA ###
        # logger.info('Running the FBA (fraction of reaction)...')
        # rpsbml.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
        cobra_results = runCobra(
            rpsbml = rpsbml,
            objective_id = src_obj_id,
            ignore_met = ignore_met,
            logger = logger
        )
        if cobra_results is None:
            return None, rpsbml

        write_results(
            rpsbml = rpsbml,
            objective_id = src_obj_id,
            cobra_results = cobra_results,
            pathway_id = pathway_id,
            logger = logger
        )

        # cobra_results.objective_value
        fbc_obj_annot = get_annot_objective(rpsbml, src_obj_id)
        if fbc_obj_annot is None:
            logger.error('No annotation available for: '+str(src_obj_id))

    source_flux = float(
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

    # TODO: add another to check if the objective id exists
    logger.debug('FBA source flux ('+str(src_rxn_id)+') is: '+str(source_flux))
    if not objective_id:
        objective_id = 'obj_'+str(tgt_rxn_id)+'__restricted_'+str(src_rxn_id)

    objective_id = rpsbml.find_or_create_objective(
        reactions = [tgt_rxn_id],
        coefficients = [tgt_coeff],
        is_max = is_max,
        objective_id = objective_id
    )

    logger.debug('Optimising the objective: '+str(objective_id))
    logger.debug('     Setting upper bound: '+str(source_flux*frac_of_src))
    logger.debug('     Setting lower bound: '+str(source_flux*frac_of_src))

    old_upper_bound, old_lower_bound = rpsbml.getReactionConstraints(src_rxn_id)
    rpsbml.setReactionConstraints(
        src_rxn_id,
        source_flux*frac_of_src,
        source_flux*frac_of_src
    )

    cobra_results = runCobra(
        rpsbml = rpsbml,
        objective_id = objective_id,
        ignore_met = ignore_met,
        logger = logger
    )
    if cobra_results is None:
        return None, rpsbml

    write_results(
        rpsbml = rpsbml,
        objective_id = objective_id,
        cobra_results = cobra_results,
        pathway_id = pathway_id,
        logger = logger
    )

    ##### print the biomass results ######
    logger.debug('Biomass: '+str(cobra_results.fluxes.biomass))
    logger.debug(' Target: '+str(cobra_results.fluxes.rxn_target))

    # reset the bounds to the original values for the target
    rpsbml.setReactionConstraints(
        src_rxn_id,
        old_upper_bound,
        old_lower_bound
    )

    logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))

    return cobra_results, rpsbml


def runCobra(
    rpsbml: rpSBML,
    objective_id: str,
    ignore_met: bool = True,
    logger: Logger = getLogger(__name__)
) -> cobra_solution:
    """
    Run Cobra and write results to the rpsbml object.

    Parameters
    ----------
    rpsbml: rpSBML
        rpSBML object of which reactions will be updated with results
    objective_id: str
        The id of the objective to optimise
    pathway_id: str
        The id of the pathway within reactions will be updated
    logger
        Logger object
    """

    rpsbml.activateObjective(
        objective_id = objective_id,
        plugin = 'fbc'
    )

    cobraModel = cobra_model(
        rpsbml = rpsbml,
        logger = logger
    )

    if not cobraModel:
        return None

    # print("Reactions")
    # print("---------")
    # for x in cobraModel.reactions:
    #     print("%s : %s" % (x.id, x.reaction))

    # x = cobraModel.reactions.rxn_1
    # print("%s : %s" % (x.id, x.reaction))
    # x = cobraModel.reactions.rxn_2
    # print("%s : %s" % (x.id, x.reaction))
    # x = cobraModel.reactions.rxn_3
    # print("%s : %s" % (x.id, x.reaction))
    # x = cobraModel.reactions.rxn_target
    # print("%s : %s" % (x.id, x.reaction))

    # remove isolated species from the model
    if ignore_met:
        isolated_species = [cmp.replace('__64__', '@') for cmp in rpsbml.get_isolated_species()]
        # print("")
        # print("Metabolites")
        # print("-----------")
        # for x in cobraModel.metabolites:
        #     if x.id in isolated_species:
        #         print('%9s : %s' % (x.id, x.formula))
        cobraModel.remove_metabolites([cobraModel.metabolites.get_by_id(met) for met in isolated_species])

    cobra_results = cobraModel.optimize(
        objective_sense='maximize',
        raise_error=True
    )

    logger.debug(cobra_results)

    return cobra_results


def cobra_model(
    rpsbml: rpSBML,
    logger: Logger = getLogger(__name__)
) -> cobra_model:
    """Convert the rpSBML object to cobra object

    :return: Success or failure of the function
    :rtype: bool
    """

    rpsbml.logger.info('Creating Cobra object from rpSBML...')

    with NamedTemporaryFile() as temp_f:
        rpsbml.writeToFile(temp_f.name)
        try:
            cobraModel = cobra_io.read_sbml_model(temp_f.name, use_fbc_package=True)
        except CobraSBMLError:
            logger.error('Something went wrong reading the SBML model')
            (model, errors) = validate_sbml_model(temp_f.name)
            logger.error(str(json_dumps(errors, indent=4)))
            return None

    logger.debug(cobraModel)

    return cobraModel


def write_results(
    rpsbml: rpSBML,
    objective_id: str,
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
    rpsbml.logger.debug('objective_id: ' + str(objective_id))
    rpsbml.logger.debug('cobra_results: ' + str(cobra_results))
    rpsbml.logger.debug('cobra_results.objective_value: ' + str(cobra_results.objective_value))
    rpsbml.logger.debug('pathway_id: ' + str(pathway_id))

    write_objective_to_pathway(
        cobra_results.objective_value,
        rpsbml,
        pathway_id,
        objective_id,
        logger
    )

    create_ignored_species_group(
        rpsbml,
        logger
    )

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
        objective_id,
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
        str(objective_value),
        'mmol_per_gDW_per_hr',
        False
    )

    for flux_obj in obj.getListOfFluxObjectives():

        rxn = flux_obj.getReaction()
        flux = fluxes.get(rxn)

        # sometimes flux cannot be returned
        if flux is None:
            rpsbml.logger.warning(
                'Cobra BUG: Cannot retreive ' + str(rxn) \
              + ' flux from cobrapy... setting to 0.0'
            )
            flux = 0.0

        rpsbml.logger.debug(
            'Set the reaction ' + str(rxn) \
            + ' a flux_value of ' + str(flux)
        )
        rpsbml.updateBRSynth(
            flux_obj,
            'flux_value',
            str(flux),
            'mmol_per_gDW_per_hr',
            False
        )


def write_fluxes_to_reactions(
    fluxes: np_series,
    rpsbml: rpSBML,
    pathway_id: str,
    objective_id: str,
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
          + ' a ' + str('fba_' + str(objective_id)) \
          + ' of ' + str(flux))

        rpsbml.updateBRSynth(
            rxn,
            'fba_'+str(objective_id),
            str(flux),
            'mmol_per_gDW_per_hr',
            False
        )


def write_objective_to_pathway(
    objective_value: float,
    rpsbml: rpSBML,
    pathway_id: str,
    objective_id: str,
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

    rp_pathway = rpsbml.getGroup(pathway_id)

    logger.debug(
        'Set ' + str(pathway_id) + ' with ' \
      + str('fba_'+str(objective_id)) + ' to ' \
      + str(objective_value)
    )

    rpsbml.updateBRSynth(
        rp_pathway,
        'fba_'+str(objective_id),
        str(objective_value),
        'mmol_per_gDW_per_hr',
        False
    )


def create_ignored_species_group(
    rpsbml: rpSBML,
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

    group_id = 'ignored_species_for_FBA'

    rpsbml.createGroup(
        id = group_id,
        brs_annot = False
    )

    for spe in rpsbml.get_isolated_species():
        rpsbml.addMember(
            group_id = group_id,
            idRef = spe
        )

    # for group in rpsbml.getModel().getPlugin('groups').getListOfGroups():
    #     print(group)

    # exit()


    logger.debug(
        'Create ' + str(group_id) + ' group with ' \
      + str(rpsbml.get_isolated_species())
    )

    # ignored_spe = rpsbml.getGroup(
    #     rpsbml.getPlugin('groups'),
    #     group_id
    # )

    # rpsbml.updateBRSynth(
    #     sbase_obj = ignored_spe,
    #     annot_header = 'ignored_species',
    #     value = str(rpsbml.get_isolated_species())
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


