from logging import (
    Logger,
    getLogger
)
from cobra.flux_analysis import pfba
from cobra               import io    as cobra_io
from cobra.core.model    import Model as cobra_model
from rptools.rplibs      import rpSBML
from typing import (
    List,
    Dict,
    Tuple
)
from tempfile import (
    NamedTemporaryFile
)


# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes


def rp_fba(
          rpsbml: rpSBML,
     reaction_id:    str,
     coefficient:  float = 1.0,
          is_max:   bool = True,
      pathway_id:    str = 'rp_pathway',
    objective_id:    str = None,
          logger: Logger = getLogger(__name__)
) -> Tuple[float, rpSBML]:
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
    logger.debug('reaction_id:  ' + reaction_id)
    logger.debug('coefficient:  ' + str(coefficient))
    logger.debug('is_max:       ' + str(is_max))
    logger.debug('pathway_id:   ' + pathway_id)
    logger.debug('objective_id: ' + str(objective_id))

    return _rp_fba(
              rpsbml = rpsbml,
         reaction_id = reaction_id,
            sim_type = 'fba',
         coefficient = coefficient,
              is_max = is_max,
          pathway_id = pathway_id,
        objective_id = objective_id,
              logger = logger
    )


def rp_pfba(
          rpsbml: rpSBML,
     reaction_id:    str,
     coefficient:  float = 1.0,
     frac_of_opt:  float = 0.95,
          is_max:   bool = True,
      pathway_id:    str = 'rp_pathway',
    objective_id:    str = None,
          logger: Logger = getLogger(__name__)
) -> Tuple[float, rpSBML]:
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
    logger.debug('reaction_id:  ' + reaction_id)
    logger.debug('coefficient:  ' + str(coefficient))
    logger.debug('frac_of_opt:  ' + str(frac_of_opt))
    logger.debug('is_max:       ' + str(is_max))
    logger.debug('pathway_id:   ' + pathway_id)
    logger.debug('objective_id: ' + str(objective_id))

    return _rp_fba(
              rpsbml = rpsbml,
         reaction_id = reaction_id,
            sim_type = 'pfba',
         coefficient = coefficient,
         frac_of_opt = frac_of_opt,
              is_max = is_max,
          pathway_id = pathway_id,
        objective_id = objective_id,
              logger = logger
    )


def _rp_fba(
          rpsbml: rpSBML,
     reaction_id:    str,
        sim_type:    str,
     coefficient:  float = 1.0,
     frac_of_opt:  float = 0.95,
          is_max:   bool = True,
      pathway_id:    str = 'rp_pathway',
    objective_id:    str = None,
          logger: Logger = getLogger(__name__)
) -> Tuple[float, rpSBML]:
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
    
    fbc_plugin = rpsbml.getPlugin('fbc')

    objective_id = rpsbml.FindOrCreateObjective(
        [reaction_id],
        [coefficient],
        is_max,
        objective_id
    )

    rpsbml.checklibSBML(
        fbc_plugin.setActiveObjectiveId(objective_id),
        'Setting active objective '+str(objective_id)
    )

    cobraModel = cobra(rpsbml)
    if not cobraModel:
        return -1, rpSBML
        # logger.error('cobraModel is None')
        # return 0.0, rpsbml

    # run the FBA
    if sim_type == 'fba':
        cobra_results = cobraModel.optimize()
    elif sim_type == 'pfba':
        cobra_results = pfba(cobraModel, frac_of_opt)
    else:
        logger.error('Cannot recognise sim_type: ' + str(sim_type))
        return -2, rpsbml

    addAnalysisResults(
        rpsbml,
        objective_id,
        cobra_results,
        pathway_id
    )

    return cobra_results.objective_value, rpsbml


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
          logger: Logger = getLogger(__name__)
) -> Tuple[float, rpSBML]:
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
    logger.debug('src_rxn_id:   ' + src_rxn_id)
    logger.debug('src_coeff:    ' + str(src_coeff))
    logger.debug('tgt_rxn_id:   ' + tgt_rxn_id)
    logger.debug('tgt_coeff:    ' + str(tgt_coeff))
    logger.debug('frac_of_src:  ' + str(frac_of_src))
    logger.debug('is_max:       ' + str(is_max))
    logger.debug('pathway_id:   ' + pathway_id)
    logger.debug('objective_id: ' + str(objective_id))


    # retreive the biomass objective and flux results and set as maxima
    fbc_plugin = rpsbml.getPlugin('fbc')
    source_obj_id = rpsbml.FindOrCreateObjective([src_rxn_id], [src_coeff], is_max)

    # TODO: use the rpSBML BRSynth annotation parser
    source_flux = None
    try:
        fbc_obj = fbc_plugin.getObjective(source_obj_id)
        # TODO: if this is None need to set it up
        fbc_obj_annot = fbc_obj.getAnnotation()
        if not fbc_obj_annot:
            raise ValueError
        source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        logger.debug('Already calculated flux for '+str(source_obj_id))

    except (AttributeError, ValueError) as e:
        logger.debug(e)
        logger.debug('Performing FBA to calculate the source reaction')

        ### FBA ###
        # logger.info('Running the FBA (fraction of reaction)...')
        # self.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
        rpsbml.checklibSBML(
            fbc_plugin.setActiveObjectiveId(
                source_obj_id
            ),
            'Setting active objective '+str(source_obj_id)
        )

        cobraModel = cobra(rpsbml)
        if not cobraModel:
            logger.error('Converting libSBML to CobraPy returned False')
            addAnalysisResults(
                rpsbml,
                source_obj_id,
                0.0,
                pathway_id
            )
            return 0.0, None

        cobra_results = cobraModel.optimize()
        addAnalysisResults(
            rpsbml,
            source_obj_id,
            cobra_results,
            pathway_id
        )

        # cobra_results.objective_value
        fbc_obj = fbc_plugin.getObjective(source_obj_id)
        fbc_obj_annot = fbc_obj.getAnnotation()
        if fbc_obj_annot is None:
            logger.error('No annotation available for: '+str(source_obj_id))
        source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))

    # TODO: add another to check if the objective id exists
    logger.debug('FBA source flux ('+str(src_rxn_id)+') is: '+str(source_flux))
    if not objective_id:
        objective_id = 'obj_'+str(tgt_rxn_id)+'__restricted_'+str(src_rxn_id)

    objective_id = rpsbml.FindOrCreateObjective(
        [tgt_rxn_id],
        [tgt_coeff],
        is_max,
        objective_id
    )

    logger.debug('Optimising the objective: '+str(objective_id))
    logger.debug('Setting upper bound: '+str(source_flux*frac_of_src))
    logger.debug('Setting loer bound: '+str(source_flux*frac_of_src))
    old_upper_bound, old_lower_bound = rpsbml.setReactionConstraints(
        src_rxn_id,
        source_flux*frac_of_src,
        source_flux*frac_of_src
    )
    rpsbml.checklibSBML(
        fbc_plugin.setActiveObjectiveId(objective_id),
        'Setting active objective '+str(objective_id)
    )
    cobraModel = cobra(rpsbml)
    if not cobraModel:
        # although this may not be the greatest idea, set flux to 0.0 when cobrapy error
        addAnalysisResults(
            rpsbml,
            objective_id,
            0.0,
            pathway_id
        )
        return 0.0, None
    cobra_results = cobraModel.optimize()
    addAnalysisResults(
        rpsbml,
        objective_id,
        cobra_results,
        pathway_id
    )
    ##### print the biomass results ######
    logger.debug('Biomass: '+str(cobra_results.fluxes.biomass))
    logger.debug('Target: '+str(cobra_results.fluxes.Rxn_sink))
    # reset the bounds to the original values for the target
    old_upper_bound, old_lower_bound = rpsbml.setReactionConstraints(
        src_rxn_id,
        old_upper_bound,
        old_lower_bound
    )
    logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))
    return cobra_results.objective_value, rpsbml


# TODO: do not use the species_group_id and the sink_species_group_id. Loop through all the groups (and if the same) and overwrite the annotation instead
def runFBA(
              rpsbml_path: str,
            gem_sbml_path: str,
                  outFile: str,
                 sim_type: str,
               src_rxn_id: str,
                src_coeff: float,
               tgt_rxn_id: str,
                tgt_coeff: float,
                   is_max: bool = True,
              frac_of_src: float = 0.75,
               dont_merge: bool = True,
               pathway_id: str = 'rp_pathway',
             objective_id: str = None,
           compartment_id: str = 'MNXC3',
    # fill_orphan_species=False,
         species_group_id: str = 'central_species',
    sink_species_group_id: str = 'rp_sink_species',
                   logger: Logger = getLogger(__name__)
) -> Tuple[float, rpSBML]:
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
    :param dont_merge: Output the merged model (Default: True)
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
    :type dont_merge: bool
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

    logger.debug('rpsbml_path:           ' + str(rpsbml_path))
    logger.debug('gem_sbml_path:         ' + str(gem_sbml_path))
    logger.debug('outFile:               ' + str(outFile))
    logger.debug('sim_type:              ' + str(sim_type))
    logger.debug('src_rxn_id:            ' + src_rxn_id)
    logger.debug('src_coeff:             ' + str(src_coeff))
    logger.debug('tgt_rxn_id:            ' + tgt_rxn_id)
    logger.debug('tgt_coeff:             ' + str(tgt_coeff))
    logger.debug('frac_of_src:           ' + str(frac_of_src))
    logger.debug('is_max:                ' + str(is_max))
    logger.debug('dont_merge:            ' + str(dont_merge))
    logger.debug('pathway_id:            ' + pathway_id)
    logger.debug('objective_id:          ' + str(objective_id))
    logger.debug('compartment_id:        ' + str(compartment_id))
    logger.debug('species_group_id:      ' + str(species_group_id))
    logger.debug('sink_species_group_id: ' + str(sink_species_group_id))

    rpsbml = rpSBML(rpsbml_path, logger=logger)
    logger.debug('input_sbml: ' + str(rpsbml))

    # Save the central species
    # groups = rpsbml.getPlugin('groups')
    # central = groups.getGroup(species_group_id)
    # sink_group = groups.getGroup(sink_species_group_id)
    # rp_group = groups.getGroup(pathway_id)
    # cent_spe = [str(i.getIdRef()) for i in central.getListOfMembers()]
    # sink_spe = [str(i.getIdRef()) for i in sink_group.getListOfMembers()]
    # rp_reac  = [str(i.getIdRef()) for i in rp_group.getListOfMembers()]
    # logger.debug('old central species: '+str(cent_spe))
    # logger.debug('old sink species: '+str(sink_spe))
    # logger.debug('old rp reactions: '+str(rp_reac))

    rpsbml_gem = rpSBML(gem_sbml_path, logger=logger)
    logger.debug('rpsbml_gem: ' + str(rpsbml_gem))

    logger.info('Merging rpSBML models: ' + rpsbml.getName() + ' and ' + rpsbml_gem.getName() + '...')
    species, reactions = rpSBML.mergeModels(
        source_rpsbml = rpsbml,
        target_rpsbml = rpsbml_gem,
        logger = logger
    )
    logger.debug('input_sbml: ' + str(rpsbml))
    logger.debug('rpsbml_gem: ' + str(rpsbml_gem))

    # NOTE: reactions is organised with key being the rpsbml reaction and value being the rpsbml_gem value`
    # BUG: when merging the Rxn_sink (very rare cases) can be recognised if another reaction contains the same species as a reactant
    ## under such as scenario the algorithm will consider that they are the same -- TODO: overwrite it
    if tgt_rxn_id in reactions:
        logger.warning('The target_reaction ('+str(tgt_rxn_id)+') has been detected in model '+str(outFile)+', ignoring this model...')
        return False
    # rev_reactions = {v: k for k, v in reactions.items()}
    # logger.debug('species: ' + str(species))
    # logger.debug('reactions:     ' + str(reactions))
    # logger.debug('rev_reactions: ' + str(rev_reactions))

    ######## FBA ########
    if sim_type == 'fraction':
        obj_val, rpsbml_gem = rp_fraction(
                  rpsbml = rpsbml_gem,
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
    elif sim_type == 'fba':
        obj_val, rpsbml_gem = rp_fba(
                  rpsbml = rpsbml_gem,
             reaction_id = tgt_rxn_id,
             coefficient = tgt_coeff,
                  is_max = is_max,
              pathway_id = pathway_id,
            objective_id = objective_id,
                  logger = logger
        )
    ####### pFBA #######
    elif sim_type == 'pfba':
        obj_val, rpsbml_gem = rp_pfba(
                  rpsbml = rpsbml_gem,
             reaction_id = tgt_rxn_id,
             coefficient = tgt_coeff,
             frac_of_opt = frac_of_src,
                  is_max = is_max,
              pathway_id = pathway_id,
            objective_id = objective_id,
                  logger = logger
        )
    else:
        logger.error('Cannot recognise sim_type: ' + str(sim_type))
        return 0.0, None

    if obj_val < 0:
        return obj_val, rpsbml

    '''
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    '''
    if dont_merge:
        logger.info('Building model with heterologous pathway only')
        groups = rpsbml_gem.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        logger.debug('---- Reactions ----')
        for member in rp_pathway.getListOfMembers():
            #### reaction annotation
            logger.debug(member.getIdRef())
            reacFBA = rpsbml_gem.getModel().getReaction(member.getIdRef())
            logger.debug(reacFBA)
            try:
                #reacIN = rpsbml.model.getReaction(reactions_convert[member.getIdRef()])
                reacIN = rpsbml.getModel().getReaction(rev_reactions_convert[member.getIdRef()])
            except KeyError:
                reacIN = rpsbml.getModel().getReaction(member.getIdRef())
            logger.debug(reacIN)
            logger.debug(reacFBA.getAnnotation())
            reacIN.setAnnotation(reacFBA.getAnnotation())
            #### species TODO: only for shadow price
        #### add groups ####
        source_groups = rpsbml_gem.getPlugin('groups')
        target_groups = rpsbml.getPlugin('groups')
        target_groupsID = [i.getId() for i in target_groups.getListOfGroups()]
        for source_group in source_groups.getListOfGroups():
            logger.info('Replacing group id: '+str(source_group.getId()))
            if source_group.getId()==species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                # TODO: #### replace the new potentially incorect central species with the normal ones #####
                # delete all the previous members
                logger.info('Removing central_species')
                for i in range(target_group.getNumMembers()):
                    logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                # add the new ones
                for cs in cent_spe:
                    logger.info('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId()==sink_species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                logger.info('Removing sink species')
                for i in range(target_group.getNumMembers()):
                    logger.info('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in sink_spe:
                    logger.info('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId() in target_groupsID:
                target_group = target_groups.getGroup(source_group.getId())
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
        #### add objectives ####
        source_fbc = rpsbml_gem.getPlugin('fbc')
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
        #rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
        target_fbc.setActiveObjectiveId(source_obj_id) #tmp random assigenement of objective
        logger.info('Writing model with heterologous pathway only into file..')
        rpsbml.writeSBML(outFile)
    else:
        logger.info('Writing the full model into file..')
        rpsbml_gem.writeSBML(outFile)

    logger.info('  |--> written in ' + outFile)

    return True


def cobra(
    rpsbml: rpSBML
) -> cobra_model:
    """Convert the rpSBML object to cobra object

    :return: Success or failure of the function
    :rtype: bool
    """

    rpsbml.logger.info('Creating Cobra object from rpSBML...')

    try:
        with NamedTemporaryFile() as temp_f:
            rpsbml.writeSBML(temp_f.name)
            #self.logger.info(glob(tmpOutputFolder+'/*'))
            #self.logger.info(cobra.io.validate_sbml_model(glob(tmpOutputFolder+'/*')[0]))
            cobraModel = cobra_io.read_sbml_model(temp_f.name, use_fbc_package=True)
        #self.cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(), use_fbc_package=True)
        #use CPLEX
        # self.cobraModel.solver = 'cplex'

    except cobra_io.sbml.CobraSBMLError as e:
        rpsbml.logger.error(e)
        return None

    return cobraModel


def addAnalysisResults(
    rpsbml: rpSBML,
    objective_id: str,
    cobra_results,
    pathway_id: str = 'rp_pathway'
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

    rp_pathway = rpsbml.getGroup(
        rpsbml.getPlugin('groups'),
        pathway_id
    )

    # write the results to the rp_pathway
    rpsbml.logger.debug('Set '+str(pathway_id)+' with '+str('fba_'+str(objective_id))+' to '+str(cobra_results.objective_value))
    rpsbml.updateBRSynth(
        rp_pathway,
        'fba_'+str(objective_id),
        str(cobra_results.objective_value),
        'mmol_per_gDW_per_hr',
        False
    )

    # get the objective
    obj = getObjective(
        rpsbml,
        objective_id,
        cobra_results.objective_value
    )

    rpsbml.logger.debug('Set the objective '+str(objective_id)+' a flux_value of '+str(cobra_results.objective_value))
    for flux_obj in obj.getListOfFluxObjectives():
        # sometimes flux cannot be returned
        if cobra_results.fluxes.get(flux_obj.getReaction()) is None:
            rpsbml.logger.warning('Cobra BUG: Cannot retreive '+str(flux_obj.getReaction())+' flux from cobrapy... setting to 0.0')
            rpsbml.updateBRSynth(
                flux_obj,
                'flux_value',
                str(0.0),
                'mmol_per_gDW_per_hr',
                False
            )
            rpsbml.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(0.0))
        else:
            rpsbml.updateBRSynth(
                flux_obj,
                'flux_value',
                str(cobra_results.fluxes.get(flux_obj.getReaction())),
                'mmol_per_gDW_per_hr',
                False
            )
            rpsbml.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(cobra_results.fluxes.get(flux_obj.getReaction())))

    # write all the results to the reactions of pathway_id
    for member in rp_pathway.getListOfMembers():
        reac = rpsbml.getModel().getReaction(member.getIdRef())
        if reac is None:
            rpsbml.logger.error('Cannot retreive the following reaction: '+str(member.getIdRef()))
            #return False
            continue
        rpsbml.logger.debug('Set the reaction '+str(member.getIdRef())+' a '+str('fba_'+str(objective_id))+' of '+str(cobra_results.fluxes.get(reac.getId())))
        rpsbml.updateBRSynth(
            reac,
            'fba_'+str(objective_id),
            str(cobra_results.fluxes.get(reac.getId())),
            'mmol_per_gDW_per_hr',
            False
        )


def getObjective(
    rpsbml: rpSBML,
    objective_id: str,
    objective_value: float 
) -> None:

    fbc_plugin = rpsbml.getPlugin('fbc')

    objective = fbc_plugin.getObjective(objective_id)

    rpsbml.checklibSBML(
        objective,
        'Getting objective '+str(objective_id)
    )

    rpsbml.updateBRSynth(
        objective,
        'flux_value',
        str(objective_value),
        'mmol_per_gDW_per_hr',
        False
    )

    return objective




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
#     rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
#     return rpsbml


