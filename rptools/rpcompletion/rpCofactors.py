import logging
from copy import deepcopy
from typing import (
    Dict,
    List,
    Tuple
)



################################################################
######################### PUBLIC FUNCTIONS #####################
################################################################

## Given a dictionnary describing a monocomponent reaction, add the cofactors by comparing it with the original reaction
#
# @param step Dictionnary describing the reaction
# @param reac_side String 'right' or 'left' describing the direction of the monocomponent reaction compared with the original reaction
# @param rr_reac Dictionnary describing the monocomponent reaction from RetroRules
# @param f_reac Dictionnary describing the full original reaction
# @param pathway_cmp Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway
#
def complete_reac(cache, rxn_side, rr_reac, full_reac, mono_side, reac_smiles_side, pathway_cmp, logger=logging.getLogger(__name__)):

    if mono_side:
        ## add the unknown species to pathway_cmp for the next steps
        rr_mono_cmp   = list(rr_reac.keys())
        step_mono_cmp = list(rxn_side.keys())
        if (len(rr_mono_cmp) == len(step_mono_cmp) == 1):
            # this is purposely overwitten since the main cmp between reactions can change
            pathway_cmp[step_mono_cmp[0]] = rr_mono_cmp[0]
        else:
            logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
            return False, {}, {}

    ## add the side species
    smiles_side, side = add_side_species(cache, full_reac, rr_reac, logger=logger)

    ## Update the stochio
    smiles_side, side_stochio = update_stochio(cache, side, full_reac, reac_smiles_side+smiles_side, pathway_cmp, logger=logger)
    side.update(side_stochio)
    return True, smiles_side, side


def add_side_species(cache, full_reac, rr_reac, logger=logging.getLogger(__name__)):
    
    smiles = ''
    side = {}

    for toAdd in full_reac.keys()-rr_reac.keys():

        side.update({toAdd: full_reac[toAdd]})

        ### update the reaction rule string
        try:
            smi = cache.get('cid_strc')[toAdd]['smiles']
            if smi is not None:
                for sto_add in range(int(full_reac[toAdd])):
                    smiles += '.'+str(smi)
        except KeyError:
            logger.debug('Cannot find smiles structure for '+str(toAdd))

    return smiles, side


def update_stochio(cache, rxn_side, full_reac, reac_smiles_side, pathway_cmp, logger=logging.getLogger(__name__)):
    
    side = deepcopy(rxn_side)

    for step_spe in side:

        if step_spe in full_reac:

            if not side[step_spe] == full_reac[step_spe]:
                stochio_diff = full_reac[step_spe]-side[step_spe]
                side[step_spe] = full_reac[step_spe]
                if stochio_diff<0:
                    logger.warning('full_reac stochio should never be smaller than step')
                    continue
                for i in range(stochio_diff):
                    ### update the reaction rule string
                    try:
                        smi = cache.get('cid_strc')[step_spe]['smiles']
                        if not smi==None:
                            reac_smiles_side += '.'+str(smi)
                    except KeyError:
                        #@Mel toAdd -> step_spe
                        logger.warning('Cannot find smiles structure for '+str(step_spe))

        elif step_spe in pathway_cmp:
            if pathway_cmp[step_spe] in full_reac:
                if side[step_spe] != full_reac[pathway_cmp[step_spe]]:
                    side[step_spe] = full_reac[pathway_cmp[step_spe]]

    return reac_smiles_side, side


def complete_reaction(cache, rxn, reaction_from_rr,
                      full_reaction_from_rr_1, full_reaction_from_rr_2,
                      mono_side, pathway_cmp,
                      logger=logging.getLogger(__name__)):

    reac_smiles = {}
    rxn_annot = rxn['brsynth']
    reac_smiles['left'], reac_smiles['right'] = rxn_annot['smiles'].split('>>')

    try:
        # LEFT SIDE
        isSuccess, \
        reac_smiles['left'], \
        rxn_right = complete_reac(cache,
                                  rxn['right'],
                                  reaction_from_rr['left'],
                                  full_reaction_from_rr_1,
                                  True,
                                  reac_smiles['left'],
                                  pathway_cmp,
                                  logger=logger)
        if not isSuccess:
            logger.warning('Could not recognise reaction rule for step: '+str(rxn))
            return {}, {}, reac_smiles['left'], reac_smiles['right']

        # RIGHT SIDE
        isSuccess, \
        reac_smiles['right'], \
        rxn_left = complete_reac(cache,
                                 rxn['left'],
                                 reaction_from_rr['right'],
                                 full_reaction_from_rr_2,
                                 False,
                                 reac_smiles['right'],
                                 pathway_cmp,
                                 logger=logger)
        if not isSuccess:
            logger.warning('Could not recognise reaction rule for step (2): '+str(rxn_annot['rxn_idx']))
            return {}, rxn_right, reac_smiles['left'], reac_smiles['right']

    except KeyError:
        logger.warning('Could not find the full reaction for reaction: '+str(rxn))
        return {}, {}, reac_smiles['left'], reac_smiles['right']

    return rxn_left, rxn_right, reac_smiles['left'], reac_smiles['right']


## Get the cofactors to monocomponent reactions
#
# @param rxn Reaction in the pathway
# @param pathway_cmp Dictionnary of intermediate compounds with their public ID's
# @return Boolean determine if the step is to be added
def get_cofactors_rxn(
    cache: 'rrCache',
    rxn: Dict,
    logger=logging.getLogger(__name__)
) -> Dict:

    # rxn_annot = rxn['brsynth']

    # print(dumps(rxn, indent=4))
    # rule_id = rxn_annot['rule_id']
    # print(rule_id, cache.get('rr_reactions')[rule_id])
    # rxn_id = cache._checkRIDdeprecated(rxn_annot['rule_ori_reac'], cache.get('deprecatedRID_rid'))
    # print(rxn_id, cache.get('template_reactions')[rxn_id])
    # print(pathway_cmp)

    from rxn_rebuild import rebuild_rxn

    rxn_annot = rxn['brsynth']
    rule_id = rxn_annot['rule_id']
    transfo = rxn_annot['smiles']
    # rxn_id = cache._checkRIDdeprecated(rxn_annot['rule_ori_reac'], cache.get('deprecatedRID_rid'))
    rxn_id = rxn_annot['rule_ori_reac']

    # rule_id = 'RR-02-f85f00f767901186-16-F'
    # transfo = '[H]OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])[H].O=P(O)(O)OP(=O)(O)O>>[H]OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C([H])([H])C(=C([H])[H])C([H])([H])[H].[H]OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])[H]'
    # rxn_id = 'MNXR100137'

    # if rxn_id == 'MNXR95743':
    #     print(rule_id)
    #     from json import dumps
    #     print(dumps(cache.get('rr_reactions')[rule_id], indent=4))
    #     print(dumps(rxn, indent=4))

    return rebuild_rxn(
        cache=cache,
        rxn_rule_id=rule_id,
        transfo=transfo,
        tmpl_rxn_id=rxn_id,
        logger=logger
    )

    # exit()
    reaction_from_rr      = cache.get('rr_reactions')[rule_id][rxn_id]
    full_reaction_from_rr = cache.get('template_reactions')[rxn_id]

    # Prepare arguments
    if reaction_from_rr['rel_direction'] == -1:
        full_reaction_from_rr_1 = full_reaction_from_rr['right']
        full_reaction_from_rr_2 = full_reaction_from_rr['left']

    elif reaction_from_rr['rel_direction'] == 1:
        full_reaction_from_rr_1 = full_reaction_from_rr['left']
        full_reaction_from_rr_2 = full_reaction_from_rr['right']

    else:
        logger.error('Relative direction can only be 1 or -1: '+str(reaction_from_rr['rel_direction']))
        return {}, {}, {}

    # Getting elements to complete the reaction
    rxn_left, rxn_right, \
    smiles_left, smiles_right = complete_reaction(cache,
                                                  rxn,
                                                  reaction_from_rr,
                                                  full_reaction_from_rr_1,
                                                  full_reaction_from_rr_2,
                                                  True,
                                                  pathway_cmp,
                                                  logger=logger)

    if 'full_transfo' in completed_transfos[rxn_id]:
        print(completed_transfos[rxn_id]['full_transfo'])
    else:
        print(transfo)
        print(dumps(completed_transfos[rxn_id], indent=4))
    print(rxn_left)
    print(rxn_right)
    print(smiles_left+'>>'+smiles_right)
    exit()
    return rxn_left, rxn_right, smiles_left+'>>'+smiles_right


def retrieve_infos(
    cache: 'rrCache',
    species: str,
    tmp_species: str,
    logger=logging.getLogger(__name__)
) -> Dict:

    xref      = {}
    inchi     = None
    inchikey  = None
    smiles    = None
    chem_name = None

    ###### Try to retreive the InChI ############
    try:
        inchi = cache.get('cid_strc')[tmp_species]['inchi']
    except KeyError:
        logger.debug('Cannot find the inchi for this species: '+str(tmp_species))

    ###### Try to retreive the InChIKey ############
    try:
        inchikey = cache.get('cid_strc')[tmp_species]['inchikey']
        logger.debug('Found the inchikey: '+str(inchikey))
        # #### TODO: find a better way to check if two species are the same ####
        # isfound = False
        # for rpsbml_species in tmp_species:
        #     # TODO add a comparison by xref as well
        #     logger.debug(str(species[rpsbml_species]['brsynth']['inchikey'])+' <--> '+str(inchikey))
        #     if str(species[rpsbml_species]['brsynth']['inchikey'])==str(inchikey):
        #         spe_conv[tmp_species] = rpsbml_species
        #         logger.debug('The species '+str(tmp_species)+' is the same as '+str(rpsbml_species))
        #         isfound = True
        #         break
        # # if isfound:
        # #     continue
    except KeyError:
        logger.debug('Cannot find the inchikey for this species: '+str(species))

    ##### Try to retreive the SMILES ############
    try:
        smiles = cache.get('cid_strc')[tmp_species]['smiles']
    except KeyError:
        logger.debug('Cannot find the smiles for this species: '+str(species))

    ###### Try to retreive the xref, using the inchikey if the cid fails #######
    try:
        xref = cache.get('cid_xref')[tmp_species]
    except KeyError:
        try:
            xref = cache.get('cid_xref')[tmp_species]
        except KeyError:
            # if you cannot find using cid, try to retreive it using its inchikey
            try:
                if inchikey:
                    # @Joan: Can you think of a better way of doing that?
                    # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                    tmp_cids = [i for i in cache.get('inchikey_cid')[inchikey] if i[:3]=='MNX']
                    # TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                    if tmp_cids:
                        xref = cache.get('cid_xref')[cache._checkCIDdeprecated(tmp_cids[0], cache.get('deprecatedCID_cid'))]
            except KeyError:
                logger.debug('Cannot find the xref for this species: '+str(species))
                xref = {}

    #### Common Name ####
    try:
        chem_name = cache.get('cid_name')[cache._checkCIDdeprecated(tmp_species, cache.get('deprecatedCID_cid'))]
    except KeyError:
        # if you cannot find using cid, try to retreive it using its inchikey
        try:
            if inchikey:
                # @Joan: Same question as above
                tmp_cids = [i for i in cache.get('inchikey_cid')[inchikey] if i[:3]=='MNX']
                if tmp_cids:
                    chem_name = cache.get('cid_name')[cache._checkCIDdeprecated(tmp_cids[0], cache.get('deprecatedCID_cid'))]
        except KeyError:
            logger.debug('Cannot find the name for this species: '+str(species))

    return {
        'xref'     : xref,
        'inchi'    : inchi,
        'inchikey' : inchikey,
        'smiles'   : smiles,
        'chem_name': chem_name
    }

## Function to reconstruct the heterologous pathway
#
#  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
#
#  @param self Object pointer
#  @param rpsbml rpSBML object with a single model
#  @return Boolean if True then you keep that model for the next step, if not then ignore it
def add_cofactors(
    cache: 'rrCache',
    rpsbml: 'rpSBML',
    compartment_id: str='MNXC3',
    pathway_id: str='rp_pathway',
    logger=logging.getLogger(__name__)
) -> 'rpSBML':
    
    # This keeps the IDs conversions to the pathway
    pathway_cmp = {}
    spe_conv = {}
    rpsbml_dict = rpsbml.toDict(pathway_id)
    # rp_path = rpsbml.convert_pathway_to_dict(pathway_id)
    ori_rpsbml_dict = deepcopy(rpsbml_dict)
    # #We reverse the loop to ID the intermediate CMP to their original ones
    # for stepNum in sorted(list(rp_path), reverse=True):

    # For each reaction id taken in forward direction
    for rxn_id in sorted(list(rpsbml_dict['reactions']), key=lambda x:rpsbml_dict['reactions'][x]['brsynth']['rxn_idx']):

        rxn     =     rpsbml_dict['reactions'][rxn_id]
        ori_rxn = ori_rpsbml_dict['reactions'][rxn_id]

        from json import dumps

        # print('COFACTORS:', dumps(cofactors, indent=4))
        # print('RXN [LEFT]:', dumps(rxn['left'], indent=4))
        # print('RXN [RIGHT]:', dumps(rxn['right'], indent=4))
        # exit()


        # rxn['left'].update(left)
        # rxn['right'].update(right)
        # rxn['smiles'] = smiles

        # if rxn['left'] or rxn['right'] or rxn['smiles']:

        ### add the new cofactors to the SBML
        # remove the original species from the monocomponent reaction
        # reactants = set( set(rxn['left'].keys())  - set(ori_rxn['left'].keys())  )
        # products  = set( set(rxn['right'].keys()) - set(ori_rxn['right'].keys()) )

        for side in ['left', 'right']:
            for species in rxn[side]:
                tmp_species = cache._checkCIDdeprecated(
                    species,
                    cache.get('deprecatedCID_cid')
                )
                # check to make sure that they do not yet exist and if not create a new one
                # TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                if not rpsbml.speciesExists(tmp_species, compartment_id):
                    spe_infos = retrieve_infos(
                        cache,
                        species,
                        tmp_species,
                        logger=logger
                    )
                    #### Finally create the species in the SBML file ######
                    rpsbml.createSpecies(
                        tmp_species,
                        compartment_id,
                        spe_infos['chem_name'],
                        spe_infos['xref'],
                        spe_infos['inchi'],
                        spe_infos['inchikey'],
                        spe_infos['smiles']
                    )

        added_species = get_cofactors_rxn(
            cache,
            rxn,
            logger=logger
        )

        for side in ['left', 'right']:
            # Compounds with structure
            for spe_name, spe_infos in added_species[rxn['brsynth']['rule_ori_reac']]['added_cmpds'][side].items():
                tmp_species = cache._checkCIDdeprecated(
                    spe_name,
                    cache.get('deprecatedCID_cid')
                )
                # check to make sure that they do not yet exist and if not create a new one
                # TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                if not rpsbml.speciesExists(tmp_species, compartment_id):
                    #### Finally create the species in the SBML file ######
                    rpsbml.createSpecies(
                        species_id=tmp_species,
                        compartment_id=compartment_id,
                        species_name=spe_infos['name'],
                        inchi=spe_infos['inchi'],
                        inchikey=spe_infos['inchikey'],
                        smiles=spe_infos['smiles']
                    )
            # Compounds with no structure
            for spe_name, spe_infos in added_species[rxn['brsynth']['rule_ori_reac']]['added_cmpds'][side+'_nostruct'].items():
                tmp_species = cache._checkCIDdeprecated(
                    spe_name,
                    cache.get('deprecatedCID_cid')
                )
                # check to make sure that they do not yet exist and if not create a new one
                # TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                if not rpsbml.speciesExists(tmp_species, compartment_id):
                    #### Finally create the species in the SBML file ######
                    rpsbml.createSpecies(
                        species_id=tmp_species,
                        compartment_id=compartment_id
                    )

        # add the new species to the RP reactions
        rxn_rpsbml = rpsbml.getModel().getReaction(rxn_id)
        for side in ['left', 'right']:
            for spe_name, spe_infos in {
                **added_species[rxn['brsynth']['rule_ori_reac']]['added_cmpds'][side],
                **added_species[rxn['brsynth']['rule_ori_reac']]['added_cmpds'][side+'_nostruct']
            }.items():
                species = cache._checkCIDdeprecated(
                    spe_name,
                    cache.get('deprecatedCID_cid')
                )
                # if not species.endswith('__64__'+str(compartment_id)):
                #     species += '__64__'+str(compartment_id)
                if side == 'left':
                    spe = rxn_rpsbml.createProduct()
                else:
                    spe = rxn_rpsbml.createReactant()
                spe.setSpecies(species)
                spe.setConstant(True)
                spe.setStoichiometry(spe_infos['stoichio'])

        # replace the reaction rule with new one
        rpsbml.updateBRSynth(rxn_rpsbml, 'smiles', added_species[rxn['brsynth']['rule_ori_reac']]['full_transfo'], None, True)

        # else:
        #     # if the cofactors cannot be found delete it from the list
        #     logger.warning('Cannot find cofactors... skipping')
        #     return False

    return rpsbml