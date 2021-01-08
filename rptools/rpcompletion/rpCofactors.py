import logging
from copy import deepcopy



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
def completeReac(cache, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp, logger=None):

    logger = logger or logging.getLogger(__name__)

    if mono_side:
        ## add the unknown species to pathway_cmp for the next steps
        rr_mono_cmp = list(rr_reac.keys())
        step_mono_cmp = list(step.keys())
        if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
            #this is purposely overwitten since the main cmp between reactions can change
            pathway_cmp[step_mono_cmp[0]] = rr_mono_cmp[0]
        else:
            logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
            return False

    ## add the side species
    rr_string += add_side_species(cache, step, full_reac, rr_reac, logger=logger)

    ## Update the stochio
    return True, update_stochio(cache, step, full_reac, rr_string, pathway_cmp, logger=logger)


def add_side_species(cache, step, full_reac, rr_reac, logger=None):
    logger = logger or logging.getLogger(__name__)
    rr_string = ''
    for toAdd in full_reac.keys()-rr_reac.keys():
        step.update({toAdd: full_reac[toAdd]})
        ### update the reaction rule string
        try:
            smi = cache.cid_strc[toAdd]['smiles']
            if not smi==None:
                for sto_add in range(int(full_reac[toAdd])):
                    rr_string += '.'+str(smi)
        except KeyError:
            logger.warning('Cannot find smiles structure for '+str(toAdd))
    return rr_string


def update_stochio(cache, step, full_reac, rr_string, pathway_cmp, logger=None):
    logger = logger or logging.getLogger(__name__)
    for step_spe in step:
        if step_spe in full_reac:
            if not step[step_spe]==full_reac[step_spe]:
                stochio_diff = full_reac[step_spe]-step[step_spe]
                step[step_spe] = full_reac[step_spe]
                if stochio_diff<0:
                    logger.warning('full_reac stochio should never be smaller than step')
                    continue
                for i in range(stochio_diff):
                    ### update the reaction rule string
                    try:
                        smi = cache.cid_strc[step_spe]['smiles']
                        if not smi==None:
                            rr_string += '.'+str(smi)
                    except KeyError:
                        #@Mel toAdd -> step_spe
                        logger.warning('Cannot find smiles structure for '+str(step_spe))
        elif step_spe in pathway_cmp:
            if pathway_cmp[step_spe] in full_reac:
                if not step[step_spe]==full_reac[pathway_cmp[step_spe]]:
                    step[step_spe] = full_reac[pathway_cmp[step_spe]]
        #Its fine if the stochio is not updated, better than ignoring a whole pathway
            #else:
            #    logger.warning('Cannot find '+str(step_spe)+' in full reaction')
            #    return False
        #else:
        #    logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp')
        #    return False
        return rr_string


## Add the cofactors to monocomponent reactions
#
# @param step Step in a pathway
# @param pathway_cmp Dictionnary of intermediate compounds with their public ID's
# @return Boolean determine if the step is to be added
def addCofactors_step(cache, step, pathway_cmp, logger=None):
    logger = logger or logging.getLogger(__name__)
    reac_smiles_left = step['reaction_rule'].split('>>')[0]
    reac_smiles_right = step['reaction_rule'].split('>>')[1]
    if cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==-1:
        try:
            isSuccess, reac_smiles_left = completeReac(cache, step['right'],
                    cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                    cache.rr_full_reactions[cache._checkRIDdeprecated(step['rule_ori_reac'], cache.deprecatedRID_rid)]['right'],
                    True,
                    reac_smiles_left,
                    pathway_cmp,
                    logger=logger)
            if not isSuccess:
                logger.warning('Could not recognise reaction rule for step (1): '+str(step))
                return False
        except KeyError:
            logger.warning('Could not find the full reaction for reaction (1): '+str(step))
            return False
        try:
            isSuccess, reac_smiles_right = completeReac(cache, step['left'],
                    cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                    cache.rr_full_reactions[cache._checkRIDdeprecated(step['rule_ori_reac'], cache.deprecatedRID_rid)]['left'],
                    False,
                    reac_smiles_right,
                    pathway_cmp,
                    logger=logger)
            if not isSuccess:
                logger.warning('Could not recognise reaction rule for step (2): '+str(step))
                return False
        except KeyError:
            logger.warning('Could not find the full reaction for reaction (2): '+str(step))
            return False
    elif cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==1:
        try:
            isSuccess, reac_smiles_left = completeReac(cache, step['right'],
                    cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                    cache.rr_full_reactions[cache._checkRIDdeprecated(step['rule_ori_reac'], cache.deprecatedRID_rid)]['left'],
                    True,
                    reac_smiles_left,
                    pathway_cmp,
                    logger=logger)
            if not isSuccess:
                logger.error('Could not recognise reaction rule for step (3): '+str(step))
                return False
        except KeyError:
            logger.warning('Could not find the full reaction for reaction (3): '+str(step))
            return False
        try:
            isSuccess, reac_smiles_right = completeReac(cache, step['left'],
                    cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                    cache.rr_full_reactions[cache._checkRIDdeprecated(step['rule_ori_reac'], cache.deprecatedRID_rid)]['right'],
                    False,
                    reac_smiles_right,
                    pathway_cmp,
                    logger=logger)
            if not isSuccess:
                logger.error('Could not recognise reaction rule for step (4): '+str(step))
                return False
        except KeyError:
            logger.warning('Could not find the full reaction for reaction (4): '+str(step))
            return False
    else:
        logger.error('Relative direction can only be 1 or -1: '+str(cache.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']))
        return False
    step['reaction_rule'] = reac_smiles_left+'>>'+reac_smiles_right
    return True


## Function to reconstruct the heterologous pathway
#
#  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
#
#  @param self Object pointer
#  @param rpsbml rpSBML object with a single model
#  @return Boolean if True then you keep that model for the next step, if not then ignore it
def addCofactors(cache, rpsbml, compartment_id='MNXC3', pathway_id='rp_pathway', logger=None):
    logger = logger or logging.getLogger(__name__)
    #This keeps the IDs conversions to the pathway
    pathway_cmp = {}
    spe_conv = {}
    rpsbml_json = rpsbml.genJSON(pathway_id)
    rp_path = rpsbml.convert_pathways_to_dict(pathway_id)
    ori_rp_path = deepcopy(rp_path)
    #We reverse the loop to ID the intermediate CMP to their original ones
    for stepNum in sorted(list(rp_path), reverse=True):
    #for stepNum in sorted(list(rp_path)):
        if addCofactors_step(cache, rp_path[stepNum], pathway_cmp, logger=logger):
            ###add the new cofactors to the SBML
            #remove the original species from the monocomponent reaction
            reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
            products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
            for species in reactants|products:
                tmp_species = cache._checkCIDdeprecated(species, cache.deprecatedCID_cid)
                #check to make sure that they do not yet exist and if not create a new one
                #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                if not rpsbml.speciesExists(tmp_species, compartment_id):
                    xref = {}
                    inchi = None
                    inchikey = None
                    smiles = None
                    chem_name = None
                    ###### Try to retreive the InChI ############
                    try:
                        inchi = cache.cid_strc[tmp_species]['inchi']
                    except KeyError:
                        logger.warning('Cannot find the inchi for this species: '+str(tmp_species))
                    try:
                        inchikey = cache.cid_strc[tmp_species]['inchikey']
                        #logger.debug('Found the inchikey: '+str(inchikey))
                        #### TODO: find a better way to check if two species are the same ####
                        isfound = False
                        for rpsbml_species in rpsbml_json['species']:
                            #TODO add a comparison by xref as well
                            #logger.debug(str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])+' <--> '+str(inchikey))
                            if str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])==str(inchikey):
                                spe_conv[tmp_species] = rpsbml_species
                                logger.debug('The species '+str(tmp_species)+' is the same as '+str(rpsbml_species))
                                isfound = True
                                break
                        if isfound:
                            continue
                    except KeyError:
                        logger.warning('Cannot find the inchikey for this species: '+str(species))
                    ##### Try to retreive the SMILES ############
                    try:
                        smiles = cache.cid_strc[tmp_species]['smiles']
                    except KeyError:
                        logger.warning('Cannot find the smiles for this species: '+str(species))
                    ###### Try to retreive the xref, using the inchikey if the cid fails #######
                    try:
                        xref = cache.cid_xref[tmp_species]
                    except KeyError:
                        try:
                            xref = cache.cid_xref[tmp_species]
                        except KeyError:
                            #if you cannot find using cid, try to retreive it using its inchikey
                            try:
                                if inchikey:
                                    #@Joan: Can you think of a better way of doing that?
                                    # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                                    tmp_cids = [i for i in cache.inchikey_cid[inchikey] if i[:3]=='MNX']
                                    #TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                                    if tmp_cids:
                                        xref = cache.cid_xref[cache._checkCIDdeprecated(tmp_cids[0], cache.deprecatedCID_cid)]
                            except KeyError:
                                logger.warning('Cannot find the xref for this species: '+str(species))
                                xref = {}
                    #### Common Name ####
                    try:
                        chem_name = cache.cid_name[cache._checkCIDdeprecated(tmp_species, cache.deprecatedCID_cid)]
                    except KeyError:
                        #if you cannot find using cid, try to retreive it using its inchikey
                        try:
                            if inchikey:
                                #@Joan: Same question as above
                                tmp_cids = [i for i in cache.inchikey_cid[inchikey] if i[:3]=='MNX']
                                if tmp_cids:
                                    chem_name = cache.cid_name[cache._checkCIDdeprecated(tmp_cids[0], cache.deprecatedCID_cid)]
                        except KeyError:
                            logger.warning('Cannot find the name for this species: '+str(species))
                    #### Finally create the species in the SBML file ######
                    rpsbml.createSpecies(tmp_species,
                            compartment_id,
                            chem_name,
                            xref,
                            inchi,
                            inchikey,
                            smiles)
            #add the new species to the RP reactions
            reac = rpsbml.getModel().getReaction(rp_path[stepNum]['reaction_id'])
            pre_reactants = [i.species for i in reac.getListOfReactants()]
            pre_products = [i.species for i in reac.getListOfProducts()]
            for pro in products:
                if cache._checkCIDdeprecated(pro, cache.deprecatedCID_cid) in spe_conv:
                    toadd = spe_conv[cache._checkCIDdeprecated(pro, cache.deprecatedCID_cid)]
                else:
                    toadd = str(cache._checkCIDdeprecated(pro, cache.deprecatedCID_cid))+'__64__'+str(compartment_id)
                #prod.setSpecies(str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id))
                if toadd in pre_products:
                    continue
                prod = reac.createProduct()
                prod.setSpecies(toadd)
                prod.setConstant(True)
                prod.setStoichiometry(rp_path[stepNum]['right'][pro])
            for sub in reactants:
                if cache._checkCIDdeprecated(sub, cache.deprecatedCID_cid) in spe_conv:
                    toadd = spe_conv[cache._checkCIDdeprecated(sub, cache.deprecatedCID_cid)]
                else:
                    toadd = str(cache._checkCIDdeprecated(sub, cache.deprecatedCID_cid))+'__64__'+str(compartment_id)
                #prod.setSpecies(str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id))
                if toadd in pre_reactants:
                    continue
                subs = reac.createReactant()
                subs.setSpecies(toadd)
                subs.setConstant(True)
                subs.setStoichiometry(rp_path[stepNum]['left'][sub])
            #replace the reaction rule with new one
            rpsbml.updateBRSynth(reac, 'smiles', rp_path[stepNum]['reaction_rule'], None, True)
        else:
            #if the cofactors cannot be found delete it from the list
            logger.warning('Cannot find cofactors... skipping')
            return False
    return True
