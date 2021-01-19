from csv       import reader         as csv_reader
from requests  import post           as r_post
from requests  import get            as r_get
from itertools import product        as itertools_product
from time      import time           as time_time
from time      import sleep          as time_sleep
from argparse  import ArgumentParser as argparse_ArgumentParser
from os        import path           as os_path
from os        import mkdir          as os_mkdir
from json      import decoder        as json_decoder
from io        import StringIO
from copy      import deepcopy
from brs_utils import insert_and_or_replace_in_sorted_list
from rptools.rplibs                   import rpSBML
from rptools.rpcompletion.rpCofactors import add_cofactors
import logging
import pandas as pd


class Species:
 
 
    def __init__(self, inchi, inchikey, smiles, xref):
        self.inchi = inchi
        self.inchikey = inchikey
        self.smiles = smiles
        self.xref = xref

class SBML_Item:


    def __init__(self, score, index, rpsbml_obj):
        self.score = score
        self.index = index
        self.rpsbml_obj = rpsbml_obj


    def __eq__(self, sbml_item):
        return self.rpsbml_obj == sbml_item.rpsbml_obj


    def __lt__(self, sbml_item):
        return self.score < sbml_item.score


    def __gt__(self, sbml_item):
        return self.score > sbml_item.score


    def __str__(self):
        return 'SBML_Item' + '\n' \
             + '\t' + 'score:      ' + str(self.score)      + '\n' \
             + '\t' + 'index:      ' + str(self.index)      + '\n' \
             + '\t' + 'rpsbml_obj: ' + str(self.rpsbml_obj) + '\n' \


#####################
pubchem_min_count = 0
pubchem_min_start = 0.0


#######################################################################
############################# PRIVATE FUNCTIONS #######################
#######################################################################

def _pubChemLimit(logger=logging.getLogger(__name__)):
    global pubchem_min_count, pubchem_min_start

    if pubchem_min_start==0.0:
        pubchem_min_start = time_time()
    pubchem_min_count += 1
    #### requests per minute ####
    if pubchem_min_count>=500 and time_time()-pubchem_min_start<=60.0:
        logger.warning('Reached 500 requests per minute for pubchem... waiting a minute')
        time_sleep(60.0)
        pubchem_min_start = time_time()
        pubchem_min_count = 0
    elif time_time()-pubchem_min_start>60.0:
        pubchem_min_start = time_time()
        pubchem_min_count = 0


## Try to retreive the xref from an inchi structure using pubchem
#
# No more than 5 requests per second. No more than 400 requests per minute. No longer than 300 second running time per minute. Requests exceeding limits are rejected (HTTP 503 error)
# @param self The object pointer
# @param strct Strucutre string of the molecule
# @param itype Input type of the structure. Accepted values are: inchi, inchikey, smiles
# @return Dict with the following structure: {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}
def _pubchemStrctSearch(strct, itype='inchi', logger=logging.getLogger(__name__)):

    _pubChemLimit(logger=logger)
    try:
        r = r_post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
        res_list = r.json()
    except json_decoder.JSONDecodeError:
        logger.warning('JSON decode error')
        return {}
    try:
        res_list = res_list['InformationList']['Information']
    except KeyError:
        logger.warning('pubchem JSON keyerror: '+str(res_list))
        return {}
    xref = {}
    if len(res_list)==1:
        _pubChemLimit(logger=logger)
        try:
            prop = r_get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON')
            prop_list = prop.json()
        except json_decoder.JSONDecodeError:
            logger.warning('JSON decode error')
            return {}
        try:
            name = prop_list['PropertyTable']['Properties'][0]['IUPACName']
            inchi = prop_list['PropertyTable']['Properties'][0]['InChI']
            inchikey = prop_list['PropertyTable']['Properties'][0]['InChIKey']
            smiles = prop_list['PropertyTable']['Properties'][0]['CanonicalSMILES']
        except KeyError:
            logger.warning('pubchem JSON keyerror: '+str(prop_list))
            return {}
        #TODO: need to determine how long cobra cannot handle this
        #TODO: determine if names that are too long is the problem and if not remove this part
        if len(name)>30:
            _pubChemLimit(logger=logger)
            try:
                syn = r_get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/synonyms/JSON')
                syn_lst = syn.json()
            except json_decoder.JSONDecodeError:
                logger.warning('pubchem JSON decode error')
                return {}
            try:
                syn_lst = syn_lst['InformationList']['Information'][0]['Synonym']
                syn_lst = [x for x in syn_lst if not 'CHEBI' in x and not x.isupper()]
                name = syn_lst[0] #need a better way instead of just the firs tone
            except KeyError:
                logger.warning('pubchem JSON keyerror: '+str(syn.json()))
                return {}
            except IndexError:
                name = ''
        xref['pubchem'] = [str(res_list[0]['CID'])]
        for url in res_list[0]['SBURL']:
            if 'https://biocyc.org/compound?orgid=META&id=' in url:
                if 'biocyc' not in xref:
                    xref['biocyc'] = []
                xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
            if 'http://www.hmdb.ca/metabolites/' in url:
                if 'hmdb' not in xref:
                    xref['hmdb'] = []
                xref['hmdb'].append(url.replace('http://www.hmdb.ca/metabolites/', ''))
            if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                if 'kegg_c' not in xref:
                    xref['kegg_c'] = []
                xref['kegg_c'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
            if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                if 'chebi' not in xref:
                    xref['chebi'] = []
                xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
    elif len(res_list)==0:
        logger.warning('Could not find results for: '+str(strct))
        return {}
    else:
        logger.warning('There are more than one result for '+str(strct)+'... Ignoring')
        return {}
    return {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}


###############################################################
############################ RP2paths entry functions #########
###############################################################

## Function to group all the functions for parsing RP2 output to rpSBML files
#
# Takes RP2paths's compounds.txt and out_paths.csv and RetroPaths's *_scope.csv files and generates rpSBML files, where each individual file contains a single heterologous file
#
# @param rp2_pathways Path (string) to the output file of RetroPath2.0
# @param rp2paths_compounds Path (string) to the RP2paths compounds file
# @param rp2paths_pathways Path (string) to the RP2paths pathways file
# @param outdir folder where to write files
# @param upper_flux_bound Upper flux bound for all reactions that will be created (default: 999999)
# @param lower_flux_bound Lower flux bound for all reactions that will be created (default: 0)
# @param max_subpaths_filter The maximal number of subpaths per path (default: 10)
# @param pathway_id Groups id that will contain all the heterologous reactions
# @param compartment_id The Groups id of the SBML's model compartment where to add the reactions
# @param species_group_id The Groups id of the central species of the heterologous pathway
# @param species_group_id The Groups id of the sink species of the heterologous pathway
# @return Boolean The success or failure of the function
def rp_completion(cache,
                  rp2_pathways,
                  rp2paths_compounds,
                  rp2paths_pathways,
                  outdir,
                  upper_flux_bound=999999,
                  lower_flux_bound=0,
                  max_subpaths_filter=10,
                  pathway_id='rp_pathway',
                  compartment_id='MNXC3',
                  species_group_id='central_species',
                  sink_species_group_id='rp_sink_species',
                  pubchem_search=False,
                  logger=logging.getLogger(__name__)):

    if max_subpaths_filter<0:
        raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

    if not os_path.exists(outdir):
        os_mkdir(outdir)
    elif os_path.isfile(outdir):
        logger.error('Outdir name '+outdir+' already exists and is actually file. Stopping the process...')
        exit()

    rp_strc = _compounds(cache, rp2paths_compounds, logger=logger)
    rp_transformation, sink_molecules = _transformation(rp2_pathways, logger=logger)

    return write_rp2paths_to_rpSBML(cache,
                                    rp_strc,
                                    rp_transformation,
                                    sink_molecules,
                                    rp2paths_pathways,
                                    outdir,
                                    upper_flux_bound,
                                    lower_flux_bound,
                                    max_subpaths_filter,
                                    pathway_id,
                                    compartment_id,
                                    species_group_id,
                                    sink_species_group_id,
                                    pubchem_search,
                                    logger=logger)


## Function to parse the compounds.txt file
#
#  Extract the smile and the structure of each compounds of RP2Path output
#  Method to parse all the RP output compounds.
#
#  @param path The compounds.txt file path
#  @return rp_compounds Dictionnary of smile and structure for each compound
def _compounds(cache, path, logger=logging.getLogger(__name__)):

    #self.rp_strc = {}
    rp_strc = {}
    try:
        if isinstance(path, bytes):
            reader = csv_reader(StringIO(path.decode('utf-8')), delimiter='\t')
        else:
            reader = csv_reader(open(path, 'r', encoding='utf-8'), delimiter='\t')
        next(reader)
        for row in reader:
            rp_strc[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
            try:
                rp_strc[row[0]]['inchi'] = cache.cid_strc[row[0]]['inchi']
            except KeyError:
                #try to generate them yourself by converting them directly
                try:
                    resConv = cache._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                    rp_strc[row[0]]['inchi'] = resConv['inchi']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
            try:
                rp_strc[row[0]]['inchikey'] = cache.cid_strc[row[0]]['inchikey']
                #try to generate them yourself by converting them directly
                #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
            except KeyError:
                try:
                    resConv = cache._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                    rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
    except (TypeError, FileNotFoundError) as e:
        logger.error('Could not read the compounds file ('+str(path)+')')
        raise RuntimeError
    return rp_strc


## Function to parse the scope.csv file
#
# Extract the reaction rules from the retroPath2.0 output using the scope.csv file
#
# @param path The scope.csv file path
# @return tuple with dictionnary of all the reactions rules and the list unique molecules that these apply them to
def _transformation(path, logger=logging.getLogger(__name__)):

    rp_transformation = {}
    sink_molecules = []
    #### we might pass binary in the REST version
    reader = None
    if isinstance(path, bytes):
        reader = csv_reader(StringIO(path.decode('utf-8')), delimiter=',')
    else:
        try:
            reader = csv_reader(open(path, 'r'), delimiter=',')
        except FileNotFoundError:
            logger.error('Could not read the compounds file: '+str(path))
            return {}
    next(reader)
    for row in reader:
        if not row[1] in rp_transformation:
            rp_transformation[row[1]] = {}
            rp_transformation[row[1]]['rule'] = row[2]
            rp_transformation[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
        if row[7]=='1':
            for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                sink_molecules.append(i)

    # logger.info(rp_transformation)
    # logger.info(sink_molecules)
    return rp_transformation, list(set(sink_molecules))


## Function that converts pathways from output of rp2paths to a dictionnary with the following structure:
# {
#     path_base_idx: {
#         path_step: {
#             'left': { components },
#             'right': { products },
#             'transformation_id': transformation_id,
#             'reactions': {
#                 rxn_id: {
#                     'rule_id': rule_id,
#                     'rule_score': rule_score
#                 }
#             }
#         }
#     }
# }
#
# @param infile rp2_pathways file
# @param rr_reactions RetroRules reactions
# @param deprecatedCID_cid Deprecated cid
# @return Dictionnary of the pathways
def rp2paths_to_dict(infile, rr_reactions, deprecatedCID_cid, logger=logging.getLogger(__name__)):

    df = pd.read_csv(infile)

    check = check_pathways(df)
    if not check:
        logger.error(check)
        exit()

    rp2paths_pathways = {}

    # store new Path ID
    path_base_idx = None

    for index, row in df.iterrows():

        ### INIT ###
        # if 'Path ID' has changed
        if row['Path ID'] != path_base_idx:
            # store new Path ID
            path_base_idx = row['Path ID']
            # create new entry in rp2paths_pathways
            rp2paths_pathways[path_base_idx] = {}
            # re-init path step (rxn number)
            path_step = 1
        else: # else increment path step (rxn number)
            path_step += 1

        ### REACTION RULE COMMON PART ###
        # create new entry in rp2paths_pathways
        rp2paths_pathways[path_base_idx][path_step] = {}
        # Common fields to all reactions for current step (all reactions from all rules at a given step)
        rp2paths_pathways[path_base_idx][path_step] = {
                                   # remove all illegal characters in SBML ids
            'left': build_side_rxn(row['Left'].replace("'", "").replace('-', '_').replace('+', ''),
                                   deprecatedCID_cid),
            'right': build_side_rxn(row['Right'], deprecatedCID_cid),
            'transformation_id': row['Unique ID'][:-2],
            'reactions': {}
        }

        ### REACTION SPECIFIC PART ###
        ## Generate all reactions from each rule list
        # For each rule
        for rule_id in row['Rule ID'].split(','):
            # For each reaction of a rule
            for rxn_id in rr_reactions[rule_id]:
                # Add reaction to rp2paths_pathways
                rp2paths_pathways[path_base_idx][path_step]['reactions'][rxn_id] = {
                    'rule_id'       : rule_id,
                    'rule_score'    : rr_reactions[rule_id][rxn_id]['rule_score']
                }

    return rp2paths_pathways

        # # if 'Path ID' has changed
        # if row['Path ID'] != current_path_base_idx:
        #     # re-init index of pathway variant
        #     path_variant_idx = 1
        #     # store new Path ID
        #     current_path_base_idx = row['Path ID']
        # else: # else increment index of pathway variant
        #     path_variant_idx += 1


        # rule_ids_dict = rxns_from_rules(row['Rule ID'].split(','), rr_reactions)
        # path_variant_idx = 1

        # for rule_id in rule_ids_dict:
        #     ## Build reaction
        #     rxn = {
        #             'rule_id'           : rule_id.split('__')[0],
        #             'rule_ori_reac'     : rule_id.split('__')[1],
        #             'rule_score'        : rr_reactions[rule_id.split('__')[0]][rule_id.split('__')[1]]['rule_score'],
        #                                                  # remove all illegal characters in SBML ids
        #             'left'              : build_side_rxn(row['Left'].replace("'", "").replace('-', '_').replace('+', ''),
        #                                                  deprecatedCID_cid),
        #             'right'             : build_side_rxn(row['Right'], deprecatedCID_cid),
        #             'rxn_idx'              : path_step,
        #             'transformation_id' : row['Unique ID'][:-2]
        #             }
        #     rp2paths_pathways = add_rxn_to_rppaths(rp2paths_pathways, row['Path ID'], path_variant_idx, rxn)
        #     path_variant_idx += 1


## Function that checks pathways data
#
# @param df pathways read with pandas
# @return True or a message if something went wrong
def check_pathways(df):
    if len(df) == 0:
        return 'infile is empty'

    if df['Path ID'].dtypes != 'int64':
        return '\'Path ID\' column contain non integer value(s)'

    return True


## Function that returns the left/right part of reaction as a dict
#
# @param side is left or right part of the reaction read in pathways
# @param deprecatedCID_cid
# @return Dictionnary of the left part of the reaction
def build_side_rxn(side, deprecatedCID_cid, logger=logging.getLogger(__name__)):

    logger=logging.getLogger(__name__)
    rxn_side = {}

    for s in side.split(':'):
 
        tmp = s.split('.')
 
        try:
            #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
            cid = '' # TODO: change this
            if tmp[1] in deprecatedCID_cid:
                cid = deprecatedCID_cid[tmp[1]]
            else:
                cid = tmp[1]
            rxn_side[cid] = int(tmp[0])
 
        except ValueError:
            logger.error('Cannot convert tmp[0] to int ('+str(tmp[0])+')')
            #return {}
            raise

    return rxn_side


## Function to parse the out_paths.csv file
#
#  Reading the RP2path output and extract all the information for each pathway
#  RP2path Metabolic pathways from out_paths.csv
#  create all the different values for heterologous paths from the RP2path out_paths.csv file
#  Note that path_step are in reverse order here
#
#  @param self Object pointer
#  @param path The out_path.csv file path
#  @max_subpaths_filter maximal numer of subpaths per paths
#  @outFolder folder where to write files
#  @return Boolean The success or failure of the function
def write_rp2paths_to_rpSBML(cache,
                             rp_strc, rp_transformation, sink_molecules,
                             rp2paths_pathways,
                             outFolder,
                             upper_flux_bound=999999,
                             lower_flux_bound=0,
                             max_subpaths_filter=10,
                             pathway_id='rp_pathway',
                             compartment_id='MNXC3',
                             species_group_id='central_species',
                             sink_species_group_id='rp_sink_species',
                             pubchem_search=False,
                             logger=logging.getLogger(__name__)):
    # TODO: make sure that you account for the fact that each reaction may have multiple associated reactions

    rp2paths_pathways = rp2paths_to_dict(rp2paths_pathways, cache.rr_reactions, cache.deprecatedCID_cid, logger=logger)

    # for each line or rp2paths_pathways:
    #     generate comb
    #     for each combinant:
    #         rank
    #         process
    #         add cofactors
    #         dedup

    for path_base_idx in rp2paths_pathways:

        # topX pathways
        best_rpsbml = []

        # Build all pathways from reaction list over step
        pathways_comb = list(map(list,list(itertools_product(*[list(rp2paths_pathways[path_base_idx][step]['reactions'].keys()) for step in rp2paths_pathways[path_base_idx]]))))

        path_variant_idx = 1
        for path_variant in pathways_comb:

            path_id_idx = str(path_base_idx)+'_'+str(path_variant_idx)
            path_id     = 'rp_'+path_id_idx

            # 1) Create an rpSBML object with species
            rpsbml, species = create_rpSBML(pathway_id, path_base_idx, path_variant_idx,
                                            path_id, path_id_idx,
                                            rp2paths_pathways[path_base_idx],
                                            cache, compartment_id, upper_flux_bound, lower_flux_bound,
                                            rp_strc, rp_transformation, sink_molecules,
                                            species_group_id, sink_species_group_id,
                                            pubchem_search)

            # 2) Add complete reactions
            rpsbml = complete_reactions(rpsbml, species, rp2paths_pathways[path_base_idx],
                                        pathway_id, path_variant,
                                        compartment_id, upper_flux_bound, lower_flux_bound,
                                        rp_transformation)

            # 3) Add the cofactors
            rpsbml = add_cofactors(cache, rpsbml, logger=logger)

            # 4) Apply to best rpsbml list
            best_rpsbml = apply_best_rpsbml(best_rpsbml, max_subpaths_filter,
                                            rpsbml, path_id)

            path_variant_idx += 1

        # Write results to files
        for rpsbml_item in best_rpsbml:
            rpsbml_item.rpsbml_obj.writeSBML(os_path.join(outFolder, str(rpsbml_item.rpsbml_obj.modelName))+'_sbml.xml')

    return 0


def apply_best_rpsbml(best_rpsbml, max_subpaths_filter,
                      rpsbml, path_id):

    # Insert the new rpsbml object in sorted list
    sbml_item = SBML_Item(rpsbml.compute_score(), path_id, rpsbml)
    best_rpsbml = insert_and_or_replace_in_sorted_list(sbml_item, best_rpsbml)

    # Keep only topX
    best_rpsbml = best_rpsbml[-max_subpaths_filter:]

    return best_rpsbml


def create_rpSBML(pathway_id, path_base_idx, path_variant_idx,
                  path_id, path_id_idx,
                  steps,
                  cache, compartment_id, upper_flux_bound, lower_flux_bound,
                  rp_strc, rp_transformation, sink_molecules,
                  species_group_id, sink_species_group_id,
                  pubchem_search,
                  logger=logging.getLogger(__name__)):

    rpsbml = rpSBML(name=path_id, logger=logger)

    # 1) Create a generic Model, ie the structure and unit definitions that we will use the most
    ##### TODO: give the user more control over a generic model creation:
    # -> special attention to the compartment
    rpsbml.genericModel(
            'RetroPath_Pathway_'+path_id_idx,
            'RP_model_'+path_id_idx,
            cache.comp_xref[cache.deprecatedCompID_compid[compartment_id]],
            compartment_id,
            upper_flux_bound,
            lower_flux_bound)

    # 2) Create the pathway, species, sink species (groups)
    rpsbml = create_groups(rpsbml,
                           pathway_id, path_id, path_base_idx, path_variant_idx,
                           species_group_id, sink_species_group_id,
                           logger)

    # 3) List species and add them to the model
    # List of species of the current pathway variant
    species = get_species(steps)
    rpsbml = add_unique_species(rpsbml, species,
                                rp_strc, sink_molecules, compartment_id, species_group_id, sink_species_group_id, pubchem_search,
                                cache,
                                logger=logger)
    
    return rpsbml, species


def complete_reactions(rpsbml, species, steps,
                       pathway_id, path_variant,
                       compartment_id, upper_flux_bound, lower_flux_bound,
                       rp_transformation):

    # 4) Add the complete reactions of the pathway (all steps) and their annotations
    rpsbml = add_reactions(path_variant,
                           steps,
                           upper_flux_bound, lower_flux_bound,
                           rp_transformation,
                           compartment_id, pathway_id,
                           rpsbml)

    # 5) Adding the consumption of the target
    rxn_target = build_rxn(left={[i for i in species if i.startswith('TARGET')][0]: 1},
                           right=[])
    rpsbml.createReaction('rxn_target',
                            upper_flux_bound, lower_flux_bound,
                            rxn_target,
                            compartment_id)
    
    return rpsbml


## Function that returns rpSBML object with groups (pathway, species, sink species) added
#
#  @param rpsbml rpSBML The in-building rpSBML object
#  @param pathway_id Str The id of the pathway to add reactions to
#  @param path_id Str The name of the pathway to add reactions to
#  @param path_base_idx Int The index of base pathway (rules based)
#  @param path_variant_idx Int The index of variant pathway (reactions based)
#  @param species_group_id Str The Groups id to add the species
#  @param sink_group_id Str The Groups id sink species to add the species
#  @return rpSBML The updated rpSBML object
def create_groups(rpsbml,
                  pathway_id, path_id, path_base_idx, path_variant_idx,
                  species_group_id, sink_species_group_id,
                  logger=logging.getLogger(__name__)):

    logger.debug('Create pathway group: '+pathway_id)

    # Create new group
    rpsbml.createGroup(pathway_id)

    # Add pathway id
    rpsbml_dict = rpsbml.toDict(pathway_id)
    # path id
    rpsbml_dict['pathway']['brsynth']['path_id']          = path_id
    # base path (from reaction rules) index
    rpsbml_dict['pathway']['brsynth']['path_base_idx']    = path_base_idx
    # variant path (from real reactions) index
    rpsbml_dict['pathway']['brsynth']['path_variant_idx'] = path_variant_idx
    # list of reaction ids in the pathway, in forward direction
    rpsbml_dict['pathway']['brsynth']['reactions']        = []

    # Update rpSBML object
    rpsbml.updateBRSynthPathway(rpsbml_dict, pathway_id)
    logger.debug('Create species group: '+species_group_id)
    rpsbml.createGroup(species_group_id)
    logger.debug('Create sink species group: '+sink_species_group_id)
    rpsbml.createGroup(sink_species_group_id)

    return rpsbml


## Function that returns the dictionnary of elements passed
#
#  @return Dict
def build_rxn(rule_id=None, rule_ori_reac=None, rule_score=None, left=None, right=None, index=None):
    return {
        'rule_id':       rule_id,
        'rule_ori_reac': rule_ori_reac,
        'rule_score':    rule_score,
        'left':          left,
        'right':         right,
        'rxn_idx':       index
        }


## Function that returns the list of species in rp2paths_pathways
#
#  @param path_variant List The list of reactions in the pathway variant
#  @param steps Dict The sub-dictionnary of steps of a base pathway
#  @param upper_flux_bound Int The reaction fbc upper bound
#  @param lower_flux_bound Int The reaction fbc lower bound
#  @param lower_flux_bound List Unique molecules that reactions apply them to
#  @param compartment_id Str The id of the compartment to add the reaction
#  @param pathway_id Str The id of the pathway to add reactions to
#  @param rpsbml rpSBML The in-building rpSBML object
#  @return rpSBML The updated rpSBML object
def add_reactions(path_variant, steps,
                  upper_flux_bound, lower_flux_bound,
                  rp_transformation,
                  compartment_id, pathway_id,
                  rpsbml):

    for step in range(len(path_variant)):

        rule_ori_reac = path_variant[step]
        index         = len(path_variant)-int(step) # switch rxn number order from reverse to forward
        rxn_id        = 'rxn_'+str(index)

        # Build the complete reaction
        rxn = build_rxn(
            steps[step+1]['reactions'][rule_ori_reac]['rule_id'],
            rule_ori_reac,
            steps[step+1]['reactions'][rule_ori_reac]['rule_score'],
            steps[step+1]['left'],
            steps[step+1]['right'],
            index
        )

        transformation_id = steps[step+1]['transformation_id']
        # Create the reaction in the model
        rpsbml.createReaction(
                rxn_id, # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['rxn_idx']),
                upper_flux_bound, lower_flux_bound,
                rxn,
                compartment_id,
                rp_transformation[transformation_id]['rule'],
                {'ec': rp_transformation[transformation_id]['ec']},
                pathway_id)

    return rpsbml


## Function that returns the list of species in rp2paths_pathways
#
#  @param steps Dictionnary of steps in a given pathway
#  @return List The list of species
def get_species(steps):

    species = set()

    for step in steps:
        species.update(list(steps[step]['left'].keys()))
        species.update(list(steps[step]['right'].keys()))
    
    return species


def unique_species(cache, meta, rp_strc, pubchem_search, logger=logging.getLogger(__name__)):

    try:
        chemName = cache.cid_strc[meta]['name']
    except KeyError:
        chemName = None

    # compile as much info as you can

    # xref
    try: xref = cache.cid_xref[meta]
    except KeyError: xref = {}
    spe = Species(None, None, None, xref)

    ###### Try to recover the structures ####
    pubchem = Species(None, None, None, {})

    # inchi
    try:
        # @Joan: Do you make sure here that the spe.inchi is not None?
        spe.inchi = rp_strc[meta]['inchi']
        if not spe.xref and pubchem_search:
            try:
                # print()
                # print("*************")
                # print("pubchem_species READ (INCHI)", spe.inchi)
                # print("*************")
                # print()
                pubchem.inchi    = cache._pubchem_species[spe.inchi]['inchi']
                pubchem.inchikey = cache._pubchem_species[spe.inchi]['inchikey']
                pubchem.smiles   = cache._pubchem_species[spe.inchi]['smiles']
                pubchem.xref     = cache._pubchem_species[spe.inchi]['xref']
            except KeyError:
                # print()
                # print("*************")
                # print("pubchem_species SEARCH (INCHI)", spe.inchi)
                # print("*************")
                # print()
                pubres = _pubchemStrctSearch(spe.inchi, 'inchi', logger=logger)
                if not chemName:
                    chemName = pubres['name']
                if 'chebi' in pubres['xref']:
                    try:
                        spe.xref = cache.cid_xref[cache.chebi_cid[pubres['xref']['chebi'][0]]]
                    except KeyError:
                        pass
                # pubchem.fill_missing(pubres)
                if not pubchem.inchi:
                    pubchem.inchi = pubres['inchi']
                if not pubchem.xref:
                    pubchem.xref = pubres['xref']
                if not pubchem.inchikey:
                    pubchem.inchikey = pubres['inchikey']
                if not pubchem.smiles:
                    pubchem.smiles = pubres['smiles']
    except KeyError:
        pass
    # inchikey
    try:
        # @Joan: The same question with InchI. Are you making sure that inchikey is not None
        spe.inchikey = rp_strc[meta]['inchikey']
        if not spe.xref and pubchem_search:
            # print("*************")
            # print("pubchem_species SEARCH (INCHIKEY)", spe.inchi)
            # print("*************")
            # print()
            pubres = _pubchemStrctSearch(spe.inchikey, 'inchikey', logger=logger)
            if not chemName:
                chemName = pubres['name']
            if 'chebi' in pubres['xref']:
                try:
                    spe.xref = cache.cid_xref[cache.chebi_cid[pubres['xref']['chebi'][0]]]
                except KeyError:
                    pass
            if not pubchem.xref:
                pubchem.xref = pubres['xref']
            if not pubchem.inchi:
                pubchem.inchi = pubres['inchi']
            if not pubchem.smiles:
                pubchem.smiles = pubres['smiles']
    except KeyError:
        pass
    #smiles
    try:
        #@Joan: The same question with InchI. Are you making sure that SMILES is not None
        spe.smiles = rp_strc[meta]['smiles']
        if not spe.xref and pubchem_search:
            # print()
            # print("*************")
            # print("pubchem_species SEARCH (SMILES)", spe.inchi)
            # print("*************")
            # print()
            pubres = _pubchemStrctSearch(spe.smiles, 'smiles', logger=logger)
            #print(pubres)
            if not chemName:
                chemName = pubres['name']
            if 'chebi' in pubres['xref']:
                try:
                    spe.xref = cache.cid_xref[cache.chebi_cid[pubres['xref']['chebi'][0]]]
                except KeyError:
                    pass
            if not pubchem.xref:
                pubchem.xref = pubres['xref']
            if not pubchem.inchi:
                pubchem.inchi = pubres['inchi']
            if not pubchem.inchikey:
                pubchem.inchikey = pubres['inchikey']
    except KeyError:
        pass

    if not spe.inchi:
        spe.inchi = pubchem.inchi
    if not spe.inchikey:
        spe.inchikey = pubchem.inchikey
    if not spe.smiles:
        spe.smiles = pubchem.smiles
    if not spe.xref:
        spe.xref = pubchem.xref
    if pubchem.inchi:
        # print()
        # print("*************")
        # print("pubchem_species WRITTEN")
        # print("*************")
        # print()
        cache._pubchem_species[pubchem.inchi] = {'inchi': pubchem.inchi, 'smiles': pubchem.smiles, 'inchikey': pubchem.inchikey, 'xref': pubchem.xref}

    return (chemName, spe)


def add_unique_species(rpsbml, species,
                       rp_strc, sink_molecules, compartment_id, species_group_id, sink_species_group_id, pubchem_search,
                       cache, logger=logging.getLogger(__name__)):
    for specie in species:
        (chemName, spe) = unique_species(cache, specie, rp_strc, pubchem_search, logger=logger)
        if chemName:
            chemName = chemName.replace("'", "")
        # pass the information to create the species
        rpsbml = add_species(rpsbml, specie, sink_molecules, compartment_id, chemName, spe, species_group_id, sink_species_group_id, logger=logger)
    return rpsbml


def add_species(rpsbml, meta, sink_molecules, compartment_id, chemName, spe, species_group_id, sink_species_group_id, logger=logging.getLogger(__name__)):
    if meta in sink_molecules:
        rpsbml.createSpecies(meta,
                             compartment_id,
                             chemName,
                             spe.xref,
                             spe.inchi,
                             spe.inchikey,
                             spe.smiles,
                             species_group_id,
                             sink_species_group_id)
    else:
        rpsbml.createSpecies(meta,
                             compartment_id,
                             chemName,
                             spe.xref,
                             spe.inchi,
                             spe.inchikey,
                             spe.smiles,
                             species_group_id)

    return rpsbml


#############################################################################################
############################### TSV data tsv ################################################
#############################################################################################

## Function to parse the TSV of measured heterologous pathways to SBML
#
# TODO: update this to the new compartements and others
# Given the TSV of measured pathways, parse them to a dictionnary, readable to next be parsed
# to SBML
#
# @param self object pointer
# @param inFile The input JSON file
# @param mnxHeader Reorganise the results around the target MNX products
# @return Dictionnary of SBML
# def _parseTSV(inFile, remove_inchi_4p=False, mnxHeader=False, logger=logging.getLogger(__name__)):
#     #     data = {}
#     try:
#         for row in csv_DictReader(open(inFile), delimiter='\t'):
#             ######## path_base_idxx ######
#             try:
#                 pathID = int(row['pathway_ID'])
#             except ValueError:
#                 logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
#                 continue
#             if not pathID in data:
#                 data[pathID] = {}
#                 data[pathID]['isValid'] = True
#                 data[pathID]['steps'] = {}
#             ####### target #########
#             if not 'target' in data[pathID]:
#                 data[pathID]['target'] = {}
#                 data[pathID]['target']['name'] = row['target_name']
#                 if remove_inchi_4p:
#                     data[pathID]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
#                 else:
#                     data[pathID]['target']['inchi'] = row['target_structure']
#             ####### step #########
#             try:
#                 stepID = int(row['rxn_idx'])
#             except ValueError:
#                 logger.error('Cannot convert step ID: '+str(row['rxn_idx']))
#                 data[pathID]['isValid'] = False
#                 continue
#             if stepID==0:
#                 continue
#             elif stepID==1:
#                 data[pathID]['organism'] = row['organism'].replace(' ', '')
#                 data[pathID]['reference'] = row['reference'].replace(' ', '')
#             data[pathID]['steps'][stepID] = {}
#             ##### substrates #########
#             data[pathID]['steps'][stepID]['substrates'] = []
#             lenDBref = len(row['substrate_dbref'].split(';'))
#             for i in row['substrate_dbref'].split(';'):
#                 if i=='':
#                     lenDBref -= 1
#             lenStrc = len(row['substrate_structure'].split('_'))
#             for i in row['substrate_structure'].split('_'):
#                 if i=='':
#                     lenStrc -= 1
#             lenSub = len(row['substrate_name'].split(';'))
#             for i in row['substrate_name'].split(';'):
#                 if i=='':
#                     lenSub -= 1
#             if lenSub==lenStrc==lenSub:
#                 for name, inchi, dbrefs in zip(row['substrate_name'].split(';'),
#                         row['substrate_structure'].split('_'),
#                         row['substrate_dbref'].split(';')):
#                     tmp = {}
#                     if remove_inchi_4p:
#                         tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
#                     else:
#                         tmp['inchi'] = inchi.replace(' ', '')
#                     tmp['name'] = name
#                     tmp['dbref'] = {}
#                     for dbref in dbrefs.split('|'):
#                         if len(dbref.split(':'))==2:
#                             db_name = dbref.split(':')[0].replace(' ', '').lower()
#                             db_cid = dbref.split(':')[1].replace(' ', '')
#                             if not db_name in tmp['dbref']:
#                                 tmp['dbref'][db_name] = []
#                             tmp['dbref'][db_name].append(db_cid)
#                         else:
#                             logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
#                             data[pathID]['isValid'] = False
#                     data[pathID]['steps'][stepID]['substrates'].append(tmp)
#             else:
#                 logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
#                 data[pathID]['isValid'] = False
#                 continue
#             ##### products #########
#             data[pathID]['steps'][stepID]['products'] = []
#             lenDBref = len(row['product_dbref'].split(';'))
#             for i in row['product_dbref'].split(';'):
#                 if i=='':
#                     lenDBref -= 1
#             lenStrc = len(row['product_structure'].split('_'))
#             for i in row['product_structure'].split('_'):
#                 if i=='':
#                     lenStrc -= 1
#             lenSub = len(row['product_name'].split(';'))
#             for i in row['product_name'].split(';'):
#                 if i=='':
#                     lenSub -= 1
#             if lenSub==lenStrc==lenDBref:
#                 for name, inchi, dbrefs in zip(row['product_name'].split(';'),
#                         row['product_structure'].split('_'),
#                         row['product_dbref'].split(';')):
#                     tmp = {}
#                     if remove_inchi_4p:
#                         tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
#                     else:
#                         tmp['inchi'] = inchi.replace(' ', '')
#                     tmp['name'] = name
#                     tmp['dbref'] = {}
#                     for dbref in dbrefs.split('|'):
#                         if len(dbref.split(':'))==2:
#                             db_name = dbref.split(':')[0].replace(' ', '').lower()
#                             db_cid = dbref.split(':')[1].replace(' ', '')
#                             if not db_name in tmp['dbref']:
#                                 tmp['dbref'][db_name] = []
#                             tmp['dbref'][db_name].append(db_cid)
#                         else:
#                             data[pathID]['isValid'] = False
#                             logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
#                     data[pathID]['steps'][stepID]['products'].append(tmp)
#             else:
#                 logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
#                 data[pathID]['isValid'] = False
#             if not row['uniprot']=='':
#                 data[pathID]['steps'][stepID]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
#             if not row['EC_number']=='':
#                 data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
#             data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
#             data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
#     except FileNotFoundError:
#         logger.error('Cannot open the file: '+str(inFile))
#     #now loop through all of them and remove the invalid paths
#     toRet = deepcopy(data)
#     for path_base_idxx in data.keys():
#         if toRet[path_base_idxx]['isValid']==False:
#             del toRet[path_base_idxx]
#         else:
#             del toRet[path_base_idxx]['isValid']
#     #reorganise the results around the target products mnx
#     if not mnxHeader:
#         return toRet
#     else:
#         toRetTwo = {}
#         for path_base_idxx in toRet:
#             try:
#                 final_pro_mnx = toRet[path_base_idx]['steps'][max(toRet[path_base_idx]['steps'])]['products'][0]['dbref']['mnx'][0]
#             except KeyError:
#                 logger.error('The species '+str(toRet[path_base_idx]['steps'][max(toRet[path_base_idx]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_base_idx))
#                 #continue
#             if not final_pro_mnx in toRetTwo:
#                 toRetTwo[final_pro_mnx] = {}
#             toRetTwo[final_pro_mnx][path_base_idx] = toRet[path_base_idx]
#         return toRetTwo


## Parse the validation TSV to SBML
#
# Parse the TSV file to SBML format and adds them to the sbml_paths
#
# @param self Object pointer
# @param inFile Input file
# @param compartment_id compartment of the
# TODO: update this with the new SBML groups
# def TSVtoSBML(cache,
#               inFile,
#               tmpOutputFolder=None,
#               upper_flux_bound=99999,
#               lower_flux_bound=0,
#               compartment_id='MNXC3',
#               pathway_id='rp_pathway',
#               species_group_id='central_species',
#               header_name='',
#               logger=logging.getLogger(__name__)):
#     #     data = _parseTSV(inFile, logger=logger)
#     sbml_paths = {}
#     if header_name=='':
#         header_name = inFile.split('/')[-1].replace('.tsv', '').replace('.csv', '')
#     # TODO: need to exit at this loop
#     for path_base_idx in data:
#         try:
#             mnxc = cache.deprecatedCompID_compid[compartment_id]
#         except KeyError:
#             logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
#             return False
#         rpsbml = rpSBML.rpSBML(name=header_name+'_'+str(path_base_idx), logger=logger)
#         # 1) create a generic Model, ie the structure and unit definitions that we will use the most
#         ##### TODO: give the user more control over a generic model creation:
#         # -> special attention to the compartment
#         rpsbml.genericModel(header_name+'_Path'+str(path_base_idx),
#                             header_name+'_Path'+str(path_base_idx),
#                             cache.comp_xref[mnxc],
#                             compartment_id,
#                             upper_flux_bound,
#                             lower_flux_bound)
#         # 2) create the pathway (groups)
#         rpsbml.createGroup(pathway_id)
#         rpsbml.createGroup(species_group_id)
#         # 3) find all the unique species and add them to the model
#         allChem = []
#         for stepNum in data[path_base_idx]['steps']:
#             # because of the nature of the input we need to remove duplicates
#             for i in data[path_base_idx]['steps'][stepNum]['substrates']+data[path_base_idx]['steps'][stepNum]['products']:
#                 if not i in allChem:
#                     allChem.append(i)
#         # add them to the SBML
#         for chem in allChem:
#             # PROBLEM: as it stands one expects the meta to be MNX
#             if 'mnx' in chem['dbref']:
#                 # must list the different models
#                 meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
#             else:
#                 # TODO: add the species with other types of xref in annotation
#                 logger.warning('Some species are not referenced by a MNX id and will be ignored')
#                 # try CHEBI
#                 try:
#                     meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
#                     meta = 'CHEBI_'+str(meta)
#                 except KeyError:
#                     # TODO: need to find a better way
#                     logger.warning('Cannot determine MNX or CHEBI entry, using random')
#                     tmpDB_name = list(chem['dbref'].keys())[0]
#                     meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
#                     meta = str(tmpDB_name)+'_'+str(meta)
#                 # break
#             # try to conver the inchi into the other structures
#             smiles = None
#             inchikey = None
#             try:
#                 resConv = cache._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
#                 smiles = resConv['smiles']
#                 inchikey = resConv['inchikey']
#             except NotImplementedError as e:
#                 logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
#             # create a new species
#             # here we want to gather the info from rpReader's rp_strc and cid_strc
#             try:
#                 chem_name = cache.cid_strc[meta]['name']
#             except KeyError:
#                 chem_name = meta
#             # compile as much info as you can
#             # xref
#             try:
#                 # TODO: add the xref from the document
#                 spe_xref = cache.cid_xref[meta]
#             except KeyError:
#                 #spe_xref = {}
#                 spe_xref = chem['dbref']
#             # inchi
#             try:
#                 spe_inchi = cache.cid_strc[meta]['inchi']
#             except KeyError:
#                 spe_inchi = chem['inchi']
#             # inchikey
#             try:
#                 spe_inchikey = cache.cid_strc[meta]['inchikey']
#             except KeyError:
#                 spe_inchikey =  resConv['inchikey']
#             # smiles
#             try:
#                 spe_smiles = cache.cid_strc[meta]['smiles']
#             except KeyError:
#                 spe_smiles = resConv['smiles']
#             # pass the information to create the species
#             rpsbml.createSpecies(meta,
#                                  compartment_id,
#                                  chem_name,
#                                  spe_xref,
#                                  spe_inchi,
#                                  spe_inchikey,
#                                  spe_smiles,
#                                  species_group_id)
#         # 4) add the complete reactions and their annotations
#         # create a new group for the measured pathway
#         # need to convert the validation to step for reactions
#         for stepNum in data[path_base_idx]['steps']:
#             toSend = {
#                 'left': {},
#                 'right': {},
#                 'rule_id': None,
#                 'rule_ori_reac': None,
#                 'rule_score': None,
#                 # 'path_base_idx': path_base_idx,
#                 'rxn_idx': stepNum,
#                 # 'sub_step': None
#                 }
#             for chem in data[path_base_idx]['steps'][stepNum]['substrates']:
#                 if 'mnx' in chem['dbref']:
#                     meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
#                     # try CHEBI
#                 else:
#                     logger.warning('Not all the species to have a MNX ID')
#                     # break
#                     try:
#                         meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
#                         meta = 'CHEBI_'+str(meta)
#                     except KeyError:
#                         # TODO: need to find a better way
#                         logger.warning('Cannot determine MNX or CHEBI entry, using random')
#                         tmpDB_name = list(chem['dbref'].keys())[0]
#                         meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
#                         meta = str(tmpDB_name)+'_'+str(meta)
#                 toSend['left'][meta] = 1
#             for chem in data[path_base_idx]['steps'][stepNum]['products']:
#                 if 'mnx' in chem['dbref']:
#                     meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
#                     # try CHEBI
#                 else:
#                     logger.warning('Need all the species to have a MNX ID')
#                     try:
#                         meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
#                         meta = 'CHEBI_'+str(meta)
#                     except KeyError:
#                         # TODO: need to find a better way
#                         logger.warning('Cannot determine MNX or CHEBI entry, using random')
#                         tmpDB_name = list(chem['dbref'].keys())[0]
#                         meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
#                         meta = str(tmpDB_name)+'_'+str(meta)
#                 toSend['right'][meta] = 1
#                     # break
#             # if all are full add it
#             reac_xref = {}
#             if 'ec_numbers' in data[path_base_idx]['steps'][stepNum]:
#                 reac_xref['ec'] = data[path_base_idx]['steps'][stepNum]['ec_numbers']
#             if 'uniprot' in data[path_base_idx]['steps'][stepNum]:
#                 reac_xref['uniprot'] = data[path_base_idx]['steps'][stepNum]['uniprot']
#             logger.debug('#########################################')
#             logger.debug(toSend)
#             logger.debug('#########################################')
#             rpsbml.createReaction(header_name+'_Step'+str(stepNum),
#                                   upper_flux_bound,
#                                   lower_flux_bound,
#                                   toSend,
#                                   compartment_id,
#                                   None,
#                                   reac_xref,
#                                   pathway_id)
#             if stepNum==1:
#                 # adding the consumption of the target
#                 targetStep = {'rule_id': None,
#                               'left': {},
#                               'right': {},
#                               'rxn_idx': None,
#                             #   'sub_step': None,
#                             #   'path_base_idx': None,
#                               'transformation_id': None,
#                               'rule_score': None,
#                               'rule_ori_reac': None}
#                 for chem in data[path_base_idx]['steps'][stepNum]['products']:
#                     try:
#                         # smallest MNX
#                         meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
#                     except KeyError:
#                         # try CHEBI
#                         try:
#                             meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
#                             meta = 'CHEBI_'+str(meta)
#                         except KeyError:
#                             logger.warning('Cannot determine MNX or CHEBI entry, using random')
#                             tmpDB_name = list(chem['dbref'].keys())[0]
#                             meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
#                             meta = str(tmpDB_name)+'_'+str(meta)
#                     targetStep['left'][meta] = 1
#                 rpsbml.createReaction(header_name+'_Step1_sink',
#                                       upper_flux_bound,
#                                       lower_flux_bound,
#                                       targetStep,
#                                       compartment_id)
#                 rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
#         if tmpOutputFolder:
#             rpsbml.writeSBML(tmpOutputFolder)
#         else:
#             sbml_paths[header_name+'_Path'+str(path_base_idx)] = rpsbml
#     if tmpOutputFolder:
#         return {}
#     else:
#         return sbml_paths
