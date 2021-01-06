from csv       import reader         as csv_reader
from csv       import DictReader     as csv_DictReader
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
from rptools.rpcompletion.rpCofactors import addCofactors
import logging

#import rpCofactors

## @package rpCompletion
#
# Collection of functions that convert the outputs from various sources to the SBML format (rpSBML) for further analyses

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




# ## Class to read all the input files
# #
# # Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
# class rpCompletion():
#     ## InputReader constructor
#     #
#     # @param self The object pointer
#     # @param db Type of cache database to use
#     def __init__(db='file'):
#         super().__init__(db)
        # self.rpcofactors = rpCofactors(db, self.print)

#####################
pubchem_min_count = 0
pubchem_min_start = 0.0


#######################################################################
############################# PRIVATE FUNCTIONS #######################
#######################################################################

def _pubChemLimit(logger=None):
    global pubchem_min_count, pubchem_min_start
    logger = logger or logging.getLogger(__name__)
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
def _pubchemStrctSearch(strct, itype='inchi', logger=None):
    logger = logger or logging.getLogger(__name__)
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
                  logger=None):

    logger = logger or logging.getLogger(__name__)

    if max_subpaths_filter<0:
        raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

    if not os_path.exists(outdir):
        os_mkdir(outdir)

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
#  @param self Object pointer
#  @param path The compounds.txt file path
#  @return rp_compounds Dictionnary of smile and structure for each compound
def _compounds(cache, path, logger=None):
    logger = logger or logging.getLogger(__name__)
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
# @param self Object pointer
# @param path The scope.csv file path
# @return tuple with dictionnary of all the reactions rules and the list unique molecules that these apply them to
def _transformation(path, logger=None):
    logger = logger or logging.getLogger(__name__)
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


## Function that reads the pathway output of rp2paths
#
# @param self Object pointer
# @param rp2paths_pathways file parsing. Can either be a string of the file path or a bytes object
# @return Dictionnary of the pathways
def _read_paths(cache, rp2paths_pathways, logger=None):

    logger = logger or logging.getLogger(__name__)

    #### we might pass binary in the REST version
    if isinstance(rp2paths_pathways, bytes):
        reader = csv_reader(StringIO(rp2paths_pathways.decode('utf-8')))
    else:
        reader = csv_reader(open(rp2paths_pathways, 'r'))
    next(reader)
    current_path_id = 0
    path_step = 1

    rp_paths = {}

    for row in reader:
        try:
            #Remove all illegal characters in SBML ids
            row[3] = row[3].replace("'", "").replace('-', '_').replace('+', '')
            if not int(row[0])==current_path_id:
                path_step = 1
            else:
                path_step += 1
            #important to leave them in order
            current_path_id = int(row[0])
        except ValueError:
            logger.error('Cannot convert path_id to int ('+str(row[0])+')')
            raise
        #################################
        ruleIds = row[2].split(',')
        if ruleIds==None:
            logger.warning('The rulesIds is None')
            #pass # or continue
            continue
        ###WARNING: This is the part where we select some rules over others
        # we do it by sorting the list according to their score and taking the topx
        tmp_rr_reactions = {}
        for r_id in ruleIds:
            for rea_id in cache.rr_reactions[r_id]:
                tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = cache.rr_reactions[r_id][rea_id]
        # if len(ruleIds)>int(maxRuleIds):
        #     logger.warning('There are too many rules, limiting the number to random top '+str(maxRuleIds))
        #     try:
        #         ruleIds = [y for y,_ in sorted([(i, tmp_rr_reactions[i]['rule_score']) for i in tmp_rr_reactions])][:int(maxRuleIds)]
        #     except KeyError:
        #         logger.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
        #         ruleIds = random.sample(tmp_rr_reactions, int(maxRuleIds))
        # else:
        ruleIds = tmp_rr_reactions
        sub_path_step = 1
        for singleRule in ruleIds:
            tmpReac = {'rule_id': singleRule.split('__')[0],
                       'rule_ori_reac': singleRule.split('__')[1],
                       'rule_score': cache.rr_reactions[singleRule.split('__')[0]][singleRule.split('__')[1]]['rule_score'],
                       'right': {},
                       'left': {},
                       'path_id': int(row[0]),
                       'step': path_step,
                       'transformation_id': row[1][:-2]}
            ############ LEFT ##############
            for l in row[3].split(':'):
                tmp_l = l.split('.')
                try:
                    #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                    cid = '' #TODO: change this
                    if tmp_l[1] in cache.deprecatedCID_cid:
                        cid = cache.deprecatedCID_cid[tmp_l[1]]
                    else:
                        cid = tmp_l[1]
                    tmpReac['left'][cid] = int(tmp_l[0])
                except ValueError:
                    logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                    #return {}
                    return False
            ############## RIGHT ###########
            for r in row[4].split(':'):
                tmp_r = r.split('.')
                try:
                    #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                    cid = '' #TODO change this
                    if tmp_r[1] in cache.deprecatedCID_cid:
                        cid = cache.deprecatedCID_cid[tmp_r[1]]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                    else:
                        cid = tmp_r[1]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                    tmpReac['right'][cid] = int(tmp_r[0])
                except ValueError:
                    logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                    return {}
            #################################
            if not int(row[0]) in rp_paths:
                rp_paths[int(row[0])] = {}
            if not int(path_step) in rp_paths[int(row[0])]:
                rp_paths[int(row[0])][int(path_step)] = {}
            rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
            sub_path_step += 1

    return rp_paths


def _unique_species(cache, meta, rp_strc, pubchem_search, logger=None):

    logger = logger or logging.getLogger(__name__)

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
                             rp_strc, rp_transformation,
                             sink_molecules,
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
                             logger=None):
    # TODO: make sure that you account for the fact that each reaction may have multiple associated reactions

    logger = logger or logging.getLogger(__name__)

    rp_paths = _read_paths(cache, rp2paths_pathways, logger=logger)
    sink_species = []

    # for each line or rp2paths_pathways:
    #     generate comb
    #     for each combinant:
    #         rank
    #         process
    #         add cofactors
    #         dedup

    #### pathToSBML ####
    try:
        compid = cache.deprecatedCompID_compid[compartment_id]
    except KeyError:
        logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
        return 1


    for pathNum in rp_paths:

        # first level is the list of lists of sub_steps
        # second is itertools all possible combinations using product
        altPathNum = 1
        # topX subpaths of the current rp2path pathway
        local_SBMLItems = []

        for comb_path in list(itertools_product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
            steps = []
            for i, y in comb_path:
                steps.append(rp_paths[pathNum][i][y])
            path_id = steps[0]['path_id']
            rpsbml = rpSBML(name='rp_'+str(path_id)+'_'+str(altPathNum))

            # 1) Create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            # -> special attention to the compartment
            rpsbml.genericModel(
                    'RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                    'RP_model_'+str(path_id)+'_'+str(altPathNum),
                    cache.comp_xref[compid],
                    compartment_id,
                    upper_flux_bound,
                    lower_flux_bound)

            # 2) Create the pathway (groups)
            rpsbml.createPathway(pathway_id)
            rpsbml.createPathway(species_group_id)
            rpsbml.createPathway(sink_species_group_id)

            # 3) Find all unique species and add them to the model
            all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
            for meta in all_meta:
                (chemName, spe) = _unique_species(cache, meta, rp_strc, pubchem_search, logger=logger)
                if chemName:
                    chemName = chemName.replace("'", "")
                # pass the information to create the species
                rpsbml = add_species(rpsbml, meta, sink_molecules, compartment_id, chemName, spe, species_group_id, sink_species_group_id, logger=logger)

            # 4) Add the complete reactions and their annotations
            for step in steps:
                # add the substep to the model
                step['sub_step'] = altPathNum
                rpsbml.createReaction(
                        'RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                        upper_flux_bound, lower_flux_bound, step, compartment_id,
                        rp_transformation[step['transformation_id']]['rule'],
                        {'ec': rp_transformation[step['transformation_id']]['ec']},
                        pathway_id)

            # 5) Adding the consumption of the target
            targetStep = {
                    'rule_id': None,
                    'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1},
                    'right': [],
                    'step': None,
                    'sub_step': None,
                    'path_id': None,
                    'transformation_id': None,
                    'rule_score': None,
                    'rule_ori_reac': None
                    }
            rpsbml.createReaction('RP1_sink',
                                  upper_flux_bound, lower_flux_bound,
                                  targetStep,
                                  compartment_id)

            # 6) Adding the cofactors
            addCofactors(cache, rpsbml, logger=logger)

            # 7) Insert the new rpsbml object in sorted rpsbml_items list
            sbml_item = SBML_Item(rpsbml.compute_score(),
                                  'rp_'+str(path_id)+'_'+str(altPathNum),
                                  rpsbml)
            local_SBMLItems = insert_and_or_replace_in_sorted_list(sbml_item, local_SBMLItems)

            # 8) Keep only topX
            local_SBMLItems = local_SBMLItems[-max_subpaths_filter:]

            altPathNum += 1

        # Write results to files
        for rpsbml_item in local_SBMLItems:
            rpsbml_item.rpsbml_obj.writeSBML(os_path.join(outFolder, str(rpsbml_item.rpsbml_obj.modelName))+'_sbml.xml')

    return 0


def add_species(rpsbml, meta, sink_molecules, compartment_id, chemName, spe, species_group_id, sink_species_group_id, logger=None):
    logger = logger or logging.getLogger(__name__)
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
def _parseTSV(inFile, remove_inchi_4p=False, mnxHeader=False, logger=None):
    logger = logger or logging.getLogger(__name__)
    data = {}
    try:
        for row in csv_DictReader(open(inFile), delimiter='\t'):
            ######## path_id ######
            try:
                pathID = int(row['pathway_ID'])
            except ValueError:
                logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                continue
            if not pathID in data:
                data[pathID] = {}
                data[pathID]['isValid'] = True
                data[pathID]['steps'] = {}
            ####### target #########
            if not 'target' in data[pathID]:
                data[pathID]['target'] = {}
                data[pathID]['target']['name'] = row['target_name']
                if remove_inchi_4p:
                    data[pathID]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
                else:
                    data[pathID]['target']['inchi'] = row['target_structure']
            ####### step #########
            try:
                stepID = int(row['step'])
            except ValueError:
                logger.error('Cannot convert step ID: '+str(row['step']))
                data[pathID]['isValid'] = False
                continue
            if stepID==0:
                continue
            elif stepID==1:
                data[pathID]['organism'] = row['organism'].replace(' ', '')
                data[pathID]['reference'] = row['reference'].replace(' ', '')
            data[pathID]['steps'][stepID] = {}
            ##### substrates #########
            data[pathID]['steps'][stepID]['substrates'] = []
            lenDBref = len(row['substrate_dbref'].split(';'))
            for i in row['substrate_dbref'].split(';'):
                if i=='':
                    lenDBref -= 1
            lenStrc = len(row['substrate_structure'].split('_'))
            for i in row['substrate_structure'].split('_'):
                if i=='':
                    lenStrc -= 1
            lenSub = len(row['substrate_name'].split(';'))
            for i in row['substrate_name'].split(';'):
                if i=='':
                    lenSub -= 1
            if lenSub==lenStrc==lenSub:
                for name, inchi, dbrefs in zip(row['substrate_name'].split(';'),
                        row['substrate_structure'].split('_'),
                        row['substrate_dbref'].split(';')):
                    tmp = {}
                    if remove_inchi_4p:
                        tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                    else:
                        tmp['inchi'] = inchi.replace(' ', '')
                    tmp['name'] = name
                    tmp['dbref'] = {}
                    for dbref in dbrefs.split('|'):
                        if len(dbref.split(':'))==2:
                            db_name = dbref.split(':')[0].replace(' ', '').lower()
                            db_cid = dbref.split(':')[1].replace(' ', '')
                            if not db_name in tmp['dbref']:
                                tmp['dbref'][db_name] = []
                            tmp['dbref'][db_name].append(db_cid)
                        else:
                            logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                            data[pathID]['isValid'] = False
                    data[pathID]['steps'][stepID]['substrates'].append(tmp)
            else:
                logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                data[pathID]['isValid'] = False
                continue
            ##### products #########
            data[pathID]['steps'][stepID]['products'] = []
            lenDBref = len(row['product_dbref'].split(';'))
            for i in row['product_dbref'].split(';'):
                if i=='':
                    lenDBref -= 1
            lenStrc = len(row['product_structure'].split('_'))
            for i in row['product_structure'].split('_'):
                if i=='':
                    lenStrc -= 1
            lenSub = len(row['product_name'].split(';'))
            for i in row['product_name'].split(';'):
                if i=='':
                    lenSub -= 1
            if lenSub==lenStrc==lenDBref:
                for name, inchi, dbrefs in zip(row['product_name'].split(';'),
                        row['product_structure'].split('_'),
                        row['product_dbref'].split(';')):
                    tmp = {}
                    if remove_inchi_4p:
                        tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                    else:
                        tmp['inchi'] = inchi.replace(' ', '')
                    tmp['name'] = name
                    tmp['dbref'] = {}
                    for dbref in dbrefs.split('|'):
                        if len(dbref.split(':'))==2:
                            db_name = dbref.split(':')[0].replace(' ', '').lower()
                            db_cid = dbref.split(':')[1].replace(' ', '')
                            if not db_name in tmp['dbref']:
                                tmp['dbref'][db_name] = []
                            tmp['dbref'][db_name].append(db_cid)
                        else:
                            data[pathID]['isValid'] = False
                            logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                    data[pathID]['steps'][stepID]['products'].append(tmp)
            else:
                logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                data[pathID]['isValid'] = False
            if not row['uniprot']=='':
                data[pathID]['steps'][stepID]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
            if not row['EC_number']=='':
                data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
            data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
            data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
    except FileNotFoundError:
        logger.error('Cannot open the file: '+str(inFile))
    #now loop through all of them and remove the invalid paths
    toRet = deepcopy(data)
    for path_id in data.keys():
        if toRet[path_id]['isValid']==False:
            del toRet[path_id]
        else:
            del toRet[path_id]['isValid']
    #reorganise the results around the target products mnx
    if not mnxHeader:
        return toRet
    else:
        toRetTwo = {}
        for path_id in toRet:
            try:
                final_pro_mnx = toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['dbref']['mnx'][0]
            except KeyError:
                logger.error('The species '+str(toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_id))
                #continue
            if not final_pro_mnx in toRetTwo:
                toRetTwo[final_pro_mnx] = {}
            toRetTwo[final_pro_mnx][path_id] = toRet[path_id]
        return toRetTwo


## Parse the validation TSV to SBML
#
# Parse the TSV file to SBML format and adds them to the sbml_paths
#
# @param self Object pointer
# @param inFile Input file
# @param compartment_id compartment of the
# TODO: update this with the new SBML groups
def TSVtoSBML(cache,
              inFile,
              tmpOutputFolder=None,
              upper_flux_bound=99999,
              lower_flux_bound=0,
              compartment_id='MNXC3',
              pathway_id='rp_pathway',
              species_group_id='central_species',
              header_name='',
              logger=None):
    logger = logger or logging.getLogger(__name__)
    data = _parseTSV(inFile, logger=logger)
    sbml_paths = {}
    if header_name=='':
        header_name = inFile.split('/')[-1].replace('.tsv', '').replace('.csv', '')
    # TODO: need to exit at this loop
    for path_id in data:
        try:
            mnxc = cache.deprecatedCompID_compid[compartment_id]
        except KeyError:
            logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        rpsbml = rpSBML.rpSBML(name=header_name+'_'+str(path_id), logger=logger)
        # 1) create a generic Model, ie the structure and unit definitions that we will use the most
        ##### TODO: give the user more control over a generic model creation:
        # -> special attention to the compartment
        rpsbml.genericModel(header_name+'_Path'+str(path_id),
                            header_name+'_Path'+str(path_id),
                            cache.comp_xref[mnxc],
                            compartment_id,
                            upper_flux_bound,
                            lower_flux_bound)
        # 2) create the pathway (groups)
        rpsbml.createPathway(pathway_id)
        rpsbml.createPathway(species_group_id)
        # 3) find all the unique species and add them to the model
        allChem = []
        for stepNum in data[path_id]['steps']:
            # because of the nature of the input we need to remove duplicates
            for i in data[path_id]['steps'][stepNum]['substrates']+data[path_id]['steps'][stepNum]['products']:
                if not i in allChem:
                    allChem.append(i)
        # add them to the SBML
        for chem in allChem:
            # PROBLEM: as it stands one expects the meta to be MNX
            if 'mnx' in chem['dbref']:
                # must list the different models
                meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
            else:
                # TODO: add the species with other types of xref in annotation
                logger.warning('Some species are not referenced by a MNX id and will be ignored')
                # try CHEBI
                try:
                    meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                    meta = 'CHEBI_'+str(meta)
                except KeyError:
                    # TODO: need to find a better way
                    logger.warning('Cannot determine MNX or CHEBI entry, using random')
                    tmpDB_name = list(chem['dbref'].keys())[0]
                    meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                    meta = str(tmpDB_name)+'_'+str(meta)
                # break
            # try to conver the inchi into the other structures
            smiles = None
            inchikey = None
            try:
                resConv = cache._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                smiles = resConv['smiles']
                inchikey = resConv['inchikey']
            except NotImplementedError as e:
                logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
            # create a new species
            # here we want to gather the info from rpReader's rp_strc and cid_strc
            try:
                chem_name = cache.cid_strc[meta]['name']
            except KeyError:
                chem_name = meta
            # compile as much info as you can
            # xref
            try:
                # TODO: add the xref from the document
                spe_xref = cache.cid_xref[meta]
            except KeyError:
                #spe_xref = {}
                spe_xref = chem['dbref']
            # inchi
            try:
                spe_inchi = cache.cid_strc[meta]['inchi']
            except KeyError:
                spe_inchi = chem['inchi']
            # inchikey
            try:
                spe_inchikey = cache.cid_strc[meta]['inchikey']
            except KeyError:
                spe_inchikey =  resConv['inchikey']
            # smiles
            try:
                spe_smiles = cache.cid_strc[meta]['smiles']
            except KeyError:
                spe_smiles = resConv['smiles']
            # pass the information to create the species
            rpsbml.createSpecies(meta,
                                 compartment_id,
                                 chem_name,
                                 spe_xref,
                                 spe_inchi,
                                 spe_inchikey,
                                 spe_smiles,
                                 species_group_id)
        # 4) add the complete reactions and their annotations
        # create a new group for the measured pathway
        # need to convert the validation to step for reactions
        for stepNum in data[path_id]['steps']:
            toSend = {'left': {}, 'right': {}, 'rule_id': None, 'rule_ori_reac': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
            for chem in data[path_id]['steps'][stepNum]['substrates']:
                if 'mnx' in chem['dbref']:
                    meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                    # try CHEBI
                else:
                    logger.warning('Not all the species to have a MNX ID')
                    # break
                    try:
                        meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        meta = 'CHEBI_'+str(meta)
                    except KeyError:
                        # TODO: need to find a better way
                        logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmpDB_name = list(chem['dbref'].keys())[0]
                        meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        meta = str(tmpDB_name)+'_'+str(meta)
                toSend['left'][meta] = 1
            for chem in data[path_id]['steps'][stepNum]['products']:
                if 'mnx' in chem['dbref']:
                    meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                    # try CHEBI
                else:
                    logger.warning('Need all the species to have a MNX ID')
                    try:
                        meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        meta = 'CHEBI_'+str(meta)
                    except KeyError:
                        # TODO: need to find a better way
                        logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmpDB_name = list(chem['dbref'].keys())[0]
                        meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        meta = str(tmpDB_name)+'_'+str(meta)
                toSend['right'][meta] = 1
                    # break
            # if all are full add it
            reac_xref = {}
            if 'ec_numbers' in data[path_id]['steps'][stepNum]:
                reac_xref['ec'] = data[path_id]['steps'][stepNum]['ec_numbers']
            if 'uniprot' in data[path_id]['steps'][stepNum]:
                reac_xref['uniprot'] = data[path_id]['steps'][stepNum]['uniprot']
            logger.debug('#########################################')
            logger.debug(toSend)
            logger.debug('#########################################')
            rpsbml.createReaction(header_name+'_Step'+str(stepNum),
                                  upper_flux_bound,
                                  lower_flux_bound,
                                  toSend,
                                  compartment_id,
                                  None,
                                  reac_xref,
                                  pathway_id)
            if stepNum==1:
                # adding the consumption of the target
                targetStep = {'rule_id': None,
                              'left': {},
                              'right': {},
                              'step': None,
                              'sub_step': None,
                              'path_id': None,
                              'transformation_id': None,
                              'rule_score': None,
                              'rule_ori_reac': None}
                for chem in data[path_id]['steps'][stepNum]['products']:
                    try:
                        # smallest MNX
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                    except KeyError:
                        # try CHEBI
                        try:
                            meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            meta = 'CHEBI_'+str(meta)
                        except KeyError:
                            logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            meta = str(tmpDB_name)+'_'+str(meta)
                    targetStep['left'][meta] = 1
                rpsbml.createReaction(header_name+'_Step1_sink',
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      targetStep,
                                      compartment_id)
                rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
        if tmpOutputFolder:
            rpsbml.writeSBML(tmpOutputFolder)
        else:
            sbml_paths[header_name+'_Path'+str(path_id)] = rpsbml
    if tmpOutputFolder:
        return {}
    else:
        return sbml_paths
