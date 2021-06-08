from csv import reader as csv_reader
from requests import (
    post as r_post,
    get as r_get
)
from itertools import product as itertools_product
from time import (
    time as time_time,
    sleep as time_sleep
)
from argparse import ArgumentParser as argparse_ArgumentParser
from json import decoder as json_decoder
from io import StringIO
import pandas as pd
from typing import (
    List,
    Dict,
    Tuple,
    TypeVar
)
from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from hashlib import md5
from copy import deepcopy
from colored import fg, bg, attr
from brs_utils import (
    insert_and_or_replace_in_sorted_list,
    Item,
    Cache
)
from rr_cache import rrCache
from rxn_rebuild import rebuild_rxn
from chemlite import (
    Reaction,
    Compound
)
from rptools.rplibs import rpPathway

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
def rp_completion(
    rp2_metnet,
    rp2paths_compounds,
    rp2paths_pathways,
    cache: rrCache = None,
    upper_flux_bound=999999,
    lower_flux_bound=0,
    pathway_id='rp_pathway',
    compartment_id='MNXC3',
    species_group_id='central_species',
    sink_species_group_id='sink',
    # pubchem_search=False,
    logger: Logger = getLogger(__name__)
):

    if cache is None:
        cache = rrCache(
            db='file',
            attrs=[
                'rr_reactions',
                'template_reactions',
                'cid_strc',
                # 'deprecatedCID_cid',
                'comp_xref',
                'deprecatedCompID_compid',
                # 'cid_xref',
                # 'cid_name',
                # 'deprecatedRID_rid'
            ]
            # logger=logger
        )

    ## READ
    read_compounds(
        rp2paths_compounds,
        cache,
        logger=logger
    )
    pathways, transfos = read_pathways(
        rp2paths_pathways,
        cache.get('rr_reactions'),
        logger=logger
    )
    ec_numbers, sink_molecules = read_rp2_metnet(
        rp2_metnet,
        logger=logger
    )

    # COMPLETE TRANSFORMATIONS
    full_transfos = complete_transformations(
        transfos=transfos,
        ec_numbers=ec_numbers,
        cache=cache,
        logger=logger
    )

    # GENERATE THE COMBINATORY OF SUB-PATHWAYS
    # Build pathways over:
    #   - multiple reaction rules per transformation (TRS) and
    #   - multiple template reactions per reaction rule
    pathway_combinatorics = build_pathway_combinatorics(
        full_transfos,
        pathways,
        logger=logger
    )

    # BUILD + RANK SUB-PATHWAYS 
    all_pathways = build_all_pathways(
        pathways=pathway_combinatorics,
        transfos=full_transfos,
        sink_molecules=sink_molecules,
        rr_reactions=cache.get('rr_reactions'),
        compounds_cache=cache.get('cid_strc'),
        logger=logger
    )

    # for sub_pathways in all_pathways.values():
    #     for sub_pathway in sub_pathways:
    #         print(sub_pathway)

    return all_pathways


def complete_transformations(
    transfos: Dict,
    ec_numbers: Dict,
    cache: rrCache,
    logger: Logger = getLogger(__name__)
) -> Dict:

    full_transfos = {}

    # For each transformation
    for transfo_id, transfo in transfos.items():

        full_transfos[transfo_id] = {}
        full_transfos[transfo_id]['ec'] = ec_numbers[transfo_id]['ec']
        # Convert transformation into SMILES
        transfo_smi = '{left}>>{right}'.format(
                left=build_smiles(transfo['left']),
                right=build_smiles(transfo['right'])
            )

        # Add compounds of the current transformation
        full_transfos[transfo_id]['left'] = dict(transfos[transfo_id]['left'])
        full_transfos[transfo_id]['right'] = dict(transfos[transfo_id]['right'])
        full_transfos[transfo_id]['complement'] = {}

        # MULTIPLE RR FOR ONE TRANSFO
        for rule_id in transfo['rule_id']:

            # MULTIPLE TEMPLATE REACTIONS FOR ONE RR
            # If 'tmpl_rxn_id' is not given,
            # the transformation will be completed
            # for each template reaction from reaction rule was built from
            full_transfos[transfo_id]['complement'][rule_id] = rebuild_rxn(
                cache=cache,
                rxn_rule_id=rule_id,
                transfo=transfo_smi,
                direction='forward',
                # tmpl_rxn_id=tmpl_rxn_id,
                logger=logger
            )

    return full_transfos


def build_smiles(
    side: Dict,
    logger: Logger = getLogger(__name__)
) -> str:
    return '.'.join([Cache.get(spe_id).get_smiles() for spe_id in side.keys()])

## Function to parse the compounds.txt file
#
#  Extract the smile and the structure of each compounds of RP2Path output
#  Method to parse all the RP output compounds.
#
#  @param path The compounds.txt file path
#  @return rp_compounds Dictionnary of smile and structure for each compound
def read_compounds(path, cache, logger=getLogger(__name__)):

    try:
        if isinstance(path, bytes):
            reader = csv_reader(StringIO(path.decode('utf-8')), delimiter='\t')
        else:
            reader = csv_reader(open(path, 'r', encoding='utf-8'), delimiter='\t')
        next(reader)
        for row in reader:
            smiles = row[1]  #, 'structure':row[1].replace('[','').replace(']','')
            try:
                inchi = cache.get('cid_strc')[row[0]]['inchi']
            except KeyError:
                # try to generate them yourself by converting them directly
                try:
                    resConv = cache._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                    inchi = resConv['inchi']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
            try:
                inchikey = cache.get('cid_strc')[row[0]]['inchikey']
                # try to generate them yourself by converting them directly
                # TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
            except KeyError:
                try:
                    resConv = cache._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                    inchikey = resConv['inchikey']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
            compound = Compound(
                id=row[0],
                smiles=smiles,
                inchi=inchi,
                inchikey=inchikey
            )
            Cache.add(compound)

    except (TypeError, FileNotFoundError) as e:
        logger.error('Could not read the compounds file ('+str(path)+')')
        raise RuntimeError


## Function to parse the scope.csv file
#
# Extract the reaction rules from the retroPath2.0 output using the scope.csv file
#
# @param path The scope.csv file path
# @return tuple with dictionnary of all the reactions rules and the list unique molecules that these apply them to
def read_rp2_metnet(path, logger=getLogger(__name__)):

    ec_numbers = {}
    sink_molecules = set()
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
        if row[1] not in ec_numbers:
            ec_numbers[row[1]] = {
                # 'source': row[0],
                # 'transfo_smiles': row[2],
                # 'substrates': {
                #     'smiles': row[3],
                #     'inchi': row[4],
                # },
                # 'products': {
                #     'smiles': row[5],
                #     'inchi': row[6],
                # },
                # 'in_sink': row[7],
                # 'sink_name': [i.replace(' ', '') for i in row[8][1:-1].split(',')],
                # 'diameter': row[9],
                # 'rule_ids': [i.replace(' ', '') for i in row[10][1:-1].split(',')],
                'ec': [i.replace(' ', '') for i in row[11][1:-1].split(',') if i.replace(' ', '')!='NOEC'],
                # 'score': row[12],
                # 'start_src_smiles': row[13],
                # 'iteration': row[14]
            }
        if row[7]=='1':
            for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                sink_molecules.add(i)

    logger.debug(ec_numbers)
    logger.debug(list(sink_molecules))
    return ec_numbers, list(set(sink_molecules))
    # return list(sink_molecules)


# @param infile rp2_pathways file
# @param rr_reactions RetroRules reactions
# @param deprecatedCID_cid Deprecated cid
# @return Dictionnary of the pathways
def read_pathways(infile, rr_reactions, logger=getLogger(__name__)):

    df = pd.read_csv(infile)

    check = check_pathways(df)
    if not check:
        logger.error(check)
        exit()

    pathways = {}
    transfos = {}

    for index, row in df.iterrows():

        path_id = row['Path ID']
        transfo_id = row['Unique ID'][:-2]

        if path_id in pathways:
            pathways[path_id] += [transfo_id]
        else:
            pathways[path_id] = [transfo_id]

        if transfo_id not in transfos:
            transfos[transfo_id] = {}
            transfos[transfo_id]['rule_id'] = row['Rule ID'].split(',')
            for side in ['left', 'right']:
                transfos[transfo_id][side] = {}
                # split compounds
                compounds = row[side[0].upper()+side[1:]].split(':')
                # read compound and its stochio 
                for compound in compounds:
                    sto, spe = compound.split('.')
                    transfos[transfo_id][side][spe] = int(sto)

    return pathways, transfos


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


def build_pathway_combinatorics(
    full_transfos: Dict,
    pathways: Dict,
    logger=getLogger(__name__)
) -> Dict:
    '''Build all combinations of sub-pathways based on these facts:
          - one single transformation can have been formed from multiple reaction rules, and
          - one single reaction rule can have been generated from multiple template reactions
       To each such a combination corresponds a different complete transformation, i.e. reaction,
       then a different sub-pathway.

    Parameters
    ----------
    full_transfos : Dict
        Set of completed transformations
    pathways : Dict
        Set of pathways (pathway: list of transformation IDs)
    logger : Logger

    Returns
    -------
    Set of sub-pathways. Each transfomation IDs in 'pathways'
    has been replaced by a list of triplet(transfo_id, rule_id, tmpl_rxn_id)

    '''

    # BUILD PATHWAYS WITH ALL REACTIONS (OVER REACTION RULES * TEMPLATE REACTIONS)
    pathways_all_reactions = {}

    ## ITERATE OVER PATHWAYS
    for pathway, transfos_lst in pathways.items():

        # index of chemical reaction within the pathway
        transfo_idx = 0

        ## ITERATE OVER TRANSFORMATIONS
        # For each transformation of the current pathway
        # Iterate in retrosynthesis order (reverse)
        #   to better combine all sub-pathways
        for transfo_id in transfos_lst:

            transfo_idx += 1
            # Compounds from original transformation
            compounds = {
                'right': dict(full_transfos[transfo_id]['right']),
                'left': dict(full_transfos[transfo_id]['left'])
            }
            # Build list of transformations
            # where each transfo can correspond to multiple reactions
            # due to multiple reaction rules and/or multiple template reactions
            if pathway not in pathways_all_reactions:
                pathways_all_reactions[pathway] = []

            ## ITERATE OVER REACTION RULES
            # Multiple reaction rules for the current transformation?
            for rule_id, tmpl_rxns in full_transfos[transfo_id]['complement'].items():

                pathways_all_reactions[pathway].append([])
                ## ITERATE OVER TEMPLATE REACTIONS
                # Current reaction rule generated from multiple template reactions?
                for tmpl_rxn_id, tmpl_rxn in tmpl_rxns.items():

                    # Add template reaction compounds
                    compounds = add_compounds(
                        compounds,
                        tmpl_rxn['added_cmpds']
                    )

                    # Add the triplet ID to identify the sub_pathway
                    pathways_all_reactions[pathway][-1].append(
                        {
                            'transfo_id': transfo_id,
                            'rule_id': rule_id,
                            'tmpl_rxn_id': tmpl_rxn_id
                        }
                    )

    return pathways_all_reactions


def build_all_pathways(
    pathways: Dict,
    transfos: Dict,
    sink_molecules: List,
    rr_reactions: Dict,
    compounds_cache: Dict,
    logger: Logger = getLogger(__name__)
) -> Dict:

    res_pathways = {}

    ## PATHWAYS
    for path_idx, transfos_lst in pathways.items():

        # Combine over multiple template reactions
        sub_pathways = list(itertools_product(*transfos_lst))

        ## SUB-PATHWAYS
        # # Keep only topX best sub_pathways
        # # within a same master pathway
        res_pathways[path_idx] = []
        for sub_path_idx in range(len(sub_pathways)):

            pathway = rpPathway(
                id=str(path_idx).zfill(3)+'_'+str(sub_path_idx+1).zfill(4),
                logger=logger
            )

            ## ITERATE OVER REACTIONS
            nb_reactions = len(sub_pathways[sub_path_idx])
            for rxn_idx in range(nb_reactions):

                rxn = sub_pathways[sub_path_idx][rxn_idx]
                transfo_id = rxn['transfo_id']
                transfo = transfos[transfo_id]
                rule_id = rxn['rule_id']
                tmpl_rxn_id = rxn['tmpl_rxn_id']

                ## COMPOUNDS
                # Compounds from original transformation
                compounds = {
                    'right': deepcopy(transfo['right']),
                    'left': deepcopy(transfo['left'])
                }
                # Add template reaction compounds
                compounds = add_compounds(
                    compounds,
                    transfo['complement'][rule_id][tmpl_rxn_id]['added_cmpds']
                )

                ## REACTION
                # revert reaction index (forward)
                rxn_idx_forward = nb_reactions - rxn_idx
                rxn = Reaction(
                    id='rxn_'+str(rxn_idx_forward),
                    ec_numbers=transfo['ec'],
                    # reactants=dict(compounds['left']),
                    # products=dict(compounds['right'])
                )
                rxn.set_infos(
                    {
                        **sub_pathways[sub_path_idx][rxn_idx],
                        **{
                            'rule_score': rr_reactions[rule_id][tmpl_rxn_id]['rule_score'],
                            'idx_in_path': rxn_idx_forward
                        }
                    }
                )
                # Reactants
                for spe_id, spe_sto in compounds['left'].items():
                    compound = build_compound(
                        spe_id,
                        spe_sto,
                        compounds_cache
                    )
                    rxn.add_reactant(
                        compound=compound,
                        stoichio=spe_sto
                    )
                # Products
                for spe_id, spe_sto in compounds['right'].items():
                    compound = build_compound(
                        spe_id,
                        spe_sto,
                        compounds_cache
                    )
                    rxn.add_product(
                        compound=compound,
                        stoichio=spe_sto
                    )
                # Add at the beginning of the pathway
                # to have the pathway in forward direction
                # Search for the target in the current reaction
                target_id = [spe_id for spe_id in rxn.get_products_ids() if 'TARGET' in spe_id]
                if target_id != []:
                    target_id = target_id[0]
                else:
                    target_id = None
                logger.debug(rxn.to_dict())
                pathway.add_reaction(
                    rxn=rxn,
                    target_id=target_id
                )

            ## SINK
            pathway.set_sink(
                list(
                    set(pathway.get_species_ids()) & set(sink_molecules)
                )
            )

            # RANK AMONG ALL SUB-PATHWAYS OF THE CURRENT MASTER PATHWAY
            res_pathways[path_idx] = apply_to_best_pathways(
                res_pathways[path_idx],
                pathway,
                logger
            )

    # Transform the list of Item into a list of Pathway
    results = {}
    for res_pathway_idx, res_pathway in res_pathways.items():
        results[res_pathway_idx] = [pathway.object for pathway in res_pathway]

    return results

def build_compound(
    spe_id: str,
    spe_sto: int,
    compounds_cache: Dict,
    logger=getLogger(__name__)
) -> Compound:
    # Look if compound has been read from rp2paths compounds
    if spe_id in Cache.get_objects():
        compound = Cache.get(spe_id)
        # Add infos from global cache
        for key in ['name', 'formula']:
            if spe_id in compounds_cache:
                setattr(
                    compound,
                    'set_'+key,
                    compounds_cache[spe_id][key]
                )
    # Else get it from the cache
    else:
        compound = Compound(
            id=spe_id,
            smiles=compounds_cache[spe_id]['smiles'],
            inchi=compounds_cache[spe_id]['inchi'],
            inchikey=compounds_cache[spe_id]['inchikey'],
            formula=compounds_cache[spe_id]['formula'],
            name=compounds_cache[spe_id]['name']
        )
    return compound


def add_compounds(
    compounds: Dict,
    compounds_to_add: Dict,
    logger=getLogger(__name__)
) -> Dict:
    _compounds = dict(compounds)
    for side in ['right', 'left']:
        # added compounds with struct
        for cmpd_id, cmpd in compounds_to_add[side].items():
            if cmpd_id in _compounds[side]:
                _compounds[side][cmpd_id] += cmpd['stoichio']
            else:
                _compounds[side][cmpd_id] = cmpd['stoichio']
        # added compounds with no struct
        for cmpd_id, cmpd in compounds_to_add[side+'_nostruct'].items():
            if cmpd_id in _compounds[side]:
                _compounds[side][cmpd_id] += cmpd['stoichio']
            else:
                _compounds[side][cmpd_id] = cmpd['stoichio']
    return _compounds


def apply_to_best_pathways(
    pathways: List[Dict],
    # max_subpaths_filter: int,
    pathway: rpPathway,
    logger=getLogger(__name__)
) -> List[Dict]:
    '''
    Given a pathway object, looks if an equivalent pathway (cf rpSBML::__eq__ method)
    is present in the given list. If found, then compare scores and keep the highest.
    Otherwise, insert it in the list.

    Parameters
    ----------
    pathways: List[Dict]
        List of pathways sorted by increasing scores
    max_subpaths_filter: int
        Number of top elements to return
    pathway: Dict
        Pathway to insert
    logger : Logger
        The logger object.

    Returns
    -------
    best_pathways: List[Dict]
        List of pathways with highest scores
    '''

    logger.debug('Best pathways:       ' + str([item for item in pathways]))
    # logger.debug('max_subpaths_filter: ' + str(max_subpaths_filter))
    logger.debug('pathway:             ' + str(pathway))

    # Compute the score of applicant pathway
    # sum of rule_scores over reactions
    score = sum(rxn.get_info('rule_score') for rxn in pathway.get_reactions())
    # normalize
    score /= pathway.get_nb_reactions()

    # from bisect import insort as bisect_insort
    # bisect_insort(best_rpsbml, sbml_item)

    # Insert pathway in best_pathways list by increasing score
    pathways = insert_and_or_replace_in_sorted_list(
        Item(pathway, score),
        pathways
    )

    # for item in best_rpsbml:
    #     print(item.rpsbml_obj._get_reactions_with_species_keys())
    # logger.debug(str([item.rpsbml_obj._get_reactions_with_species_keys() for item in best_rpsbml]))

    # Keep only topX
    # best_rpsbml = best_rpsbml[-max_subpaths_filter:]

    return pathways

