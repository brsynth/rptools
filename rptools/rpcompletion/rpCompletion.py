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
from json import dump as json_dump
from colored import fg, bg, attr
from brs_utils import (
    insert_and_or_replace_in_sorted_list,
    Item
)
from rr_cache import rrCache
from rxn_rebuild import rebuild_rxn
from rptools.rplibs import rpSBML

# Allowed ouput formats
FORMATS = {
    'RPSBML': '_sbml.xml',
    'JSON': '.json'
}


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
    outdir,
    out_format='RPSBML',
    upper_flux_bound=999999,
    lower_flux_bound=0,
    max_subpaths_filter=10,
    pathway_id='rp_pathway',
    compartment_id='MNXC3',
    species_group_id='central_species',
    sink_species_group_id='rp_sink_species',
    pubchem_search=False,
    logger=getLogger(__name__)
):

    out_format = out_format.upper()
    if out_format not in FORMATS.keys():
        raise ValueError(
            'Output format {format} is not recognized (choices: {formats})'.format(
                format='\''+out_format+'\'',
                formats=', '.join(['\''+format+'\'' for format in FORMATS.keys()])
            )
        )
        exit(-1)

    if max_subpaths_filter < 0:
        raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

    if os_path.exists(outdir) and os_path.isfile(outdir):
        logger.error('Outdir name '+outdir+' already exists and is actually file. Stopping the process...')
        exit(-1)

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
    compounds_strc = read_compounds(
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
        compounds_strc=compounds_strc,
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
    rp_infos = {
        'pathway_id': 'rp_pathway',
        'compartment_id': 'MNXC3',
        'species_group_id': 'central_species',
        'sink_species_group_id': 'rp_sink_species',
        'upper_flux_bound': upper_flux_bound,
        'lower_flux_bound': lower_flux_bound
    }
    all_pathways = build_all_pathways(
        pathways=pathway_combinatorics,
        transfos=full_transfos,
        sink_molecules=sink_molecules,
        rp_infos=rp_infos,
        rr_reactions=cache.get('rr_reactions'),
        compounds_strc=compounds_strc,
        compounds_cache=cache.get('cid_strc'),
        logger=logger
    )

    # WRITE OUT
    if not os_path.exists(outdir):
        os_mkdir(outdir)
    # Write out only topX sub-pathways per master pathway
    for pathway_id, sub_pathways in all_pathways.items():
        for sub_pathway in sub_pathways[-max_subpaths_filter:]:
            if out_format == 'RPSBML':
                write_to_RPSBML(
                    pathway=sub_pathway.object,
                    cache=cache,
                    outdir=outdir,
                    logger=logger
                )
            elif out_format == 'JSON':
                write_to_JSON(
                    pathway=sub_pathway.object,
                    outdir=outdir,
                    logger=logger
                )


def complete_transformations(
    transfos: Dict,
    ec_numbers: Dict,
    compounds_strc: Dict,
    cache: rrCache,
    logger: Logger = getLogger(__name__)
) -> Dict:

    full_transfos = {}

    # For each transformations read
    for transfo_id, transfo in transfos.items():

        full_transfos[transfo_id] = {}
        full_transfos[transfo_id]['ec'] = ec_numbers[transfo_id]['ec']
        # Convert transformation into SMILES
        transfo_smi = '{left}>>{right}'.format(
                left=build_smiles(transfo['left'], compounds_strc),
                right=build_smiles(transfo['right'], compounds_strc)
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
    compd_strc: Dict,
    logger: Logger = getLogger(__name__)
) -> str:
    return '.'.join([compd_strc[spe]['smiles'] for spe in side.keys()])


## Function to parse the compounds.txt file
#
#  Extract the smile and the structure of each compounds of RP2Path output
#  Method to parse all the RP output compounds.
#
#  @param path The compounds.txt file path
#  @return rp_compounds Dictionnary of smile and structure for each compound
def read_compounds(path, cache, logger=getLogger(__name__)):

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
                rp_strc[row[0]]['inchi'] = cache.get('cid_strc')[row[0]]['inchi']
            except KeyError:
                #try to generate them yourself by converting them directly
                try:
                    resConv = cache._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                    rp_strc[row[0]]['inchi'] = resConv['inchi']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
            try:
                rp_strc[row[0]]['inchikey'] = cache.get('cid_strc')[row[0]]['inchikey']
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
    """Build all combinations of sub-pathways based on these facts:
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

    """

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
    rp_infos: Dict,
    rr_reactions: Dict,
    compounds_strc: Dict,
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
        # best_sub_pathways = []
        res_pathways[path_idx] = []
        for sub_path_idx in range(len(sub_pathways)):

            pathway = {
                'id': str(path_idx).zfill(3)+'_'+str(sub_path_idx+1).zfill(4),
                'species': [],
                'reactions': [],
                'sink': [],
                'global_infos': rp_infos
            }

            ## ITERATE OVER REACTIONS
            for rxn in sub_pathways[sub_path_idx]:

                transfo_id = rxn['transfo_id']
                transfo = transfos[transfo_id]
                rule_id = rxn['rule_id']
                tmpl_rxn_id = rxn['tmpl_rxn_id']

                ## COMPOUNDS
                # Compounds from original transformation
                compounds = {
                    'right': dict(transfo['right']),
                    'left': dict(transfo['left'])
                }
                # Add template reaction compounds
                compounds = add_compounds(
                    compounds,
                    transfo['complement'][rule_id][tmpl_rxn_id]['added_cmpds']
                )

                ## REACTION
                reaction = {
                    'transfo_id': rxn['transfo_id'],
                    'ec': transfo['ec'],
                    'rule_id': rule_id,
                    'rule_score': rr_reactions[rule_id][tmpl_rxn_id]['rule_score'],
                    'tmpl_rxn_id': tmpl_rxn_id,
                    'full_transfo': transfo['complement'][rule_id][tmpl_rxn_id]['full_transfo'],
                    'reactants': dict(compounds['left']),
                    'products': dict(compounds['right'])
                }
                # Add at the beginning of the pathway
                # to have the pathway in forward direction
                pathway['reactions'].insert(0, reaction)

                ## 1/2 SPECIES
                # Handle duplicates later on
                pathway['species'] += (
                    list(reaction['reactants'].keys())
                  + list(reaction['products'].keys())
                )

            ## 2/2 SPECIES
            # Remove duplicates with sets
            species = set(pathway['species'])
            pathway['species'] = {}
            for specie in species:
                # Look if compound has been read from rp2paths compounds
                if specie in compounds_strc:
                    compound = compounds_strc[specie]
                    compound['name'] = specie
                # Else get it from the cache
                else:
                    compound = compounds_cache[specie]
                pathway['species'][specie] = compound

            ## SINK
            pathway['sink'] = list(
                species & set(sink_molecules)
            )

            # RANK AMONG ALL SUB-PATHWAYS OF THE CURRENT MASTER PATHWAY
            res_pathways[path_idx] = apply_to_best_pathways(
                res_pathways[path_idx],
                pathway,
                logger
            )

    return res_pathways


def write_to_JSON(
    pathway: Dict,
    outdir: str,
    logger: Logger = getLogger(__name__)
) -> None:

    out_filename = os_path.join(
        outdir,
        'rp_'+pathway['id']
    ) + FORMATS['JSON']

    with open(out_filename, 'w') as fp:
        json_dump(pathway, fp, indent=4)


def write_to_RPSBML(
    pathway: Dict,
    cache: rrCache,
    outdir: str,
    logger: Logger = getLogger(__name__)
) -> None:

    rpsbml = rpSBML(name='rp_'+pathway['id'], logger=logger)

    ## Create a generic Model, ie the structure and unit definitions that we will use the most
    rpsbml.genericModel(
        'RetroPath_Pathway_'+pathway['id'],
        'RP_model_'+pathway['id'],
        cache.get('comp_xref')[cache.get('deprecatedCompID_compid')[pathway['global_infos']['compartment_id']]],
        pathway['global_infos']['compartment_id'],
        pathway['global_infos']['upper_flux_bound'],
        pathway['global_infos']['lower_flux_bound']
    )

    ## Create the groups (pathway, species, sink species)
    rpsbml = create_rpSBML_groups(
        rpsbml=rpsbml,
        pathway=pathway,
        logger=logger
    )

    ## Add species to the model
    rpsbml, target_id = create_rpSBML_species(
        rpsbml=rpsbml,
        pathway=pathway,
        logger=logger
    )

    ## Add reactions to the model
    rpsbml = create_rpSBML_reactions(
        rpsbml=rpsbml,
        pathway=pathway,
        target_id=target_id,
        logger=logger
    )

    ## Write to file
    rpsbml.writeToFile(
        os_path.join(
            outdir,
            str(rpsbml.modelName)
        ) + FORMATS['RPSBML']
    )


def create_rpSBML_reactions(
    rpsbml: rpSBML,
    pathway: Dict,
    target_id: str,
    logger: Logger = getLogger(__name__)
) -> rpSBML:

    rxn_idx = 0

    for rxn in pathway['reactions']:

        rxn_idx += 1
        rxn_id = 'rxn_' + str(rxn_idx)

        # Format according to what rpSBML method expects
        rxn_rpSBML = {
            'rule_id': rxn['rule_id'],
            'rule_ori_reac': rxn['tmpl_rxn_id'],
            'rule_score': rxn['rule_score'],
            'left': rxn['reactants'],
            'right': rxn['products'],
            'rxn_idx': rxn_idx
            }
        # Add the reaction in the model
        rpsbml.createReaction(
            reac_id='rxn_' + str(rxn_idx),
            fluxUpperBound=pathway['global_infos']['upper_flux_bound'],
            fluxLowerBound=pathway['global_infos']['lower_flux_bound'],
            rxn=rxn_rpSBML,
            compartment_id=pathway['global_infos']['compartment_id'],
            reaction_smiles= rxn['full_transfo'],
            reacXref={'ec': rxn['ec']},
            pathway_id=pathway['global_infos']['pathway_id']
        )

    # Add the consumption of the target
    # Build the target reaction
    # according to what rpSBML method expects
    rxn_target_rpSBML = {
        'rule_id': None,
        'rule_ori_reac': None,
        'rule_score': None,
        'left': { target_id: 1 },
        'right': {},
        'rxn_idx': None
        }
    # Create the reaction in the model
    rpsbml.createReaction(
        reac_id='rxn_target',
        fluxUpperBound=pathway['global_infos']['upper_flux_bound'],
        fluxLowerBound=pathway['global_infos']['lower_flux_bound'],
        rxn=rxn_target_rpSBML,
        compartment_id=pathway['global_infos']['compartment_id']
    )

    return rpsbml


def create_rpSBML_species(
    rpsbml: rpSBML,
    pathway: Dict,
    logger: Logger = getLogger(__name__)
) -> rpSBML:

    target_id = None

    for specie_id, specie in pathway['species'].items():

        if specie_id.startswith('TARGET'):
            target_id = specie_id

        # Handle the sink
        if specie_id in pathway['sink']:
            sink_species_group_id = pathway['global_infos']['sink_species_group_id']
        else:
            sink_species_group_id = None

        rpsbml.createSpecies(
            species_id=specie_id,
            species_name=specie['name'],
            compartment_id=pathway['global_infos']['compartment_id'],
            inchi=specie['inchi'],
            inchikey=specie['inchikey'],
            smiles=specie['smiles'],
            species_group_id=pathway['global_infos']['species_group_id'],
            in_sink_group_id=sink_species_group_id
        )

    return rpsbml, target_id


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
def create_rpSBML_groups(
    rpsbml: rpSBML,
    pathway: Dict,
    logger: Logger = getLogger(__name__)
) -> rpSBML:

    path_id = rpsbml.getName()
    path_base_idx, path_variant_idx = pathway['id'].split('_')
    pathway_id = pathway['global_infos']['pathway_id']

    logger.debug('Create pathway group: '+pathway_id)

    # Create groups
    rpsbml.createGroup(pathway_id)
    rpsbml.createGroup(pathway['global_infos']['species_group_id'])
    rpsbml.createGroup(pathway['global_infos']['sink_species_group_id'])

    # Add pathway id
    rpsbml_dict = rpsbml.toDict(pathway_id)
    # path id
    rpsbml_dict['pathway']['brsynth']['path_id']          = path_id
    # base path (from reaction rules) index
    rpsbml_dict['pathway']['brsynth']['path_base_idx']    = path_base_idx
    # variant path (from real reactions) index
    rpsbml_dict['pathway']['brsynth']['path_variant_idx'] = path_variant_idx

    # Update rpSBML object
    rpsbml.updateBRSynthPathway(rpsbml_dict, pathway_id)

    return rpsbml

    # exit()
    #     pathways_comb = list(map(list,list(itertools_product(*[list(rp2paths_pathways[path_base_idx][step]['reactions'].keys()) for step in rp2paths_pathways[path_base_idx]]))))

    #     path_variant_idx = 1
    #     for path_variant in pathways_comb:

    #         path_id_idx = str(path_base_idx).zfill(3)+'_'+str(path_variant_idx).zfill(4)
    #         path_id     = 'rp_'+path_id_idx

    #         logger.debug(path_id+', '+path_id_idx)

    #         StreamHandler.terminator = ""
    #         logger.info(
    #             '{color}{typo}Complete reactions for{rst} '.format(
    #                 color=fg('white'),
    #                 typo=attr('bold'),
    #                 rst=attr('reset')
    #             )
    #         )
    #         StreamHandler.terminator = "\r"
    #         logger.info(
    #             '{color}pathway {path_id}'.format(
    #                 color=fg('grey_70'),
    #                 path_id=path_id
    #             )
    #         )
    #         # logger.info('pathway '+path_id)
    #         StreamHandler.terminator = "\n"

    #         # 1) Create an rpSBML object with species
    #         rpsbml, species = create_rpSBML(
    #             pathway_id, path_base_idx, path_variant_idx,
    #             path_id, path_id_idx,
    #             rp2paths_pathways[path_base_idx],
    #             cache, compartment_id, upper_flux_bound, lower_flux_bound,
    #             rp_strc, sink_molecules,
    #             species_group_id, sink_species_group_id,
    #             pubchem_search,
    #             logger
    #         )

    #         # 2) Add complete reactions
    #         rpsbml = complete_reactions(
    #             cache,
    #             rpsbml, species, rp2paths_pathways[path_base_idx],
    #             pathway_id, path_variant,
    #             compartment_id, upper_flux_bound, lower_flux_bound,
    #             rp_transfos,
    #             logger
    #         )

    #         # 3) Get the cofactors
    #         rpsbml = add_cofactors(cache, rpsbml, logger=logger)

    #         # for rxn_id in rpsbml.readGroupMembers():
    #         #     reaction = rpsbml.getModel().getReaction(rxn_id)
    #         #     print(rpsbml.getName())
    #         #     print(rpsbml.getScore())
    #         #     print(
    #         #         rpSBML.readBRSYNTHAnnotation(
    #         #             reaction.getAnnotation(),
    #         #             logger=rpsbml.logger
    #         #         )
    #         #     )

    #         # 4) Apply to best rpsbml list
    #         best_rpsbml = apply_best_rpsbml(
    #             best_rpsbml,
    #             max_subpaths_filter,
    #             rpsbml,
    #             path_id,
    #             logger
    #         )

    #         path_variant_idx += 1

    # # Erase the end of line
    # StreamHandler.terminator = ""
    # logger.info('\x1b[2K\r')

    # logger.info(
    #     '{color}{typo}Complete pathway reactions{rst}'.format(
    #         color=fg('white'),
    #         typo=attr('bold'),
    #         rst=attr('reset')
    #     )
    # )
    # StreamHandler.terminator = "\n"
    # logger.info(
    #     '{color}{typo} OK{rst}'.format(
    #         color=fg('green'),
    #         typo=attr('bold'),
    #         rst=attr('reset')
    #     )
    # )

    # # Write topX results to files
    # for rpsbml_item in best_rpsbml[-max_subpaths_filter:]:
    #     rpsbml_item.rpsbml_obj.writeToFile(
    #         os_path.join(
    #             outFolder,
    #             str(rpsbml_item.rpsbml_obj.modelName)
    #         ) + '_sbml.xml'
    #     )

    # return 0


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
    pathway: Dict,
    logger=getLogger(__name__)
) -> List[Dict]:
    """
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
    """

    logger.debug('Best pathways:       ' + str([item for item in pathways]))
    # logger.debug('max_subpaths_filter: ' + str(max_subpaths_filter))
    logger.debug('pathway:             ' + str(pathway))

    # Compute the score of applicant pathway
    # sum of rule_scores over reactions
    score = sum(rxn['rule_score'] for rxn in pathway['reactions'])
    # normalize
    score /= len(pathway['reactions'])

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

