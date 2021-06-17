# import libsbml
# import argparse
# import sys #exit using sys exit if any error is encountered
# import os

# import io
# #import zipfile
# import tarfile
# import glob
# import tempfile
# import shutil

from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from typing import (
    Dict,
    List,
    Tuple
)
from json.decoder import JSONDecodeError
from time import sleep
from equilibrator_api import (
    ComponentContribution,
    Q_
)
from numpy import (
    zeros as np_zeros,
    where as np_where,
    add as np_add,
    arange as np_arange,
    ndarray as np_ndarray
)
from scipy.optimize import linprog
from colored import fg, bg, attr
from brs_utils import (
    print_OK_adv as print_OK,
    print_title_adv as print_title
)
from chemlite import(
    Reaction,
    Compound
)
from rptools.rplibs import rpPathway

# Name of sides with sign
SIDES = [
    {
        'name': 'reactants',
        'sign': -1
    },
    {
        'name': 'products',
        'sign': 1
    },
]


def thermo(
    pathway: rpPathway,
    cc: ComponentContribution = None,
    ph: float = None,
    ionic_strength: float = None,
    pMg: float = None,
    temp_k: float = None,
    logger: Logger = getLogger(__name__)
) -> Dict:
    """Given a tar input file, perform thermodynamics analysis for each rpSBML file.

    :param inFile: The path to the input file
    :param outFile: The path to the output file
    :param pathway_id: The id of the heterologous pathway of interest
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type pathway: Dict
    :type pathway_id: str
    :type ph: float
    :type ionic_strength: float
    :type pMg: float
    :type temp_k: float
    :type logger: Logger

    :rtype: Dict
    :return: Pathway updated with thermodynalics values
    """

    ## REACTIONS
    # reactions = rpsbml.read_reactions(
    #     pathway_id='rp_pathway',
    #     logger=logger
    # )

    # Store thermo values for the net reactions
    # and for each of the reactions within the pathway
    results = {
        'net_reaction': {},
        'reactions': {},
        'species': {}
    }

    print_title(
        txt='Pathway Reactions',
        logger=logger,
        waiting=False
    )
    for rxn in pathway.get_reactions():
        print_reaction(
            rxn=rxn,
            logger=logger
        )
    # logger = getLogger(__name__)

    ## eQuilibrator
    print_title(
        txt='Initialising eQuilibrator...',
        logger=logger,
        waiting=True
    )
    if cc is None:
        cc = initThermo(
            ph,
            ionic_strength,
            pMg,
            temp_k,
            logger
        )
    print_OK(logger)

    # ## COMPOUNDS
    # Get Compounds objects from eQuilibrator cache
    # compound = None if ID does not exist in the cache
    print_title(
        txt='Identifying compounds with cache (eQuilibrator)...',
        logger=logger,
        waiting=True
    )
    species, unk_compounds = get_compounds_from_cache(
        compounds=pathway.get_species(),
        cc=cc,
        logger=logger
    )
    print_OK(logger)

    # Get the formation energy for each compound
    for spe_id, spe in species.items():
        results['species'][spe_id] = {
            'standard_dg_formation': {
                'value': cc.standard_dg_formation(cc.get_compound(spe.get_accession()))[0],  # get .mu
                'units': 'kilojoule / mole'
            }
        }

    ## REACTIONS
    # Compute thermo for each reaction
    for rxn in pathway.get_reactions():
        results['reactions'][rxn.get_id()] = eQuilibrator(
            reaction=rxn,
            species=species,
            cc=cc,
            logger=logger
        )

    ## UNKNOWN COMPOUNDS
    # Remove unknown compounds
    reactions = remove_unknown_compounds(
        unk_compounds=unk_compounds,
        reactions=pathway.get_reactions(),
        rxn_target_id=pathway.get_target_rxn_id(),
        logger=logger
    )

    ## THERMO
    print_title(
        txt='Computing thermodynamics (eQuilibrator)...',
        logger=logger,
        waiting=True
    )
    results['net_reaction'] = eQuilibrator(
        reaction=pathway.net_reaction(),
        species=species,
        cc=cc,
        logger=logger
    )
    print_OK(logger)

    return results

def eQuilibrator(
    reaction: Reaction,
    species: Dict,
    cc: 'ComponentContribution',
    logger: Logger=getLogger(__name__)
) -> Dict:

    measures = {
        'dG0_prime': 'standard_dg_prime',
        'dGm_prime': 'physiological_dg_prime',
        'dG_prime': 'dg_prime',
        'dG': 'standard_dg',
    }

    thermo = {}
    for key in measures.keys():
        thermo[key] = {
            'value': 'NA',
            'error': 'NA',
            'units': 'kilojoule / mole',
        }

    ## Format reaction to what eQuilibrator expects
    compounds = {SIDES[0]['name']: [], SIDES[1]['name']: []}

    try:
        # For both sides left and right
        for cmpd_id, cmpd_sto in reaction.get_left().items():
            compounds[SIDES[0]['name']] += [
                f'{cmpd_sto} {species[cmpd_id].inchi_key}'
            ]
        for cmpd_id, cmpd_sto in reaction.get_right().items():
            compounds[SIDES[1]['name']] += [
                f'{cmpd_sto} {species[cmpd_id].inchi_key}'
            ]
    except KeyError:  # the compound is unknown
        return thermo

    # Join both sides
    rxn_str = '{left} = {right}'.format(
        left=' + '.join(compounds[SIDES[0]['name']]),
        right=' + '.join(compounds[SIDES[1]['name']])
    )

    # Parse formula by eQuilibrator
    rxn = cc.parse_reaction_formula(rxn_str)

    try:
        # Apply each CC method to each required measure
        for key in measures.keys():
            thermo[key] = getattr(cc, measures[key])(rxn)
        thermo = {
            key:{
                'value': float(str(thermo[key].value).split()[0]),
                'error': float(str(thermo[key].error).split()[0]),
                'units': str(thermo[key].units),
            } for key in thermo.keys()
        }
    except Exception as e:
        logger.error(e)

    return thermo
    # {
    #     'dG0_prime': {
    #         'value': float(str(dG0_prime.value).split()[0]),
    #         'error': float(str(dG0_prime.error).split()[0]),
    #         'units': str(dG0_prime.units),
    #     },
    #     'dGm_prime': {
    #         'value': float(str(dGm_prime.value).split()[0]),
    #         'error': float(str(dGm_prime.error).split()[0]),
    #         'units': str(dGm_prime.units),
    #     },
    #     'dG_prime': {
    #         'value': float(str(dG_prime.value).split()[0]),
    #         'error': float(str(dG_prime.error).split()[0]),
    #         'units': str(dG_prime.units),
    #     }
    # }


# def net_reaction(
#     reactions: List[Reaction],
#     logger: Logger=getLogger(__name__)
# ) -> Dict:
#     '''
#     '''

#     # SUM ALL SPECIES
#     species = {}
#     for rxn in reactions:
#         # For both sides left and right
#         for spe_id, spe_sto in rxn.get_left().items():
#             if spe_id in species:
#                 species[spe_id] -= spe_sto
#             else:
#                 species[spe_id] = -spe_sto
#         for spe_id, spe_sto in rxn.get_right().items():
#             if spe_id in species:
#                 species[spe_id] += spe_sto
#             else:
#                 species[spe_id] = spe_sto
#         # for side in SIDES:
#         #     for spe_id, spe_sto in getattr(rxn, 'get_'+side['name']).items():
#         #         if spe_id in species:
#         #             species[spe_id] += side['sign']*spe_sto
#         #         else:
#         #             species[spe_id] = side['sign']*spe_sto
#     # WRITE INTO REACTIONS
#     net_reaction = {
#         SIDES[0]['name']: {},
#         SIDES[1]['name']: {}
#     }
#     for spe_id, spe_sto in species.items():
#         # Ignore compounds with stochio = 0
#         if spe_sto < 0:
#             net_reaction[SIDES[0]['name']][spe_id] = SIDES[0]['sign']*spe_sto
#         elif spe_sto > 0:
#             net_reaction[SIDES[1]['name']][spe_id] = SIDES[1]['sign']*spe_sto

#     return net_reaction


def remove_unknown_compounds(
    reactions: List[Reaction],
    rxn_target_id: str,
    unk_compounds: List = [],
    logger: Logger=getLogger(__name__)
) -> Dict:
    '''Try to remove compounds that are unknown in the eQuilibrator cache from reaction set (pathway).
    For this purpose, a stoichiometric matrix is built with only coefficients of these compounds.
    Then, try to solve as a linear equations system and apply new coeffs to reactions
    '''

    # # From compounds from cache, get those which are None (unknown)
    # unk_compounds = list(
    #     dict(
    #         filter(
    #             lambda elem: elem[1] is None,
    #             compounds.items()
    #         )
    #     ).keys()
    # )

    # unk_compounds = ['CMPD_0000000003', 'CMPD_0000000010', 'CMPD_0000000025']

    if unk_compounds == []:
        return reactions

    # # If the target is unknown in eQuilibrator cache,
    # # then stop the thermo
    # result = filter(lambda x: x.startswith('TARGET'), unk_compounds)
    # print(list(result))

    # Build the stoichio matrix for unknown compounds
    ## S
    sto_mat = build_stoichio_matrix(
        reactions=reactions,
        compounds=unk_compounds,
        logger=logger
    )
    print(sto_mat)

    # Get the target reaction to maximize
    rxn_target_idx = get_target_rxn_idx(
        reactions,
        rxn_target_id,
        logger
    )

    coeffs = minimize(
        sto_mat,
        rxn_target_idx,
        logger
    )

    print(coeffs)

    ## Impact coeff to reactions
    for rxn_idx in range(len(reactions)):
        if coeffs[rxn_idx] != 0:
            reactions[rxn_idx].mult_stoichio_coeff(coeffs[rxn_idx])

    return reactions


def get_target_rxn_idx(
    reactions: List[Reaction],
    rxn_target_id: str,
    logger: Logger=getLogger(__name__)
) -> int:
    for rxn_idx in range(len(reactions)):
        rxn = reactions[rxn_idx]
        if rxn_target_id == rxn.get_id():
            return rxn_idx
    return None


def minimize(
    S: 'np_ndarray',
    rxn_target_idx: int,
    logger: Logger=getLogger(__name__)
) -> 'np_ndarray':

    nb_compounds, nb_reactions = S.shape

    # S[0][2] = 2

    ## S.x = 0
    # Init with zeros
    b_eq = np_zeros(nb_compounds)

    # 0 <= x <= 1
    # Init with (0, 1)
    # bounds = [(0, 1)]*len(reactions)
    # the same
    # lower bound (lb) and upper bound (ub) will be applied
    # to all variables.
    bounds = (0, 1)

    ## Max rxn_x
    c = np_zeros(nb_reactions)
    # Search the rxn ID that
    # maximizes the production of the target

    # translate it in the stoichio matrix index
    c[rxn_target_idx] = 0

    ## Solve
    res = linprog(
        c,
        A_eq=S, 
        b_eq=b_eq, 
        bounds=bounds, 
        method='simplex'
    )

    return res.x


def build_stoichio_matrix(
    reactions: List[Reaction],
    compounds: List[str] = [],
    logger: Logger = getLogger(__name__)
) -> 'np_ndarray':
    '''Build the stoichiometric matrix of reactions.
       If compounds is not None, then fill the matrix only for these
    '''

    ## Build list of compounds to put in the matrix
    # If compounds not passed in arg,
    # then detect them from reactions
    # else only put in the matrix compounds contained in 'compounds'
    if compounds == []:
        species = [rxn.get_species_ids() for rxn in reactions]
        _compounds = list(
            set(
                [
                    # spe_id for rxn in reactions for side in SIDES for spe_id in list(rxn[side['name']].keys()) 
                    spe_id for sub_species in species for spe_id in sub_species
                ]
            )
        )
    else:
        _compounds = list(compounds)

    ## Init the stoichio matrix
    sto_mat = np_zeros(
        (
            len(_compounds),  # rows
            len(reactions)  # columns
        )
    )

    ## Fill up the matrix
    # For each compound...
    for cmpd_idx in range(len(_compounds)):
        cmpd_id = _compounds[cmpd_idx]
        # look where it appears in each reaction
        # For each reaction
        for rxn_idx in range(len(reactions)):
            rxn = reactions[rxn_idx]
            if cmpd_id in rxn.get_reactants_ids():
                # Put the stochio value in the stochio matrix
                sto_mat[cmpd_idx][rxn_idx] = -rxn.get_reactants_stoichio()[cmpd_id]
            elif cmpd_id in rxn.get_products_ids():
                # Put the stochio value in the stochio matrix
                sto_mat[cmpd_idx][rxn_idx] = rxn.get_products_stoichio()[cmpd_id]

    return sto_mat


def get_compounds_from_cache(
    compounds: List[Compound],
    cc: ComponentContribution,
    logger: Logger=getLogger(__name__)
) -> Dict:
    """Get compounds accession from the cache
    Compounds are None if not found

    Parameters
    ----------
    compounds : Dict
        Compounds to check the existence in the cache
    cc : ComponentContribution
        ComponentContribution
    logger : Logger

    Returns
    -------
    Dictionary of compounds

    """
    compounds_dict = {}
    unknown_compounds = []

    for cmpd in compounds:
        compound = cc.get_compound(cmpd.get_id())
        # If ID not found,
        # then search with inchikey
        try:
            if compound is None:
                compound = cc.search_compound(cmpd.get_inchikey())
            # If the first level of inchikeys are the same,
            # then substitute
            if compound.inchi_key.split('-')[0] == cmpd.get_inchikey().split('-')[0]:
                compounds_dict[cmpd.get_id()] = compound
            else:
                # search by inchikey
                try:
                    compounds_dict[cmpd.get_id()] = cc.search_compound_by_inchi_key(cmpd.get_inchikey())[0]
                except IndexError:  # the compound is considered as unknown
                    unknown_compounds += [cmpd.get_id()]
        except (AttributeError, ValueError):
            try:
                compound = cc.get_compound_by_inchi(cmpd.get_inchi())
            except ValueError:
                unknown_compounds += [cmpd.get_id()]

        if compound is None:
            unknown_compounds += [cmpd.get_id()]

    return compounds_dict, unknown_compounds


# def net_reaction(
#     rpsbml: "rpSBML",
#     pathway_id: str='rp_pathway',
#     logger: Logger=getLogger(__name__)
# ):

#     rp_pathway = rpsbml.getModel().getPlugin('groups').getGroup(pathway_id)

#     reactions = list(
#         reversed(
#             [rpsbml.getModel().getReaction(i.getIdRef()) for i in rpsbml.getModel().getPlugin('groups').getGroup(pathway_id).getListOfMembers()]
#         )
#     )

#     # Build the stoichiometric matrix
#     stoichio_matrix, species_table = build_stoichio_matrix(
#         reactions,
#         logger
#     )

#     # print(stoichio_matrix)
#     # print(species_table)
#     # exit()

#     # Print out reactions
#     print_reactions(
#         reactions,
#         stoichio_matrix,
#         species_table,
#         logger
#     )

#     # Scan all species to find those that we want to remove from left part,
#     # i.e. start with 'CMPD_'
#     species_to_rm = []
#     for spe in species_table.keys():
#         if spe.startswith('CMPD_'):
#             species_to_rm += [species_table[spe]]
#         #     prio_rm_spe.insert(0, species_table[spe])
#         # else:
#         #     prio_rm_spe.append(species_table[spe])
#     # print(species_table)
#     # # prio_rm_spe = []
#     # stoichio_matrix[5][1] = -2
#     # stoichio_matrix[4][1] = -0.5
#     # from numpy import array
#     # stoichio_matrix = array([
#     #     [-2., 1., 2.],
#     #     [0., -1., 1.],
#     # ])
#     # # CONFLICT!
#     # stoichio_matrix = array([
#     #     [0., 1., -2.],
#     #     [0., -1., 1.],
#     # ])
#     # stoichio_matrix = array([
#     #     [2., -1., -2.],
#     #     [0., -1., 1.],
#     # ])
#     # prio_rm_spe = [1, 0]
#     # # stoichio_matrix[4][0] = -2
#     # # Print out reactions
#     # for i_rxn in range(len(reactions)):
#     #     rxn_id = reactions[i_rxn].getIdAttribute()
#     #     reactants, products = build_reaction(
#     #         rxn_id = rxn_id,
#     #         sto_mat = stoichio_matrix.transpose()[i_rxn],
#     #         spe_tbl = species_table,
#     #         logger = logger
#     #     )
#     #     print_reaction(
#     #         rxn_id = rxn_id,
#     #         reactants = reactants,
#     #         products = products,
#     #         logger = logger
#     #     )
#     # Remove intermediate species
#     stoichio_matrix = remove_species(
#         sto_mat = stoichio_matrix,
#         species = species_to_rm,
#         logger = logger
#     )

#     # Print out the net reaction
#     # For each species:
#     #   sum neg coeff between them (reactants)
#     #   sum pos coeff between them (products)
#     sto_rxn = np_add(
#         [coeff[coeff<0].sum() for coeff in stoichio_matrix],
#         [coeff[coeff>0].sum() for coeff in stoichio_matrix]
#     )
#     rxn_id = 'rxn_net'
#     reactants, products = build_reaction(
#         rxn_id = rxn_id,
#         sto_mat = sto_rxn,
#         spe_tbl = species_table,
#         logger = logger
#     )
#     logger.info(
#         '{color}{typo}Net Reaction{rst}'.format(
#             color=c_fg('white'),
#             typo=c_attr('bold'),
#             rst=c_attr('reset')
#         )
#     )
#     print_reaction(
#         rxn_id = rxn_id,
#         reactants = reactants,
#         products = products,
#         logger = logger
#     )


#     exit()

#     # Detect intermediate species
#     intermediate_species = detect_intermediate_species(pathway_reactions, logger)

#     # Adjust stoichiometric coeff of intermediate species so that we can remove them
#     pathway_reactions = adjust_stoichio(pathway_reactions)

#     logger.info('intermediate_species = ' + str(intermediate_species))


# def print_reactions(
#     reactions: List,
#     sto_mat: 'numpy.ndarray',
#     spe_tbl = Dict,
#     logger: Logger=getLogger(__name__)
# ) -> None:
#     logger.info(
#         '{color}{typo}Reactions{rst}'.format(
#             color=c_fg('white'),
#             typo=c_attr('bold'),
#             rst=c_attr('reset')
#         )
#     )
#     for i_rxn in range(len(reactions)):
#         rxn_id = reactions[i_rxn].getIdAttribute()
#         reactants, products = build_reaction(
#             rxn_id = rxn_id,
#             sto_mat = sto_mat.transpose()[i_rxn],
#             spe_tbl = spe_tbl,
#             logger = logger
#         )
#         print_reaction(
#             rxn_id = rxn_id,
#             reactants = reactants,
#             products = products,
#             logger = logger
#         )


# def build_reaction(
#     rxn_id = str,
#     sto_mat = 'numpy.ndarray',
#     spe_tbl = Dict,
#     logger: Logger=getLogger(__name__)
# ) -> Tuple[Dict, Dict]:
#     """
#     Returns reactants and products ids.

#     Parameters
#     ----------
#     rxn_id: str
#         Id of the reaction
#     sto_mat: 'numpy.ndarray'
#         Stoichiometrix matrix of the reaction
#     spe_tbl: Dict
#         Correspondance between species id and and its index within sto_mat

#     Returns
#     -------
#     reactants: Dict
#         Dictionary of reactants: {spe_id: sto_coeff}
#     products: Dict
#         Dictionary of products: {spe_id: sto_coeff}
#     """

#     reactants = {}
#     products = {}

#     for i_sto_coeff in range(len(sto_mat)):
#         sto_coeff = sto_mat[i_sto_coeff]
#         # Retieve the species id corresponding to the current coeff index
#         spe_id = list(spe_tbl.keys())[list(spe_tbl.values()).index(i_sto_coeff)]
#         if sto_coeff < 0: # reactants
#             reactants[spe_id] = sto_coeff
#         elif sto_coeff > 0: # products
#             products[spe_id] = sto_coeff

#     return reactants, products


def print_reaction(
    rxn: Reaction,
    logger: Logger=getLogger(__name__)
) -> None:
    """
    Print out the reaction.

    Parameters
    ----------
    rxn_id: str
        Id of the reaction
    reactants: Dict
        Dictionary of reactants: {spe_id: sto_coeff}
    products: Dict
        Dictionary of products: {spe_id: sto_coeff}
    """
    logger.info(
        '{color}{typo}   |- {rxn_id}: {rst}{reactants} --> {products}'.format(
            rxn_id=rxn.get_id(),
            reactants=' + '.join([str(sto)+' '+str(spe) for spe,sto in rxn.get_reactants_stoichio().items()]),
            products=' + '.join([str(sto)+' '+str(spe) for spe,sto in rxn.get_products_stoichio().items()]),
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )


# def build_stoichio_matrix(
#     reactions: List,
#     logger: Logger=getLogger(__name__)
# ) -> Tuple['numpy.ndarray', Dict]:
#     """
#     Returns the stoichiometric matrix of the equations system.

#     Parameters
#     ----------
#     reactions: List
#         List of reactions (SBML Reaction)

#     Returns
#     -------
#     sto_mat: 'numpy.ndarray'
#         Stoichiometrix matrix
#     species_table: Dict
#         Correspondance between species ids and row number in sto_mat
#     """

#     # Put in a dictionary all species of the system to have:
#     #   1. the dimension of stoichiometrix matrix
#     #   2. the correspondance between species ids and their indices in the stoichiometrix matrix
#     species_table = {}
#     for react in reactions:
#         for spe in list(react.getListOfReactants()) + list(react.getListOfProducts()):
#             if spe.getSpecies() not in species_table.keys():
#                 species_table[spe.getSpecies()] = len(species_table.keys())

#     # Init the stoichio matrix
#     sto_mat = np_zeros(
#         (
#             len(species_table), # rows
#             len(reactions) # columns
#         )
#     )

#     for i_react in range(len(reactions)):

#         rxn = reactions[i_react]

#         # Add stoichio coeff for reactants
#         for spe in rxn.getListOfReactants():
#             sto_mat[species_table[spe.getSpecies()]][i_react] = -spe.getStoichiometry()

#         # Add stoichio coeff for products
#         for spe in rxn.getListOfProducts():
#             sto_mat[species_table[spe.getSpecies()]][i_react] = spe.getStoichiometry()
    
#     logger.debug(sto_mat)

#     return sto_mat, species_table


# def remove_species(
#     sto_mat: 'numpy.ndarray',
#     species: List,
#     logger: Logger=getLogger(__name__)
# ) -> 'numpy.ndarray':
#     """
#     Returns stoichio matrix after having removed species.

#     Parameters
#     ----------
#     sto_mat: 'numpy.ndarray'
#         Stoichiometrix matrix
#     prio_rm_spe: List
#         Species to remove in priority

#     Returns
#     -------
#     sto_mat: 'numpy.ndarray'
#         Reduced stoichiometrix matrix
#     """

#     # if prio_rm_spe == []:
#     #     prio_rm_spe = range(sto_mat.shape[0])
#     # else:
#     #     # Species to be removed in priority have to be processed last
#     #     # so that coeff mult will not be changed thereafter
#     #     prio_rm_spe.reverse()

#     # Try to balance species that appears both in left and right sides of the equations system
#     print(sto_mat)
#     for i_spe in species:
#         spe = sto_mat[i_spe]
#         # Check if there are + and - coeff, otherwise it is not possible to balance
#         pos_coeff = spe[spe > 0]
#         neg_coeff = spe[spe < 0]
#         if len(pos_coeff) > 0 and len(neg_coeff) > 0:
#             sum_pos_coeff = pos_coeff.sum()
#             sum_neg_coeff = abs(neg_coeff.sum())
#             # Choose to multiply all neg coeffs by the following ratio
#             # (we could do the inverse, i.e. mult all pos coeff by sum_neg / sum_pos)
#             mult = sum_pos_coeff / sum_neg_coeff
#             # Try to reduce reactants coeff
#             if mult < 1: # more neg coeff (reactants) than pos (products)
#                 # then reduce neg coeff
#                 rxn_ids = np_where(spe > 0)
#                 mult = 1. / mult
#             else: # more or equal pos coeff (products) than neg (reactants)
#                 # then reduce pos coeff
#                 rxn_ids = np_where(spe < 0)
#             # Multiply coeffs of all species in the reaction where the species to balanced has a negative coeff
#             for rxn_id in rxn_ids:
#                 sto_mat.transpose()[rxn_id] *= mult
#         print(species)
#         print(sto_mat)

#     # Remove species that are overall balanced
#     # Has to be done after the balancing because balancing has to be done on all species
#     for i_spe in range(sto_mat.shape[0]):
#         spe = sto_mat[i_spe]
#         # If sum of coeffs is 0
#         if spe.sum() == 0:
#             # then remove species from overall equation by setting coeffs to 0
#             spe.fill(0)

#     return sto_mat


# def detect_intermediate_species(
#     reactions: Dict,
#     logger: Logger=getLogger(__name__)
# ) -> List:
#     """
#     Returns list of species that appear both in left (reactants) and right (products) sides.

#     Parameters
#     ----------
#     reactions: Dict
#         Dictionary of reactions

#     Returns
#     -------
#     intermediate_species: List[str]
#         List of intermediate species names
#     """

#     logger.debug('reactions: ' + str(reactions))

#     reactants = products = set()

#     for rxn in reactions.keys():
#         reactants = reactants.union(set(reactions[rxn]['reactants']))
#         products = products.union(set(reactions[rxn]['products']))
    
#     return list(reactants & products)


# used to initialise and download the data for equilibrator
def initThermo(
    ph: float,
    ionic_strength: float,
    pMg: float,
    temp_k: float,
    logger: Logger=getLogger(__name__)
) -> ComponentContribution:

    # cc = None
    # while cc is None:
        # try:
    cc = ComponentContribution()
        # except JSONDecodeError:
        #     logger.warning('Waiting for zenodo.org... Retrying in 5s')
        #     sleep(5)

    if ph is not None:
        cc.p_h = Q_(ph)
    if ionic_strength is not None:
        cc.ionic_strength = Q_(str(ionic_strength)+' mM')
    if pMg is not None:
        cc.p_mg = Q_(pMg)
    if temp_k is not None:
        cc.temperature = Q_(str(temp_k)+' K')

    return cc


# def runMDF_hdd(inputTar, outputTar, pathway_id='rp_pathway', thermo_id='dfG_prime_o', fba_id='fba_obj_fraction', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
#     """Given a tar input file, perform MDF analysis for each rpSBML file.

#     :param inputTar: The path to the input TAR file
#     :param outputTar: The path to the output TAR file
#     :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
#     :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
#     :param fba_id: The id of the FBA value (Default: fba_obj_fraction)
#     :param ph: The pH of the host organism (Default: 7.0)
#     :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
#     :param pMg: The pMg of the host organism (Default: 10.0)
#     :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
#     :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

#     :type inputTar: str
#     :type outputTar: str
#     :type pathway_id: str
#     :type tmpOutputFolder: str
#     :type ph: float
#     :type ionic_strength: float
#     :type pMg: float
#     :type temp_k: float
#     :type stdev_factor: float

#     :rtype: bool
#     :return: Success or failure of the function
#     """
#     rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
#     with tempfile.TemporaryDirectory() as tmpInputFolder:
#         with tempfile.TemporaryDirectory() as tmpOutputFolder:
#             tar = tarfile.open(inputTar, mode='r')
#             tar.extractall(path=tmpInputFolder)
#             tar.close()
#             if len(glob.glob(tmpInputFolder+'/*'))==0:
#                 logging.error('Input file is empty')
#                 return False
#             for sbml_path in glob.glob(tmpInputFolder+'/*'):
#                 logging.debug('=========== '+str(sbml_path)+' ============')
#                 fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
#                 rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
#                 rpequilibrator.rpsbml = rpsbml
#                 res = rpequilibrator.MDF(pathway_id, thermo_id, fba_id, stdev_factor, True) #ignore the results since written to SBML file
#                 #ignore res since we are passing write to SBML
#                 rpsbml.writeSBML(tmpOutputFolder)
#                 rpsbml = None
#             if len(glob.glob(tmpOutputFolder+'/*'))==0:
#                 logging.error('rpThermo has not produced any results')
#                 return False
#             with tarfile.open(outputTar, mode='w:gz') as ot:
#                 for sbml_path in glob.glob(tmpOutputFolder+'/*'):
#                     fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', ''))
#                     fileName += '.sbml.xml'
#                     info = tarfile.TarInfo(fileName)
#                     info.size = os.path.getsize(sbml_path)
#                     ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
#     return True


# def runEqSBtab_hdd(inputTar, outputTar, pathway_id='rp_pathway', fba_id=None, thermo_id='dfG_prime_o', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
#     """Given a tar input file, perform MDF analysis for each rpSBML file.

#     :param inputTar: The path to the input TAR file
#     :param outputTar: The path to the output TAR file
#     :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
#     :param fba_id: The id of the FBA value. Default sets all FBA values to 1.0 and if specified (Default: None)
#     :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
#     :param ph: The pH of the host organism (Default: 7.0)
#     :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
#     :param pMg: The pMg of the host organism (Default: 10.0)
#     :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
#     :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

#     :type inputTar: str
#     :type outputTar: str
#     :type pathway_id: str
#     :type tmpOutputFolder: str
#     :type ph: float
#     :type ionic_strength: float
#     :type pMg: float
#     :type temp_k: float
#     :type stdev_factor: float

#     :rtype: bool
#     :return: Success or failure of the function
#     """
#     rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
#     with tempfile.TemporaryDirectory() as tmpInputFolder:
#         with tempfile.TemporaryDirectory() as tmpOutputFolder:
#             tar = tarfile.open(inputTar, mode='r')
#             tar.extractall(path=tmpInputFolder)
#             tar.close()
#             if len(glob.glob(tmpInputFolder+'/*'))==0:
#                 logging.error('Input file is empty')
#                 return False
#             for sbml_path in glob.glob(tmpInputFolder+'/*'):
#                 logging.debug('=========== '+str(sbml_path)+' ============')
#                 fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
#                 rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
#                 rpequilibrator.rpsbml = rpsbml
#                 status = rpequilibrator.toNetworkSBtab(os.path.join(tmpOutputFolder, fileName+'.tsv'), pathway_id, thermo_id, fba_id, stdev_factor)
#                 rpsbml = None
#             if len(glob.glob(tmpOutputFolder+'/*'))==0:
#                 logging.error('rpThermo has not produced any results')
#                 return False
#             with tarfile.open(outputTar, mode='w:gz') as ot:
#                 for sbml_path in glob.glob(tmpOutputFolder+'/*'):
#                     fileName = str(sbml_path.split('/')[-1])
#                     info = tarfile.TarInfo(fileName)
#                     info.size = os.path.getsize(sbml_path)
#                     ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
#     return True
