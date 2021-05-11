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
    getLogger
)
from typing import (
    Dict,
    List,
    Tuple
)
from equilibrator_api import ComponentContribution, Q_
from rptools.rplibs import rpSBML
from rptools.rpthermo.rpEquilibrator import pathway, MDF
from numpy import (
    zeros as np_zeros,
    where as np_where,
    add as np_add,
    arange as np_arange,
)


def runThermo(
    inFile: str,
    outFile: str,
    pathway_id: str = 'rp_pathway',
    ph: float = 7.0,
    ionic_strength: float = 200.0,
    pMg: float = 10.0,
    temp_k: float = 298.15,
    logger: Logger = getLogger(__name__)
) -> None:
    """Given a tar input file, perform thermodynamics analysis for each rpSBML file.

    :param inFile: The path to the input file
    :param outFile: The path to the output file
    :param pathway_id: The id of the heterologous pathway of interest
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type inFile: str
    :type outFile: str
    :type pathway_id: str
    :type tmpOutputFolder: str
    :type ph: float
    :type ionic_strength: float
    :type pMg: float
    :type temp_k: float
    :type stdev_factor: float

    :rtype: bool
    :return: Success or failure of the function
    """

    # cc = initThermo(ph, ionic_strength, pMg, temp_k)

# get_compound('MNXM736240@MNXD1')

    # print(cc.ccache.get_compound_by_inchi('1S/H4O7P2/c1-8(2,3)7-9(4,5)6/h(H2,1,2,3)(H2,4,5,6)'))
    # exit()

    rpsbml = rpSBML(inFile)

    net_reaction(
        rpsbml=rpsbml,
        pathway_id=pathway_id,
        logger=logger
    )
    exit()

    pathway(
        rpsbml=rpsbml,
        cc=cc,
        pathway_id=pathway_id,
        update_rpsbml=True,
        logger=logger
    ) # ignore the results since written to SBML file
 
    MDF(
        rpsbml=rpsbml,
        cc=cc,
        logger=logger
    )

    rpsbml.writeToFile(outFile)

    # #mnx_default_conc = json.load(open('data/mnx_default_conc.json', 'r'))
    # mnx_default_conc = json.load(open(os.path.join(os.path.dirname(os.path.abspath( __file__ )), 'data', 'mnx_default_conc.json'), 'r'))


def net_reaction(
    rpsbml: "rpSBML",
    pathway_id: str='rp_pathway',
    logger: Logger=getLogger(__name__)
):

    rp_pathway = rpsbml.getModel().getPlugin('groups').getGroup(pathway_id)

    reactions = [rpsbml.getModel().getReaction(i.getIdRef()) for i in rpsbml.getModel().getPlugin('groups').getGroup(pathway_id).getListOfMembers()]

    # Build the stoichiometric matrix
    stoichio_matrix, species_table = build_stoichio_matrix(
        reactions,
        logger
    )

    # Print out reactions
    for i_rxn in range(len(reactions)):
        rxn_id = reactions[i_rxn].getIdAttribute()
        reactants, products = build_reaction(
            rxn_id = rxn_id,
            sto_mat = stoichio_matrix.transpose()[i_rxn],
            spe_tbl = species_table,
            logger = logger
        )
        print_reaction(
            rxn_id = rxn_id,
            reactants = reactants,
            products = products,
            logger = logger
        )

    # Scan all species to find those that we want to remove from left part,
    # i.e. start with 'CMPD_'
    prio_rm_spe = []
    for spe in species_table.keys():
        if spe.startswith('CMPD_'):
            prio_rm_spe.insert(0, species_table[spe])
        else:
            prio_rm_spe.append(species_table[spe])

    print(species_table)
    # prio_rm_spe = []
    stoichio_matrix[5][1] = -2
    stoichio_matrix[4][1] = -0.5
    from numpy import array
    stoichio_matrix = array([
        [-2., 1., 2.],
        [0., -1., 1.],
    ])
    # CONFLICT!
    stoichio_matrix = array([
        [0., 1., -2.],
        [0., -1., 1.],
    ])
    stoichio_matrix = array([
        [2., -1., -2.],
        [0., -1., 1.],
    ])
    prio_rm_spe = [1, 0]
    # stoichio_matrix[4][0] = -2
    # Print out reactions
    for i_rxn in range(len(reactions)):
        rxn_id = reactions[i_rxn].getIdAttribute()
        reactants, products = build_reaction(
            rxn_id = rxn_id,
            sto_mat = stoichio_matrix.transpose()[i_rxn],
            spe_tbl = species_table,
            logger = logger
        )
        print_reaction(
            rxn_id = rxn_id,
            reactants = reactants,
            products = products,
            logger = logger
        )
    # Remove intermediate species
    stoichio_matrix = remove_intermediate_species(
        sto_mat = stoichio_matrix,
        prio_rm_spe = prio_rm_spe,
        logger = logger
    )

    # Print out the net reaction
    # For each species:
    #   sum neg coeff between them (reactants)
    #   sum pos coeff between them (products)
    sto_rxn = np_add(
        [coeff[coeff<0].sum() for coeff in stoichio_matrix],
        [coeff[coeff>0].sum() for coeff in stoichio_matrix]
    )
    rxn_id = 'Net Reaction'
    reactants, products = build_reaction(
        rxn_id = rxn_id,
        sto_mat = sto_rxn,
        spe_tbl = species_table,
        logger = logger
    )
    print_reaction(
        rxn_id = rxn_id,
        reactants = reactants,
        products = products,
        logger = logger
    )


    exit()

    # Detect intermediate species
    intermediate_species = detect_intermediate_species(pathway_reactions, logger)

    # Adjust stoichiometric coeff ofr intermediate species so that we can remove them
    pathway_reactions = adjust_stoichio(pathway_reactions)

    logger.info('intermediate_species = ' + str(intermediate_species))


def build_reaction(
    rxn_id = str,
    sto_mat = 'numpy.ndarray',
    spe_tbl = Dict,
    logger: Logger=getLogger(__name__)
) -> Tuple[Dict, Dict]:
    """
    Returns reactants and products ids.

    Parameters
    ----------
    rxn_id: str
        Id of the reaction
    sto_mat: 'numpy.ndarray'
        Stoichiometrix matrix of the reaction
    spe_tbl: Dict
        Correspondance between species id and and its index within sto_mat

    Returns
    -------
    reactants: Dict
        Dictionary of reactants: {spe_id: sto_coeff}
    products: Dict
        Dictionary of products: {spe_id: sto_coeff}
    """

    reactants = {}
    products = {}

    for i_sto_coeff in range(len(sto_mat)):
        sto_coeff = sto_mat[i_sto_coeff]
        # Retieve the species id corresponding to the current coeff index
        spe_id = list(spe_tbl.keys())[list(spe_tbl.values()).index(i_sto_coeff)]
        if sto_coeff < 0: # reactants
            reactants[spe_id] = sto_coeff
        elif sto_coeff > 0: # products
            products[spe_id] = sto_coeff

    return reactants, products


def print_reaction(
    rxn_id: str,
    reactants: Dict,
    products: Dict,
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
            rxn_id + ': ' \
        + ' + '.join([str(-sto)+' '+str(spe) for spe,sto in reactants.items()]) \
        + ' --> ' \
        + ' + '.join([str(sto)+' '+str(spe) for spe,sto in products.items()]) \
    )


def build_stoichio_matrix(
    reactions: List,
    logger: Logger=getLogger(__name__)
) -> Tuple['numpy.ndarray', Dict]:
    """
    Returns the stoichiometric matrix of the equations system.

    Parameters
    ----------
    reactions: List
        List of reactions (SBML Reaction)

    Returns
    -------
    sto_mat: 'numpy.ndarray'
        Stoichiometrix matrix
    species_table: Dict
        Correspondance between species ids and row number in sto_mat
    """

    # Put in a dictionary all species of the system to have:
    #   1. the dimension of stoichiometrix matrix
    #   2. the correspondance between species ids and their indices in the stoichiometrix matrix
    species_table = {}
    for react in reactions:
        for spe in list(react.getListOfReactants()) + list(react.getListOfProducts()):
            if spe.getSpecies() not in species_table.keys():
                species_table[spe.getSpecies()] = len(species_table.keys())

    # Init the stoichio matrix
    sto_mat = np_zeros(
        (
            len(species_table), # rows
            len(reactions) # columns
        )
    )

    for i_react in range(len(reactions)):

        rxn = reactions[i_react]

        # Add stoichio coeff for reactants
        for spe in rxn.getListOfReactants():
            sto_mat[species_table[spe.getSpecies()]][i_react] = -spe.getStoichiometry()

        # Add stoichio coeff for products
        for spe in rxn.getListOfProducts():
            sto_mat[species_table[spe.getSpecies()]][i_react] = spe.getStoichiometry()
    
    logger.debug(sto_mat)

    return sto_mat, species_table


def remove_intermediate_species(
    sto_mat: 'numpy.ndarray',
    prio_rm_spe: List = [],
    logger: Logger=getLogger(__name__)
) -> 'numpy.ndarray':
    """
    Returns stoichio matrix after having remove intermediate species.

    Parameters
    ----------
    sto_mat: 'numpy.ndarray'
        Stoichiometrix matrix
    prio_rm_spe: List
        Species to remove in priority

    Returns
    -------
    sto_mat: 'numpy.ndarray'
        Reduced stoichiometrix matrix
    """

    if prio_rm_spe == []:
        prio_rm_spe = range(sto_mat.shape[0])
    else:
        # Species to be removed in priority have to be processed last
        # so that coeff mult will not be changed thereafter
        prio_rm_spe.reverse()

    # Try to balance species that appears both in left and right sides of the equations system
    for i_spe in prio_rm_spe:
        print(sto_mat)

        print(i_spe)
        spe = sto_mat[i_spe]
        # Check if there are + and - coeff, otherwise it is not possible to balance
        pos_coeff = spe[spe > 0]
        neg_coeff = spe[spe < 0]
        if len(pos_coeff) > 0 and len(neg_coeff) > 0:
            sum_pos_coeff = pos_coeff.sum()
            sum_neg_coeff = abs(neg_coeff.sum())
            # Choose to multiply all neg coeffs by the following ratio
            # (we could do the inverse, i.e. mult all pos coeff by sum_neg / sum_pos)
            mult = sum_pos_coeff / sum_neg_coeff
            # Try to reduce reactants coeff
            if mult > 1: # more pos coeff (products) than neg (reactants)
                # then reduce pos coeff
                rxn_ids = np_where(spe > 0)
                mult = 1. / mult
            else: # more or equal neg coeff (reactants) than pos (products)
                # then reduce neg coeff
                rxn_ids = np_where(spe < 0)
            # Multiply coeffs of all species in the reaction where the species to balanced has a negative coeff
            for rxn_id in rxn_ids:
                sto_mat.transpose()[rxn_id] *= mult
    print(sto_mat)


    # Remove species that are overall balanced
    # Has to be done after the balancing because balancing has to be done on all species
    for i_spe in range(sto_mat.shape[0]):
        spe = sto_mat[i_spe]
        # If sum of coeffs is 0
        if spe.sum() == 0:
            # then remove species from overall equation by setting coeffs to 0
            spe.fill(0)

    return sto_mat


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
    temp_k: float
) -> ComponentContribution:

    cc = ComponentContribution()
    cc.p_h = Q_(ph)
    cc.ionic_strength = Q_(str(ionic_strength)+' mM')
    cc.p_mg = Q_(pMg)
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
