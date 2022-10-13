from os import path as os_path
from logging import (
    Logger,
    getLogger,
    StreamHandler
)
from typing import (
    Dict,
    List
)
from csv import reader as csv_reader
from copy import deepcopy
from equilibrator_api import (
    ComponentContribution,
    Q_
)
from numpy import (
    zeros as np_zeros,
    ndarray as np_ndarray
)
from scipy.optimize import linprog
from colored import fg, bg, attr
from brs_utils import (
    print_OK_adv as print_OK,
    print_title_adv as print_title
)
from chemlite import Reaction
from rptools.rplibs import rpPathway
from .Args import (
    DEFAULT_pH,
    DEFAULT_pMg,
    DEFAULT_ionic_strength,
)

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


def runThermo(
    pathway: rpPathway,
    cc: ComponentContribution=None,
    ph: float=DEFAULT_pH,
    ionic_strength: float=DEFAULT_ionic_strength,
    pMg: float=DEFAULT_pMg,
    compound_substitutes: Dict = None,
    logger: Logger = getLogger(__name__)
) -> Dict:
    """Given a tar input file, perform thermodynamics analysis for each rpSBML file.

    :param inFile: The path to the input file
    :param outFile: The path to the output file
    :param pathway_id: The id of the heterologous pathway of interest
    :param ph: The pH of the host organism (Default: 7.5)
    :param ionic_strength: Ionic strenght of the host organism (Default: 0.25M)
    :param pMg: The pMg of the host organism (Default: 3.0)
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

    print_title(
        txt='Pathway Reactions',
        logger=logger,
        waiting=False
    )
    for rxn in pathway.get_list_of_reactions():
        print_reaction(
            rxn=rxn,
            logger=logger
        )

    ## INTERMEDIATE COMPOUNDS
    # Optimise the production of target
    # and remove (if possible) intermediate compounds
    reactions = remove_compounds(
        compounds=pathway.get_intermediate_species(),
        reactions=pathway.get_list_of_reactions(),
        rxn_target_id=pathway.get_target_rxn_id(),
        logger=logger
    )

    ## eQuilibrator
    if cc is None:
        cc = initThermo(
            ph,
            ionic_strength,
            pMg,
            logger
        )

    # Search for the key ID known by eQuilibrator
    cc_species = {}
    substituted_species = {}
    sep = '__64__'
    if compound_substitutes is None:
        compound_substitutes = read_compound_substitutes(
            os_path.join(
                os_path.dirname(os_path.realpath(__file__)),
                'data',
                'compound_substitutes.csv'
            )
        )
    for spe in pathway.get_species():
        spe_split = spe.get_id().split(sep)
        if len(spe_split) > 1:
            _compound_substitutes = {k+sep+spe_split[1]: v for k, v in compound_substitutes.items()}
        else:
            _compound_substitutes = deepcopy(compound_substitutes)
        # If the specie is listed in substitutes file, then take search values from it
        # Check if starts with in case of compound names are like CMPD_NAME__64__COMPID
        if spe.get_id() in _compound_substitutes:
            cc_species[spe.get_id()] = search_equilibrator_compound(
                cc=cc,
                id=_compound_substitutes[spe.get_id()]['id'],
                inchikey=_compound_substitutes[spe.get_id()]['inchikey'],
                inchi=_compound_substitutes[spe.get_id()]['inchi'],
                logger=logger
            )
        # Else, take search values from rpCompound
        else:
            cc_species[spe.get_id()] = search_equilibrator_compound(
                cc=cc,
                id=spe.get_id(),
                inchikey=spe.get_inchikey(),
                inchi=spe.get_inchi(),
                smiles=spe.get_smiles(),
                logger=logger
            )
        if cc_species[spe.get_id()] != {}:
            if spe.get_id() != cc_species[spe.get_id()]['id']:
                substituted_species[spe.get_id()] = cc_species[spe.get_id()][cc_species[spe.get_id()]['cc_key']]
        else:
            logger.warning(f'Compound {spe.get_id()} has not been found within eQuilibrator cache')

    # Store thermo values for the net reactions
    # and for each of the reactions within the pathway
    results = {
        'net_reaction': {},
        'optimized_net_reaction': Reaction.sum_stoichio(reactions),
        'reactions': {},
        'optimized_reactions': {
            rxn.get_id(): rxn
            for rxn in reactions
        },
        'species': {},
        'substituted_species': substituted_species
    }

    # Get the formation energy for each compound
    for spe_id, cc_spe in cc_species.items():
        try:
            value = cc.standard_dg_formation(
                cc.get_compound(
                    cc_spe[cc_spe['cc_key']]
                )
            )[0]  # get .mu
        except Exception as e:
            value = None
            logger.debug(e)
        if value is None:
            value = 'NaN'
        results['species'][spe_id] = {
            'standard_dg_formation': {
                'value': value,
                'units': 'kilojoule / mole'
            }
        }

    # Build the list of IDs known by eQuilibrator
    species_cc_ids = {}
    for spe_id, cc_spe in cc_species.items():
        if cc_spe == {}:
            species_cc_ids[spe_id] = spe_id
        else:
            species_cc_ids[spe_id] = cc_spe[cc_spe['cc_key']]

    ## REACTIONS
    # Compute thermo for each reaction
    for rxn in pathway.get_list_of_reactions():
        results['reactions'][rxn.get_id()] = eQuilibrator(
            species_stoichio=rxn.get_species(),
            species_ids=species_cc_ids,
            cc=cc,
            logger=logger
        )

    ## THERMO
    print_title(
        txt='Computing thermodynamics (eQuilibrator)...',
        logger=logger,
        waiting=True
    )

    results['net_reaction'] = eQuilibrator(
        species_stoichio=Reaction.sum_stoichio(reactions),
        species_ids=species_cc_ids,
        cc=cc,
        logger=logger
    )

    print_OK(logger)

    # Write results into the pathway
    write_results_to_pathway(pathway, results, logger)

    return results


def read_compound_substitutes(filename: str) -> Dict[str, str]:
    comp_sub = {}
    with open(filename, 'r') as csv_file:
        reader = csv_reader(csv_file, delimiter=';')
        for row in reader:
            if (
                not row[0].startswith('#')
                and not row[6] == row[7] == row[8] == ''
            ):
                comp_sub[row[0]] = {
                    'id': row[6],
                    'inchi': row[7],
                    'inchikey': row[8]
                }
    return comp_sub


def search_equilibrator_compound(
    cc: ComponentContribution,
    id: str = None,
    inchikey: str = None,
    inchi: str = None,
    smiles: str = None,
    logger: Logger = getLogger(__name__)
) -> Dict[str, str]:

    def copy_data(
        compound,
        data: Dict,
        overwrite: bool = False,
    ) -> Dict:
        # ...copy initial data into result compound
        _compound = deepcopy(data)
        # fill with eQuilibrator data for empty fields
        for k, v in _compound.items():
            if (
                overwrite
                or v is None
                or v == ''
            ): _compound[k] = getattr(compound, k)
        # keep the key known by eQuilibrator
        _compound['cc_key'] = key
        # keep the value known by eQuilibrator
        _compound[key] = val
        return _compound

    data = {
        'id': id,
        'inchi_key': inchikey,
        'inchi': inchi,
        'smiles': smiles
    }
    for key, val in data.items():
        if val:
            compound = cc.get_compound(val)
            # If compound is found in eQuilibrator, then...
            if compound is not None:
                # ...copy initial data into result compound
                _compound = copy_data(compound, data)
                return _compound

    if inchikey:
        # In last resort, try to search only with the first part of inchikey
        compounds = cc.search_compound_by_inchi_key(
            # first part of inchikey
            inchikey.split('-')[0]
        )
        # eQuilibrator returns a list of compounds
        if compounds:
            # first compound in the list, hope it is sorted by decrease relevance
            _compound = copy_data(compounds[0], data, overwrite=True)
            # make inchi_key the ID key
            _compound['cc_key'] = 'inchi_key'
            return _compound

    return {}


def write_results_to_pathway(
  pathway: rpPathway,
  results: Dict,
  logger: Logger = getLogger(__name__)
) -> None:
    # Write species results
    for spe_id, score in results['species'].items():
        for k, v in score.items():
            pathway.get_specie(spe_id).add_thermo_info(
                key=k,
                value=v
            )
    # Write reactions results
    for rxn_id, score in results['reactions'].items():
        for k, v in score.items():
            pathway.get_reaction(rxn_id).add_thermo_info(
                key=k,
                value=v
            )
    # Write pathway result
    for k, v in results['net_reaction'].items():
        pathway.add_thermo_info(
        key=k,
        value=v
        )
    pathway.set_thermo_substituted_species(
        results['substituted_species']
    )

def eQuilibrator(
    species_stoichio: Dict[str, float],
    species_ids: Dict,
    cc: 'ComponentContribution',
    logger: Logger=getLogger(__name__)
) -> Dict:

    measures = {
        'dG0_prime': 'standard_dg_prime',
        'dGm_prime': 'physiological_dg_prime',
        'dG_prime': 'dg_prime',
        'dG': 'standard_dg',
    }

    ## Format reaction to what eQuilibrator expects
    compounds = {SIDES[0]['name']: [], SIDES[1]['name']: []}
    reactants = {spe_id: -spe_sto for (spe_id, spe_sto) in species_stoichio.items() if spe_sto < 0}
    products = {spe_id: spe_sto for (spe_id, spe_sto) in species_stoichio.items() if spe_sto > 0}

    # For both sides left and right
    for cmpd_id, cmpd_sto in reactants.items():
        compounds[SIDES[0]['name']] += [
            f'{cmpd_sto} {species_ids[cmpd_id]}'
        ]
    for cmpd_id, cmpd_sto in products.items():
        compounds[SIDES[1]['name']] += [
            f'{cmpd_sto} {species_ids[cmpd_id]}'
        ]

    # Join both sides
    rxn_str = '{left} = {right}'.format(
        left=' + '.join(compounds[SIDES[0]['name']]),
        right=' + '.join(compounds[SIDES[1]['name']])
    )

    logger.debug(rxn_str)

    results = {}
    try:
        # Parse formula by eQuilibrator
        rxn = cc.parse_reaction_formula(rxn_str)
        # Apply each CC method to each required measure
        thermo = {}
        for key in measures.keys():
            thermo[key] = getattr(cc, measures[key])(rxn)
        results = {
            key:{
                'value': float(str(thermo[key].value).split()[0]),
                'error': float(str(thermo[key].error).split()[0]),
                'units': str(thermo[key].units),
            } for key in thermo.keys()
        }

    except Exception as e:
        for key in measures.keys():
            results[key] = {
                'value': 'NaN',
                'error': 'NaN',
                'units': 'kilojoule / mole',
            }
        logger.debug(e)

    return results


def remove_compounds(
    reactions: List[Reaction],
    rxn_target_id: str,
    compounds: List = [],
    logger: Logger=getLogger(__name__)
) -> Dict:
    '''Try to remove compounds that are unknown in the eQuilibrator cache from reaction set (pathway).
    For this purpose, a stoichiometric matrix is built with only coefficients of these compounds.
    Then, try to solve as a linear equations system and apply new coeffs to reactions
    '''

    # unk_compounds = ['CMPD_0000000003', 'CMPD_0000000010', 'CMPD_0000000025']

    if not compounds:
        return reactions

    logger.debug(compounds)
    # Build the stoichio matrix for unknown compounds
    ## S
    sto_mat = build_stoichio_matrix(
        reactions=reactions,
        compounds=compounds,
        logger=logger
    )

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

    # print(sto_mat)
    # print(f'rxn_target_idx: {rxn_target_idx}')
    # print(coeffs)
    # exit()
    _reactions = []
    from copy import deepcopy
    ## Impact coeff to reactions
    for rxn_idx in range(len(reactions)):
        _reactions += [deepcopy(reactions[rxn_idx])]
        if (
            coeffs[rxn_idx] != 0
            and coeffs[rxn_idx] != abs(float("inf"))
        ):
            _reactions[rxn_idx].mult_stoichio_coeff(coeffs[rxn_idx])

    return _reactions


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

    ## S.x = 0
    # Init with zeros
    b_eq = np_zeros(nb_compounds)

    # 0 <= x <= 1
    # Init with (0, 1)
    # bounds = [(0, 1)]*len(reactions)
    # the same
    # lower bound (lb) and upper bound (ub) will be applied
    # to all variables.
    # bounds = (0, 1)
    bounds = (1, None)

    ## Max rxn_x
    c = np_zeros(nb_reactions)
    # Search the rxn ID that
    # maximizes the production of the target
    # translate it in the stoichio matrix index
    c[rxn_target_idx] = -1

    ## Solve
    res = linprog(
        c,
        A_eq=S,
        b_eq=b_eq,
        bounds=bounds,
        method='revised simplex'
    )

    logger.debug(f'COEFFS: {res}')

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
    if not compounds:
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
                sto_mat[cmpd_idx][rxn_idx] = -rxn.get_reactants()[cmpd_id]
            elif cmpd_id in rxn.get_products_ids():
                # Put the stochio value in the stochio matrix
                sto_mat[cmpd_idx][rxn_idx] = rxn.get_products()[cmpd_id]

    return sto_mat


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
            reactants=' + '.join([str(sto)+' '+str(spe) for spe,sto in rxn.get_reactants().items()]),
            products=' + '.join([str(sto)+' '+str(spe) for spe,sto in rxn.get_products().items()]),
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )


# used to initialise and download the data for equilibrator
def initThermo(
    ph: float=DEFAULT_pH,
    ionic_strength: float=DEFAULT_ionic_strength,
    pMg: float=DEFAULT_pMg,
    logger: Logger=getLogger(__name__)
) -> ComponentContribution:

    print_title(
        txt='Initialising eQuilibrator...',
        logger=logger,
        waiting=True
    )

    cc = ComponentContribution()

    cc.p_h = Q_(ph)
    cc.ionic_strength = Q_(f'{ionic_strength}M')
    cc.p_mg = Q_(pMg)
    # if temp_k is not None:
    #     cc.temperature = Q_(str(temp_k)+' K')

    print_OK(logger)

    return cc
