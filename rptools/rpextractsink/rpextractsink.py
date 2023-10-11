# from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
from logging import (
    Logger,
    getLogger
)
from cobra          import io            as cobra_io
from cobra          import flux_analysis as cobra_flux_analysis
from tempfile       import TemporaryDirectory
from rptools.rplibs import rpSBML
from .Args import default_comp
from os             import path          as os_path
from requests import get as r_get
from re import search as re_search

from rr_cache import rrCache
from rptools.rpfba import cobra_format

# because cobrapy is terrible
from brs_utils import timeout
__TIMEOUT = 5

## Taken from Thomas Duigou's code
#
# @param input Cobra model object
#
def _reduce_model(
    cobraModel,
    logger: Logger = getLogger(__name__)
):
    """Reduces the model by removing reaction that cannot carry any flux and orphan metabolites

    :param model: cobra model object
    :return: reduced cobra model object
    """
    
    lof_zero_flux_rxn = cobra_flux_analysis.find_blocked_reactions(cobraModel, open_exchanges=True)
    # For assert and logger: Backup the list of metabolites and reactions
    # nb_metabolite_model_ids = set([m.id for m in cobraModel.metabolites])
    nb_reaction_model_ids   = set([m.id for m in cobraModel.reactions])
    # Remove unwanted reactions and metabolites
    cobraModel.remove_reactions(lof_zero_flux_rxn, remove_orphans=True)
    # # Assert the number are expected numbers
    # assert len(set([m.id for m in cobraModel.reactions])) == len(nb_reaction_model_ids) - len(lof_zero_flux_rxn)
    if len(set([m.id for m in cobraModel.reactions])) != len(nb_reaction_model_ids) - len(lof_zero_flux_rxn):
        logger.error(" *** Nb of reactions incorrect. Exiting...")
        exit()
    return cobraModel

##
#
#
@timeout(__TIMEOUT*60.0)
def _removeDeadEnd(sbml_path) -> rpSBML:
    cobraModel = cobra_io.read_sbml_model(sbml_path, use_fbc_package=True)
    cobraModel = _reduce_model(cobraModel)
    with TemporaryDirectory() as tmpOutputFolder:
        cobra_io.write_sbml_model(cobraModel, os_path.join(tmpOutputFolder, 'tmp.xml'))
        rpsbml = rpSBML(os_path.join(tmpOutputFolder, 'tmp.xml'))
        return rpsbml


@timeout(__TIMEOUT*60.0)
def _get_dead_end_metabolites(
    sbml_path: str,
    logger: Logger = getLogger(__name__)
) -> list:
    """Search for dead end metabolites

    Metabolites are iteratively tested for production
    by adding a demand bound on the metabolite and 
    trying to maximize the flux passing through it.

    Parameters
    ----------
    model : rpSBML
        model to investigate
    logger : Logger, optional
        logger object, by default getLogger(__name__)

    Returns
    -------
    list
        dead end metabolite IDs
    """
    logger.debug("Searching for dead end metabolites...")
    model = cobra_io.read_sbml_model(sbml_path, use_fbc_package=True)
    dead_ends = []
    for met in model.metabolites:
        with model:  # The context manager is unecessary
            try:
                logger.debug(f"Creating demand bound for metabolite {met.id}")
                rxn = model.add_boundary(met, type="demand")
            except ValueError as e:
                logger.debug(f"Cannot create a demand bound for metabolite {met.id}... {e}")
                if f"DM_{met}" in model.demands:
                    logger.debug(f"Using existing demand bound for metabolite {met.id}... {e}")
                    rxn = model.reactions.get_by_id(f"DM_{met}")
                else:
                    logger.warning(
                        "Cannot create a demand bound for metabolite {} "
                        "while searching for dead end metabolites... "
                        "This metabolite will be considered as producible"
                        )
            model.objective = rxn
            value = model.slim_optimize(error_value=0.0)
            if value == 0.0:
                dead_ends.append(met.id)
    return dead_ends

#######################################################################
############################# PUBLIC FUNCTIONS ########################
#######################################################################


def get_inchi_from_crossid(
    id: str,
    logger: Logger = getLogger(__name__)
) -> str:
    '''
    Get the InChI from a given ID on MetaNetX

    :param id: ID to retrieve the InChI from
    :type id: str
    :param logger: logger object, by default getLogger(__name__)
    :type logger: Logger, optional
    :return: InChI
    :rtype: str
    '''
    # Get the InChI structure from MetaNetX
    logger.debug(f'Retrieving InChI from MetaNetX for {id}...')
    url_mnx = 'https://www.metanetx.org'
    url_search = f'{url_mnx}/cgi-bin/mnxweb/search'
    try:
        page = r_get(f'{url_search}?query={id}')
    except Exception as e:
        logger.warning(f'Connection lost from {url_search}')
        return ''
    url_crossid = re_search(r'/chem_info/\w+', page.text).group()
    return get_inchi_from_url(f'{url_mnx}{url_crossid}', logger)


def get_inchi_from_url(
    url: str,
    logger: Logger = getLogger(__name__)
) -> str:
    '''
    Get the InChI from a given URL
    
    :param url: URL
    :type url: str
    :param logger: logger object, by default getLogger(__name__)
    :type logger: Logger, optional
    :return: InChI
    :rtype: str
    '''
    logger.debug(f'Retrieving InChI from {url}...')
    try:
        page = r_get(url)
    except Exception as e:
        logger.warning(f'Connection lost from {url}')
        return ''
    x = re_search(r'InChI=[^<"\n-]+', page.text)
    if x:
        inchi = x.group()
        return inchi
    else:
        return ''


def genSink(
    cache: rrCache,
    input_sbml: str,
    remove_dead_end: bool = False,
    compartment_id: str = default_comp,
    standalone: bool = False,
    logger: Logger = getLogger(__name__)
) -> dict:
    '''
    Generate the sink dictionary from a given SBML file
    
    :param cache: cache object
    :type cache: rrCache
    :param input_sbml: path to the SBML file
    :type input_sbml: str
    :param remove_dead_end: remove dead end metabolites, by default False
    :type remove_dead_end: bool, optional
    :param compartment_id: compartment ID, by default default_comp
    :type compartment_id: str, optional
    :param standalone: do not use MetaNetX, by default False
    :type standalone: bool, optional
    :param logger: logger object, by default getLogger(__name__)
    :type logger: Logger, optional
    :return: sink dictionary
    :rtype: dict
    '''
    logger.debug('Extracting the sink from: '+str(input_sbml))

    sink = {}

    dead_ends = _get_dead_end_metabolites(input_sbml, logger) if remove_dead_end else []
    species = []
    rpsbml = rpSBML(input_sbml)
    compartments = [comp.getId() for comp in rpsbml.getModel().getListOfCompartments()]

    # Check if given compartment is in the model
    if compartment_id not in compartments:
        logger.error(f'Unable to find the compartment \'{compartment_id}\' in the model.')
        logger.error(f'Available compartments are {compartments}.')
        return sink

    logger.debug(f'List of Species: {list(rpsbml.getModel().getListOfSpecies())}')
    logger.debug(f'List of Compartments: {compartments}')
    logger.debug(f'Retrieving the species of the compartment {compartment_id}...')

    for i in rpsbml.getModel().getListOfSpecies():
        if (i.getCompartment() == compartment_id and
            cobra_format.to_cobra(i.id) not in dead_ends):
            species.append(i)
    if not species:
        logger.warning(f'Could not retrieve any species in the compartment: {compartment_id}')

    for spe in species:

        logger.debug(f'Processing species: {spe.getId()}')
        miriam = rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
        logger.debug(f'MIRIAM: {miriam}')
        inchi = ''

        # Find MNX ID from MIRIAM
        mnx_id = find_mnx_id(miriam)

        if mnx_id:
            # Get InChI from cache
            if mnx_id in cache.get('cid_strc'):
                inchi = cache.get('cid_strc')[mnx_id]['inchi']
            else:
                logger.debug(f'Could not retrieve InChI for {mnx_id} from cache')
                if not standalone:
                    # Get InChI from MetaNetX
                    inchi = get_inchi_from_url(
                        f'https://www.metanetx.org/chem_info/{mnx_id}',
                        logger
                    )
        elif not standalone:
            # Search on MetaNetX with MIRIAM cross-references
            i = 0
            while not inchi and i < len(miriam):
                url = miriam[i]
                id = url.split('/')[-1]
                inchi = get_inchi_from_crossid(id, logger)
                i += 1

        if not inchi:
            logger.warning(f'Could not retrieve any InChI for {spe.getId()}')
        else:
            logger.debug(f'InChI: {inchi}')
            if spe.getId() in sink:
                logger.warning(f'MetaNetX ID {spe.getId()} already in sink')
            sink[spe.getId()] = inchi

    return sink


def find_mnx_id(
    miriam: list,
    logger: Logger = getLogger(__name__)
) -> str:
    '''
    Find the MetaNetX ID from the MIRIAM annotation
    
    :param miriam: list of MIRIAM annotation
    :type miriam: list
    :param logger: logger object, by default getLogger(__name__)
    :type logger: Logger, optional
    :return: MetaNetX ID
    :rtype: str
    '''
    mnx_id = ''
    for url in miriam:
        if 'metanetx' in url:
            mnx_id = url.split('/')[-1].split(':')[-1]
            logger.debug(f'MetaNetX ID: {mnx_id}')
            break
    return mnx_id
