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


def get_inchi_from_url(
    url: str,
    logger: Logger = getLogger(__name__)
) -> str:
    '''
    Get the InChI from a given URL

    :param url: URL to retrieve the InChI from
    :type url: str
    :param logger: logger object, by default getLogger(__name__)
    :type logger: Logger, optional
    :return: InChI
    :rtype: str
    '''
    # print()
    # print(res[0])
    # page = r_get(res[0])
    # soup = BeautifulSoup(page.text, features='html.parser')
    # data = json_loads(soup.find('script', type='application/ld+json').text)
    # print(data)
    # Get the InChI structure from MetaNetX
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


## Generate the sink from a given model and the
#
# NOTE: this only works for MNX models, since we are parsing the id
# TODO: change this to read the annotations and extract the MNX id's
#
def genSink(
    cache,
    input_sbml,
    output_sink,
    remove_dead_end=False,
    compartment_id=default_comp,
    logger: Logger = getLogger(__name__)
):
    logger.debug('Extracting the sink from: '+str(input_sbml))
    dead_ends = _get_dead_end_metabolites(input_sbml, logger) if remove_dead_end else []
    species = []
    rpsbml = rpSBML(input_sbml)
    compartments = [comp.getId() for comp in rpsbml.getModel().getListOfCompartments()]
    # Check if given compartment is in the model
    if compartment_id not in compartments:
        logger.error(f'Unable to find the compartment \'{compartment_id}\' in the model.')
        logger.error(f'Available compartments are {compartments}.')
        return False
    logger.debug(f'List of Species: {list(rpsbml.getModel().getListOfSpecies())}')
    logger.debug(f'List of Compartments: {compartments}')
    logger.debug(f'Retrieving the species of the compartment {compartment_id}...')
    for i in rpsbml.getModel().getListOfSpecies():
        if (i.getCompartment() == compartment_id and
            cobra_format.to_cobra(i.id) not in dead_ends):
            species.append(i)
    if not species:
        logger.warning(f'Could not retrieve any species in the compartment: {compartment_id}')

    sink = {}
    for spe in species:
        logger.debug(f'Processing species: {spe.getId()}')
        res = rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
        logger.debug(f'MIRIAM: {res}')
        # search for metanetx URL
        logger.debug(f'Searching for MetaNetX ID...')
        mnx_url = ''
        for url in res:
            if 'metanetx' in url:
                mnx_url = url
                break
        inchi = ''
        mnx_id = ''
        if mnx_url:
            # Get InChI from MetaNetX cache
            mnx_id = mnx_url.split('/')[-1].split(':')[-1]
            logger.debug(f'MetaNetX ID: {mnx_id}')
            if mnx_id in cache.get('cid_strc'):
                inchi = cache.get('cid_strc')[mnx_id]['inchi']
        if not inchi:
            # Try to get InChI from MetaNetX
            if mnx_url:
                inchi = get_inchi_from_url(mnx_url, logger)
            else:
                # Try to get InChI from each other URL
                for url in res:
                    inchi = get_inchi_from_url(url, logger)
                    if inchi:
                        break
        if not inchi:
            logger.warning(f'Could not retrieve any InChI for {spe.getId()}')
        else:
            logger.debug(f'InChI: {inchi}')
            if not mnx_id:
                mnx_id = cobra_format.uncobraize(spe.getId())
            if mnx_id in sink:
                logger.warning(f'MetaNetX ID {mnx_id} already in sink')
            sink[mnx_id] = inchi

    logger.debug(f'Writing sink to {output_sink}...')
    # Write the sink file
    with open(output_sink, 'w', encoding='utf-8') as outS:
        write(outS, ['Name', 'InChI'])
        for _mnx_id, _inchi in sink.items():
            write(outS, [_mnx_id, _inchi])


def write(outFile, elts, delimiter=',', quotechar='"'):
    """
    Write elements of elts list into file as 'csv' would do

    :param file: file to write into
    :param elts: list of elements to write
    :param delimiter: character to insert between each element
    :param quotechar: character to put around each element
    """
    if elts:
        outFile.write(quotechar+elts[0]+quotechar)
        for elt in elts[1:]:
            outFile.write(delimiter+quotechar+elt+quotechar)
    outFile.write('\n')
