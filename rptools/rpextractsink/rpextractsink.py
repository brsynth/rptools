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
# because cobrapy is terrible
from brs_utils import timeout
# from timeout_decorator import timeout           as timeout_decorator_timeout
# from timeout_decorator import timeout_decorator as timeout_decorator_timeout_decorator
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
    model = cobra_io.read_sbml_model(sbml_path, use_fbc_package=True)
    dead_ends = []
    for met in model.metabolites:
        with model:  # The context manager is unecessary
            try:
                rxn = model.add_boundary(met, type="demand")
            except ValueError as e:
                if f"DM_{met}" in model.demands:
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
    dead_ends = _get_dead_end_metabolites(input_sbml) if remove_dead_end else []
    print(len(dead_ends))
    species = []
    rpsbml = rpSBML(input_sbml)
    for i in rpsbml.getModel().getListOfSpecies():
        if (
            i.getCompartment() == compartment_id and
            i.id not in dead_ends
        ):
            species.append(i)
    print("================== DEAD ENDS")
    print(dead_ends)
    print("================== SPECIES")
    print([s.id for s in species])
    if not species:
        logger.error('Could not retreive any species in the compartment: '+str(compartment_id))
        logger.error('Is the right compartment set?')
        return False
    with open(output_sink, 'w', encoding='utf-8') as outS:
        write(outS, ['Name', 'InChI'])
        for i in species:
            res = rpsbml.readMIRIAMAnnotation(i.getAnnotation())
            # extract the MNX id's
            try:
                mnx = res['metanetx'][0]
            except KeyError:
                logger.warning('Cannot find MetaNetX ID for '+str(i.getId()))
                continue
            try:
                inchi = cache.get('cid_strc')[mnx]['inchi']
            except KeyError:
                inchi = None
            if inchi:
                write(outS, [mnx, inchi])


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
