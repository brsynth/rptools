# from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import logging
from cobra          import io            as cobra_io
from cobra          import flux_analysis as cobra_flux_analysis
from tempfile       import TemporaryDirectory
from rptools.rplibs import rpSBML
from os             import path          as os_path
# because cobrapy is terrible
from timeout_decorator import timeout           as timeout_decorator_timeout
from timeout_decorator import timeout_decorator as timeout_decorator_timeout_decorator
TIMEOUT = 5

## Taken from Thomas Duigou's code
#
# @param input Cobra model object
#
def _reduce_model(cobraModel, logger=logging.getLogger(__name__)):
    """
    Reduce the model by removing reaction that cannot carry any flux and orphan metabolites

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
@timeout_decorator_timeout(TIMEOUT*60.0)
def _removeDeadEnd(sbml_path):
    cobraModel = cobra_io.read_sbml_model(sbml_path, use_fbc_package=True)
    cobraModel = _reduce_model(cobraModel)
    with TemporaryDirectory() as tmpOutputFolder:
        cobra_io.write_sbml_model(cobraModel, os_path.join(tmpOutputFolder, 'tmp.xml'))
        rpsbml = rpSBML(os_path.join(tmpOutputFolder, 'tmp.xml'))
        return rpsbml


#######################################################################
############################# PUBLIC FUNCTIONS ########################
#######################################################################


## Generate the sink from a given model and the
#
# NOTE: this only works for MNX models, since we are parsing the id
# TODO: change this to read the annotations and extract the MNX id's
#
def genSink(cache, input_sbml, output_sink, remove_dead_end=False, compartment_id='MNXC3', logger=logging.getLogger(__name__)):
    
    ### because cobrapy can be terrible and cause infinite loop depending on the input SBML model
    if remove_dead_end:
        try:
            rpsbml = _removeDeadEnd(input_sbml)
        except timeout_decorator_timeout_decorator.TimeoutError:
            logger.warning('removeDeadEnd reached its timeout... parsing the whole model')
            rpsbml = rpSBML(input_sbml)
    else:
        rpsbml = rpSBML(input_sbml)
    ### open the cache ###
    cytoplasm_species = []
    for i in rpsbml.getModel().getListOfSpecies():
        if i.getCompartment() == compartment_id:
            cytoplasm_species.append(i)
    if not cytoplasm_species:
        logger.error('Could not retreive any species in the compartment: '+str(compartment_id))
        logger.error('Is the right compartment set?')
        return False
    with open(output_sink, 'w', encoding='utf-8') as outS:
        # writer = csv_writer(outS, delimiter=',', quotechar='"', quoting=QUOTE_NONNUMERIC)
        # writer.writerow(['Name','InChI'])
        write(outS, ['Name', 'InChI'])
        for i in cytoplasm_species:
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
                # writer.writerow([mnx,inchi])


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
