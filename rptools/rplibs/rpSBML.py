import libsbml
import re
import numpy as np

from copy import deepcopy
from filetype import guess
from hashlib import sha256
from math import isnan
from pandas import DataFrame  as pd_DataFrame

from os import (
        path as os_path,
        remove
)
from json import (
    load as json_load,
    dump as json_dump,
    dumps as json_dumps
)
from inspect import (
    getmembers as inspect_getmembers,
    ismethod as inspect_ismethod
)
from logging import (
    Logger,
    getLogger
)
from typing import (
    List,
    Dict,
    Tuple,
    TypeVar,
    Union
)
from tempfile import (
    NamedTemporaryFile,
    TemporaryDirectory,
    gettempdir,
)
import cobra
from cobra import (
    io as cobra_io
)
from cobra.io.sbml import (
    CobraSBMLError,
    validate_sbml_model
)
from cobra.medium.annotations import (
    compartment_shortlist,
    excludes,
    sbo_terms
)

from brs_utils import(
    extract_gz,
)
from .rpGraph import rpGraph

## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call BRSYNTH nodes.
# The object holds an SBML object and a series of methods to write and access BRSYNTH related annotations


def err_code(code: int) -> str:
    codes = {
        libsbml.LIBSBML_LEVEL_MISMATCH: 'libsbml.LIBSBML_LEVEL_MISMATCH',
        libsbml.LIBSBML_VERSION_MISMATCH: 'libsbml.LIBSBML_VERSION_MISMATCH',
        libsbml.LIBSBML_DUPLICATE_OBJECT_ID: 'libsbml.LIBSBML_DUPLICATE_OBJECT_ID',
        libsbml.LIBSBML_INVALID_OBJECT: 'libsbml.LIBSBML_INVALID_OBJECT',
        libsbml.LIBSBML_OPERATION_FAILED: 'libsbml.LIBLIBSBML_OPERATION_FAILEDSBML_LEVEL_MISMATCH',
    }
    return codes.get(code, f'Unknown code {code}')


##################################################################
############################### rpSBML ###########################
##################################################################


class rpSBML:

    """This class uses the libSBML object and handles it by adding BRSynth annotation
    """
    def __init__(
        self,
        inFile: str = None,
        rpsbml: 'rpSBML' = None,
        name: str = None,
        logger: Logger = getLogger(__name__)
    ) -> 'rpSBML':
        """Constructor for the rpSBML class

        Note that the user can pass either a document libSBML object or a path to a SBML file. If a path is passed it overwrite the passed document object.

        :param modelName: The Name of the model
        :param document: The libSBML document class (Default: None)
        :param inFile: The path of a SBML file (Default: '')

        :type modelName: str
        :type path: str
        :type document: libsbml.SBMLDocument
        """

        # logger
        self.logger = logger

        self.logger.debug('New instance of rpSBML')
        self.logger.debug(f'inFile: {inFile}')
        self.logger.debug(f'rpsbml: {rpsbml}')
        self.logger.debug(f'name:   {name}')

        # document
        # if an sbml file is given, then read it (self.document will be created)
        if inFile is not None:
            infile = inFile
            try:
                kind = guess(infile)
            except FileNotFoundError as e:
                self.logger.error(f'{e}')
                self.logger.error('Exiting...')
                exit()
            except TypeError as e:
                self.logger.error(f'{e}')
                self.logger.error('Exiting...')
                exit()
            with TemporaryDirectory() as temp_d:
                if kind:
                    self.logger.debug(f'inFile is detected as {kind}')
                    if kind.mime == 'application/gzip':
                        infile = extract_gz(inFile, temp_d)
                self.logger.debug(f'Reading {infile} file...')
                self.readSBML(infile)
                self.logger.debug('File READ')
        else:
            if rpsbml is None:
                self.document = None
            else:
                self.logger.debug(f'Cloning SBMLDocument...')
                self.document = rpsbml.getDocument().clone()
                self.logger.debug(f'SBMLDocument CLONED')

        if rpsbml is not None and not self.checkSBML():
            self.logger.error('SBML document not valid')
            self.logger.error('Exiting...')
            exit()

        # model name
        self.logger.debug('Setting name...')
        self.setName(name if name else self.getName())
        if self.getModel():
            if self.getModel().getName() == '':
                self.getModel().setName(self.getName())

        # headers
        self.miriam_header = {
            'compartment': {
                'mnx': 'metanetx.compartment/',
                'bigg': 'bigg.compartment/',
                'seed': 'seed/',
                'name': 'name/'
            },
            'reaction': {
                'mnx': 'metanetx.reaction/',
                'rhea': 'rhea/',
                'reactome': 'reactome/',
                'bigg': 'bigg.reaction/',
                'sabiork': 'sabiork.reaction/',
                'ec-code': 'ec-code/',
                'biocyc': 'biocyc/',
                'lipidmaps': 'lipidmaps/',
                'uniprot': 'uniprot/'
            },
            'species': {
                'inchikey': 'inchikey/',
                'pubchem': 'pubchem.compound/',
                'mnx': 'metanetx.chemical/',
                'chebi': 'chebi/',
                'bigg': 'bigg.metabolite/',
                'hmdb': 'hmdb/',
                'kegg_c': 'kegg.compound/',
                'kegg_d': 'kegg.drug/',
                'biocyc': 'biocyc/',
                'seed': 'seed.compound/',
                'metacyc': 'metacyc.compound/',
                'sabiork': 'sabiork.compound/',
                'reactome': 'reactome/'
            }
        }
        self.header_miriam = {
            'compartment': {
                'metanetx.compartment': 'mnx',
                'bigg.compartment': 'bigg',
                'seed': 'seed',
                'name': 'name'
                },
            'reaction': {
                'metanetx.reaction': 'mnx',
                'rhea': 'rhea',
                'reactome': 'reactome',
                'bigg.reaction': 'bigg',
                'sabiork.reaction': 'sabiork',
                'ec-code': 'ec-code',
                'biocyc': 'biocyc',
                'lipidmaps': 'lipidmaps',
                'uniprot': 'uniprot'
            },
            'species': {
                'inchikey': 'inchikey',
                'pubchem.compound': 'pubchem',
                'metanetx.chemical': 'mnx',
                'chebi': 'chebi',
                'bigg.metabolite': 'bigg',
                'hmdb': 'hmdb',
                'kegg.compound': 'kegg_c',
                'kegg.drug': 'kegg_d',
                'biocyc': 'biocyc',
                'seed.compound': 'seed',
                'metacyc.compound': 'metacyc',
                'sabiork.compound': 'sabiork',
                'reactome': 'reactome'
            }
        }

        self.logger.debug('rpSBML instance created')


    def checkSBML(self) -> bool:
        self.logger.debug('Checking SBML format...')
        try:
            return libsbml.SBMLValidator().validate(self.getDocument()) == 0
        except ValueError:
            return False


    def getModel(self):
        if self.getDocument():
            return self.getDocument().getModel()
        else:
            return None


    def getDocument(self):
        return self.document


    def getName(self):

        name = ''

        # try if modelName already exists
        try:
            name = self.modelName
        except AttributeError:
            pass

        # self.modelName exists, returns it
        if name:
            return name
        else: # else try to get name from model
            return self.getModel().getName() if self.getModel() else None


    def getObjective(
        self,
        objective_id: str,
        objective_value: float 
    ) -> None:

        objective = self.getPlugin(
            'fbc'
        ).getObjective(objective_id)

        self.checklibSBML(
            objective,
            'Getting objective '+str(objective_id)
        )

        return objective


    def setName(self, name):
        if not name:
            if self.getModel():
                self.modelName = self.getModel().getId()
            else:
                self.modelName = 'dummy'
        else:
            self.modelName = name
        # check if name starts with a letter
        # if not, add rp_ at the beginning
        self.modelName = self.modelName if self.modelName[0].isalpha() else 'rp_'+self.modelName
        self.logger.debug(f'Name set to {self.modelName}')


    def setLogger(self, logger):
        self.logger = logger


    ##########################################################################
    ############################ QUERY #######################################
    ##########################################################################
    def search_compartment(
        self,
        compartment: str
    ) -> libsbml.Compartment:
        """Search in the model if a compartment id exists.

        :param compartment: An id or a name to search in the model

        :type compartment: str

        :return: Return a compartment if the specie is in the model, None otherwise
        :rtype: libsbml.Compartment
        """
        self.logger.debug(f'compartment: {compartment}')

        # Search in model.
        if self.getModel() is None:
            return None
        # Build data.
        comp_models = self.getModel().getListOfCompartments()
        # Find synonyms
        comp_synonyms = []
        for c_short, c_long in compartment_shortlist.items():
            c_long.append(c_short)
            c_long = [x.lower() for x in c_long]
            if compartment.lower() in c_long:
                comp_synonyms = c_long
                break
        if len(comp_synonyms) == 0:
            comp_synonyms.append(compartment.lower())
        # Not strict
        for comp_synonym in comp_synonyms:
            for comp_model in comp_models:
                if comp_synonym == comp_model.getId().lower() \
                    or comp_synonym == comp_model.getName().lower():
                    return comp_model
        return None

    def search_specie(
        self,
        specie: str
    ) -> libsbml.Species:
        """Search in model if a specie exists.

        :param specie: A specie id to search in the model

        :type specie: str

        :return: Return a specie if the specie is in the model, None otherwise
        :rtype: libsbml.Specie
        """
        # Search in model.
        if self.getModel() is None:
            return None
        for spe in self.getModel().getListOfSpecies():
            if re.search(specie, spe.getId(), re.IGNORECASE) \
                or re.search(specie, spe.getName(), re.IGNORECASE):
                return spe
        return None

    def search_reaction(
        self,
        reaction: str
    ) -> libsbml.Reaction:
        """Search in model if a reaction exists.

        :param reaction: A reaction id to search in the model

        :type reaction: str

        :return: Return a reaction if the reaction is in the model, None otherwise
        :rtype: libsbml.Reaction
        """
        # Search in model.
        if self.getModel() is None:
            return None
        for rxn in self.getModel().getListOfReactions():
            if re.search(reaction, rxn.getId(), re.IGNORECASE) \
                or re.search(reaction, rxn.getName(), re.IGNORECASE):
                return rxn
        return None

    def has_compartment(
        self,
        compartment: str,
        strict: bool=False
    ) -> bool:
        """Check in the model if a compartment id exists.

        :param compartment: An id (or name if not strict mode) \
            to search in the model
        :param strict: Perform research with the exact id provided \
            without mapping (optional: false)

        :type compartment: str
        :type strict: bool

        :return: Success or Failure if the compartment is in the model
        :rtype: bool
        """
        if self.getModel() is None:
            return None
        if strict:
            return self.getModel().getCompartment(compartment)
        else:
            return self.search_compartment(compartment)

    def has_specie(
        self,
        specie: str,
        strict: bool=False
    ) -> bool:
        """Check in the model if a specie exists.

        :param specie: A specie to extract its id
        :param strict: Perform research with the exact id provided \
            without mapping (optional: false)

        :type specie: str
        :type strict: bool

        :return: Success or Failure if the specie is in the model
        :rtype: bool
        """
        # Check
        if self.getModel() is None:
            return False
        if strict:
            if self.getModel().getSpecies(specie) is None:
                return False
            return True
        else:
            if self.search_specie(specie) is None:
                return False
            return True

    def has_reaction(
        self,
        reaction: str,
        strict: bool=False
    ) -> bool:
        """Check in the model if a reaction exists.

        :param reaction: A reaction to extract its id
        :param strict: Perform research with the exact id provided \
            without mapping (optional: false)

        :type reaction: str
        :type strict: bool

        :return: Success or Failure if the reaction is in the model
        :rtype: bool
        """
        # Check
        if self.getModel() is None:
            return False
        if strict:
            if self.getModel().getReactions(reaction) is None:
                return False
            return True
        else:
            if self.search_reaction(reaction) is None:
                return False
            return True

    #############################################################################################################
    ############################################ MERGE ##########################################################
    #############################################################################################################

    @staticmethod
    def mergeFiles(
        input_sbml: str,
        input_target: str,
        output_merged: str,
        compartment_id: str,
        logger: Logger = getLogger(__name__)
    ) -> bool:
        """Public function that merges two SBML files together

        :param path_source: Path of the source SBML file
        :param path_target: Path of the target SBML file
        :param path_merge: Path of the output SBML file

        :type path_source: str
        :type path_target: str
        :type path_merge: str

        :return: Success or failure of the function
        :rtype: bool
        """

        if not os_path.exists(input_sbml):
            logger.error('Source SBML file is invalid: '+str(input_sbml))
            return False

        if not os_path.exists(input_target):
            logger.error('Target SBML file is invalid: '+str(input_target))
            return False

        source_rpsbml = rpSBML(
            inFile = input_sbml,
            name = 'source',
            logger = logger
        )
        target_rpsbml = rpSBML(
            inFile = input_target,
            name = 'target',
            logger = logger
        )

        merged_rpsbml, reactions_in_both, missing_species = rpSBML.merge(
            source_rpsbml=source_rpsbml,
            target_rpsbml=target_rpsbml,
            compartment_id=compartment_id,
            logger=logger
        )

        merged_rpsbml.writeToFile(output_merged)

        return True


    @staticmethod
    # TODO: add a confidence in the merge using the score in
    def merge(
        pathway: 'rpSBML',
        model: 'rpSBML',
        compartment_id: str,
        logger: Logger = getLogger(__name__)
    ) -> Tuple['rpSBML', Dict, List, str]:
        # # Copy target rpSBML object into a new one so that
        # # it can be modified and returned
        # merged_rpsbml = rpSBML(
        #     rpsbml = model,
        #     logger = logger
        # )
        # merged_rpsbml.getModel().setId(
        #     merged_rpsbml.getModel().getId()+'__'+pathway.getModel().getId()
        # )
        # merged_rpsbml.getModel().setName(
        #     merged_rpsbml.getModel().getName()+' merged with '+pathway.getModel().getId()
        # )
        # merged_rpsbml.getModel().appendFrom(pathway.getModel())
        # if merged_rpsbml.checkSBML():
        #     return (
        #         merged_rpsbml,
        #         None,
        #         None,
        #         compartment_id
        #     )
        # else:
        #     logger.error('Merging rpSBML objects results in a invalid SBML format')
        #     return None, None, None, None
        """
        Merge two models species and reactions using the annotations to recognise the same species and reactions

        The source model has to have both the GROUPS and FBC packages enabled in its SBML.
        The course must have a group called rp_pathway. If not, use the readSBML() function to create a model
        We add the reactions and species from the rpsbml to the target_model

        Parameters
        ----------
        source_rpsbml: rpSBML
            The source rpSBML object
        target_rpsbml: rpSBML
            The target rpSBML object
        logger : Logger
            The logger object.

        Returns
        -------
        None or str If outfile is None, returns the resulting csv format as a string. Otherwise returns None..

            :param source_rpsbml: The source rpSBML object
            :param target_rpsbml: The target rpSBML object

            :type source_rpsbml: rpSBML
            :type target_rpsbml: rpSBML

            :return: Tuple of dict where the first entry is the species source to target conversion and the second is the reaction source to target conversion
            :rtype: tuple
        """
        logger.debug('pathway: ' + str(pathway))
        logger.debug('model: ' + str(model))
        logger.info(
            f"Merging : {model.getName()} model and {pathway.getName()} pathway..."
        )

        # Copy target rpSBML object into a new one so that
        # it can be modified and returned
        merged_rpsbml = rpSBML(
            rpsbml = model,
            logger = logger
        )

        ## MODEL FBC ###################################
        # Find the ID's of the similar target_rpsbml.model species
        pkg = 'fbc'
        url = 'http://www.sbml.org/sbml/level3/version1/fbc/version2'
        pathway.enable_package(pkg, url)
        merged_rpsbml.enable_package(pkg, url)

        # pkg = 'groups'
        # url = 'http://www.sbml.org/sbml/level3/version1/groups/version1'
        # source_fbc = source_rpsbml.enable_package(pkg, url)
        # target_fbc = merged_rpsbml.enable_package(pkg, url)

        # target_fbc, source_fbc = rpSBML.getFBCModel(
        #     source_sbml_doc = source_rpsbml.getDocument(),
        #     target_sbml_doc = merged_rpsbml.getDocument(),
        #     logger = logger
        # )

        ## UNIT DEFINITIONS ############################
        # return the list of unit definitions id's for the target to avoid overwritting
        # WARNING: this means that the original unit definitions will be prefered over the new one
        merged_rpsbml.copyUnitDefinitions(pathway)

        ## COMPARTMENTS ################################
        # Compare by MIRIAM annotations
        # Note that key is source and value is target conversion
        merged_rpsbml.copyCompartments(pathway)

        ## PARAMETERS ##################################
        # WARNING: here we compare by ID
        # target_paramsID = rpSBML.getParameters(
        merged_rpsbml.copyParameters(pathway)

        ## FBC GENE PRODUCTS ###########################
        # WARNING: here we compare by ID
        merged_rpsbml.copyFBCGeneProducts(pathway)

        ## FBC OBJECTIVES ##############################
        # WARNING: here we compare by ID
        merged_rpsbml.copyFBCObjectives(pathway)

        ## SPECIES #####################################
        species_pathway_model_transl, missing_species = merged_rpsbml.copySpecies(
            source_sbml=pathway,
            compartment_id=compartment_id
        )

        ## REACTIONS ###################################
        reactions_in_both = merged_rpsbml.copyReactions(
            source_sbml=pathway,
            species=species_pathway_model_transl
        )

        ## GROUPS ######################################
        merged_rpsbml.copyGroups(
            source_sbml=pathway,
            species=species_pathway_model_transl,
            reactions=reactions_in_both
        )

        ## TITLES ######################################
        merged_rpsbml.copyTitles(
            source_sbml=pathway
        )

        # merged_rpsbml._checkSingleParent()


        # merged_rpsbml.write_to_file('joan.xml')

        if merged_rpsbml.checkSBML():
            return (
                merged_rpsbml,
                reactions_in_both,
                missing_species,
                compartment_id
            )
        else:
            logger.error('Merging rpSBML objects results in a invalid SBML format')
            return None, None, None, None


    def check_SBML_compartment(self, compartment_id: str) -> Tuple[str, str]:
        # Check model compartment ID
        # Set new compartment ID in case it exists under another ID
        compartment = self.has_compartment(compartment_id)
        if compartment is not None:
            self.logger.debug(
                f"Compartment '{compartment_id}' found in the model '{self.getName()}' as {compartment.getId()}"
            )
            if compartment_id != compartment.getId():
                self.logger.warning(
                    f"Compartment '{compartment_id}' has been replaced by {compartment.getId()}"
                )
            return compartment.getId()
        else:
            self.logger.error(
                f"Compartment '{compartment_id}' not found in the model '{self.getName()}'"
            )
            self.logger.error(
                "Available compartments: {comp_ids}".format(
                    comp_ids=", ".join(
                        [
                            comp.getId()
                            for comp in list(self.getModel().getListOfCompartments())
                        ]
                    )
                )
            )
            return None


    def check_SBML_rxnid(self, rxn_id: str) -> Tuple[str, str]:
        self.logger.debug(f"rpsbml: {self.getName()}")
        self.logger.debug(f"rxn_id: {rxn_id}")

        # Check model reaction ID
        # Get reaction from the rpSBML
        _rxn = self.getModel().getReaction(rxn_id)
        if _rxn is None:
            self.logger.error(
                f"Reaction ID '{rxn_id}' not found in the model '{self.getName()}'"
            )
            possible_rxn_ids = [
                rxn.getId()
                for rxn in list(self.getModel().getListOfReactions())
                if rxn_id in rxn.getId().lower()
            ]
            if possible_rxn_ids != []:
                self.logger.error(
                    "Possible reactions: {rxn_ids}".format(
                        rxn_ids=", ".join(possible_rxn_ids)
                    )
                )
            return None

        self.logger.debug(f"rxn_id found as {rxn_id}")

        return _rxn.getId()


    @staticmethod
    def renameSpecies(
        rpsbml: 'rpSBML',
        transl_dict_species: Dict[str, str]
    ) -> None:
        with NamedTemporaryFile(delete=False) as tempf:
            rpsbml.write_to_file(tempf.name)
            tempf.close()
            with open(tempf.name, 'r') as in_f:
                text = in_f.read()
                for speID, speID_transl in transl_dict_species.items():
                    text = text.replace(speID+'"', speID_transl+'"')
                with NamedTemporaryFile(mode='w', delete=False) as out_tempf:
                    out_tempf.write(text)
                    out_tempf.close()
                    rpsbml = rpSBML(inFile=out_tempf.name)
        remove(tempf.name)
        remove(out_tempf.name)
        return rpsbml

    @staticmethod
    def cobraize(rpsbml: 'rpSBML') -> 'rpSBML':
        from rptools.rpfba.cobra_format import (
            cobraize
        )
        return rpSBML.renameSpecies(
            rpsbml=rpsbml,
            transl_dict_species={
                specie.getId(): cobraize(
                    specie.getId(),
                    specie.getCompartment()
                ) for specie in list(rpsbml.getModel().getListOfSpecies())
            }
        )

    def enable_package(
        self,
        pkg: str,
        url: str
    ) -> None:
    # ) -> libsbml.SBasePlugin:
        if not self.getModel().isPackageEnabled(pkg):
            rpSBML.checklibSBML(
                self.getModel().enablePackage(url, pkg, True),
                'Enabling the ' + pkg + 'package'
            )
        # # note sure why one needs to set this as False
        # rpSBML.checklibSBML(
        #     self.getDocument().setPackageRequired(pkg, False),
        #     pkg + ' package not required'
        # )
        # return self.getModel().getPlugin(pkg)


    def copyUnitDefinitions(
        self,
        source_sbml: 'rpSBML'
    ) -> None:

        source_sbml_doc = source_sbml.getDocument()
        self.logger.debug(f'source_sbml_doc: {source_sbml_doc}')

        target_unitDefID = [i.getId() for i in self.getModel().getListOfUnitDefinitions()]

        for source_unitDef in source_sbml_doc.getModel().getListOfUnitDefinitions():

            if not source_unitDef.getId() in target_unitDefID: # have to compare by ID since no annotation
                # create a new unitDef in the target
                target_unitDef = self.getModel().createUnitDefinition()
                rpSBML.checklibSBML(
                    target_unitDef,
                    'fetching target unit definition'
                )
                # copy unitDef info to the target
                rpSBML.checklibSBML(
                    target_unitDef.setId(source_unitDef.getId()),
                    'setting target unit definition ID'
                )
                rpSBML.checklibSBML(
                    target_unitDef.setAnnotation(source_unitDef.getAnnotation()),
                    'setting target unit definition Annotation'
                )

                for source_unit in source_unitDef.getListOfUnits():
                    # copy unit info to the target unitDef
                    target_unit = target_unitDef.createUnit()
                    rpSBML.checklibSBML(
                        target_unit,
                        'creating target unit'
                    )
                    rpSBML.checklibSBML(
                        target_unit.setKind(
                            source_unit.getKind()
                        ),
                        'setting target unit kind'
                    )
                    rpSBML.checklibSBML(
                        target_unit.setExponent(source_unit.getExponent()),
                        'setting target unit exponent'
                    )
                    rpSBML.checklibSBML(
                        target_unit.setScale(source_unit.getScale()),
                        'setting target unit scale'
                    )
                    rpSBML.checklibSBML(
                        target_unit.setMultiplier(source_unit.getMultiplier()),
                        'setting target unit multiplier'
                    )

                # add to the list to make sure its not added twice
                target_unitDefID.append(source_unitDef.getId())


    def copyCompartments(
        self,
        rpsbml: 'rpSBML'
    ) -> None:

        for source_compartment in rpsbml.getModel().getListOfCompartments():

            found = False

            source_annotation = source_compartment.getAnnotation()

            if not source_annotation:
                self.logger.warning(f'No annotation for the source of compartment {source_compartment.getId()}')

            for target_compartment in self.getModel().getListOfCompartments():
                self.logger.debug(f'source_compartment id: {source_compartment.getId()}; target_compartment id: {target_compartment.getId()}')
                self.logger.debug(f'source_compartment name: {source_compartment.getName()}; target_compartment name: {target_compartment.getName()}')
                # compare by ID or name first
                if (
                    source_compartment.getId() == target_compartment.getId()
                    or source_compartment.getName() == target_compartment.getName()
                ):
                    found = True
                    break
                # then, compare by MIRIAM
                else:
                    target_annotation = target_compartment.getAnnotation()
                    self.logger.debug(f'source_annotation: {source_annotation}')
                    self.logger.debug(f'target_annotation: {target_annotation}')
                    if not target_annotation:
                        self.logger.warning(f'No annotation for the target of compartment: {target_compartment.getId()}')
                    elif rpSBML.compareMIRIAMAnnotations(
                        source_annotation,
                        target_annotation,
                        self.logger
                    ):
                        found = True
                        break

            # If not found, add the compartment to the model
            if not found:
                target_compartment = self.getModel().createCompartment()
                rpSBML.checklibSBML(
                    target_compartment,
                    'Creating target compartment'
                )
                rpSBML.checklibSBML(
                    target_compartment.setMetaId(
                        source_compartment.getMetaId()
                    ),
                    'setting target metaId'
                )
                # make sure that the ID is different
                if source_compartment.getId() == target_compartment.getId():
                    rpSBML.checklibSBML(
                        target_compartment.setId(
                            source_compartment.getId()+'_sourceModel'
                        ),
                        'setting target id'
                    )
                else:
                    rpSBML.checklibSBML(
                        target_compartment.setId(
                            source_compartment.getId()
                        ),
                        'setting target id'
                    )
                rpSBML.checklibSBML(
                    target_compartment.setName(
                        source_compartment.getName()
                    ),
                    'setting target name'
                )
                rpSBML.checklibSBML(
                    target_compartment.setConstant(
                        source_compartment.getConstant()
                    ),
                    'setting target constant'
                )
                rpSBML.checklibSBML(
                    target_compartment.setAnnotation(
                        source_compartment.getAnnotation()
                    ),
                    'setting target annotation'
                )
                rpSBML.checklibSBML(
                    target_compartment.setSBOTerm(
                        source_compartment.getSBOTerm()
                    ),
                    'setting target annotation'
                    )

        # self.logger.debug('comp_source_target: '+str(comp_source_target))


    def copyParameters(
        self,
        source_sbml: 'rpSBML'
    ) -> None:

        source_sbml_doc = source_sbml.getDocument()

        self.logger.debug(f'source_sbml_doc: {source_sbml_doc}')

        target_paramsID = [i.getId() for i in self.getModel().getListOfParameters()]

        for source_parameter in source_sbml_doc.getModel().getListOfParameters():

            if source_parameter.getId() not in target_paramsID:
                target_parameter = self.getModel().createParameter()
                rpSBML.checklibSBML(
                    target_parameter,
                    'creating target parameter'
                )
                rpSBML.checklibSBML(
                    target_parameter.setId(
                        source_parameter.getId()
                    ),
                        'setting target parameter ID'
                )
                rpSBML.checklibSBML(
                    target_parameter.setSBOTerm(
                        source_parameter.getSBOTerm()
                    ),
                    'setting target parameter SBO'
                )
                rpSBML.checklibSBML(
                    target_parameter.setUnits(
                        source_parameter.getUnits()
                    ),
                    'setting target parameter Units'
                )
                rpSBML.checklibSBML(
                    target_parameter.setValue(
                        source_parameter.getValue()
                    ),
                    'setting target parameter Value'
                )
                rpSBML.checklibSBML(
                    target_parameter.setConstant(
                        source_parameter.getConstant()
                    ),
                    'setting target parameter ID'
                )

        # return target_paramsID


    def copyFBCGeneProducts(
        self,
        source_rpsbml: 'rpSBML'
    ) -> None:

        source_fbc = source_rpsbml.getModel().getPlugin('fbc')
        target_fbc = self.getModel().getPlugin('fbc')

        self.logger.debug(f'source_fbc: {source_fbc}')
        self.logger.debug(f'target_fbc: {target_fbc}')

        targetGenProductID = [i.getId() for i in target_fbc.getListOfGeneProducts()]

        for source_geneProduct in source_fbc.getListOfGeneProducts():

            if not source_geneProduct.getId() in targetGenProductID:
                target_geneProduct = target_fbc.createGeneProduct()
                rpSBML.checklibSBML(
                    target_geneProduct, 'creating target gene product')
                rpSBML.checklibSBML(
                    target_geneProduct.setId(
                        source_geneProduct.getId()
                    ),
                    'setting target gene product id'
                )
                rpSBML.checklibSBML(
                    target_geneProduct.setLabel(
                        source_geneProduct.getLabel()
                    ),
                    'setting target gene product label'
                )
                rpSBML.checklibSBML(
                    target_geneProduct.setName(
                        source_geneProduct.getName()
                    ),
                    'setting target gene product name'
                )
                rpSBML.checklibSBML(
                    target_geneProduct.setMetaId(
                        source_geneProduct.getMetaId()
                    ),
                    'setting target gene product meta_id'
                )

        # return targetGenProductID


    def copyFBCObjectives(
        self,
        source_rpsbml: 'rpSBML'
    ) -> None:
    # ) -> Tuple[List[str], List[str]]:

        # TODO: if overlapping id's need to replace the id with modified, as for the species

        source_fbc = source_rpsbml.getModel().getPlugin('fbc')
        target_fbc = self.getModel().getPlugin('fbc')

        self.logger.debug(f'source_fbc: {source_fbc}')
        self.logger.debug(f'target_fbc: {target_fbc}')

        targetObjectiveID = [i.getId() for i in target_fbc.getListOfObjectives()]
        # sourceObjectiveID = [i.getId() for i in source_fbc.getListOfObjectives()]

        for source_objective in source_fbc.getListOfObjectives():

            if not source_objective.getId() in targetObjectiveID:
                target_objective = target_fbc.createObjective()
                rpSBML.checklibSBML(
                    target_objective,
                    'creating target objective'
                )
                rpSBML.checklibSBML(
                    target_objective.setId(
                        source_objective.getId()
                    ),
                    'setting target objective'
                )
                rpSBML.checklibSBML(
                    target_objective.setName(
                        source_objective.getName()
                    ),
                    'setting target objective'
                )
                rpSBML.checklibSBML(
                    target_objective.setType(
                        source_objective.getType()
                    ),
                    'setting target objective type'
                )
                for source_fluxObjective in source_objective.getListOfFluxObjectives():
                    target_fluxObjective = target_objective.createFluxObjective()
                    rpSBML.checklibSBML(
                        target_fluxObjective,
                        'creating target flux objective'
                    )
                    rpSBML.checklibSBML(
                        target_fluxObjective.setName(
                            source_fluxObjective.getName()
                        ),
                        'setting target flux objective name'
                    )
                    rpSBML.checklibSBML(
                        target_fluxObjective.setCoefficient(
                            source_fluxObjective.getCoefficient()
                        ),
                        'setting target flux objective coefficient'
                    )
                    rpSBML.checklibSBML(
                        target_fluxObjective.setReaction(
                            source_fluxObjective.getReaction()
                        ),
                        'setting target flux objective reaction'
                    )
                    rpSBML.checklibSBML(
                        target_fluxObjective.setAnnotation(
                            source_fluxObjective.getAnnotation()
                        ),
                        'setting target flux obj annotation from source flux obj'
                    )
                rpSBML.checklibSBML(
                    target_objective.setAnnotation(
                        source_objective.getAnnotation()
                    ),
                    'setting target obj annotation from source obj'
                )
        # self.logger.debug('targetObjectiveID: '+str(targetObjectiveID))
        # self.logger.debug('sourceObjectiveID: '+str(sourceObjectiveID))

        # return sourceObjectiveID, targetObjectiveID

    def copySpecies(
        self,
        source_sbml: 'rpSBML',
        compartment_id: str,
    ) -> Tuple[List[str], List[str]]:

        self.logger.debug(f'compartment_id: {compartment_id}')

        # 'corr_species' is the dictionary of correspondace
        # between species both in pathway and model.
        species_ids = [
            spe.getId()
            for spe in list(
                source_sbml.getModel().getListOfSpecies()
            )
        ]
        corr_species, miss_species = self.speciesMatchWith(
            species_ids=species_ids,
            compartment_id=compartment_id
        )
        self.logger.debug(f'Species found in the model: {list(corr_species.keys())}')
        self.logger.debug(f'Species not found in the model: {miss_species}')

        source_species_ids = [spe.id for spe in source_sbml.getModel().getListOfSpecies()]
        # target_species_ids = [spe.id for spe in self.getModel().getListOfSpecies()]

        for source_spe_id in source_species_ids:

            # list_target = [i for i in corr_species[source_species]]
            # if source_species in list_target:
            #     self.logger.debug('The source ('+str(source_species)+') and target species ids ('+str(list_target)+') are the same')

            # If match, copy the BRSynth annotation from the source to the target
            if source_spe_id in corr_species:
                self.logger.debug(f'{source_spe_id} IN MODEL')
                # list_species = [i for i in corr_species[source_species]]
                # # self.logger.debug('list_species: '+str(list_species))
                # if len(list_species)==0:
                #     continue
                #     # self.logger.warning('Source species '+str(member.getIdRef())+' has been created in the target model')
                # elif len(list_species)>1:
                #     self.logger.debug('There are multiple matches to the species '+str(source_species)+'... taking the first one: '+str(list_species))
                # # TODO: loop throught the annotations and replace the non-overlapping information
                target_member = self.getModel().getSpecies(corr_species[source_spe_id])
                source_member = source_sbml.getModel().getSpecies(source_spe_id)
                # print(target_member.toXMLNode().toXMLString())
                # print('-------------------------------------')
                # print(source_member.toXMLNode().toXMLString())
                rpSBML.checklibSBML(
                    target_member,
                    f'Retrieving the target species: {corr_species[source_spe_id]}'
                )
                rpSBML.checklibSBML(
                    source_member,
                    f'Retrieving the source species: {source_spe_id}'
                )
                rpSBML.checklibSBML(
                    target_member.getAnnotation(
                        ).getChild(
                            'RDF'
                        ).addChild(
                            source_member.getAnnotation(
                            ).getChild(
                                'RDF'
                            ).getChild(
                                'BRSynth'
                            )
                    ),
                    'Replacing the annotations'
                )
                source_member.setId(corr_species[source_spe_id])
                # print('-------------------------------------')
                # print(target_member.toXMLNode().toXMLString())
                # print('-------------------------------------')
                # print(source_member.toXMLNode().toXMLString())
                # print('=====================================')

            # if no match then add it to the target model
            # with the corresponding compartment:
            else:
                self.logger.debug(f'{source_spe_id} OUT OF MODEL')
                # self.logger.debug('Creating source species '+str(source_species)+' in target rpsbml')
                source_spe = source_sbml.getModel().getSpecies(source_spe_id)
                if not source_spe:
                    self.logger.error(f'Cannot retreive model species: {source_spe_id}')
                else:
                    rpSBML.checklibSBML(
                        source_spe,
                        'fetching source species'
                    )
                    target_spe = self.getModel().createSpecies()
                    rpSBML.checklibSBML(
                        target_spe,
                        'creating species'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setMetaId(
                            source_spe.getMetaId()
                        ),
                        'setting target metaId'
                    )
                    # ## need to check if the id of the source species does not already exist in the target model
                    # if source_spe.getId() in target_species_ids:
                    #     target_species_id = source_sbml.getModel().id+'__'+str(source_spe.getId())
                    #     if not source_spe.getId() in corr_species:
                    #         corr_species[source_spe.getId()] = {}
                    #     corr_species[source_spe.getId()][source_sbml.getModel().id+'__'+str(source_spe.getId())] = 1.0
                    # else:
                    #     target_species_id = source_spe.getId()
                    rpSBML.checklibSBML(
                        target_spe.setId(
                            source_spe.getId()
                        ),
                        'setting target id'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setCompartment(
                            compartment_id
                        ),
                        'setting target compartment'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setInitialConcentration(
                            source_spe.getInitialConcentration()
                        ),
                        'setting target initial concentration'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setBoundaryCondition(
                            source_spe.getBoundaryCondition()
                        ),
                        'setting target boundary concentration'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setHasOnlySubstanceUnits(
                            source_spe.getHasOnlySubstanceUnits()
                        ),
                        'setting target has only substance unit'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setBoundaryCondition(
                            source_spe.getBoundaryCondition()
                        ),
                        'setting target boundary condition'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setConstant(
                            source_spe.getConstant()
                        ),
                        'setting target constant'
                    )
                    rpSBML.checklibSBML(
                        target_spe.setAnnotation(
                            source_spe.getAnnotation()
                        ),
                        'setting target annotation'
                    )

        return corr_species, miss_species


    @staticmethod
    def reactionsInBoth(
        rpsbml_1: 'rpSBML',
        rpsbml_2: 'rpSBML',
        species: Dict,
        logger: Logger = getLogger(__name__)
    ) -> Dict:

        reactions_in_both = {}

        for rxn_1 in rpsbml_1.getModel().getListOfReactions():

            logger.debug('rxn_1: ' + str(rxn_1))

            for rxn_2 in rpsbml_2.getModel().getListOfReactions():
                match = rpSBML.reactionsAreEqual(
                    species,
                    rxn_1,
                    rxn_2,
                    logger = logger
                )
                if match:
                    reactions_in_both[rxn_1.getId()] = rxn_2.getId()
        
        return reactions_in_both


    def copyReactions(
        self,
        source_sbml: 'rpSBML',
        species: Dict
    ) -> Dict:

        self.logger.debug(f'source_sbml: {source_sbml}')
        self.logger.debug(f'species: {species}')

        def copySpecies(
            source_reaction: libsbml.Reaction,
            target_reaction: libsbml.Reaction,
            species: Dict,
            spe_type: str,
            logger: Logger = getLogger(__name__)
        ) -> None:

            logger.debug('source_reaction: ' + str(source_reaction))
            logger.debug('target_reaction: ' + str(target_reaction))
            logger.debug('species: ' + str(species))
            logger.debug('spe_type: ' + str(spe_type))

            source_species_l = getattr(source_reaction, f'getListOf{spe_type.capitalize()}s')()
            for source_species in source_species_l:
                src_spe_id = source_species.species
                target_species = getattr(target_reaction, f'create{spe_type.capitalize()}')()
                rpSBML.checklibSBML(
                    target_species,
                    f'create target {spe_type}'
                )
                if src_spe_id in species:
                    rpSBML.checklibSBML(
                        target_species.setSpecies(
                            species[src_spe_id]
                        ),
                        f'assign {spe_type} species'
                    )
                else:
                    rpSBML.checklibSBML(
                        target_species.setSpecies(
                            src_spe_id
                        ),
                        f'assign {spe_type} species'
                    )
                rpSBML.checklibSBML(
                    source_species,
                    f'fetch source {spe_type}'
                )
                rpSBML.checklibSBML(
                    target_species.setConstant(
                        source_species.getConstant()
                    ),
                    'set "constant" on species '+str(source_species.getConstant())
                )
                rpSBML.checklibSBML(
                    target_species.setStoichiometry(
                        source_species.getStoichiometry()
                    ),
                    'set stoichiometry ('+str(source_species.getStoichiometry)+')'
                )

        def copyProducts(
            source_reaction: libsbml.Reaction,
            target_reaction: libsbml.Reaction,
            species: Dict,
            logger: Logger = getLogger(__name__)
        ) -> None:

            logger.debug('source_reaction: ' + str(source_reaction))
            logger.debug('target_reaction: ' + str(target_reaction))
            logger.debug('species: ' + str(species))

            for source_reaction_productID in [i.species for i in source_reaction.getListOfProducts()]:
                # self.logger.debug('\tAdding '+str(source_reaction_productID))
                target_product = target_reaction.createProduct()
                rpSBML.checklibSBML(
                    target_product,
                    'create target reactant'
                )
                if source_reaction_productID in species:
                    if species[source_reaction_productID]:
                        if len(species[source_reaction_productID]) > 1:
                            logger.warning('Multiple matches for '+str(source_reaction_productID)+': '+str(species[source_reaction_productID]))
                            logger.warning('Taking one arbitrarely')
                        # WARNING: taking the first one arbitrarely
                        rpSBML.checklibSBML(
                            target_product.setSpecies(
                                [i for i in species[source_reaction_productID]][0]
                            ),
                            'assign reactant product'
                        )
                    else:
                        rpSBML.checklibSBML(
                            target_product.setSpecies(
                                source_reaction_productID
                            ),
                            'assign reactant product'
                        )
                else:
                    rpSBML.checklibSBML(
                        target_product.setSpecies(
                            source_reaction_productID
                        ),
                        'assign reactant product'
                    )
                source_product = source_reaction.getProduct(source_reaction_productID)
                rpSBML.checklibSBML(
                    source_product,
                    'fetch source reactant'
                )
                rpSBML.checklibSBML(
                    target_product.setConstant(
                        source_product.getConstant()
                    ),
                    'set "constant" on product '+str(source_product.getConstant())
                )
                rpSBML.checklibSBML(
                    target_product.setStoichiometry(
                        source_product.getStoichiometry()
                    ),
                    'set stoichiometry ('+str(source_product.getStoichiometry)+')'
                )

            # return target_product

        def setParameters(
            source_reaction: libsbml.Reaction,
            target_reaction: libsbml.Reaction,
            logger: Logger = getLogger(__name__)
        ) -> None:
            # self.logger.debug('Cannot find source reaction: '+str(source_reaction.getId()))
            rpSBML.checklibSBML(
                source_reaction,
                'fetching source reaction'
            )
            rpSBML.checklibSBML(
                target_reaction,
                'create reaction'
            )
            target_fbc = target_reaction.getPlugin('fbc')
            rpSBML.checklibSBML(
                target_fbc,
                'fetching target FBC package'
            )
            source_fbc = source_reaction.getPlugin('fbc')
            rpSBML.checklibSBML(
                source_fbc,
                'fetching source FBC package'
            )
            source_upperFluxBound = source_fbc.getUpperFluxBound()
            rpSBML.checklibSBML(
                source_upperFluxBound,
                'fetching upper flux bound'
            )
            rpSBML.checklibSBML(
                target_fbc.setUpperFluxBound(
                    source_upperFluxBound
                ),
                'setting upper flux bound'
            )
            source_lowerFluxBound = source_fbc.getLowerFluxBound()
            rpSBML.checklibSBML(
                source_lowerFluxBound,
                'fetching lower flux bound'
            )
            rpSBML.checklibSBML(
                target_fbc.setLowerFluxBound(
                    source_lowerFluxBound
                ),
                'setting lower flux bound'
            )
            rpSBML.checklibSBML(
                target_reaction.setId(
                    source_reaction.getId()
                ),
                'set reaction id'
            )
            rpSBML.checklibSBML(
                target_reaction.setName(
                    source_reaction.getName()
                ),
                'set name'
            )
            rpSBML.checklibSBML(
                target_reaction.setSBOTerm(
                    source_reaction.getSBOTerm()
                ),
                'setting the reaction system biology ontology (SBO)'
            ) # set as process
            # TODO: consider having the two parameters as input to the function
            rpSBML.checklibSBML(
                target_reaction.setReversible(
                    source_reaction.getReversible()
                ),
                'set reaction reversibility flag'
            )
            rpSBML.checklibSBML(
                target_reaction.setFast(
                    source_reaction.getFast()
                ),
                'set reaction "fast" attribute'
            )
            rpSBML.checklibSBML(
                target_reaction.setMetaId(
                    source_reaction.getMetaId()
                ),
                'setting species meta_id'
            )
            rpSBML.checklibSBML(
                target_reaction.setAnnotation(
                    source_reaction.getAnnotation()
                ),
                'setting annotation for source reaction'
            )

        reactions_in_both = {}

        # Convert species from pathway to model naming
        # source_reactions = {}
        for reaction in source_sbml.getModel().getListOfReactions():
            # source_reactions[reaction.getId()] = {
            #     'reactants': [],
            #     'products': []
            # }
            # Reactants
            for spe in reaction.getListOfReactants():
                if spe.getSpecies() in species:
                    # source_reactions[reaction.getId()]['reactants'].append(
                    #     species[spe.species]
                    # )
                    # rename species in the source (pathway)
                    spe.setSpecies(species[spe.getSpecies()])
                # else:
                #     source_reactions[reaction.getId()]['reactants'].append(spe.species)
            # Products
            for spe in reaction.getListOfProducts():
                if spe.getSpecies() in species:
                    spe.setSpecies(species[spe.getSpecies()])
                #     source_reactions[reaction.getId()]['products'].append(
                #         species[spe.species]
                #     )
                # else:
                #     source_reactions[reaction.getId()]['products'].append(spe.species)

        for source_reaction in source_sbml.getModel().getListOfReactions():

            self.logger.debug(f'source_reaction: {source_reaction.getId()}')

            is_found = False

            for target_reaction in self.getModel().getListOfReactions():
                # _target_reaction = {
                #     'reactants': [spe.species for spe in target_reaction.getListOfReactants()],
                #     'products': [spe.species for spe in target_reaction.getListOfProducts()]
                # }
                # match = (
                #     len(source_reaction_infos['reactants']) == len(_target_reaction['reactants'])
                #     and len(source_reaction_infos['products']) == len(_target_reaction['products'])
                #     and sorted(source_reaction_infos['reactants']) == sorted(_target_reaction['reactants'])
                #     and sorted(source_reaction_infos['products']) == sorted(_target_reaction['products'])
                # )
                if rpSBML.reactionsAreEqual(
                    source_reaction,
                    target_reaction
                ):
                # match = rpSBML.reactionsAreEqual(
                #     species,
                #     source_reaction,
                #     ,
                #     logger = self.logger
                # )
                # if match:
                    # self.logger.debug('Source reaction '+str(source_reaction)+' matches with target reaction '+str(target_reaction))
                    # source_reaction[source_reaction.getId()] = target_reaction.getId()
                    reactions_in_both[source_reaction.getId()] = target_reaction.getId()
                    is_found = True
                    break

            if not is_found:
                self.getModel().addReaction(
                    source_sbml.getModel().getReaction(source_reaction.getId())
                )

                # source_reaction = source_sbml.getModel().getReaction(source_reaction_id)
                # target_reaction = self.getModel().createReaction()
                # setParameters(
                #     source_reaction=source_reaction,
                #     target_reaction=target_reaction,
                #     logger=self.logger
                #     )
                # # Reactants
                # copySpecies(
                #     source_reaction=source_reaction,
                #     target_reaction=target_reaction,
                #     species=species,
                #     spe_type='reactant',
                #     logger=self.logger
                # )
                # # Products
                # # NOTE: source_reaction_reactantID returned is the last one of a loop
                # #       in rpSBML.getReactants() but used below in copyProducts()
                # copySpecies(
                #     source_reaction=source_reaction,
                #     target_reaction=target_reaction,
                #     species=species,
                #     spe_type='product',
                #     logger=self.logger
                # )
                # print(target_reaction.toXMLNode().toXMLString())
                # print('--------------------------------------------')
                # print(source_reaction.clone().toXMLNode().toXMLString())
                # print('============================================')

        return reactions_in_both


    def copyGroups(
        self,
        source_sbml: 'rpSBML',
        species: Dict,
        reactions: Dict
    ) -> None:

        source_sbml_doc = source_sbml.getDocument()

        self.logger.debug(f'source_sbml_doc: {source_sbml_doc}')
        self.logger.debug(f'species: {species}')
        self.logger.debug(f'reactions: {reactions}')

        # TODO loop through the groups to add them

        if not self.getModel().isPackageEnabled('groups'):
            rpSBML.checklibSBML(
                self.getModel().enablePackage(
                    'http://www.sbml.org/sbml/level3/version1/groups/version1',
                    'groups',
                    True
                ),
                'Enabling the GROUPS package'
            )
        # !!!! must be set to false for no apparent reason
        rpSBML.checklibSBML(
            source_sbml_doc.setPackageRequired(
                'groups',
                False
            ),
            'enabling groups package'
        )
        source_groups = source_sbml_doc.getModel().getPlugin('groups')
        rpSBML.checklibSBML(
            source_groups,
            'fetching the source model groups'
        )
        target_groups = self.getModel().getPlugin('groups')
        rpSBML.checklibSBML(
            target_groups,
            'fetching the target model groups'
        )

        self.logger.debug(f'species: {species}')
        self.logger.debug(f'reactions: {reactions}')

        source_groups_ids = [i.id for i in source_groups.getListOfGroups()]
        target_groups_ids = [i.id for i in target_groups.getListOfGroups()]
        # NOTE: only need to update the source species since these are the ones that are replaced with their equivalent
        for source_group in source_groups.getListOfGroups():
            # overwrite in the group the reaction members that have been replaced
            for member in source_group.getListOfMembers():
                if member.getIdRef() in reactions:
                    if reactions[member.getIdRef()]:
                        member.setIdRef(reactions[member.getIdRef()])
            # overwrite in the group the species members that have been replaced
            for member in source_group.getListOfMembers():
                if member.getIdRef() in species:
                    # if species[member.getIdRef()]:
                    # list_species = [i for i in species[member.getIdRef()]]
                    # self.logger.debug('species: '+str(species))
                    # self.logger.debug('list_species: '+str(list_species))
                    # if len(list_species)==0:
                    #     continue
                    #     # self.logger.warning('Source species '+str(member.getIdRef())+' has been created in the target model')
                    # elif len(list_species)>1:
                    #     self.logger.warning('There are multiple matches to the species '+str(member.getIdRef())+'... taking the first one: '+str(list_species))
                    rpSBML.checklibSBML(
                        member.setIdRef(species[member.getIdRef()]),
                        'Setting name to the groups member'
                    )
            # create and add the groups if a source group does not exist in the target
            if not source_group.id in target_groups_ids:
                rpSBML.checklibSBML(
                    target_groups.addGroup(source_group),
                    'copy the source groups to the target groups'
                )
            # if the group already exists in the target then need to add new members
            else:
                target_group = target_groups.getGroup(source_group.id)
                target_group_ids = [i.getIdRef() for i in target_group.getListOfMembers()]
                for member in source_group.getListOfMembers():
                    if member.getIdRef() not in target_group_ids:
                        new_member = target_group.createMember()
                        rpSBML.checklibSBML(
                            new_member,
                            'Creating a new groups member'
                        )
                        rpSBML.checklibSBML(
                            new_member.setIdRef(member.getIdRef()),
                            'Setting name to the groups member'
                        )

        # return source_groups, target_groups


    def getListOfGroups(self) -> List[libsbml.Group]:

        return list(self.getPlugin('groups').getListOfGroups())


    def getGroup(
        self,
        group_id: str
    ) -> libsbml.Group:

        # model_plugin = self.getPlugin('groups')

        # print(list(model_plugin.getListOfGroups()))
        # print(model_plugin.getElementBySId())
        # exit()

        return self.getPlugin('groups').getElementBySId(group_id)

        # if group is None:
        #     self.logger.warning('The group '+str(group_id)+' does not exist... creating it')
        #     group = self.createGroup(group_id)
        
        # return group

        # rp_pathway = groups.getGroup(pathway_id)
    
        # if rp_pathway is None:
        #     self.logger.warning('The group '+str(pathway_id)+' does not exist... creating it')
        #     self.createGroup(pathway_id)
        #     rp_pathway = groups.getGroup(pathway_id)
    
        # self.checklibSBML(rp_pathway, 'Getting RP pathway')
    
        # return rp_pathway


    def copyTitles(
        self,
        source_sbml: 'rpSBML',
        logger: Logger = getLogger(__name__)
    ) -> None:

        source_sbml_doc = source_sbml.getDocument()

        logger.debug('source_sbml_doc: ' + str(source_sbml_doc))

        self.getModel().setId(
            self.getModel().getId()+'__'+source_sbml_doc.getModel().getId()
        )
        self.getModel().setName(
            self.getModel().getName()+' merged with '+source_sbml_doc.getModel().getId()
        )


    def get_isolated_species(self) -> List[str]:
        """
        Return speices that are isolated, i.e. are only consumed or only produced.

        Returns
        -------
            List of isolated species
        """
        try:
            return self.isolated_species
        except AttributeError:
            return []


    def set_isolated_species(self, species: List[str]) -> None:
        """
        Set isolated species.

        Parameters
        ----------
            species: List [str]
        """
        self.isolated_species = deepcopy(species)

    # @staticmethod
    def search_isolated_species(
    # def completeHeterologousPathway(
        self,
        species: List[str] = [],
        pathway_id: str = 'rp_pathway',
        central_species_group_id: str = 'rp_trunk_species',
        sink_species_group_id: str = 'rp_sink_species'
    ) -> bool:
        """Check if there are any single parent species in a heterologous pathways and if there are, either delete them or add reaction to complete the heterologous pathway

        :param rpsbml: The rpSBML object
        :param upper_flux_bound: The upper flux bounds unit definitions default when adding new reaction (Default: 999999.0)
        :param lower_flux_bound: The lower flux bounds unit definitions default when adding new reaction (Defaul: 0.0)
        :param compartment_id: The id of the model compartment
        :param pathway_id: The pathway ID (Default: rp_pathway)
        :param central_species_group_id: The central species Groups id (Default: central_species)
        :param sink_species_group_id: The sink specues Groups id (Default: sink_species_group_id)

        :type rpsbml: rpSBML
        :type upper_flux_bound: float
        :type lower_flux_bound: float
        :type compartment_id: str
        :type pathway_id: str
        :type central_species_group_id: str
        :type sink_species_group_id: str

        :rtype: bool
        :return: Success of failure of the function
        """

        # from cobra.io.sbml       import validate_sbml_model
        # from json import dumps as json_dumps
        # merged_rpsbml_path = '/tmp/merged.sbml'
        # from inspect import currentframe, getframeinfo
        # self.writeToFile(merged_rpsbml_path)
        # (model, errors) = validate_sbml_model(merged_rpsbml_path)
        # if model is None:
        #     frameinfo = getframeinfo(currentframe())
        #     print(frameinfo.filename, frameinfo.lineno)
        #     self.logger.error('Something went wrong reading the SBML model')
        #     self.logger.error(str(json_dumps(errors, indent=4)))
        #     exit()

        rpgraph = rpGraph(
            self,
            True,
            pathway_id,
            central_species_group_id,
            sink_species_group_id,
            logger = self.logger
        )

        self.isolated_species = list(
            set(
                rpgraph.onlyConsumedSpecies(species) +
                rpgraph.onlyProducedSpecies(species)
            )
        )

        # return self.get_isolated_species()

        # for pro in produced_species_nid:

        #     step = {
        #         'rule_id': None,
        #         'left': {pro.split('__')[0]: 1},
        #         'right': {},
        #         'transformation_id': None,
        #         'rule_score': None,
        #         'tmpl_rxn_id': None
        #     }
        #     # note that here the pathways are passed as NOT being part of the heterologous pathways and
        #     # thus will be ignored when/if we extract the rp_pathway from the full GEM model
        #     self.createReaction(
        #         pro+'__consumption',
        #         upper_flux_bound,
        #         lower_flux_bound,
        #         step,
        #         compartment_id
        #     )
        #     # print(pro)
        #     # self.writeToFile(merged_rpsbml_path)
        #     # (model, errors) = validate_sbml_model(merged_rpsbml_path)
        #     # if model is None:
        #     #     frameinfo = getframeinfo(currentframe())
        #     #     print(frameinfo.filename, frameinfo.lineno)
        #     #     self.logger.error('Something went wrong reading the SBML model')
        #     #     self.logger.error(str(json_dumps(errors, indent=4)))
        #     #     exit()

        # for react in consumed_species_nid:

        #     step = {
        #         'rule_id': None,
        #         'left': {},
        #         'right': {react.split('__')[0]: 1},
        #         'transformation_id': None,
        #         'rule_score': None,
        #         'tmpl_rxn_id': None}
        #     #note that here the pathwats are passed as NOT being part of the heterologous pathways and
        #     #thus will be ignored when/if we extract the rp_pathway from the full GEM model
        #     self.createReaction(
        #         react+'__production',
        #         upper_flux_bound,
        #         lower_flux_bound,
        #         step,
        #         compartment_id
        #     )

        # # self.writeToFile(merged_rpsbml_path)
        # # (model, errors) = validate_sbml_model(merged_rpsbml_path)
        # # if model is None:
        # #     frameinfo = getframeinfo(currentframe())
        # #     print(frameinfo.filename, frameinfo.lineno)
        # #     self.logger.error('Something went wrong reading the SBML model')
        # #     self.logger.error(str(json_dumps(errors, indent=4)))
        # #     return None

        # return True


    @staticmethod
    def _findUniqueRowColumn(pd_matrix, logger=getLogger(__name__)):
        """Private function that takes the matrix of similarity scores between the reactions or species of two models and finds the unqiue matches

        pd_matrix is organised such that the rows are the simulated species and the columns are the measured ones

        :param pd_matrix: Matrix of reactions or species of two models

        :type pd_matrix: np.array

        :return: Dictionary of matches
        :rtype: dict
        """
        
        # self.logger.debug(pd_matrix)
        to_ret = {}
        ######################## filter by the global top values ################
        # self.logger.debug('################ Filter best #############')
        # transform to np.array
        x = pd_matrix.values
        # resolve the rouding issues to find the max
        x = np.around(x, decimals=5)
        # first round involves finding the highest values and if found set to 0.0 the rows and columns (if unique)
        top = np.where(x == np.max(x))
        # as long as its unique keep looping
        if np.count_nonzero(x)==0:
            return to_ret
        while len(top[0])==1 and len(top[1])==1:
            if np.count_nonzero(x)==0:
                return to_ret
            pd_entry = pd_matrix.iloc[[top[0][0]],[top[1][0]]]
            row_name = str(pd_entry.index[0])
            col_name = str(pd_entry.columns[0])
            # if col_name in to_ret:
                # self.logger.debug('Overwriting (1): '+str(col_name))
                # self.logger.debug(x)
            to_ret[col_name] = [row_name]
            # delete the rows and the columns
            # self.logger.debug('==================')
            # self.logger.debug('Column: '+str(col_name))
            # self.logger.debug('Row: '+str(row_name))
            pd_matrix.loc[:, col_name] = 0.0
            pd_matrix.loc[row_name, :] = 0.0
            x = pd_matrix.values
            x = np.around(x, decimals=5)
            top = np.where(x == np.max(x))
            # self.logger.debug(pd_matrix)
            # self.logger.debug(top)
            # self.logger.debug('==================')
        #################### filter by columns (measured) top values ##############
        # self.logger.debug('################ Filter by column best ############')
        x = pd_matrix.values
        x = np.around(x, decimals=5)
        if np.count_nonzero(x)==0:
            return to_ret
        reloop = True
        while reloop:
            if np.count_nonzero(x)==0:
                return to_ret
            reloop = False
            for col in range(len(x[0])):
                if np.count_nonzero(x[:,col])==0:
                    continue
                top_row = np.where(x[:,col]==np.max(x[:,col]))[0]
                if len(top_row)==1:
                    top_row = top_row[0]
                    # if top_row == 0.0:
                    #    continue
                    # check to see if any other measured pathways have the same or larger score (accross)
                    row = list(x[top_row, :])
                    # remove current score consideration
                    row.pop(col)
                    if max(row)>=x[top_row, col]:
                        logger.warning('For col '+str(col)+' there are either better or equal values: '+str(row))
                        logger.warning(x)
                        continue
                    # if you perform any changes on the rows and columns, then you can perform the loop again
                    reloop = True
                    pd_entry = pd_matrix.iloc[[top_row],[col]]
                    # self.logger.debug('==================')
                    row_name = pd_entry.index[0]
                    col_name = pd_entry.columns[0]
                    # self.logger.debug('Column: '+str(col_name))
                    # self.logger.debug('Row: '+str(row_name))
                    # if col_name in to_ret:
                        # self.logger.debug('Overwriting (2): '+str(col_name))
                        # self.logger.debug(pd_matrix.values)
                    to_ret[col_name] = [row_name]
                    # delete the rows and the columns
                    pd_matrix.loc[:, col_name] = 0.0
                    pd_matrix.loc[row_name, :] = 0.0
                    x = pd_matrix.values
                    x = np.around(x, decimals=5)
                    # self.logger.debug(pd_matrix)
                    # self.logger.debug('==================')
        ################## laslty if there are multiple values that are not 0.0 then account for that ######
        # self.logger.debug('################# get the rest ##########')
        x = pd_matrix.values
        x = np.around(x, decimals=5)
        if np.count_nonzero(x)==0:
            return to_ret
        for col in range(len(x[0])):
            if not np.count_nonzero(x[:,col])==0:
                top_rows = np.where(x[:,col]==np.max(x[:,col]))[0]
                if len(top_rows)==1:
                    top_row = top_rows[0]
                    pd_entry = pd_matrix.iloc[[top_row],[col]]
                    row_name = pd_entry.index[0]
                    col_name = pd_entry.columns[0]
                    if col_name not in to_ret:
                        to_ret[col_name] = [row_name]
                    else:
                        logger.warning('At this point should never have only one: '+str(x[:,col]))
                        logger.warning(x)
                else:
                    for top_row in top_rows:
                        pd_entry = pd_matrix.iloc[[top_row],[col]]
                        row_name = pd_entry.index[0]
                        col_name = pd_entry.columns[0]
                        if col_name not in to_ret:
                            to_ret[col_name] = []
                        to_ret[col_name].append(row_name)
        # self.logger.debug(pd_matrix)
        # self.logger.debug('###################')
        return to_ret


    ##########################################################################################
    #################################### REACTION ############################################
    ##########################################################################################


    @staticmethod
    def reactionsAreEqual(
        source_reaction: libsbml.Reaction,
        target_reaction: libsbml.Reaction,
        logger: Logger = getLogger(__name__)
    ) -> bool:
        """Compare two reactions and elect that they are the same
        if they have exactly the same reactants and products

        species_source_target: {'MNXM4__64__MNXC3': {'M_o2_c': 1.0}, 'MNXM10__64__MNXC3': {'M_nadh_c': 1.0}, 'CMPD_0000000003__64__MNXC3': {}, 'TARGET_0000000001__64__MNXC3': {}, 'MNXM188__64__MNXC3': {'M_anth_c': 1.0}, 'BC_32877__64__MNXC3': {'M_nh4_c': 0.8}, 'BC_32401__64__MNXC3': {'M_nad_c': 0.2}, 'BC_26705__64__MNXC3': {'M_h_c': 1.0}, 'BC_20662__64__MNXC3': {'M_co2_c': 1.0}}
        the first keys are the source compartment ids
        the second key is the source species id
        the value is the target species id
        Note that we assure that the match is 1:1 between species using the species match

        :param species_source_target: The comparison dictionary between the species of two SBML files
        :param source_reaction: The target reaction
        :param target_reaction: The source reaction

        :type species_source_target: dict
        :type source_reaction: libsbml.Reaction
        :type target_reaction: libsbml.Reaction

        :return: The score of the match and boolean if its a match or not
        :rtype: tuple
        """
        source_reactants = [spe.species for spe in source_reaction.getListOfReactants()]
        target_reactants = [spe.species for spe in target_reaction.getListOfReactants()]
        source_products = [spe.species for spe in source_reaction.getListOfProducts()]
        target_products = [spe.species for spe in target_reaction.getListOfProducts()]
        return (
            source_reaction.getNumReactants() == target_reaction.getNumReactants()
            and source_reaction.getNumProducts() == target_reaction.getNumProducts()
            and sorted(source_reactants) == sorted(target_reactants)
            and sorted(source_products) == sorted(target_products)
        )
        # print(
        #     target_reaction.getId(),
        #     [spe.species for spe in target_reaction.getListOfReactants()],
        #     [spe.species for spe in target_reaction.getListOfProducts()],
        #     )

        # print(
        #     ' + '.join([spe.species for spe in source_reaction.getListOfReactants()])
        #     + ' = '
        #     + ' + '.join([spe.species for spe in source_reaction.getListOfProducts()])
        # )
        # print(
        #     ' + '.join([spe.species for spe in target_reaction.getListOfReactants()])
        #     + ' = '
        #     + ' + '.join([spe.species for spe in target_reaction.getListOfProducts()])
        # )
        # print()
        # print()
        # def fill_targets(
        #     target_reaction: libsbml.Reaction,
        #     species_source_target: Dict,
        #     logger: Logger = getLogger(__name__)
        # ) -> List[libsbml.SpeciesReference]:
        #     targets = []
        #     for i in target_reaction.getListOfReactants():
        #         if i.species in species_source_target:
        #             if species_source_target[i.species]:
        #                 # WARNING: Taking the first one arbitrarely
        #                 conv_spe = [y for y in species_source_target[i.species]][0]
        #                 targets.append(conv_spe)
        #                 # scores.append(species_source_target[i.species][conv_spe])
        #             else:
        #                 targets.append(i.species)
        #                 # scores.append(1.0)
        #         else:
        #             targets.append(i.species)
        #             # scores.append(1.0)
        #     return targets

        # # REACTANTS
        # source_reactants = [i.species for i in source_reaction.getListOfReactants()]
        # target_reactants = fill_targets(
        #     target_reaction = target_reaction,
        #     species_source_target = species_source_target,
        #     logger = logger
        # )

        # ## PRODUCTS
        # source_products = [i.species for i in source_reaction.getListOfProducts()]
        # target_products = fill_targets(
        #     target_reaction = target_reaction,
        #     species_source_target = species_source_target,
        #     logger = logger
        # )


        # print(
        #     'REACTANTS',
        #     sorted(source_reaction['reactants']),
        #     sorted([i.species for i in target_reaction.getListOfReactants()])
        # )
        # print(
        #     'PRODUCTS',
        #     sorted(source_reaction['products']),
        #     sorted([i.species for i in target_reaction.getListOfProducts()])
        # )

        # return False

        # reactants = set(source_reactants) - set(target_reactants)
        # products  = set(source_products)  - set(target_products)

        # logger.debug('source_reactants: '+str(source_reactants))
        # logger.debug('target_reactants: '+str(target_reactants))
        # logger.debug('source_products: '+str(source_products))
        # logger.debug('target_products: '+str(target_products))
        # logger.debug('reactants: '+str(reactants))
        # logger.debug('products: '+str(products))

        # return reactants is products is None

        # if reactants is products is None:
        #     return np.mean(scores)
        # else:
        #     return None


    ##########################################################################################
    ##################################### SPECIES ############################################
    ##########################################################################################


    # TODO: for all the measured species compare with the simualted one. Then find the measured
    # and simulated species that match the best and exclude the
    # simulated species from potentially matching with another
    def speciesMatchWith(
        self,
        species_ids: List[str],
        compartment_id: str
    ) -> Tuple[
        Dict[str, str],
        List[str]
    ]:
        """Match all the measured chemical species to the simulated chemical species between two SBML

        :param species_ids: Compounds to search in the target rpSBML
        :param target_rpsbml: The target rpSBML
        :param compartment_id: The id of the compartment into perform the search
        :param logger: A logging object to output information

        :type species_ids: List[str]
        :type target_rpsbml: rpSBML
        :type compartment_id: str
        :type logger: logging

        :return: A tuple corresponding to the dictionnary correspondance between species provided
        and specie in the model and a list of species not find in the model
        :rtype: Tuple[Dict[str,str], List[str]]
        """
        self.logger.debug(f'species_ids: {species_ids}')
        self.logger.debug(f'compartment_id: {compartment_id}')

        # Correspondance table for matched species
        # between source and target to rename them later
        corr_species = {}
        miss_species = set()

        if self.getModel() is None:
            return corr_species, list(miss_species)

        # Get the species IDs of the right compartment
        compartment_species = [
            spe
            for spe in
            self.getModel().getListOfSpecies()
            if spe.getCompartment() == compartment_id
        ]

        # In all reactions of the model,
        # count the occurence of each species
        # of the given compartment.
        # The score is the number of reactions
        # in which the species are involved
        species_occurence = {}
        for rxn in self.getModel().getListOfReactions():
            for spe in compartment_species:
                # Look if the species is a reactant or a product
                if spe.getId() in [
                    reactant.species
                    for reactant in rxn.getListOfReactants()
                ] or spe.getId() in [
                    product.species
                    for product in rxn.getListOfProducts()
                ]:
                    if spe.getId() in species_occurence:
                        species_occurence[spe.getId()] += 1
                    else:
                        species_occurence[spe.getId()] = 1
                else:
                    if spe.getId() not in species_occurence:
                        species_occurence[spe.getId()] = 0

        for species_id_to_match in species_ids:

            # Search for match(es) in the current sbml model
            for model_species in compartment_species:

                self.logger.debug(f'species_id_to_match/compartment: {species_id_to_match}/{compartment_id}; model_species/compartment: {model_species.getId()}/{model_species.getCompartment()}')

                # Look if the specie in source has the same ID as in the target
                if species_id_to_match == model_species.getId():
                    corr_species[species_id_to_match] = model_species.getId()
                # Else look if the specie in source has the same ID
                # as one ID among MIRIAM annotations
                else:
                    # Try to match the species with MIRIAM annotations of the model species
                    _corr_species = rpSBML.matchWithMIRIAM(
                        species_id_to_match=species_id_to_match,
                        model_species=model_species,
                        scores=species_occurence,
                        prev_matches=corr_species,
                        logger=self.logger
                    )
                    if _corr_species:
                        corr_species[species_id_to_match] = _corr_species

            # If species not found in the model, add it to missing species list
            if species_id_to_match not in corr_species:
                miss_species.add(species_id_to_match)

        return corr_species, list(miss_species)

    @staticmethod
    def matchWithMIRIAM(
        species_id_to_match: str,
        model_species: libsbml.Species,
        scores: Dict[str, float],
        prev_matches: Dict[str, str],
        logger: Logger = getLogger(__name__)
    ) -> str:
        """Match a species with a MIRIAM annotation.
        """

        miriam_annot = rpSBML.readMIRIAMAnnotation(model_species.getAnnotation(), logger)

        class MatchSpecies(Exception):
            pass

        try:
            annot_dict = rpSBML.convert_miriam_to_dict(miriam_annot)
            # For all (lists of) cross refs in MIRIAM annot
            for db, xrefs in annot_dict.items():
                # cross_ref = cross_refs.split('/')[-1]
                # # # For all single cross ref in a DB list
                # # for cross_ref in cross_refs:
                logger.debug(f'Comparing {species_id_to_match} with {xrefs} (MIRIAM)')
                # Found a match between species ID and the current cross ref
                if species_id_to_match in xrefs:
                    curr_match = model_species.getId()
                    curr_score = scores[curr_match]
                    logger.debug(f'Found match: {species_id_to_match} -> {curr_match} (involved in {curr_score} reactions)')
                    # Species ID already matched with another species
                    if species_id_to_match in prev_matches:
                        prev_match = prev_matches[species_id_to_match]
                        prev_score = scores[prev_match]
                        logger.warning(
                            f'*** A match already exists for {species_id_to_match}'
                        )
                        logger.warning(
                            f'    \__ previous match: {species_id_to_match} -> {prev_match} (involved in {prev_score} reactions)'
                        )
                        logger.warning(
                            f'    \__  current match: {species_id_to_match} -> {curr_match} (involved in {curr_score} reactions)'
                        )
                        # Compare the scores of the previous and the current match
                        if curr_score > prev_score:
                            logger.warning(
                                f'--> Keep {species_id_to_match} -> {curr_match}'
                            )
                            return curr_match
                        else:
                            logger.warning(
                                f'--> Keep {species_id_to_match} -> {prev_match}'
                            )
                    else:
                        return curr_match
                    # Stop the search in the current MIRIAM annotations
                    raise MatchSpecies
        except MatchSpecies:
            pass

        return None

    def is_boundary_type(
        self,
        reaction: libsbml.Reaction, 
        boundary_type: str, 
        external_compartment: str
    ) -> bool:
        '''Check whether a reaction is an exchange reaction.
        Adapted from "is_boundary_type" available at cobra.medium.boudary_types to fit with libsbml.Reaction.

        :param reaction: a reaction to check if it belongs to the boundary_type
        :param boundary_type: 'exchange', 'demand' or 'sink'
        :param external_compartment: id used for the external compartment in the model

        :type reaction: libsbml.Reaction
        :type boundary_type: str
        :type external_compartment: str

        :return: True if the reaction corresponds to the boundary_type filled, False otherwise
        :rtype: bool
        '''
        # Check if the reaction has an annotation. Annotations dominate everything.
        sbo_term = reaction.getSBOTerm()
        if sbo_term > -1:
            sbo_term = str(sbo_term)
            while len(sbo_term) < 7:
                sbo_term = '0' + sbo_term
            sbo_term = 'SBO:' + sbo_term

            if sbo_term == sbo_terms[boundary_type]:
                return True
            if sbo_term in [sbo_terms[k] for k in sbo_terms if k != boundary_type]:
                return False
        
        # Check if the reaction is in the correct compartment (exterior or inside)
        reaction_compartment = reaction.getCompartment()
        if reaction_compartment is None or reaction_compartment == '':
            reactants = reaction.getListOfReactants()
            if len(reactants) == 1:
                reactant = reactants[0]
                specie_id = reactants[0].getSpecies()
                reaction_compartment = self.getModel().getSpecies(specie_id).getCompartment()
        correct_compartment = external_compartment == reaction_compartment
        if boundary_type != "exchange":
            correct_compartment = not correct_compartment

        # Check if the reaction has the correct reversibility
        rev_type = True
        if boundary_type == "demand":
            rev_type = not reaction.getReversible()
        elif boundary_type == "sink":
            rev_type = reaction.getReversible()

        # Determine if reaction is "boundary"
        is_boundary = False
        if boundary_type == "exchange":
            if ((len(reaction.getListOfProducts()) == 0 and len(reaction.getListOfReactants()) == 1)
                or (len(reaction.getListOfProducts()) == 1 and len(reaction.getListOfReactants()) == 0)):
                is_boundary = True

        # In exclude fields ?
        to_exclude = not any(ex in reaction.getId() for ex in excludes[boundary_type])
        
        return (
            is_boundary 
            and to_exclude
            and correct_compartment
            and rev_type
        )

    def build_exchange_reaction(
        model: 'rpSBML',
        compartment_id: str,
        logger: Logger = getLogger(__name__)
    ) -> pd_DataFrame:
        '''Select exchange reactions in a model.

        :param model: a model into perform searching
        :param compartment_id: id used for the external compartment in the model
        :param logger: a logger object

        :type model: rpSBML
        :type compartment_id: str
        :type logger: Logger

        :return: a dataframe with two columns "model_id" supporting id of the specie implied in the reaction and "libsbml_reaction" the exchange reaction
        :rtype: pd.DataFrame
        '''

        def _specie_id_from_exchange(
            reaction: libsbml.Reaction,
            logger: Logger= getLogger(__name__)
        ) -> str:

            products = reaction.getListOfProducts()
            reactants = reaction.getListOfReactants()
            if sum([len(products), len(reactants)]) != 1:
                logger.warning('Unexpected reaction format')
                return None
            specie_id = None
            if len(products) > 0:
                specie_id = products[0].getSpecies()
            else:
                specie_id = reactants[0].getSpecies()
            return specie_id


        df = pd_DataFrame(columns=['model_id', 'libsbml_reaction'])

        # Create list of exchange reactions
        for reaction in model.getModel().getListOfReactions():
            if model.is_boundary_type(reaction, "exchange", compartment_id):
                compound = dict(
                    model_id=_specie_id_from_exchange(
                        reaction,
                        logger
                    ),
                    libsbml_reaction=reaction
                )
                df = df.append(compound, ignore_index=True)
        # Fmg
        df.sort_values('model_id', inplace=True)
        df.reset_index(inplace=True, drop=True)
        return df


    ######################################################################################################################
    ############################################### EC NUMBER ############################################################
    ######################################################################################################################


    def compareEC(self, meas_reac_miriam, sim_reac_miriam):
        """Compare two MIRIAM annotations and find the similarity of their EC number

        :param meas_reac_miriam: The annotation object of the source
        :param sim_reac_miriam: The annotation object of the target

        :type meas_reac_miriam: libsbml.XMLNode
        :type sim_reac_miriam: libsbml.XMLNode

        :return: The match score
        :rtype: float
        """
        # Warning we only match a single reaction at a time -- assume that there cannot be more than one to match at a given time
        if 'ec-code' in meas_reac_miriam and 'ec-code' in sim_reac_miriam:
            measured_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in meas_reac_miriam['ec-code']]
            sim_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in sim_reac_miriam['ec-code']]
            # complete the ec numbers with None to be length of 4
            for i in range(len(measured_frac_ec)):
                for y in range(len(measured_frac_ec[i]), 4):
                    measured_frac_ec[i].append(None)
            for i in range(len(sim_frac_ec)):
                for y in range(len(sim_frac_ec[i]), 4):
                    sim_frac_ec[i].append(None)
            # self.logger.debug('Measured: ')
            # self.logger.debug(measured_frac_ec)
            # self.logger.debug('Simulated: ')
            # self.logger.debug(sim_frac_ec)
            best_ec_compare = {'meas_ec': [], 'sim_ec': [], 'score': 0.0, 'found': False}
            for ec_m in measured_frac_ec:
                for ec_s in sim_frac_ec:
                    tmp_score = 0.0
                    for i in range(4):
                        if not ec_m[i]==None and not ec_s[i]==None:
                            if ec_m[i]==ec_s[i]:
                                tmp_score += 0.25
                                if i == 2:
                                    best_ec_compare['found'] = True
                            else:
                                break
                    if tmp_score>best_ec_compare['score']:
                        best_ec_compare['meas_ec'] = ec_m
                        best_ec_compare['sim_ec'] = ec_s
                        best_ec_compare['score'] = tmp_score
            return best_ec_compare['score']
        else:
            self.logger.warning(f'One of the two reactions does not have any EC entries.\nMeasured: {meas_reac_miriam} \nSimulated: {sim_reac_miriam}')
            return 0.0


    ## Put species in a dictionnary for further comparison
    #
    # @param pathway rpSBML object
    # @return dict object with species in it
    def _get_reactions_with_species_keys(
        self,
        chem_repr = ['inchikey', 'inchi', 'smiles'],
        group_id: str = 'rp_pathway'
    ) -> Dict:
        """
        Build a dictionary of reactions with 'Reactants' and 'Products' sub-dictionaries
        where species are referenced by one of the keys in the given list of chemical representations.

        Parameters
        ----------
        chem_repr: List[str]
            List of chemical representation names under which species will be referenced
        group_id: str
            Name of the group of whom the reactions are picked

        Returns
        -------
        reactions: Dict
            Dictionnary of Reactants and Products of all reactions in the current pathway
        """

        self.logger.debug(self)

        def _search_key(keys: List[str], dict: Dict) -> str:
            """
            From a given list of keys, returns the first one which is in the given dictionary.

            Parameters
            ----------
            keys: List[str]
                List of keys
            dict: Dict
                A dictionary

            Returns
            -------
            key: str
                The first key found in dict
            """
            for key in keys:
                if key in dict:
                    return key

        # Get Reactions
        reactions = {}
        for rxn_id in self.readGroupMembers(group_id):
            reaction = self.getModel().getReaction(rxn_id)
            reactions[rxn_id] = rpSBML.readBRSYNTHAnnotation(
                reaction.getAnnotation(),
                logger=self.logger
            )

        # Get Species
        species = {}
        for spe in self.getModel().getListOfSpecies():
            species[spe.getId()] = rpSBML.readBRSYNTHAnnotation(
                spe.getAnnotation(),
                logger=self.logger
            )

        # Pathways dict
        d_reactions = {}

        # Select reactions already loaded (w/o sink one then)
        for reaction_id in reactions:

            d_reactions[reaction_id] = {}

            # Fill the reactants in a dedicated dict
            d_reactants = {}
            for reactant in self.getModel().getReaction(reaction_id).getListOfReactants():# inchikey / inchi sinon miriam sinon IDs
                # Take the first key found in species
                key = _search_key(chem_repr, species[reactant.getSpecies()])
                key = species[reactant.getSpecies()][key] if key else reactant.getSpecies()
                d_reactants[key] = reactant.getStoichiometry()
            # Put all reactants dicts in reactions dict for which smiles notations are the keys
            d_reactions[reaction_id]['Reactants'] = d_reactants

            # Fill the products in a dedicated dict
            d_products = {}
            for product in self.getModel().getReaction(reaction_id).getListOfProducts():
                # Take the first key found in species
                key = _search_key(chem_repr, species[product.getSpecies()])
                key = species[product.getSpecies()][key] if key else product.getSpecies()
                d_products[key] = product.getStoichiometry()
            # Put all products dicts in reactions dict for which smiles notations are the keys
            d_reactions[reaction_id]['Products'] = d_products

        self.logger.debug(d_reactions)

        return d_reactions


    def __eq__(self, other):
        self.logger.debug(f'Comparing {self} and {other}')
        return \
            self._get_reactions_with_species_keys() \
            == \
            other._get_reactions_with_species_keys()


    # def __lt__(self, rpsbml):
    #     return self.getScore() < rpsbml.getScore()


    # def __gt__(self, rpsbml):
    #     return self.getScore() > rpsbml.getScore()


    def __str__(self):
        return f'''
            name: {str(self.getName())}
            document: {str(self.document)}
            model: {str(self.getModel())}
            '''
            #  + 'score: '     + str(self.getScore()) + '\n' \


    def getPlugin(self, plugin: str) -> object:
        plug = self.getModel().getPlugin(plugin)
        self.checklibSBML(
            plug,
            'Getting {} package'.format(plugin)
        )
        return plug


    #######################################################################
    ############################# PRIVATE FUNCTIONS #######################
    #######################################################################
    @staticmethod
    def checklibSBML(
        value,
        message: str,
        logger: Logger = getLogger(__name__)
    ) -> int:
        """Private function that checks the libSBML calls.

        Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html

        :param value: The libSBML command returned int
        :param message: The string that describes the call

        :type value: int
        :type message: str

        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :return: None
        :rtype: None
        """
        if value is None:
            # logger.error('LibSBML returned a null value trying to ' + message + '.')
            # return 1
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value != libsbml.LIBSBML_OPERATION_SUCCESS:
                err_msg = ''.join(
                    [
                        '\n',
                        '*** Error encountered trying to ' + message + '.\n',
                        '     \__ LibSBML returned error code ' + str(value) + ': "',
                        libsbml.OperationReturnValue_toString(value).strip() + '"'
                    ]
                )
                # logger.error(err_msg)
                # return 2
                raise SystemExit(err_msg)

    def _nameToSbmlId(self, name):
        """String to SBML id's

        Convert any String to one that is compatible with the SBML meta_id formatting requirements

        :param name: The input string

        :type name: str

        :return: SBML valid string
        :rtype: str
        """
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    def _genMetaID(self, name):
        """String to hashed id

        Hash an input string and then pass it to _nameToSbmlId()

        :param name: Input string

        :type name: str

        :return: Hashed string id
        :rtype: str
        """
        self.logger.debug(f'Generating metaID for {name}')
        return self._nameToSbmlId(
            sha256(
                str(name).encode('utf-8')
            ).hexdigest()
        )


    def _compareXref(self, current, toadd):
        """Compare two dictionaries of lists that describe the cross-reference and return the difference

        :param current: The source cross-reference dictionary
        :param toadd: The target cross-reference list

        :type current: Dict
        :type toadd: List

        :return: Difference between the two cross-reference dictionaries
        :rtype: Dict
        """
        _toadd = deepcopy(toadd)
        for database_id in current:
            try:
                list_diff = [i for i in _toadd[database_id] if i not in current[database_id]]
                if not list_diff:
                    _toadd.pop(database_id)
                else:
                    _toadd[database_id] = list_diff
            except KeyError:
                pass
        return _toadd


    ######################################################################
    ####################### Annotations ##################################
    ######################################################################


    def _defaultBothAnnot(self, meta_id):
        """Returns a default annotation string that include MIRIAM and BRSynth annotation

        :param meta_id: The meta ID to be added to the default annotation

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''
            <annotation>
                <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
                    <bqbiol:is>
                        <rdf:Bag>
                        </rdf:Bag>
                    </bqbiol:is>
                    </rdf:Description>
                    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
                    <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
                    </brsynth:brsynth>
                    </rdf:BRSynth>
                </rdf:RDF>
            </annotation>
        '''


    def _defaultBRSynthAnnot(self, meta_id):
        """Returns BRSynth default annotation string

        :param meta_id: The meta ID to be added to the annotation string

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''
            <annotation>
                <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
                    <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
                    </brsynth:brsynth>
                    </rdf:BRSynth>
                </rdf:RDF>
            </annotation>
        '''


    def _defaultMIRIAMAnnot(self, meta_id):
        """Returns MIRIAM default annotation string

        :param meta_id: The meta ID to be added to the annotation string

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''
            <annotation>
                <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
                    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
                    <bqbiol:is>
                        <rdf:Bag>
                        </rdf:Bag>
                    </bqbiol:is>
                    </rdf:Description>
                </rdf:RDF>
            </annotation>
        '''

    def updateBRSynth(
        self,
        sbase_obj,
        annot_header,
        value,
        meta_id=None
    ):
        """Append or update an entry to the BRSynth annotation of the passed libsbml.SBase object.

        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: The libSBML object to add the different
        :param annot_header: The annotation header that defines the type of entry
        :param value: The value(s) to add
        :param unit: Add a values unit to the entry
        :param isAlone: Add the entry without any unit or defined within a value child (Setting this to True will ignore any unit)
        :param isList: Define if the value entry is a list or not
        :param isSort: Sort the list that is passed (Only if the isList is True)
        :param meta_id: The meta ID to be added to the annotation string

        :type sbase_obj: libsbml.SBase
        :type annot_header: str
        :type value: Union[str, int, float, list]
        :type unit: str
        :type isAlone: bool
        :type isList: bool
        :type isSort: bool
        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """

        def __write_list(value, annot_header: str) -> str:
            annotation = '>'
            for v in value:
                annotation += f'''
                                <brsynth:{v}'''
                annotation += f'''/>'''
            annotation += f'''
                            </brsynth:{str(annot_header)}>'''
            return annotation

        def __write_dict(value, annot_header: str) -> str:
            annotation = '>'
            for k, v in value.items():
                annotation += f'''
                                <brsynth:{k}'''
                if isinstance(v, dict):
                    for _k, _v in v.items():
                        annotation += f''' {_k}="{_v}"'''
                else:
                    annotation += f''' value="{v}"'''
                annotation += f'''/>'''
            annotation += f'''
                            </brsynth:{str(annot_header)}>'''
            return annotation

        def __write_header(annot_header: str) -> str:
            return f'''
                <annotation>
                    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
                        <rdf:BRSynth rdf:about="# adding">
                            <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
                                <brsynth:{annot_header}'''

        def __write_footer() -> str:
            return '''
                            </brsynth:brsynth>
                        </rdf:BRSynth>
                    </rdf:RDF>
                </annotation>'''

        annotation = __write_header(str(annot_header))
        if isinstance(value, list):
            annotation += __write_list(value, annot_header)
        elif isinstance(value, dict):
            annotation += __write_dict(value, annot_header)
        else:
            annotation += f' value="{str(value)}"/>'
        annotation += __write_footer()

        # self.logger.debug('annotation: {0}'.format(annotation))
        annot_obj = libsbml.XMLNode.convertStringToXMLNode(annotation)
        if not annot_obj:
            self.logger.error(f'Cannot convert this string to annotation object: {annotation}')
            return False
        #### retreive the annotation object
        brsynth_annot = None
        obj_annot = sbase_obj.getAnnotation()
        if not obj_annot:
            sbase_obj.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultBRSynthAnnot(meta_id)))
            obj_annot = sbase_obj.getAnnotation()
            if not obj_annot:
                self.logger.error('Cannot update BRSynth annotation')
                return False
        brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        if not brsynth_annot:
             self.logger.error('Cannot find the BRSynth annotation')
             return False

        # try to update the annotation
        target_found = self.updateBRSynthAnnot(brsynth_annot, annot_obj, annot_header)

        # if target not found, then add the annotation as new
        if not target_found:
            source_found = self.addBRSynthAnnot(brsynth_annot, annot_obj, annot_header)
            if not source_found:
                self.logger.error(f'Cannot find {annot_header} in source annotation')
                return False

        return True


    def updateBRSynthAnnot(self, annot, annot_obj, annot_header):

        target_found = False

        for i in range(annot.getNumChildren()):

            # self.logger.debug(annot_header+' -- '+str(annot.getChild(i).getName()))

            if annot_header == annot.getChild(i).getName():

                target_found = True
                '''
                self.checklibSBML(annot.removeChild(annot.getIndex(i)),
                    'Removing annotation '+str(annot_header))
                '''
                self.checklibSBML(annot.removeChild(i), 'Removing annotation '+str(annot_header), self.logger)
                source_found = False
                source_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')

                for y in range(source_annot.getNumChildren()):
                    # self.logger.debug('\t'+annot_header+' -- '+str(source_annot.getChild(y).getName()))
                    if str(annot_header)==str(source_annot.getChild(y).getName()):
                        source_found = True
                        # self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_annot.getChild(y).toXMLString()))
                        towrite_annot = source_annot.getChild(y)
                        self.checklibSBML(annot.addChild(towrite_annot), ' 1 - Adding annotation to the brsynth annotation', self.logger)
                        break

                if not source_found:
                    self.logger.error(f'Cannot find {annot_header} in source annotation')

        return target_found


    def addBRSynthAnnot(self, brsynth_annot, annot_obj, annot_header):

        # self.logger.debug('Cannot find '+str(annot_header)+' in target annotation')

        source_found = False
        source_brsynth_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')

        for y in range(source_brsynth_annot.getNumChildren()):
            # self.logger.debug('\t'+annot_header+' -- '+str(source_brsynth_annot.getChild(y).getName()))
            if str(annot_header) == str(source_brsynth_annot.getChild(y).getName()):
                source_found = True
                # self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_brsynth_annot.getChild(y).toXMLString()))
                towrite_annot = source_brsynth_annot.getChild(y)
                self.checklibSBML(brsynth_annot.addChild(towrite_annot), '2 - Adding annotation to the brsynth annotation', self.logger)
                break

        return source_found


    def addUpdateMIRIAM(self, sbase_obj, type_param, xref, meta_id=None):
        """Append or update an entry to the MIRIAM annotation of the passed libsbml.SBase object.

        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: The libSBML object to add the different
        :param type_param: The type of parameter entered. Valid include ['compartment', 'reaction', 'species']
        :param xref: List of the cross reference
        :param meta_id: The meta ID to be added to the annotation string

        :type sbase_obj: libsbml.SBase
        :type type_param: str
        :type xref: List
        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """
        if type_param not in ['compartment', 'reaction', 'species']:
            self.logger.error(f'type_param must be [\'compartment\', \'reaction\', \'species\'] not {type_param}')
            return False
        miriam_annot = None
        isReplace = False
        try:
            miriam_annot = sbase_obj.getAnnotation().getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            miriam_elements = self.readMIRIAMAnnotation(sbase_obj.getAnnotation())
            if not miriam_elements:
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot_1 = libsbml.XMLNode.convertStringToXMLNode(self._defaultBothAnnot(meta_id))
                miriam_annot = miriam_annot_1.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            else:
                miriam_elements = None
        except AttributeError:
            try:
                # Cannot find MIRIAM annotation, create it
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                miriam_annot = miriam_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            except AttributeError:
                self.logger.error('Fatal error fetching the annotation')
                return False
        # compile the list of current species
        inside = {}
        for i in range(miriam_annot.getNumChildren()):
            single_miriam = miriam_annot.getChild(i)
            if single_miriam.getAttributes().getLength()>1:
                self.logger.error(f'MIRIAM annotations should never have more than 1: {single_miriam.toXMLString()}')
                continue
            single_miriam_attr = single_miriam.getAttributes()
            if not single_miriam_attr.isEmpty():
                try:
                    db = single_miriam_attr.getValue(0).split('/')[-2]
                    v = single_miriam_attr.getValue(0).split('/')[-1]
                    inside[self.header_miriam[type_param][db]].append(v)
                except KeyError:
                    try:
                        db = single_miriam_attr.getValue(0).split('/')[-2]
                        v = single_miriam_attr.getValue(0).split('/')[-1]
                        inside[self.header_miriam[type_param][db]] = [v]
                    except KeyError:
                        self.logger.warning(f'Cannot find the self.header_miriram entry {db}')
                        continue
            else:
                self.logger.warning('Cannot return MIRIAM attribute')
                pass
        # add or ignore
        toadd = self._compareXref(inside, xref)
        toadd = self.convert_miriam_to_dict(toadd)
        for database_id in toadd:
            for species_id in toadd[database_id]:
                # not sure how to avoid having it that way
                if database_id in self.miriam_header[type_param]:
                    try:
                        # determine if the dictionnaries
                        annotation = '''<annotation>
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
    <rdf:Description rdf:about="# tmp">
      <bqbiol:is>
        <rdf:Bag>'''
                        if type_param=='species':
                            if database_id=='kegg' and species_id[0]=='C':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_c']+str(species_id)+'''"/>'''
                            elif database_id=='kegg' and species_id[0]=='D':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_d']+str(species_id)+'''"/>'''
                            else:
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        else:
                            annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    </rdf:RDF>
    </annotation>'''
                        toPass_annot = libsbml.XMLNode.convertStringToXMLNode(annotation)
                        toWrite_annot = toPass_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag').getChild(0)
                        miriam_annot.insertChild(0, toWrite_annot)
                    except KeyError:
                        # WARNING need to check this
                        self.logger.warning(f'Cannot find {database_id} in self.miriam_header for {type_param}')
                        continue
        if isReplace:
            ori_miriam_annot = sbase_obj.getAnnotation()
            if not ori_miriam_annot:
                sbase_obj.unsetAnnotation()
                sbase_obj.setAnnotation(miriam_annot)
            else:
                rpSBML.checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').removeChild(0), 'Removing annotation "is"', self.logger)
                rpSBML.checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').addChild(miriam_annot), 'Adding annotation to the brsynth annotation', self.logger)
        return True


    def find_or_create_objective(
        self,
        rxn_id: str,
        obj_id: str,
        coeff: float = 1.0,
        is_max: bool = True,
    ) -> str:

        self.logger.debug(f'rxn_id: {rxn_id}')
        self.logger.debug(f'coeff: {coeff}')
        self.logger.debug(f'is_max: {is_max}')

        _obj_id = self.search_objective_from_rxnid(rxn_id)

        # If cannot find a valid objective create it
        if _obj_id is None:
            self.create_multiflux_objective(
                fluxobj_id = obj_id,
                reactionNames = [rxn_id],
                coefficients = [coeff],
                is_max = is_max
            )
        else:
            obj_id = _obj_id

        return obj_id


    def getListOfObjectives(self) -> List:
        return list(self.getPlugin('fbc').getListOfObjectives())

    def search_objective_from_rxnid(
        self,
        rxn_id: str
    ) -> str:
        for objective in self.getPlugin('fbc').getListOfObjectives():
            for flux_obj in objective.getListOfFluxObjectives():
                if flux_obj.getReaction() == rxn_id:
                    return objective.getId()
        return None

    def search_objective(
        self,
        objective_id: str,
        reactions: List[str] = []
    ) -> str:
        """Find the objective (with only one reaction associated) based on the reaction ID and if not found create it

        :param reactions: List of the reactions id's to set as objectives
        :param coefficients: List of the coefficients about the objectives
        :param isMax: Maximise or minimise the objective
        :param objective_id: overwite the default id if created (from obj_[reactions])

        :type reactions: list
        :type coefficients: list
        :type isMax: bool
        :type objective_id: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: str
        :return: Objective ID
        """
        fbc_plugin = self.getPlugin('fbc')

        for objective in fbc_plugin.getListOfObjectives():

            if objective.getId() == objective_id:
                self.logger.debug(f'The specified objective id {objective_id} already exists')
                return objective_id
            
            if not set([i.getReaction() for i in objective.getListOfFluxObjectives()])-set(reactions):
                # TODO: consider setting changing the name of the objective
                self.logger.debug(f'The specified objective id {objective_id} has another objective with the same reactions: {objective.getId()}')
                return objective.getId()

        return None


    def create_multiflux_objective(
        self,
        fluxobj_id: str,
        reactionNames: List[str],
        coefficients: List[float],
        is_max: bool = True,
        meta_id: str = None
    ) -> None:
        """Create libSBML flux objective

        Using the FBC package one can add the FBA flux objective directly to the model. Can add multiple reactions. This function sets a particular reaction as objective with maximization or minimization objectives

        :param fluxobj_id: The id of the flux objective
        :param reactionNames: The list of string id's of the reaction that is associated with the reaction
        :param coefficients: The list of int defining the coefficients of the flux objective
        :param isMax: Define if the objective is coefficient (Default: True)
        :param meta_id: Meta id (Default: None)

        :type fluxobj_id: str
        :type reactionNames: list
        :type coefficients: list
        :type isMax: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """

        if len(reactionNames) != len(coefficients):
            self.logger.error('The size of reactionNames is not the same as coefficients')
            exit()

        fbc_plugin = self.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        # target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)

        if is_max:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')

        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple

        for reac, coef in zip(reactionNames, coefficients):
            target_flux_obj = target_obj.createFluxObjective()
            target_flux_obj.setReaction(reac)
            target_flux_obj.setCoefficient(coef)
            # if not meta_id:
            #     meta_id = self._genMetaID(str(fluxobj_id))
            # target_flux_obj.setMetaId(meta_id)
            # target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def check_SBML_compartment(
        self, compartment_id: str, logger: Logger = getLogger(__name__)
    ) -> Tuple[str, str]:
        # Check model compartment ID
        # Set new compartment ID in case it exists under another ID
        compartment = self.has_compartment(compartment_id)
        if compartment is not None:
            logger.debug(
                f"Compartment '{compartment_id}' found in the model '{self.getName()}' as {compartment.getId()}"
            )
            if compartment_id != compartment.getId():
                logger.warning(
                    f"Compartment '{compartment_id}' has been replaced by {compartment.getId()}"
                )
            return compartment.getId()
        else:
            logger.error(
                f"Compartment '{compartment_id}' not found in the model '{self.getName()}'"
            )
            logger.error(
                "Available compartments: {comp_ids}".format(
                    comp_ids=", ".join(
                        [
                            comp.getId()
                            for comp in list(self.getModel().getListOfCompartments())
                        ]
                    )
                )
            )
            return None


    def check_SBML_rxnid(
        self, rxn_id: str, logger: Logger = getLogger(__name__)
    ) -> Tuple[str, str]:
        logger.debug(f"rpsbml: {self.getName()}")
        logger.debug(f"rxn_id: {rxn_id}")

        # Check model reaction ID
        # Get reaction from the rpSBML
        _rxn = self.getModel().getReaction(rxn_id)
        if _rxn is None:
            logger.error(
                f"Reaction ID '{rxn_id}' not found in the model '{self.getName()}'"
            )
            possible_rxn_ids = [
                rxn.getId()
                for rxn in list(self.getModel().getListOfReactions())
                if rxn_id in rxn.getId().lower()
            ]
            if possible_rxn_ids != []:
                logger.error(
                    "Possible reactions: {rxn_ids}".format(
                        rxn_ids=", ".join(possible_rxn_ids)
                    )
                )
            return None

        logger.debug(f"rxn_id found as {rxn_id}")

        return _rxn.getId()


    #####################################################################
    ########################## INPUT/OUTPUT #############################
    #####################################################################
    def readSBML(self, inFile):
        """Open an SBML file to the object

        :param inFile: Path to the input SBML file

        :type inFile: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: None
        :return: Dictionnary of the pathway annotation
        """

        self.logger.debug(f'Read SBML file from {inFile}')

        self.document = libsbml.readSBMLFromFile(inFile)
        rpSBML.checklibSBML(self.getDocument(), 'reading input file')
        errors = self.getDocument().getNumErrors()
        # display the errors in the log accordning to the severity
        for err in [self.getDocument().getError(i) for i in range(self.getDocument().getNumErrors())]:
            # TODO if the error is related to packages not enabled (like groups or fbc) activate them
            if err.isFatal:
                self.logger.error(f'libSBML reading error: {err.getShortMessage()}')
                raise FileNotFoundError
            else:
                self.logger.warning(f'libSBML reading warning: {err.getShortMessage()}')
        if not self.getModel():
            self.logger.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        # enabling the extra packages if they do not exists when reading a model
        if not self.getModel().isPackageEnabled('groups'):
            rpSBML.checklibSBML(self.getModel().enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
            rpSBML.checklibSBML(self.getDocument().setPackageRequired('groups', False), 'enabling groups package')
        if not self.getModel().isPackageEnabled('fbc'):
            rpSBML.checklibSBML(self.getModel().enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
            rpSBML.checklibSBML(self.getDocument().setPackageRequired('fbc', False), 'enabling FBC package')


    ## Export a libSBML model to file
    #
    # Export the libSBML model to an SBML file
    #
    # @param model libSBML model to be saved to file
    # @param model_id model id, note that the name of the file will be that
    # @param path Non required parameter that will define the path where the model will be saved
    def build_filename_from_name(self) -> str:
        ext = ''
        if not str(self.getName()).endswith('_sbml'):
            ext = '_sbml'
        return str(self.getName())+ext+'.xml'


    def write_to_file(
        self,
        filename: str = None
        # ,
        # outdir: str = None
    ) -> str:
        """
        Write pathway to disk (optionally into a specific folder).

        Parameters
        ----------
        filename: str
            Filename to store the file under.
        outdir: str
            Folder to store the file into.

        Returns
        -------
        filename: str
            Full path of the stored file
        """

        ext = ''

        if not str(self.getName()).endswith('_sbml'):
            ext = '_sbml'

        if filename is not None:
            out_filename = filename
        else:
            out_filename = self.build_filename_from_name()

        libsbml.writeSBMLToFile(
            self.getDocument(),
            out_filename
        )

        return out_filename


    def readReactionSpecies(self, reaction):
        """Return the products and the species associated with a reaction

        :param reaction: Reaction object of libSBML

        :type annot: libsbml.Reaction

        :rtype: dict
        :return: Dictionary of the reaction stoichiometry
        """

        # TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model

        self.logger.debug(reaction)

        toRet = {'left': {}, 'right': {}}

        # reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            toRet['left'][reactant_ref.getSpecies()] = int(reactant_ref.getStoichiometry())

        # products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            toRet['right'][product_ref.getSpecies()] = int(product_ref.getStoichiometry())

        return toRet


    #####################################################################
    ########################## READ #####################################
    #####################################################################
    # TODO: add error handling if the groups does not exist
    def readGroupMembers(self, group_id):
        """Return the members of a groups entry

        :param group_id: The pathway ID (Default: rp_pathway)

        :type group_id: str

        :rtype: list
        :return: List of member id's of a particular group
        """
        group = self.getGroup(group_id)
        try:
            rpSBML.checklibSBML(group, 'retreiving '+group_id+' group')
            return [m.getIdRef() for m in group.getListOfMembers()]
        except:
            self.logger.debug(f'Group \'{group_id}\' not found')
            return None
        # members = []
        # for member in group.getListOfMembers():
        #     members.append(member.getIdRef())
        # return members


    # TODO: merge with unique species
    # TODO: change the name of the function to read
    def readRPspecies(self, pathway_id='rp_pathway'):
        """Return the species stoichiometry of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: dict
        :return: Dictionary of the pathway species and reactions
        """
        reacMembers = {}
        for reacId in self.readGroupMembers(pathway_id):
            reacMembers[reacId] = {}
            reacMembers[reacId]['products'] = {}
            reacMembers[reacId]['reactants'] = {}
            reac = self.getModel().getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers[reacId]['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers[reacId]['reactants'][rea.getSpecies()] = rea.getStoichiometry()
        return reacMembers


    def readUniqueRPspecies(self):
        """Return the unique species of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: list
        :return: List of unique species
        """
        rpSpecies = self.readRPspecies()
        toRet = []
        for i in rpSpecies:
            for y in rpSpecies[i]:
                for z in rpSpecies[i][y]:
                    if z not in toRet:
                        toRet.append(z)
        return toRet
        # reacMembers = self.readRPspecies(pathway_id)
        # return set(set(ori_rp_path['products'].keys())|set(ori_rp_path['reactants'].keys()))


    @staticmethod
    def readMIRIAMAnnotation(
        annot: libsbml.XMLNode,
        logger: Logger = getLogger(__name__)
    ):
        """Return the MIRIAM annotations of a libSBML object

        :param annot: The annotation object of libSBML

        :type annot: libsbml.XMLNode

        :rtype: list
        :return: List of all the annotations
        """
        try:
            toRet = []
            # Get the main IDs
            bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                toRet.append(str_annot)
            # Get the secondary IDs
            bag = annot.getChild('RDF').getChild('Description').getChild('isRelatedTo').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                toRet.append(str_annot)
            return toRet
        except:
            return []

    @staticmethod
    def readBRSYNTHAnnotation(
        annot: libsbml.XMLNode,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        """Return a dictionnary of all the information in a BRSynth annotations

        :param annot: The annotation object of libSBML

        :type annot: libsbml.XMLNode

        :rtype: dict
        :return: Dictionary of all the BRSynth annotations
        """

        if not annot:
            # logger.warning('The passed annotation is None')
            return {}

        def _readBRSYNTHAnnotationToDict(
            annot: libsbml.XMLNode,
            logger: Logger = getLogger(__name__)
        ) -> Dict:
            toRet = {}
            for i_child in range(annot.getNumChildren()):
                child = annot.getChild(i_child)
                toRet[child.getName()] = {
                    child.getAttributes().getName(i_attr): child.getAttributes().getValue(i_attr)
                    for i_attr in range(child.getAttributes().getNumAttributes())
                }
                toRet[child.getName()] = eval_value(toRet[child.getName()]['value'])
            return toRet

        def _readBRSYNTHAnnotationToList(
            annot: libsbml.XMLNode,
            logger: Logger = getLogger(__name__)
        ) -> List:
            toRet = []
            for i_child in range(annot.getNumChildren()):
                child = annot.getChild(i_child)
                toRet.append(child.getName())
            return toRet

        def _readBRSYNTHAnnotationToValue(
            annot: libsbml.XMLNode,
            logger: Logger = getLogger(__name__)
        ) -> TypeVar:
            return eval_value(annot.getAttrValue('value'))

        def eval_value(value: str) -> TypeVar:
            try:
                return eval(value)
            except:
                return str(value)

        toRet = {}

        bag = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')

        for i in range(bag.getNumChildren()):

            ann = bag.getChild(i)

            if ann == '':
                logger.warning('This contains no attributes: '+str(ann.toXMLString()))
                continue

            # Read list
            if (
                ann.getName().startswith('tmpl_rxn_ids')
                or ann.getName().startswith('rule_ids')
            ):
                toRet[ann.getName()] = _readBRSYNTHAnnotationToList(
                    annot=ann,
                    logger=logger
                )

            # Read dict
            elif (
                ann.getName().startswith('thermo_')
                or ann.getName().startswith('fba_')
                or ann.getName().startswith('selenzy_')
            ):
                toRet[ann.getName()] = _readBRSYNTHAnnotationToDict(
                    annot=ann,
                    logger=logger
                )

            else:
                toRet[ann.getName()] = _readBRSYNTHAnnotationToValue(
                    annot=ann,
                    logger=logger
                )
                if ann.getName() == 'smiles':
                    toRet[ann.getName()] = toRet[ann.getName()].replace('&gt;', '>')

        return {k: v for k, v in toRet.items()}


    #####################################################################
    ######################### INQUIRE ###################################
    #####################################################################
    def speciesExists(self, species, compartment_id='MNXC3'):
        """Determine if the model already contains a species according to its ID

        :param reaction: Reaction object of libSBML

        :type annot: libsbml.Reaction

        :rtype: bool
        :return: True if exists and False if not
        """
        # if speciesName in [i.getName() for i in self.getModel().getListOfSpecies()] or speciesName+'__64__'+compartment_id in [i.getId() for i in self.getModel().getListOfSpecies()]:
        if species in [i.getId() for i in self.getModel().getListOfSpecies()]:
            return True
        return False


    def isSpeciesProduct(self, species_id, ignoreReactions=[]):
        """Function to determine if a species can be a product of any reaction.

        :param species_id: ID of the species to find
        :param ignoreReactions: List of all the reaction id's to ignore

        :type species_id: str
        :type ignoreReactions: list

        :rtype: bool
        :return: True if its a product of a reaction False if not
        """
        # return all the parameters values
        param_dict = {i.getId(): i.getValue() for i in self.getModel().parameters}
        for reaction in self.getModel().getListOfReactions():
            if reaction.getId() not in ignoreReactions:
                # check that the function is reversible by reversibility and FBC bounds
                if reaction.reversible:
                    reaction_fbc = reaction.getPlugin('fbc')
                    # strict left to right
                    if param_dict[reaction_fbc.getLowerFluxBound()]>=0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                    # can go both ways
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                        elif species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    # strict right to left
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]<=0 and param_dict[reaction_fbc.getLowerFluxBound()]<param_dict[reaction_fbc.getUpperFluxBound()]:
                        if species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    else:
                        self.logger.warning(f'isSpeciesProduct does not find the directionailty of the reaction for reaction: {species_id}')
                        return True
                else:
                    # if the reaction is not reversible then product are the only way to create it
                    if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                        return True
        return False


    def read_reactions(
        self,
        pathway_id: str = None,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        """
        Read reactions field in rpSBML file for rp_pathway.

        Parameters
        ----------
        rp_pathway: libsbml.Group
            Pathway to extract infos from

        Returns
        -------
        pathway: Dict
            Read fields
        """
        reactions = {}
        if pathway_id is None:
            rxn_l = [rxn.getId() for rxn in list(self.getModel().getListOfReactions())]
        else:
            rxn_l = self.readGroupMembers(pathway_id)

        for rxn_id in rxn_l:
            reactions[rxn_id] = self.read_reaction(rxn_id)
        return reactions


    def read_reaction(
        self,
        rxn_id: str,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        """
        Read reaction in rpSBML file.

        Parameters
        ----------
        rxn_id: str
            Reaction to read infos from

        Returns
        -------
        reaction: Dict
            Read fields
        """
        reaction = self.getModel().getReaction(rxn_id)
        if reaction is None:
            return None
        annot = reaction.getAnnotation()
        species = self.readReactionSpecies(reaction)
        fbc = reaction.getPlugin('fbc')
        return {
            # BRSynth annotations
            'brsynth': self.readBRSYNTHAnnotation(annot, self.logger),
            # species
            'left': species['left'],
            'right': species['right'],
            # MIRIAM annotations
            'miriam': self.readMIRIAMAnnotation(annot),
            # FBC values
            'fbc_lower_value': fbc.getLowerFluxBound(),
            'fbc_upper_value': fbc.getUpperFluxBound(),
            'reversible': reaction.getReversible()
        }


    def read_species(
        self,
        logger: Logger = getLogger(__name__)
    ) -> Dict:
        """
        Read species field in rpSBML file for pathway_id.

        Returns
        -------
        pathway: Dict
            Read fields
        """
        species_dict = {}
        for species in self.getModel().getListOfSpecies():
            annot = species.getAnnotation()
            species_dict[species.getId()] = {}
            species_dict[species.getId()]['object'] = species
            species_dict[species.getId()]['brsynth'] = self.readBRSYNTHAnnotation(annot, self.logger)
            species_dict[species.getId()]['miriam']  = self.readMIRIAMAnnotation(annot)

        return species_dict


    #########################################################################
    ################### CONVERT BETWEEEN FORMATS ############################
    #########################################################################
    @staticmethod
    def from_json(
        filename: str,
        logger: Logger = getLogger(__name__)
    ) -> 'rpSBML':

        with open(filename, 'r') as fp:
            pathway = json_load(fp)

        return rpSBML.from_dict(
            pathway=pathway,
            logger=logger
        )

    def create_enriched_group(
        self,
        group_id: str,
        members: List = [],
        infos: Dict = {},
        logger: Logger = getLogger(__name__)
    ) -> None:

        group = self.createGroup(group_id)
        for member_id in members:
            member = libsbml.Member()
            member.setIdRef(member_id)
            group.addMember(member)

        for key, value in infos.items():
            self.updateBRSynth(
                sbase_obj=self.getGroup('rp_pathway'),
                annot_header=key,
                value=value
            )

    def to_cobra(
        self,
        logger: Logger = getLogger(__name__)
    ) -> cobra.Model:
        """Convert rpSBML to a cobra Model

        :param logger: a logger object

        :type logger: Logger

        :return : A cobra Model
        :rtype: cobra.Model
        """
        # To handle file removing (Windows)
        cobra_model = None
        with NamedTemporaryFile(delete=False) as temp_f:
            self.write_to_file(temp_f.name)
            temp_f.close()
            try:
                cobra_model = cobra_io.read_sbml_model(temp_f.name, use_fbc_package=True)
            except CobraSBMLError:
                logger.error('Something went wrong reading the SBML model')
                (model, errors) = validate_sbml_model(temp_f.name)
                logger.error(str(json_dumps(errors, indent=4)))

        # To handle file removing (Windows)
        remove(temp_f.name)

        return cobra_model

    @staticmethod
    def from_cobra(
        model: cobra.Model,
        logger: Logger = getLogger(__name__)
    ) -> 'rpSBML':
        """Convert a cobra Model to an rpSBML object
        BE CAREFUL: some data will be lost during conversion

        :param model: a model to convert
        :param logger: a logger object

        :type model: cobra.Model
        :type logger: Logger

        :return : An rpSBML object
        :rtype: rpSBML
        """
        # To handle file removing (Windows)
        cobra_model = None
        with NamedTemporaryFile(delete=False) as temp_f:
            cobra.io.write_sbml_model(
                model,
                temp_f.name
            )
            temp_f.close()

        rpsbml = rpSBML(
            inFile=temp_f.name,
            logger=logger
        )
        # To handle file removing (Windows)
        remove(temp_f.name)

        return rpsbml


    @staticmethod
    def compareMIRIAMAnnotations(
        source_annot,
        target_annot,
        logger: Logger = getLogger(__name__)
    ) -> bool:
        """Determine if two libsbml objects have members in common in MIRIAM annotation

        Compare two dictionnaries and if any of the values of any of the same keys are the same then the function return True, and if none are found then return False

        :param source_annot: Source object of libSBML
        :param target_annot: Target object of libSBML

        :type source_annot: libsbml.Object
        :type target_annot: libsbml.Object

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = rpSBML.convert_miriam_to_dict(
            rpSBML.readMIRIAMAnnotation(source_annot)
        )
        target_dict = rpSBML.convert_miriam_to_dict(
            rpSBML.readMIRIAMAnnotation(target_annot)
        )
        common_keys = source_dict.keys() & target_dict.keys()
        for key in common_keys:
            common_values = set(source_dict[key]) & set(target_dict[key])
            if common_values:
                return True
        return False


    @staticmethod
    def convert_miriam_to_dict(miriam: List) -> Dict:
        """Convert a MIRIAM annotation to a dictionary

        :param miriam: MIRIAM annotation
        :type miriam: list

        :rtype: dict
        :return: Dictionary of the annotation
        """
        res_dict = {}
        for xref in miriam:
            try:
                m = re.match('^http(s)?://identifiers.org/(\w+\.?\w+)[:/](\w+:?\w*)', xref).groups()
                db = m[1]
                ref = m[2]
                if db in res_dict:
                    res_dict[db].append(ref)
                else:
                    res_dict[db] = [ref]
            except AttributeError:
                pass
        return res_dict


    #########################################################################
    ############################# MODEL APPEND ##############################
    #########################################################################
    def getReactionConstraints(
        self,
        rxn_id: str
    ) -> Tuple[float, float]:
        """
        Returns flux bounds.

        Parameters
        ----------
        rxn_id: str
            Reaction ID

        Returns
        -------
        bounds: Tuple[float, float]
            Tuple of lower and upper flux bounds
        """
        reac_fbc = self.getModel().getReaction(rxn_id).getPlugin('fbc')
        old_lower_value = self.getModel().getParameter(reac_fbc.getLowerFluxBound()).value
        old_upper_value = self.getModel().getParameter(reac_fbc.getUpperFluxBound()).value
        return old_lower_value, old_upper_value

    def setReactionConstraints(
        self,
        reaction_id: str,
        upper_bound: float,
        lower_bound: float,
        unit: str = 'mmol_per_gDW_per_hr',
        is_constant: bool = True
    ) -> None:
        """Set a given reaction's upper and lower bounds

        Sets the upper and lower bounds of a reaction. Note that if the numerical values passed
        are not recognised, new parameters are created for each of them

        :param reaction_id: The id of the reaction
        :param upper_bound: Reaction upper bound
        :param lower_bound: Reaction lower bound
        :param unit: Unit to the bounds (Default: mmol_per_gDW_per_hr)
        :param is_constant: Set if the reaction is constant (Default: True)

        :type reaction_id: str
        :type upper_bound: float
        :type lower_bound: float
        :type unit: str
        :type is_constant: bool

        :rtype: tuple or bool
        :return: bool if there is an error and tuple of the lower and upper bound
        """

        reaction = self.getModel().getReaction(reaction_id)

        if not reaction:
            self.logger.error(f'Cannot find the reaction: {reaction_id}')
            return False
 
        reac_fbc = reaction.getPlugin('fbc')
        rpSBML.checklibSBML(reac_fbc, 'extending reaction for FBC')
 
        ########## upper bound #############
        upper_param = self.createReturnFluxParameter(upper_bound, unit, is_constant)
        rpSBML.checklibSBML(reac_fbc.setUpperFluxBound(upper_param.getId()),
            'setting '+str(reaction_id)+' upper flux bound')
 
        ######### lower bound #############
        lower_param = self.createReturnFluxParameter(lower_bound, unit, is_constant)
        rpSBML.checklibSBML(reac_fbc.setLowerFluxBound(lower_param.getId()),
            'setting '+str(reaction_id)+' lower flux bound')


    #########################################################################
    ############################# MODEL CREATION FUNCTIONS ##################
    #########################################################################
    def createModel(self, name='dummy', model_id='dummy', meta_id=None):
        """Create libSBML model instance

        Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment

        :param name: The name of the of the model
        :param model_id: The id of the model
        :param meta_id: Meta ID of the model (Default: None)

        :type name: str
        :type model_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """

        self.logger.debug(f'Name: {name}')
        self.logger.debug(f'Model ID: {model_id}')
        self.logger.debug(f'Meta ID: {meta_id}')

        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        rpSBML.checklibSBML(self.sbmlns, 'generating model namespace')
        rpSBML.checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        rpSBML.checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        # sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        rpSBML.checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        rpSBML.checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        #!!!! must be set to false for no apparent reason
        rpSBML.checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.logger.debug('Creating SBML model...')
        self.document.createModel()
        self.logger.debug('SBML model created')
        rpSBML.checklibSBML(self.getModel(), 'generating the model')
        self.logger.debug(f'Setting the model ID to {model_id}...')
        rpSBML.checklibSBML(self.getModel().setId(model_id), 'setting the model ID')
        self.logger.debug(f'Setting the model ID to {model_id}... Done')
        model_fbc = self.getModel().getPlugin('fbc')
        model_fbc.setStrict(True)
        if not meta_id:
            meta_id = self._genMetaID(model_id)
        rpSBML.checklibSBML(self.getModel().setMetaId(meta_id), 'setting model meta_id')
        rpSBML.checklibSBML(self.getModel().setName(name), 'setting model name')
        rpSBML.checklibSBML(self.getModel().setTimeUnits('second'), 'setting model time unit')
        rpSBML.checklibSBML(self.getModel().setExtentUnits('mole'), 'setting model compartment unit')
        rpSBML.checklibSBML(self.getModel().setSubstanceUnits('mole'), 'setting model substance unit')


    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to
    #TODO: consider seperating it in another function if another compartment is to be created
    #TODO: use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, meta_id=None):
        """Create libSBML compartment

        :param size: Size of the compartment
        :param compId: Compartment id
        :param compName: Compartment Name
        :param compXref: Cross reference dictionary of the compartment
        :param meta_id: Meta id (Default: None)

        :type size: float
        :type compId: str
        :type compName: str
        :type compXref: dict
        :type meta_id: str

        :rtype: None
        :return: None
        """

        self.logger.debug(f'Size: {size}')
        self.logger.debug(f'Compartment ID: {compId}')
        self.logger.debug(f'Compartment Name: {compName}')
        self.logger.debug(f'Compartment Xref: {compXref}')
        self.logger.debug(f'Meta ID: {meta_id}')

        comp = self.getModel().createCompartment()
        rpSBML.checklibSBML(comp, 'create compartment')
        rpSBML.checklibSBML(comp.setId(compId), 'set compartment id')
        if compName:
            rpSBML.checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
        rpSBML.checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        rpSBML.checklibSBML(comp.setSize(size), 'set compartment "size"')
        rpSBML.checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if not meta_id:
            meta_id = self._genMetaID(compId)
        rpSBML.checklibSBML(comp.setMetaId(meta_id), 'set the meta_id for the compartment')
        ############################ MIRIAM ############################
        comp.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id)))
        # print(libsbml.XMLNode.convertXMLNodeToString(comp.getAnnotation()))
        self.addUpdateMIRIAM(comp, 'compartment', compXref, meta_id)
        # print(libsbml.XMLNode.convertXMLNodeToString(comp.getAnnotation()))
        # print()


    def createUnitDefinition(self, unit_id, meta_id=None):
        """Create libSBML unit definition

        Function that creates a unit definition (composed of one or more unit)

        :param unit_id: Unit id definition
        :param meta_id: Meta id (Default: None)

        :type unit_id: str
        :type meta_id: str

        :rtype: libsbml.UnitDefinition
        :return: Unit definition object created
        """
        unitDef = self.getModel().createUnitDefinition()
        rpSBML.checklibSBML(unitDef, 'creating unit definition')
        rpSBML.checklibSBML(unitDef.setId(unit_id), 'setting id')
        if not meta_id:
            meta_id = self._genMetaID(unit_id)
        rpSBML.checklibSBML(unitDef.setMetaId(meta_id), 'setting meta_id')
        # self.unitDefinitions.append(unit_id)
        return unitDef


    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        """Set or update the parameters of a libSBML unit definition

        :param unitDef: libSBML Unit
        :param libsbmlunit: String unit
        :param exponent: Exponent unit
        :param sale: Scale of the unit
        :param multiplier: Multiplier of the unit

        :type unitDef: libsbml.Unit
        :type libsbmlunit: str
        :type exponent: int
        :type sale: int
        :type multiplier: int

        :rtype: None
        :return: None
        """
        unit = unitDef.createUnit()
        rpSBML.checklibSBML(unit, 'creating unit')
        rpSBML.checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        rpSBML.checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        rpSBML.checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        rpSBML.checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    def createReturnFluxParameter(
        self,
        value,
        units='mmol_per_gDW_per_hr',
        is_constant=True,
        parameter_id=None,
        meta_id=None
    ):
        """Create libSBML flux parameters

        Parameters are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition.
        If the parameter id exists, then the function returns the libsbml.Parameter object

        :param value: Value set for the parameter
        :param unit: The unit id of the parameter (Default: mmol_per_gDW_per_hr)
        :param is_constant: Define if the parameter is constant (Default: True)
        :param parameter_id: Overwrite the default naming convention (Default: None)
        :param meta_id: Meta id (Default: None)

        :type value: float
        :type unit: str
        :type is_constant: bool
        :type parameter_id: str
        :type meta_id: str

        :rtype: libsbml.Parameter
        :return: The newly created libsbml.Parameter
        """
        if parameter_id:
            param_id = parameter_id
        else:
            if value >= 0:
                param_id = 'BRS_FBC_'+str(round(abs(value), 4)).replace('.', '_')
            else:
                param_id = 'BRS_FBC__'+str(round(abs(value), 4)).replace('.', '_')
        if param_id in [i.getId() for i in self.getModel().getListOfParameters()]:
            return self.getModel().getParameter(param_id)
        else:
            newParam = self.getModel().createParameter()
            rpSBML.checklibSBML(newParam, 'Creating a new parameter object')
            rpSBML.checklibSBML(newParam.setConstant(is_constant), 'setting as constant')
            rpSBML.checklibSBML(newParam.setId(param_id), 'setting ID')
            rpSBML.checklibSBML(newParam.setValue(value), 'setting value')
            rpSBML.checklibSBML(newParam.setUnits(units), 'setting unit')
            rpSBML.checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
            if not meta_id:
                meta_id = self._genMetaID(parameter_id)
            rpSBML.checklibSBML(newParam.setMetaId(meta_id), 'setting meta ID')
            # self.parameters.append(parameter_id)
            return newParam


    # TODO as of now not generic, works when creating a new SBML file, but no checks if modifying existing SBML file
    def createReaction(
        self,
        id: str,
        reactants: Dict[str, float],
        products: Dict[str, float],
        smiles: str,
        fbc_upper: float,
        fbc_lower: float,
        fbc_units: str,
        reversible: bool = True,
        reacXref: Dict = {},
        infos: Dict = {},
        meta_id: str = None
    ):
        """Create libSBML reaction

        Create a reaction that is added to the self.model in the input compartment id.
        fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with
        the first position that is the upper bound and the second is the lower bound.
        reactants_dict and reactants_dict are dictionnaries that hold the following
        parameters: name, compartment, stoichiometry

        :param name: Name of the reaction
        :param fluxUpperBound: The reaction fbc upper bound
        :param fluxLowerBound: The reaction fbc lower bound
        :param rxn: The id's of the reactant and products of the reactions. Example: {'left': [], 'right': []}
        :param compartment_id: The id of the compartment to add the reaction
        :param reaction_smiles: The reaction rule to add to the BRSynth annotation of the reaction (Default: None)
        :param reacXref: The dict containing the MIRIAM annotation (Default: {})
        :param pathway_id: The Groups id of the reaction to which the reacion id will be added (Default: None)
        :param meta_id: Meta id (Default: None)

        :type name: str
        :type fluxUpperBound: float
        :type fluxLowerBound: float
        :type rxn: dict
        :type compartment_id: str
        :type reaction_smiles: str
        :type reacXref: dict
        :type pathway_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """

        reac = self.getModel().createReaction()
        rpSBML.checklibSBML(reac, 'create reaction')

        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        rpSBML.checklibSBML(reac_fbc, 'extending reaction for FBC')
        ## BOUNDS
        # Upper
        upper_bound = self.createReturnFluxParameter(
            value=fbc_upper,
            units=fbc_units
        )
        rpSBML.checklibSBML(
            reac_fbc.setUpperFluxBound(
                upper_bound.getId()
            ),
            f'setting {str(id)} upper flux bound'
        )
        # Lower
        lower_bound = self.createReturnFluxParameter(
            value=fbc_lower,
            units=fbc_units
        )
        rpSBML.checklibSBML(
            reac_fbc.setLowerFluxBound(
                lower_bound.getId()
            ),
            f'setting {str(id)} lower flux bound'
        )

        #########################################
        # reactions
        rpSBML.checklibSBML(reac.setId(id), 'set reaction id') # same convention as cobrapy
        rpSBML.checklibSBML(reac.setSBOTerm(176), 'setting the system biology ontology (SBO)') # set as process
        # TODO: consider having the two parameters as input to the function
        rpSBML.checklibSBML(reac.setReversible(reversible), 'set reaction reversibility flag')
        rpSBML.checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if not meta_id:
            meta_id = self._genMetaID(id)
        rpSBML.checklibSBML(reac.setMetaId(meta_id), 'setting species meta_id')
        # TODO: check that the species exist

        for spe_id, spe_sto in reactants.items():
            spe_id = rpSBML.formatId(spe_id)
            spe = reac.createReactant()
            rpSBML.checklibSBML(
                spe,
                'create reactant'
            )
            # use the same writing convention as CobraPy
            reactant_ext = spe_id
            rpSBML.checklibSBML(
                spe.setSpecies(str(reactant_ext)),
                'assign reactant species'
            )
            # TODO: check to see the consequences of heterologous parameters not being constant
            rpSBML.checklibSBML(
                spe.setConstant(True),
                'set "constant" on species '+str(reactant_ext)
            )
            rpSBML.checklibSBML(
                spe.setStoichiometry(float(spe_sto)),
                'set stoichiometry ('+str(float(spe_sto))+')'
            )

        for spe_id, spe_sto in products.items():
            spe_id = rpSBML.formatId(spe_id)
            pro = reac.createProduct()
            rpSBML.checklibSBML(pro, 'create product')
            product_ext = spe_id
            rpSBML.checklibSBML(pro.setSpecies(str(product_ext)), 'assign product species')
            # TODO: check to see the consequences of heterologous parameters not being constant
            rpSBML.checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product_ext))
            rpSBML.checklibSBML(pro.setStoichiometry(float(spe_sto)),
                'set the stoichiometry ('+str(float(spe_sto))+')')

        ############################ MIRIAM ############################
        rpSBML.checklibSBML(reac.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        self.addUpdateMIRIAM(reac, 'reaction', reacXref, meta_id)
        ###### BRSYNTH additional information ########
        self.updateBRSynth(
            sbase_obj=reac,
            annot_header='smiles',
            value=smiles,
            meta_id=meta_id
        )
        for key, value in infos.items():
            self.updateBRSynth(
                sbase_obj=reac,
                annot_header=key,
                value=value
            )

    @staticmethod
    def formatId(id: str, logger: Logger=getLogger()) -> str:
        """Convert a string to a valid libSBML Id format.

        Parameters
        ----------
        id: str
            Id to format.
        logger : Logger, optional

        Returns
        -------
        Id well-formatted according to libSBML SID convention.
        """
        id = id.replace('-', '_')
        # if a char is neither a letter nor a digit,
        # replace it with its ascii number
        _id = id
        for i in range(len(id)):
            if not id[i].isalnum() and id[i] != '_':
                _id = f'{id[:i]}_ASCII_{str(ord(id[i]))}_ASCII_{id[i+1:]}'
        # id = comp_succ(id, '_')
        if _id[0].isdigit():
            _id = f'_{_id}'
        return _id

    def createSpecies(
        self,
        species_id: str,
        species_name: str=None,
        chemXref: List=[],
        inchi: str=None,
        inchikey: str=None,
        smiles: str=None,
        compartment: str=None,
        meta_id: str=None,
        infos: Dict={},
        is_boundary: bool=False
    ) -> None:
        """Create libSBML species

        Create a species that is added to self.model

        :param species_id: The id of the created species
        :param species_name: Overwrite the default name of the created species (Default: None)
        :param chemXref: The list containing the MIRIAM annotation (Default: {})
        :param inchi: The InChI string to be added to BRSynth annotation (Default: None)
        :param inchikey: The InChIkey string to be added to BRSynth annotation (Default: None)
        :param smiles: The SMLIES string to be added to BRSynth annotation (Default: None)
        :param compartment: The id of the compartment to add the reaction
        :param meta_id: Meta id (Default: None)
        :param infos: Supplemental informations to be added to BRSynth annotation (Default: None)
        :param is_boundary: Set if the specie is boundary (Default: False)

        :type species_id: str 
        :type species_name: str
        :type chemXref: List
        :type inchikey: str
        :type smiles: str
        :type compartment: str 
        :type meta_id: str
        :type infos: Dict
        :type is_boundary: bool 

        :rtype: None
        :return: None
        """
        spe = self.getModel().createSpecies()
        rpSBML.checklibSBML(spe, 'create species')

        # FBC.
        spe_fbc = spe.getPlugin('fbc')
        rpSBML.checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        # spe_fbc.setCharge(charge) #### These are not required for FBA
        # spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        # # if compartment_id:
        # rpSBML.checklibSBML(spe.setCompartment(compartment_id), 'set species spe compartment')
        # else:
        #    # removing this could lead to errors with xref
        #    rpSBML.checklibSBML(spe.setCompartment(self.compartment_id), 'set species spe compartment')
        # ID same structure as cobrapy
        # TODO: determine if this is always the case or it will change
        rpSBML.checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance unit')
        rpSBML.checklibSBML(spe.setBoundaryCondition(is_boundary), 'set boundary conditions')
        rpSBML.checklibSBML(spe.setConstant(False), 'set constant')
        # useless for FBA (usefull for ODE) but makes Copasi stop complaining
        rpSBML.checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        # same writting convention as COBRApy
        # rpSBML.checklibSBML(spe.setId(str(species_id)+'__64__'+str(compartment_id)), 'set species id')

        species_id = rpSBML.formatId(species_id)
        rpSBML.checklibSBML(spe.setId(str(species_id)), 'set species id')
        if not meta_id:
            meta_id = self._genMetaID(species_id)
        rpSBML.checklibSBML(spe.setMetaId(meta_id), 'setting reaction meta_id')
        if not species_name:
            rpSBML.checklibSBML(spe.setName(species_id), 'setting name for the metabolite '+str(species_id))
        else:
            rpSBML.checklibSBML(spe.setName(species_name), 'setting name for the metabolite '+str(species_name))
        if compartment:
            rpSBML.checklibSBML(spe.setCompartment(compartment), 'setting compartment')
        # this is setting MNX id as the name
        # this is setting the name as the input name
        # rpSBML.checklibSBML(spe.setAnnotation(self._defaultBRSynthAnnot(meta_id)), 'creating annotation')
        rpSBML.checklibSBML(spe.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        # Annotations.
        self.addUpdateMIRIAM(spe, 'species', chemXref, meta_id)
        if smiles is not None:
            self.updateBRSynth(
                sbase_obj=spe,
                annot_header='smiles',
                value=smiles,
                meta_id=meta_id
            )
        if inchi is not None:
            self.updateBRSynth(
                sbase_obj=spe,
                annot_header='inchi',
                value=inchi,
                meta_id=meta_id
            )
        if inchikey is not None:
            self.updateBRSynth(
                sbase_obj=spe,
                annot_header='inchikey',
                value=inchikey,
                meta_id=meta_id
            )
        for key, value in infos.items():
            self.updateBRSynth(
                sbase_obj=spe,
                annot_header=key,
                value=value,
                meta_id=meta_id
            )


    def createGroup(
        self,
        id: str,
        brs_annot: bool = True,
        meta_id: str = None
    ) -> libsbml.GroupsModelPlugin:
        """Create libSBML pathway

        Create a group that is added to self.model

        :param id: The Groups id of the pathway id
        :param meta_id: Meta id (Default: None)

        :type id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        groups_plugin = self.getModel().getPlugin('groups')
        new_group = groups_plugin.createGroup()
        new_group.setId(id)
        if meta_id is None:
            meta_id = self._genMetaID(id)
        new_group.setMetaId(meta_id)
        new_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        if brs_annot:
            new_group.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        return new_group


    def addMember(
        self,
        group_id: str,
        idRef: str
    ) -> None:
        group = self.getGroup(group_id)
        member = libsbml.Member()
        member.setIdRef(idRef)
        group.addMember(member)
        rpSBML.checklibSBML(member, 'Creating a new groups member')
        rpSBML.checklibSBML(member.setIdRef(idRef), 'Setting name to the groups member')


    def createFluxObj(self, fluxobj_id, reactionName, coefficient, isMax=True, meta_id=None):
        """Create libSBML flux objective

        WARNING DEPRECATED -- use the createMultiFluxObj() with lists of size one to define an objective function
        with a single reaction
        Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives

        :param fluxobj_id: The id of the flux objective
        :param reactionName: The id of the reaction that is associated with the reaction
        :param coefficient: The coefficient of the flux objective
        :param isMax: Define if the objective is coefficient (Default: True)
        :param meta_id: Meta id (Default: None)

        :type fluxobj_id: str
        :type reactionName: str
        :type coefficient: int
        :type isMax: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """
        fbc_plugin = self.getModel().getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        # TODO: need to define inpiut metaID
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reactionName)
        target_flux_obj.setCoefficient(coefficient)
        if not meta_id:
            meta_id = self._genMetaID(str(fluxobj_id))
        target_flux_obj.setMetaId(meta_id)
        target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def activateObjective(
        self,
        objective_id: str,
        plugin: str = 'fbc'
    ) -> None:

        self.logger.debug(f'Activate objective {objective_id}')

        plugin = self.getPlugin(plugin)
        self.checklibSBML(
            plugin.setActiveObjectiveId(objective_id),
            'Setting active objective '+str(objective_id)
        )


    ##############################################################################################
    ############################### Generic Model ################################################
    ##############################################################################################


    def genericModel(
        self,
        modelName,
        modelID,
        compartments,
        unit_def
    ):
        """Generate a generic model

        Since we will be using the same type of parameters for the RetroPath model, this function
        generates a libSBML model with parameters that will be mostly used

        :param modelName: The given name of the model
        :param modelID: The id of the model
        :param compXref: The model MIRIAM annotation
        :param compartment_id: The id of the model compartment
        :param upper_flux_bound: The upper flux bounds unit definitions default when adding new reaction (Default: 999999.0)
        :param lower_flux_bound: The lower flux bounds unit definitions default when adding new reaction (Defaul: 0.0)

        :type modelName: str
        :type modelID: str
        :type compXref: dict
        :type compartment_id: str
        :type upper_flux_bound: float
        :type lower_flux_bound: float

        :rtype: None
        :return: None
        """
        self.logger.debug(f'Model Name: {modelName}')
        self.logger.debug(f'Model ID: {modelID}')
        self.logger.debug(f'Compartments: {compartments}')
        self.logger.debug(f'Unit Definitions: {unit_def}')

        self.createModel(modelName, modelID)
        for unit_id, unit_data in unit_def.items():
            unitDef = self.createUnitDefinition(unit_id)
            for unit in unit_data:
                self.createUnit(unitDef, unit['kind'], unit['exponent'], unit['scale'], unit['multiplier'])

        # compartments
        for comp_id, comp in compartments.items():
            comp_lst = []
            for db, comps in comp['annot'].items():
                for comp_name in comps:
                    comp_lst += [f"http://identifiers.org/{self.miriam_header['compartment'][db]}{comp_name}"]
            self.createCompartment(1, comp_id, comp['name'], comp_lst)
            # self.createCompartment(1, comp_id, comp['name'], comp['annot'])

        self.logger.debug('Generic model created')