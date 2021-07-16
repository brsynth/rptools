"""A class to represent a metabolic pathway."""
# The MIT License (MIT)
#
# Copyright (c) 2018 Institute for Molecular Systems Biology, ETH Zurich.
# Copyright (c) 2019 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from typing import (
    Dict,
    List,
    Set,
    Union,
    TypeVar,
    Tuple
)
from logging import (
    Logger,
    getLogger,
    ERROR
)
from copy import deepcopy
from brs_utils import Cache
from chemlite import (
    Pathway,
    Compound
)
from numpy import isin
from .rpSBML import rpSBML
from .rpReaction import rpReaction
from .rpCompound import rpCompound
from .rpObject import rpObject


# def gen_dict_extract(key, var):
#     if hasattr(var,'items'):
#         for k, v in var.items():
#             print(k, v)
#             if k == key:
#                 yield v
#             if isinstance(v, dict):
#                 for result in gen_dict_extract(key, v):
#                     yield result
#             elif isinstance(v, list):
#                 for d in v:
#                     for result in gen_dict_extract(key, d):
#                         yield result

class rpPathway(Pathway, rpObject):

    __MNXC3 = {
        'id': 'c',
        'name': 'cytosol',
        'annot': {
            'name': ['cytosol'],
            'seed': ['cytosol', 'c0', 'c'],
            'mnx': ['MNXC3'],
            'bigg': ['c_c', 'c']
        }
    }

    def __init__(
        self,
        id: str,
        cache: Cache = None,
        logger: Logger = getLogger(__name__)
    ):
        Pathway.__init__(
            self,
            id=id,
            cache=cache,
            logger=logger
        )
        rpObject.__init__(self, logger)
        self.__species_groups = {}
        self.set_target_id(None)
        self.__unit_def = {}
        self.__compartments = {}
        self.add_compartment(
            id=rpPathway.__MNXC3['id'],
            name=rpPathway.__MNXC3['name'],
            annot=rpPathway.__MNXC3['annot'],
        )
        # Set flux bounds values/units
        self.__parameters = {}
        self.add_parameter(
            id='BRS_default_fbc_l',
            value=-10000,
            units=rpReaction.get_default_fbc_units()
        )
        self.add_parameter(
            id='BRS_default_fbc_u',
            value=10000,
            units=rpReaction.get_default_fbc_units()
        )
        # Additional names for methods
        self.get_sink = self.get_sink_species
        self.set_sink = self.set_sink_species

    ## OUT METHODS
    # def __repr__(self):
    #     return dumps(self.to_dict(), indent=4)

    def _to_dict(
        self,
        specific: bool = False
    ) -> Dict:
        if specific:
            return {
                **self.__to_dict(),
                **rpObject._to_dict(self)
            }
        else:
            return {
                **Pathway._to_dict(self),
                **rpObject._to_dict(self),
                **self.__to_dict()
            }

    def __to_dict(self) -> Dict:
        return {
            'sink': deepcopy(self.get_sink_species()),
            'target': self.get_target_id(),
            'parameters': deepcopy(self.get_parameters()),
            'unit_defs': deepcopy(self.get_unit_defs()),
            'compartments': deepcopy(self.get_compartments()),
        }

    def __eq__(self, other) -> bool:
        if not isinstance(self, other.__class__):
            return False
        # Compare with specific keys
        return all(
            (
                self._to_dict().get(key) == other._to_dict().get(key)
                or self._to_dict().get(key) is other._to_dict().get(key)
            )
            for key in [
                'reactions',
                'target',
            ]
        )

    ## READ METHODS
    def get_species_groups(self) -> Dict[str, Set]:
        return self.__species_groups
        # return rpPathway.__SPECIES_GROUPS

    def get_species_group(self, group_id: str) -> List[str]:
        return list(self.get_species_groups().get(group_id, []))

    def get_completed_species(self) -> List[str]:
        return self.get_species_group('completed')
        # return self.__completed_species

    def get_fba_ignored_species(self) -> List[str]:
        return self.get_species_group('fba_ignored')
        # return self.__fba_ignored_species

    def get_trunk_species(self) -> List[str]:
        return self.get_species_group('trunk')
        # return self.__trunk_species

    def get_thermo_substituted_species(self) -> List[str]:
        return self.get_species_group('thermo_substituted')

    def get_sink_species(self) -> List[str]:
        return self.get_species_group('sink')
        # return self.__sink

    def get_intermediate_species(self) -> List[str]:
        return self.get_species_group('intermediate')

    def get_target_rxn_id(self) -> str:
        for rxn in self.get_list_of_reactions():
            if self.get_target_id() in rxn.get_products_ids():
                return rxn.get_id()

    def get_rxn_target(self) -> rpReaction:
        return self.get_reaction(self.get_target_rxn_id())

    def get_target_id(self) -> str:
        return self.__target_id

    def get_target(self) -> Compound:
        return self.get_specie(self.get_target_id())

    def get_reactions_ids(self) -> List[str]:
        '''Returns the list of reaction IDs sorted by index within the pathway
        (forward direction).
        '''
        return [
            rxn_id for rxn_id in sorted(
                super().get_reactions_ids(),
                key=lambda x: self.get_reaction(x).get_idx_in_path()
            )
        ]

    def get_parameters(self) -> Dict:
        try:
            return self.__parameters
        except AttributeError:
            return None

    def get_parameter(self, id: str) -> Dict:
        return self.__parameters.get(id, {})

    def get_parameter_value(self, id: str) -> Dict:
        return self.get_parameter(id).get('value', 'NaN')

    def get_parameter_units(self, id: str) -> Dict:
        return self.get_parameter(id).get('units', str(''))

    def get_unit_defs(self) -> Dict:
        return self.__unit_def

    def get_unit_def(self, id: str) -> Dict:
        return self.__unit_def.get(id, {})

    def get_compartments(self) -> Dict:
        return self.__compartments


    ## WRITE METHODS
    def add_trunk_species(self, species: List[str]) -> None:
        self.set_trunk_species(
            list(
                set(
                    self.get_trunk_species()
                    + species
                )
            )
        )

    def add_completed_species(self, species: List[str]) -> None:
        self.set_completed_species(
            list(
                set(
                    self.get_completed_species()
                    + species
                )
            )
        )

    def add_species_group(
        self,
        group_id: str,
        species: List[str]
    ) -> None:
        try:
            s = set(self.__species_groups[group_id])
            s.update(species)
            # Add species to the existing group
            self.__species_groups[group_id] = deepcopy(s)
        except KeyError:
            # Create a new group
            self.set_species_group(group_id, species)

    def set_species_group(
        self,
        group_id: str,
        species: Union[List[str], Dict[str, str]]
    ) -> None:
        # Create a new group
        if isinstance(species, list):
            self.__species_groups[group_id] = list(set(deepcopy(species)))
        elif isinstance(species, dict):
            self.__species_groups[group_id] = deepcopy(species)
        else:
            self.get_logger().error(f'Wrong type {type(species)} for \'species\' argument, \'list\' or \'dict\' expected, nothing set.')

    def set_completed_species(self, species: List[str]) -> None:
        self.add_species_group('completed', species)
        # self.__completed_species = deepcopy(species)

    def set_fba_ignored_species(self, species: List[str]) -> None:
        self.add_species_group('fba_ignored', species)
        # self.__fba_ignored_species = deepcopy(species)

    def set_trunk_species(self, species: List[str]) -> None:
        self.add_species_group('trunk', species)
        # (re-)build intermediate species group
        self.__build_intermediate_species()
        # self.__trunk_species = deepcopy(species)

    def set_thermo_substituted_species(self, species: List[str]) -> None:
        self.add_species_group('thermo_substituted', species)
        # self.__trunk_species = deepcopy(species)

    def set_sink_species(self, species: List[str]) -> None:
        self.add_species_group('sink', species)
        # (re-)build intermediate species group
        self.__build_intermediate_species()
        # self.__sink = deepcopy(sink)

    def __build_intermediate_species(self) -> None:
        try:
            members = list(
                set(self.get_trunk_species())
                - set(self.get_sink_species())
                - set([self.get_target_id()])
            )
            self.set_species_group('intermediate', members)
        except TypeError:
            self.set_species_group('intermediate', [])

    @staticmethod
    def from_rpSBML(
        infile: str = None,
        rpsbml: rpSBML = None,
        logger: Logger = getLogger(__name__)
    ) -> 'rpPathway':
        def write_to(data: Dict, object: TypeVar) -> None:
            # Detect fba and thermo infos
            for key, value in data.items():
                if value == 'None':
                    value = None
                elif str(value) == 'nan':
                    value = 'NaN'
                fba_offset = len(rpObject.get_fba_prefix()) + len(rpObject.get_sep())
                thermo_offset = len(rpObject.get_thermo_prefix()) + len(rpObject.get_sep())
                if key.startswith(rpObject.get_thermo_prefix()):
                    object.set_thermo_info(
                        key[thermo_offset:],
                        value
                    )
                elif key.startswith(rpObject.get_fba_prefix()):
                    object.set_fba_info(
                        key[fba_offset:],
                        value
                    )
                else:
                    try:
                        getattr(
                            object,
                            'set_'+key.replace('rp_', '')
                        )(value)
                    except AttributeError:
                        pass

        def build_reaction(
            rxn_id: str,
            infos: Dict,
            logger: Logger = getLogger(__name__)
        ) -> Tuple[
            rpReaction,
            Union[str, None]
        ]:
            # try:
            #     ec_numbers = infos['miriam']['ec-code']
            # except KeyError:
            #     ec_numbers = []
            reaction = rpReaction(
                id=rxn_id,
                # ec_numbers=ec_numbers,
                reactants=infos['left'],
                products=infos['right'],
                lower_flux_bound=infos['fbc_lower_value'],
                upper_flux_bound=infos['fbc_upper_value'],
                flux_bound_units=infos['fbc_units'],
                reversible=infos['reversible'],
                miriam=infos['miriam'],
                logger=logger
            )
            # Add additional infos
            write_to(infos['brsynth'], reaction)
            # Detects if the current reaction produces the target
            target_id = [spe_id for spe_id in reaction.get_products_ids() if 'TARGET' in spe_id]
            if target_id != []:
                target_id = target_id[0]
            else:
                target_id = None
            return reaction, target_id

        if infile is not None:
            rpsbml = rpSBML(inFile=infile, logger=logger)

        # Create the rpPathway object
        pathway = rpPathway(
            id=rpsbml.getName(),
            logger=logger
        )

        ## COMPARTMENTS
        for compartment in rpsbml.getModel().getListOfCompartments():
            pathway.add_compartment(
                id=compartment.getId(),
                name=compartment.getName(),
                annot=rpSBML.readMIRIAMAnnotation(compartment.getAnnotation()),
            )

        ## UNIT DEFINITIONS
        for unit_defs in rpsbml.getModel().getListOfUnitDefinitions():
            for unit in unit_defs.getListOfUnits():
                pathway.add_unit_def(
                    id=unit.getId(),
                    kind=unit.getKind(),
                    exp=unit.getExponent(),
                    scale=unit.getScale(),
                    mult=unit.getMultiplier()
                )

        ## PARAMETERS
        for param in rpsbml.getModel().getListOfParameters():
            pathway.add_parameter(
                id=param.getId(),
                value=param.getValue(),
                units=param.getUnits()
            )

        ## SPECIES
        for spe_id, spe in rpsbml.read_species().items():
            infos = {}
            for key in ['smiles', 'inchi', 'inchikey']:
                try:
                    infos[key] = spe['brsynth'][key]
                except KeyError:
                    infos[key] = ''
            # Create compound to add it in the cache
            compound = rpCompound(
                id=spe_id,
                smiles=infos['smiles'],
                inchi=infos['inchi'],
                inchikey=infos['inchikey'],
                compartment_id=spe['object'].getCompartment()
            )
            write_to(spe['brsynth'], compound)

        pathway_id = 'rp_pathway'

        ## REACTIONS
        for rxn_id, rxn_infos in rpsbml.read_reactions(pathway_id).items():
            rxn_infos['fbc_lower_value'] = pathway.get_parameter_value(rxn_infos['fbc_lower_value'])
            rxn_infos['fbc_upper_value'] = pathway.get_parameter_value(rxn_infos['fbc_upper_value'])
            rxn_infos['fbc_units'] = pathway.get_parameter_units(rxn_infos['fbc_lower_value'])
            reaction, target_id = build_reaction(rxn_id, rxn_infos, logger)
            # Add the reaction to the pathway
            pathway.add_reaction(
                rxn=reaction,
                target_id=target_id
            )

        ## GROUPS
        for group in rpsbml.getPlugin('groups').getListOfGroups():
            group_id = group.getId()
            # 'rp_pathway' has no member to write into rpPathway
            if group_id == pathway_id:
                annot = rpsbml.readBRSYNTHAnnotation(
                    rpsbml.getGroup(group_id).getAnnotation(),
                    rpsbml.logger
                )
                write_to(annot, pathway)
            # 'rp_sink_species', 'rp_completed_species', 'rp_trunk_species'
            # have no annotation to write into rpPathway
            else:
                write_to(
                    {
                        group_id: rpsbml.readGroupMembers(group_id)
                    },
                    pathway
                )

        return pathway

    def to_rpSBML(self) -> rpSBML:

        rpsbml = rpSBML(name='rp_'+self.get_id(), logger=self.get_logger())

        ## Create a generic Model, ie the structure and unit definitions that we will use the most
        rpsbml.genericModel(
            self.get_id(),
            'RP_model_'+self.get_id(),
            self.get_compartments(),
            self.get_unit_defs(),
            # upper_flux_bound,
            # lower_flux_bound
        )

        ## Create the groups (pathway, species, sink species)
        rpsbml.create_enriched_group(
            group_id='rp_pathway',
            members=self.get_reactions_ids(),
            infos=rpObject._to_dict(self)
        )
        for group_id, group_members in self.get_species_groups().items():
            rpsbml.create_enriched_group(
                group_id=f'rp_{group_id}_species',
                members=group_members
            )

        ## Add species to the model
        for specie in self.get_species():
            rpsbml.createSpecies(
                species_id=specie.get_id(),
                species_name=specie.get_name(),
                inchi=specie.get_inchi(),
                inchikey=specie.get_inchikey(),
                smiles=specie.get_smiles(),
                compartment=specie.get_compartment(),
                infos=self.get_specie(specie.get_id())._to_dict(specific=True)
            )

        ## Add reactions to the model
        for rxn in self.get_list_of_reactions():
            xref = {
                'ec': rxn.get_ec_numbers(),
                **rxn.get_miriam()
            }
            # Add the reaction in the model
            rpsbml.createReaction(
                id=rxn.get_id(),
                reactants=rxn.get_reactants(),
                products=rxn.get_products(),
                smiles=rxn.get_smiles(),
                fbc_upper=rxn.get_fbc_upper(),
                fbc_lower=rxn.get_fbc_lower(),
                fbc_units=rxn.get_fbc_units(),
                reversible=rxn.reversible(),
                reacXref=xref,
                infos=rxn._to_dict(specific=True)
            )

        return rpsbml

    def add_reaction(
        self,
        rxn: rpReaction,
        rxn_id: str = None,
        target_id: str = None
    ) -> None:

        super().add_reaction(rxn, rxn_id)

        # TARGET
        if target_id is not None:
            self.set_target_id(target_id)

    def set_target_id(self, id: str) -> None:
        self.__target_id = id
        # (re-)build intermediate species group
        self.__build_intermediate_species()

    def add_unit_def(
        self,
        id: str,
        kind: int,
        exp: int,
        scale: int,
        mult: float
    ) -> None:
        if id not in self.__parameters.keys():
            self.__parameters[id] = []
        self.__parameters[id] += [
            {
                'kind': kind,
                'exponent': exp,
                'scale': scale,
                'multiplier': mult
            }
        ]

    def add_parameter(
        self,
        id: str,
        value: float,
        units: str = rpReaction.get_default_fbc_units()
    ) -> None:
        # # Check if __parameters is defined
        # if not hasattr(self, '__parameters'):
        #     self.__parameters = {}
        # Check if id does not already exist in __parameters
        if self.get_parameter(id) == {}:
            self.__parameters[id] = {
                'value': value,
                'units': units
            }
        else:
            self.get_logger.warning(f'Parameter {id} already exist in rpPathway parameters, nothing added.')

    def add_compartment(
        self,
        id: str,
        name: str,
        annot: str
    ) -> None:
        if id not in self.get_compartments():
            self.__compartments[id] = {
                'name': name,
                'annot': annot
            }

    ## MISC
    def rename_compound(self, id: str, new_id: str) -> None:
        # target
        if id == self.get_target_id():
            self.set_target_id(new_id)

        super().rename_compound(id, new_id)

        # sink
        try:
            self.__sink[self.get_sink_species().index(id)] = new_id
        except (ValueError, AttributeError):
            pass

        # trunk species
        try:
            self.__trunk_species[self.get_trunk_species().index(id)] = new_id
        except (ValueError, AttributeError):
            pass

        # completed species
        try:
            self.__completed_species[self.get_completed_species().index(id)] = new_id
        except (ValueError, AttributeError):
            pass

    def cobraize(self, compartment_id: str) -> None:
        '''Make the Pathway compliant with what Cobra expects
        Add <@compartmentID> to all compounds in species and reactions
        '''
        from rptools.rpfba.cobra_format import (
            cobra_suffix,
            cobraize,
        )
        # SPECIES
        for spe_id in self.get_species_ids():
            if not spe_id.endswith(cobra_suffix(compartment_id)):
                self.rename_compound(
                    spe_id,
                    cobraize(
                        spe_id,
                        compartment_id
                    )
                )
            # pathway.get_specie(spe_id).set_id(cobraize_string(spe_id, pathway))

        # REACTIONS
        # cobraize_reactions(pathway)

    def uncobraize(self) -> None:
        '''Make the Pathway compliant with what Cobra expects
        Remove <@compartmentID> from all compounds in species, reactions and scores
        '''
        from rptools.rpfba.cobra_format import (
            uncobraize,
        )
        for spe_id in self.get_species_ids():
            self.rename_compound(spe_id, uncobraize(spe_id))

