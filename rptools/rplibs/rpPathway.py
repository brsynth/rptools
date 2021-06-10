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
    TypeVar,
    Union,
    Callable
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
    Reaction,
    Compound
)

def gen_dict_extract(key, var):
    if hasattr(var,'items'):
        for k, v in var.items():
            print(k, v)
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result

class rpPathway(Pathway):

    def __init__(
        self,
        id: str,
        cache: Cache = None,
        logger: Logger = getLogger(__name__)
    ):
        super().__init__(
            id=id,
            cache=cache,
            logger=logger
        )
        self.set_target_id(None)
        self.set_sink(None)

    ## OUT METHODS
    # def __repr__(self):
    #     return dumps(self.to_dict(), indent=4)

    def _to_dict(self) -> Dict:
        return {
            **super()._to_dict(),
            **self.__to_dict()
        }

    def __to_dict(self) -> Dict:
        return {
            'sink': deepcopy(self.get_sink()),
            'target': self.get_target_id()
        }

    def __eq__(self, other) -> bool:
        if isinstance(self, other.__class__):
            return self._to_dict() == other._to_dict()
        return False

    ## READ METHODS
    def get_target_rxn_id(self) -> str:
        for rxn in self.get_reactions():
            if self.get_target_id() in rxn.get_products().keys():
                return rxn.get_id()

    def get_rxn_target(self) -> Reaction:
        return self.get_reaction(self.get_target_rxn_id())

    def get_target_id(self) -> str:
        return self.__target_id

    def get_target(self) -> Compound:
        return self.get_specie(self.get_target_id())

    def get_sink(self) -> List[str]:
        return self.__sink

    def get_reactions_ids(self) -> List[str]:
        '''Returns the list of reaction IDs sorted by index within the pathway
        (forward direction).
        '''
        return [
            rxn_id for rxn_id in sorted(
                super().get_reactions_ids(),
                key=lambda x: self.get_reaction(x).get_info('idx_in_path')
            )
        ]

    def get_reaction_idx_in_path(self, rxn_id: str) -> int:
        return self.get_reactions_ids().index(rxn_id) + 1

    def __get_infos(self, key: str) -> Dict[str, Dict]:
        species = {}
        for spe_id in self.get_species_ids():
            species[spe_id] = {k: v for k, v in self.get_specie(spe_id).get_infos().items() if k.startswith(key)}
        reactions = {}
        for rxn_id in self.get_reactions_ids():
            reactions[rxn_id] = {k: v for k, v in self.get_reaction(rxn_id).get_infos().items() if k.startswith(key)}
        pathway = {k: v for k, v in self.get_infos().items() if k.startswith(key)}
        return {
            'species': species,
            'reactions': reactions,
            'pathway': pathway
        }

    # def __getattr__(self, name):
    #     '''Catch get_<key>_infos() call an returns the corresponding infos'''
    #     from re import compile as re_compile
    #     def method(*args):
    #         infos = re_compile(r'get_(\w+)_infos')
    #         pattern = infos.search(name).group(1)
    #         key = pattern.split('_')[0]
    #         return self.__get_infos(key)
    #         # print(self.__get_infos(key))
    #         # print(pattern)
    #         # for i in gen_dict_extract(pattern, self.__get_infos(key)):
    #         #     print(i)
    #         return gen_dict_extract(pattern, self.__get_infos(key))
    #     return method

    ## THERMO - GLOBAL
    def __get_thermo_infos(self) -> Dict[str, Dict]:
        return self.__get_infos('thermo')

    def __get_thermo(self) -> Dict[str, Dict]:
        return self.__get_thermo_infos()['pathway']

    def __get_thermo_info(self, key: str) -> Dict:
        try:
            return self.__get_thermo()[key]
        except KeyError:
            self.get_logger().warning(f'There is no {key} value for this pathway')
            return {}

    def __get_dict_key(self, dict: Dict, key: str) -> TypeVar:
        try:
            return dict[key]
        except KeyError:
            self.get_logger().warning(f'There is no {key} key in data')
            return None

    def get_thermo_dG0_prime(self) -> Dict:
        return self.__get_thermo_info('thermo_dG0_prime')

    def get_thermo_dG0_prime_value(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG0_prime(), 'value')

    def get_thermo_dG0_prime_error(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG0_prime(), 'error')

    def get_thermo_dG0_prime_units(self) -> str:
        return self.__get_dict_key(self.get_thermo_dG0_prime(), 'units')

    def get_thermo_dGm_prime(self) -> Dict:
        return self.__get_thermo_info('thermo_dGm_prime')

    def get_thermo_dGm_prime_value(self) -> float:
        return self.__get_dict_key(self.get_thermo_dGm_prime(), 'value')

    def get_thermo_dGm_prime_error(self) -> float:
        return self.__get_dict_key(self.get_thermo_dGm_prime(), 'error')

    def get_thermo_dGm_prime_units(self) -> str:
        return self.__get_dict_key(self.get_thermo_dGm_prime(), 'units')

    def get_thermo_dG_prime(self) -> Dict:
        return self.__get_thermo_info('thermo_dG_prime')

    def get_thermo_dG_prime_value(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG_prime(), 'value')

    def get_thermo_dG_prime_error(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG_prime(), 'error')

    def get_thermo_dG_prime_units(self) -> str:
        return self.__get_dict_key(self.get_thermo_dG_prime(), 'units')

    def get_thermo_dG(self) -> Dict:
        return self.__get_dict_key(self.get_thermo_dG(), 'value')

    def get_thermo_dG_value(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG(), 'error')

    def get_thermo_dG_error(self) -> float:
        return self.__get_dict_key(self.get_thermo_dG(), 'units')

    def get_thermo_dG_units(self) -> str:
        return self.get_thermo_dG()['units']

    ## THERMO - COMPOUND
    def __get_thermo_species(self) -> Dict[str, Dict]:
        return self.get_thermo_infos()['species']

    def __get_thermo_compounds(self) -> Dict[str, Dict]:
        return self.__get_thermo_species()

    def __get_thermo_compound(self, cmpd_id: str) -> Dict:
        try:
            return self.__get_thermo_compounds()[cmpd_id]
        except KeyError:
            self.get_logger().warning(f'There is no compound {cmpd_id} in this pathway')
            return {}

    def __get_thermo_compound_info(self, cmpd_id: str, key: str) -> Dict:
        try:
            return self.__get_thermo_compound(cmpd_id)[key]
        except KeyError:
            self.get_logger().warning(f'There is no {key} value for the compound {cmpd_id}')
            return {}

    def get_thermo_compound_std_dg_form(self, cmpd_id: str) -> Dict:
        return self.__get_thermo_compound_info(cmpd_id, 'thermo_standard_dg_formation')

    def get_thermo_compound_std_dg_form_value(self, cmpd_id: str) -> float:
        return self.__get_dict_key(self._get_thermo_compound_std_dg_form(cmpd_id), 'value')

    def get_thermo_compound_std_dg_form_units(self, cmpd_id: str) -> float:
        return self.__get_dict_key(self._get_thermo_compound_std_dg_form(cmpd_id), 'units')

    ## THERMO - REACTION
    def __get_thermo_reactions(self) -> Dict[str, Dict]:
        return self.get_thermo_infos()['reactions']

    def __get_thermo_reaction(self, rxn_id: str) -> Dict:
        try:
            return self.__get_thermo_reactions()[rxn_id]
        except KeyError:
            self.get_logger().warning(f'There is no reaction {rxn_id} in this pathway')
            return {}

    def __get_thermo_reaction_info(self, rxn_id: str, key: str) -> Dict:
        try:
            return self.__get_thermo_reaction(rxn_id)[key]
        except KeyError:
            self.get_logger().warning(f'There is no {key} value for the reaction {rxn_id}')
            return {}

    def get_thermo_reaction_dG0_prime(self, rxn_id: str) -> Dict:
        return self.__get_thermo_reaction_info(rxn_id, 'thermo_dG0_prime')

    def get_thermo_reaction_dG0_prime_value(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG0_prime(rxn_id), 'value')

    def get_thermo_reaction_dG0_prime_error(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG0_prime(rxn_id), 'error')

    def get_thermo_reaction_dG0_prime_units(self, rxn_id: str) -> str:
        return self.__get_dict_key(self.get_thermo_reaction_dG0_prime(rxn_id), 'units')

    def get_thermo_reaction_dGm_prime(self, rxn_id: str) -> Dict:
        return self.__get_thermo_reaction_info(rxn_id, 'thermo_dGm_prime')

    def get_thermo_reaction_dGm_prime_value(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dGm_prime(rxn_id), 'value')

    def get_thermo_reaction_dGm_prime_error(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dGm_prime(rxn_id), 'error')

    def get_thermo_reaction_dGm_prime_units(self, rxn_id: str) -> str:
        return self.__get_dict_key(self.get_thermo_reaction_dGm_prime(rxn_id), 'units')

    def get_thermo_reaction_dG_prime(self, rxn_id: str) -> Dict:
        return self.__get_thermo_reaction_info(rxn_id, 'thermo_dG_prime')

    def get_thermo_reaction_dG_prime_value(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG_prime(rxn_id), 'value')

    def get_thermo_reaction_dG_prime_error(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG_prime(rxn_id), 'error')

    def get_thermo_reaction_dG_prime_units(self, rxn_id: str) -> str:
        return self.__get_dict_key(self.get_thermo_reaction_dG_prime(rxn_id), 'units')

    def get_thermo_reaction_dG(self, rxn_id: str) -> Dict:
        return self.__get_thermo_reaction_info(rxn_id, 'thermo_dG')

    def get_thermo_reaction_dG_value(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG(rxn_id), 'value')

    def get_thermo_reaction_dG_error(self, rxn_id: str) -> float:
        return self.__get_dict_key(self.get_thermo_reaction_dG(rxn_id), 'error')

    def get_thermo_reaction_dG_units(self, rxn_id: str) -> str:
        return self.__get_dict_key(self.get_thermo_reaction_dG(rxn_id), 'units')

    ## FBA - GLOBAL
    def __get_fba(self) -> Dict[str, Dict]:
        return self.__get_infos('fba')

    def __get_fba_biomass(self) -> Dict[str, Dict]:
        return self.__get_infos('fba_biomass')

    def __get_fba_fraction(self) -> Dict[str, Dict]:
        return self.__get_infos('fba_fraction')

    def get_fba_biomass(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_fba_biomass(),
            'pathway'
        )

    def get_fba_fraction(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_fba_fraction(),
            'pathway'
        )

    ## FBA - REACTIONS
    def __get_fba_reactions(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_infos('fba'),
            'reactions'
        )

    def __get_fba_reactions_biomass(self) -> Dict[str, Dict]:
        return self.__get_infos('fba_biomass')

    def __get_fba_reactions_fraction(self) -> Dict[str, Dict]:
        return self.__get_infos('fba_fraction')

    def get_fba_reaction_biomass(self, rxn_id: str) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_fba_reactions_biomass(),
            'reactions'
        )

    def get_fba_reaction_fraction(self, rxn_id: str) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_fba_reactions_fraction(),
            'reactions'
        )

    ## FBA - COMPOUNDS
    def __get_fba_compounds(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_infos('fba'),
            'species'
        )

    def __get_fba_compounds_biomass_shadow_prices(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_infos('fba_biomass_shadow_price'),
            'species'
        )

    def get_fba_compound_biomass_shadow_price(self, cmpd_id) -> float:
        return self.__get_dict_key(
            self.__get_fba_compounds_biomass_shadow_prices(),
            cmpd_id
        )

    def __get_fba_compounds_fraction_shadow_prices(self) -> Dict[str, Dict]:
        return self.__get_dict_key(
            self.__get_infos('fba_fraction_shadow_price'),
            'species'
        )

    def get_fba_compound_fraction_shadow_price(self, cmpd_id) -> float:
        return self.__get_dict_key(
            self.__get_fba_compounds_fraction_shadow_prices(),
            cmpd_id
        )

    def __get_rpSBML_infos(self) -> Dict[str, Dict]:
        return self.__get_infos('rpSBML')

    def __get_rpSBML_info(self, key: str) -> Dict[str, Dict]:
        try:
            return self.__get_infos('rpSBML')[key]
        except KeyError:
            self.__logger(f'There is no field {key} within rpSBML infos.')
            return {}


    ## WRITE METHODS
    def add_reaction(
        self,
        rxn: Reaction,
        rxn_id: str = None,
        target_id: str = None
    ) -> None:

        super().add_reaction(rxn, rxn_id)

        # TARGET
        if target_id is not None:
            self.set_target_id(target_id)

    def set_sink(self, sink: List) -> None:
        self.__sink = deepcopy(sink)
    
    def set_target_id(self, id: str) -> None:
        self.__target_id = id

    def rename_compound(self, id: str, new_id: str) -> None:
        # target
        if id == self.get_target_id():
            self.set_target_id(new_id)

        super().rename_compound(id, new_id)

        # sink
        try:
            self.__sink[self.__sink.index(id)] = new_id
        except ValueError:
            pass


    ## MISC
