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
        self.set_target_rxn_id(None)
        self.set_target_id(None)
        self.set_sink(None)

    ## OUT METHODS
    # def __repr__(self):
    #     return dumps(self.to_dict(), indent=4)

    def to_string(self):
        return '----------------\n' \
             + f'Pathway {self.get_id()}\n' \
             + '----------------\n' \
             + '\n'.join([rxn.__str__() for rxn in self.get_reactions()])

    def to_dict(self) -> Dict:
        return {
            **super().to_dict(),
            **{
                'sink': deepcopy(self.get_sink()),
                'target': {
                    'compd_id': deepcopy(self.get_target_id()),
                    'rxn_id': deepcopy(self.get_rxn_target_id())
                }
            }
        }

    def __eq__(self, other) -> bool:
        if isinstance(self, other.__class__):
            return self.to_dict() == other.to_dict()
        return False

    ## READ METHODS
    def get_target_rxn_id(self) -> str:
        return self.__rxn_target_id

    def get_rxn_target(self) -> Reaction:
        return self.get_reaction(self.get_target_rxn_id())

    def get_target_id(self) -> str:
        return self.__target_id

    def get_target(self) -> Compound:
        return self.get_specie(self.get_target_id())

    def get_sink(self) -> List[str]:
        return self.__sink

    def get_reactions_ids(self) -> List[str]:
        '''Returns the list of reaction IDs sorted by index within the pathway'''
        return [
            rxn_id for rxn_id in sorted(
                super().get_reactions_ids(),
                key=lambda x: self.get_reaction(x).get_info('idx_in_path')
            )
        ]

    def __get_infos(self, key: str) -> Dict[str, Dict]:
        species = {}
        for spe_id in self.get_list_of_species():
            species[spe_id] = {k:v for k,v in self.get_specie(spe_id).get_infos().items() if k.startswith(key)}
        reactions = {}
        for rxn_id in self.get_list_of_reactions():
            reactions[rxn_id] = {k:v for k,v in self.get_reaction(rxn_id).get_infos().items() if k.startswith(key)}
        pathway = {k:v for k,v in self.get_infos().items() if k.startswith(key)}
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

    def get_thermo_infos(self) -> Dict[str, Dict]:
        return self.__get_infos('thermo')

    def get_fba_infos(self) -> Dict[str, Dict]:
        return self.__get_infos('fba')

    def get_rpSBML_infos(self) -> Dict[str, Dict]:
        return self.__get_infos('rpSBML')


    ## WRITE METHODS
    def add_reaction(
        self,
        rxn: Reaction,
        rxn_id: str = None,
        target_id: str = None
    ) -> None:

        super().add_reaction(rxn, rxn_id)

        # TARGET
        if target_id != None:
            self.set_target_id(target_id)
            self.set_target_rxn_id(rxn.get_id())

    def set_sink(self, sink: List) -> None:
        self.__sink = deepcopy(sink)
    
    def set_target_id(self, id: str) -> None:
        self.__target_id = id

    def set_target_rxn_id(self, id: str) -> None:
        self.__rxn_target_id = id

    def rename_compound(self, id: str, new_id: str) -> None:
        # target
        for rxn in self.get_reactions():
            if id in rxn.get_species_ids():
                if id == self.get_target_id():
                    self.set_target_id(new_id)
                    self.set_target_rxn_id(rxn.get_id())

        super().rename_compound(id, new_id)

        # sink
        try:
            self.__sink[self.__sink.index(id)] = new_id
        except ValueError:
            pass


    ## MISC
