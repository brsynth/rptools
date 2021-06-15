"""A class to represent a chemical species."""
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
    TypeVar
)
from logging import (
    Logger,
    getLogger
)
from copy import deepcopy


def __get_dict_key(
    dict: Dict,
    key: str,
    logger: Logger = getLogger(__name__)
) -> TypeVar:
    try:
        return dict[key]
    except KeyError:
        logger.warning(f'There is no {key} key in data')
        return None


class rpObject():

    def __init__(
        self,
        logger: Logger = getLogger(__name__)
    ):
        self.__set_fba({})
        self.__set_thermo({})

    ## OUT METHODS
    # def __repr__(self):
    #     return f'Compound {self.get_id()}'

    def _to_dict(self) -> Dict:
        return {
            **super()._to_dict(),
            **self._infos_to_dict()
        }

    def _infos_to_dict(self) -> Dict:
        return {
            **self.get_fba(),
            **self.get_thermo()
        }

    def __eq__(self, other) -> bool:
        if isinstance(self, other.__class__):
            return self.__to_dict() == other.__to_dict()
        return False

    ## READ METHODS
    ### THERMO ###
    def get_thermo(self) -> Dict:
        return self.__thermo

    def get_thermo_info(self, key: str) -> Dict:
        return __get_dict_key(self.__thermo, 'thermo_'+key)

    def get_fba_info(self, key: str) -> Dict:
        return __get_dict_key(self.__fba, 'fba_'+key)

    def get_fba(self) -> Dict:
        return self.__fba

    ## WRITE METHODS
    ### THERMO ###
    def __set_thermo(self, thermo: Dict) -> None:
        self.__thermo = deepcopy(thermo)

    def set_thermo_info(self, key: str, value: TypeVar) -> None:
        prefix = 'thermo_'
        if key.startswith(prefix):
            prefix = ''
        self.__thermo[prefix+key] = deepcopy(value)

    ### FBA ###
    def __set_fba(self, fba: Dict) -> None:
        self.__fba = deepcopy(fba)

    def set_fba_info(self, key: str, value: TypeVar) -> None:
        prefix = 'fba_'
        if key.startswith(prefix):
            prefix = ''
        self.__fba[prefix+key] = deepcopy(value)

