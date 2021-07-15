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


class rpObject():

    __thermo_prefix = 'thermo_'
    __fba_prefix = 'fba_'
    __dG0_prime_str = 'dG0_prime'
    __dGm_prime_str = 'dGm_prime'
    __dG_prime_str = 'dG_prime'
    __dG_str = 'dG'
    __thermo_dG0_prime_str = __thermo_prefix+__dG0_prime_str
    __thermo_dGm_prime_str = __thermo_prefix+__dGm_prime_str
    __thermo_dG_prime_str = __thermo_prefix+__dG_prime_str
    __thermo_dG_str = __thermo_prefix+__dG_str
    __biomass_str = 'biomass'
    __fraction_str = 'fraction'
    __fba_biomass_str = __fba_prefix+__biomass_str
    __fba_fraction_str = __fba_prefix+__fraction_str
    __fba_str = 'fba'
    __pfba_str = 'pfba'
    __fba_fba_str = __fba_prefix+__fba_str
    __fba_pfba_str = __fba_prefix+__pfba_str

    @staticmethod
    def get_fba_biomass_str() -> str: return rpObject.__fba_biomass_str
    @staticmethod
    def get_fba_fraction_str() -> str: return rpObject.__fba_fraction_str

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
            # **super()._to_dict(),
            **self.__to_dict()
        }

    def __to_dict(self) -> Dict:
        return {
            **self.get_fba(),
            **self.get_thermo()
        }

    def __eq__(self, other) -> bool:
        if isinstance(self, other.__class__):
            return self._to_dict() == other._to_dict()
        return False

    ## READ METHODS
    ### THERMO ###
    def get_thermo(self) -> TypeVar:
        return self.__thermo

    def get_thermo_info(self, key: str) -> TypeVar:
        return self.__thermo.get('thermo_'+key, None)

    def get_thermo_dG0_prime(self) -> TypeVar:
        return self.__thermo.get(rpObject.__thermo_dG0_prime_str, None)

    def get_thermo_dGm_prime(self) -> TypeVar:
        return self.__thermo.get(rpObject.__thermo_dGm_prime_str, None)

    def get_thermo_dG_prime(self) -> TypeVar:
        return self.__thermo.get(rpObject.__thermo_dG_prime_str, None)

    def get_thermo_dG(self) -> TypeVar:
        return self.__thermo.get(rpObject.__thermo_dG_str, None)

    def get_fba_info(self, key: str) -> TypeVar:
        return self.__fba.get(rpObject.__fba_prefix+key, None)

    def get_fba(self) -> TypeVar:
        return self.__fba

    def __get_fba_(self, key: str) -> TypeVar:
        try:  # Returns exact key
            return self.__fba[key]
        except:  # Look for key starts with
            try:
                # Returns the value for key starts with {key}
                fba = [v for k,v in self.__fba.items() if k.startswith(key)]
                return fba[0] if len(fba) == 1 else fba
            except StopIteration:
                return None

    def get_fba_biomass(self) -> TypeVar:
        return self.__get_fba_(rpObject.__fba_biomass_str)

    def get_fba_fraction(self) -> TypeVar:
        return self.__get_fba_(rpObject.__fba_fraction_str)

    def get_fba_fba(self) -> TypeVar:
        return self.__get_fba_(rpObject.__fba_str)

    def get_fba_pfba(self) -> TypeVar:
        return self.__get_fba_(rpObject.__pfba_str)

    ## WRITE METHODS
    ### THERMO ###
    def __set_thermo(self, thermo: TypeVar) -> None:
        self.__thermo = deepcopy(thermo)

    def set_thermo_info(self, key: str, value: TypeVar) -> None:
        self.__thermo[rpObject.__thermo_prefix+key] = deepcopy(value)

    ### FBA ###
    def __set_fba(self, fba: TypeVar) -> None:
        self.__fba = deepcopy(fba)

    def set_fba_info(self, key: str, value: TypeVar) -> None:
        self.__fba[rpObject.__fba_prefix+key] = deepcopy(value)

    def set_thermo_dG0_prime(self, value: float) -> None:
        self.set_thermo_info(rpObject.__dG0_prime_str, value)

    def set_thermo_dGm_prime(self, value: float) -> None:
        self.set_thermo_info(rpObject.__dGm_prime_str, value)

    def set_thermo_dG_prime(self, value: float) -> None:
        self.set_thermo_info(rpObject.__dG_prime_str, value)

    def set_thermo_dG(self, value: float) -> None:
        self.set_thermo_info(rpObject.__dG_str, value)

    def set_fba_biomass(self, value: float) -> None:
        self.set_fba_info(rpObject.__biomass_str, value)

    def set_fba_fraction(self, value: float) -> None:
        self.set_fba_info(rpObject.__fraction_str, value)

    def set_fba_fba(self, value: float) -> None:
        self.set_fba_info(rpObject.__fba_str, value)

    def set_fba_pfba(self, value: float) -> None:
        self.set_fba_info(rpObject.__pfba_str, value)
