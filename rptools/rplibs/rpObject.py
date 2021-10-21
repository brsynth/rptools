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
    """A common class for classes which share same attributes.
    """

    __sep = '_'
    __thermo_prefix = 'thermo'
    __fba_prefix = 'fba'
    __dG0_prime_str = 'dG0_prime'
    __dGm_prime_str = 'dGm_prime'
    __dG_prime_str = 'dG_prime'
    __dG_str = 'dG'
    __biomass_str = 'biomass'
    __fraction_str = 'fraction'
    __fba_str = 'fba'
    __pfba_str = 'pfba'

    @staticmethod
    def get_sep() -> str: return rpObject.__sep
    @staticmethod
    def get_fba_prefix() -> str: return rpObject.__fba_prefix
    @staticmethod
    def get_thermo_prefix() -> str: return rpObject.__thermo_prefix

    def __init__(
        self,
        logger: Logger = getLogger(__name__)
    ):
        self.__set_fba({})
        self.__set_thermo({})
        """Create a rpObject object with default settings.

        Parameters
        ----------
        logger : Logger, optional
        """


    ## OUT METHODS
    # def __repr__(self):
    #     return f'Compound {self.get_id()}'

    def _to_dict(self) -> Dict:
        """Get attributes as a dictionary.
        Implemented for compatibility with inherited classes."""
        fba_infos = {
            rpObject.get_sep().join([
                rpObject.get_fba_prefix(),
                k
            ]):v for k,v in self.get_fba().items()
        }
        thermo_infos = {
            rpObject.get_sep().join([
                rpObject.get_thermo_prefix(),
                k
            ]):v for k,v in self.get_thermo().items()
        }
        return {
            **fba_infos,
            **thermo_infos
        }

    def __eq__(self, other) -> bool:
        """Returns the equality between two rpObject objects."""
        if isinstance(self, other.__class__):
            return self._to_dict() == other._to_dict()
        return False

    ## READ METHODS
    ### THERMO ###
    def get_thermo_infos(self) -> TypeVar:
        """Same as get_thermo()"""
        return self.get_thermo()

    def get_thermo(self) -> TypeVar:
        """Get thermodynamics values."""
        return self.__thermo

    def get_thermo_info(self, key: str) -> TypeVar:
        """Get a specific thermodynamics value.
        
        Parameters
        ----------
        key: str
            A key to return the thermodynamics value of
        """
        return self.get_thermo().get(key, None)

    def get_thermo_dG0_prime(self) -> TypeVar:
        """Get thermodynamics dG0 prime values."""
        return self.get_thermo().get(rpObject.__dG0_prime_str, None)

    def get_thermo_dGm_prime(self) -> TypeVar:
        """Get thermodynamics dGm prime values."""
        return self.get_thermo().get(rpObject.__dGm_prime_str, None)

    def get_thermo_dG_prime(self) -> TypeVar:
        """Get thermodynamics dG prime values."""
        return self.get_thermo().get(rpObject.__dG_prime_str, None)

    def get_thermo_dG(self) -> TypeVar:
        """Get thermodynamics dG prime values."""
        return self.get_thermo().get(rpObject.__dG_str, None)

    def get_fba_info(self, key: str) -> TypeVar:
        """Get a specific flux balance analysis (FBA) value.

        Parameters
        ----------
        key: str
            A key to return the FBA value of
        """
        return self.get_fba().get(key, None)

    # def get_fba_biomass_info(self) -> TypeVar:
    #     """Get flux balance analysis values concerning the production of biomass."""
    #     return self.get_fba_info(rpObject.__biomass_str)

    # def get_fba_fraction_info(self) -> TypeVar:
    #     """Get flux balance values concerning the fraction of reaction analysis."""
    #     return self.get_fba_info(rpObject.__fraction_str)

    # def get_fba_fba_info(self) -> TypeVar:
    #     """Get flux balance values concerning the fba analysis."""
    #     return self.get_fba_info(rpObject.__fba_str)

    # def get_fba_pfba_info(self) -> TypeVar:
    #     """Get flux balance values concerning the parcimonious fba analysis."""
    #     return self.get_fba_info(rpObject.__pfba_str)

    def get_fba_infos(self) -> TypeVar:
        """Same as get_fba()"""
        return self.get_fba()

    def get_fba(self) -> TypeVar:
        """Get flux balance values."""
        return self.__fba

    # def __get_fba_(self, key: str) -> TypeVar:
    #     """Get flux balance values."""
    #     try:  # Returns exact key
    #         return self.__fba[key]
    #     except:  # Look for key starts with
    #         try:
    #             # Returns the value for key starts with {key}
    #             fba = [v for k,v in self.__fba.items() if k.startswith(key)]
    #             return fba[0] if len(fba) == 1 else fba
    #         except StopIteration:
    #             return None

    def get_fba_biomass(self) -> TypeVar:
        """Get flux balance analysis values concerning the production of biomass."""
        return self.get_fba_info(rpObject.__biomass_str)

    def get_fba_fraction(self) -> TypeVar:
        """Get flux balance values concerning the fraction of reaction analysis."""
        return self.get_fba_info(rpObject.__fraction_str)

    def get_fba_fba(self) -> TypeVar:
        """Get flux balance values concerning the fba analysis."""
        return self.get_fba_info(rpObject.__fba_str)

    def get_fba_pfba(self) -> TypeVar:
        """Get flux balance values concerning the parcimonious fba analysis."""
        return self.get_fba_info(rpObject.__pfba_str)

    ## WRITE METHODS
    ### THERMO ###
    def __set_thermo(self, thermo: TypeVar) -> None:
        """Set thermodynamics values.
        
        Parameters
        ----------
        thermo: TypeVar
            Thermodynamics informations to store.
        """
        self.__thermo = deepcopy(thermo)

    def set_thermo_infos(
        self,
        infos: Dict[str, TypeVar]
    ) -> None:
        """Set the thermodynamics infos.
        
        Parameters
        ----------
        infos: Dict
            A dictionary which contains thermodynamics infos
        """
        self.__thermo = {}
        for key, value in infos:
            self.add_thermo_info(key, value)

    def add_thermo_info(
        self,
        key: str,
        value: TypeVar
    ) -> None:
        """Set a specific thermodynamics value
        within a dictionary object.
        
        Parameters
        ----------
        key: str
            Key to retrieve value later on.
        value: TypeVar
            Thermodynamics information to store.
        """
        self.__thermo[key] = deepcopy(value)

    ### FBA ###
    def __set_fba(self, fba: TypeVar) -> None:
        """Set flux balance analysis (FBA) values.
        
        Parameters
        ----------
        fba: TypeVar
            FBA informations to store.
        """
        self.__fba = deepcopy(fba)

    def set_fba_infos(
        self,
        infos: Dict[str, TypeVar]
    ) -> None:
        """Set the FBA infos.
        
        Parameters
        ----------
        infos: Dict
            A dictionary which contains FBA infos
        """
        self.__fba = {}
        for key, value in infos:
            self.add_fba_info(key, value)

    def add_fba_info(
        self,
        key: str,
        value: TypeVar
    ) -> None:
        """Set a specific flux balance analysis (FBA)
        value within a dictionary object.
        
        Parameters
        ----------
        key: str
            Key to retrieve value later on.
        value: TypeVar
            FBA information to store.
        """
        self.__fba[key] = deepcopy(value)

    def set_thermo_dG0_prime(self, value: float) -> None:
        """Set dG0 prime thermodynamics value
        within a dictionary object.
        
        Parameters
        ----------
        value: TypeVar
            dG0 prime value
        """
        self.add_thermo_info(rpObject.__dG0_prime_str, value)

    def set_thermo_dGm_prime(self, value: float) -> None:
        """Set dGm prime thermodynamics value
        within a dictionary object.
        
        Parameters
        ----------
        value: TypeVar
            dGm prime value
        """
        self.add_thermo_info(rpObject.__dGm_prime_str, value)

    def set_thermo_dG_prime(self, value: float) -> None:
        """Set dG prime thermodynamics value
        within a dictionary object.
        
        Parameters
        ----------
        value: TypeVar
            dG prime value
        """
        self.add_thermo_info(rpObject.__dG_prime_str, value)

    def set_thermo_dG(self, value: float) -> None:
        """Set dG thermodynamics value
        within a dictionary object.
        
        Parameters
        ----------
        value: TypeVar
            dG value
        """
        self.add_thermo_info(rpObject.__dG_str, value)

    def set_fba_biomass(self, value: float) -> None:
        """Set flux balance analysis value for
        biomass reaction.
        
        Parameters
        ----------
        value: float
            Flux value
        """
        self.add_fba_info(rpObject.__biomass_str, value)

    def set_fba_fraction(self, value: float) -> None:
        """Set flux balance value for
        fraction of reaction analysis.
        
        Parameters
        ----------
        value: float
            Flux value
        """
        self.add_fba_info(rpObject.__fraction_str, value)

    def set_fba_fba(self, value: float) -> None:
        """Set flux balance value for
        fba analysis.
        
        Parameters
        ----------
        value: float
            Flux value
        """
        self.add_fba_info(rpObject.__fba_str, value)

    def set_fba_pfba(self, value: float) -> None:
        """Set flux balance value for
        parcimonious fba analysis.
        
        Parameters
        ----------
        value: float
            Flux value
        """
        self.add_fba_info(rpObject.__pfba_str, value)
