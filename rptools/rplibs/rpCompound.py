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
from chemlite import Compound
from rptools.rplibs.rpObject import rpObject


class rpCompound(Compound, rpObject):
    """A class to implement a chemical species
    enriched with FBA and thermodynamics informations.
    """

    __thermo_str = 'standard_dg_formation'
    __fba_str = 'shadow_price'

    def __init__(
        self,
        id: str,
        smiles: str = '',
        inchi: str = '',
        inchikey: str = '',
        formula: str = '',
        name: str = '',
        compartment_id: str = 'c',
        logger: Logger = getLogger(__name__)
    ):
        """Create a rpCompound object with default settings.

        Parameters
        ----------
        id: str
        smiles: str, optional
        inchi: str, optional
        inchikey: str, optional
        formula: str, optional
        name: str, optional
        compartment_id: str, optional
        logger : Logger, optional
        """
        Compound.__init__(
            self,
            id=id,
            smiles=smiles,
            inchi=inchi,
            inchikey=inchikey,
            formula=formula,
            name=name,
            logger=logger
        )
        rpObject.__init__(self)
        self.set_compartment(compartment_id)

    def _to_dict(
        self,
        full: bool = True
    ) -> Dict:
        """Get attributes as a dictionary.

        Parameters
        ----------
        full: bool, optional
            If set to False, the returned dictionary will not
            contain attributes inherited from Compound class
            (default: True).
        """
        if full:
            return {
                **Compound._to_dict(self),
                **rpObject._to_dict(self),
                **self.__to_dict()
            }
        else:
            return {
                **self.__to_dict(),
                **rpObject._to_dict(self)
            }

    def __to_dict(self) -> Dict:
        """Returns a dictionary which contains attributes
        only from rpCompound class excluding inherited ones."""
        return {
            # 'compartment': self.get_compartment()
        }

    # def __eq__(self, other) -> bool:
    #     """Returns the equality between two rpCompound objects."""
    #     return super(Compound, self).__eq__(other)

    def get_thermo_standard_dg_formation(self) -> TypeVar:
        """Get thermodynamics dG formation cost."""
        return self.get_thermo_info(rpCompound.__thermo_str)

    def get_fba_biomass_shadow_price(self) -> TypeVar:
        """Get flux shadow price during biomass production."""
        return self.get_fba_info(f'biomass_{rpCompound.__fba_str}')

    def get_fba_fraction_shadow_price(self) -> TypeVar:
        """Get flux shadow price during fraction of reaction analysis."""
        return self.get_fba_info(f'fraction_{rpCompound.__fba_str}')

    def get_fba_fba_shadow_price(self) -> TypeVar:
        """Get flux shadow price during balance analysis."""
        return self.get_fba_info(f'fba_{rpCompound.__fba_str}')

    def get_fba_pfba_shadow_price(self) -> TypeVar:
        """Get flux shadow price during parcimonious balance analysis."""
        return self.get_fba_info(f'pfba_{rpCompound.__fba_str}')

    def get_compartment(self) -> str:
        """Get compound compartment ID."""
        return self.__compartment

    def set_thermo_standard_dg_formation(self, value: float) -> None:
        """Set dG formation cost.
        
        Parameters
        ----------
        value: float
        """
        self.add_thermo_info(rpCompound.__thermo_str, value)

    def set_fba_biomass_shadow_price(self, value: float) -> None:
        """Set flux shadow price during biomass production.
        
        Parameters
        ----------
        value: float
        """
        self.add_fba_info(f'biomass_{rpCompound.__fba_str}', value)

    def set_fba_fraction_shadow_price(self, value: float) -> None:
        """Set flux shadow price during fraction of reaction analysis.
        
        Parameters
        ----------
        value: float
        """
        self.add_fba_info(f'fraction_{rpCompound.__fba_str}', value)

    def set_fba_fba_shadow_price(self, value: float) -> None:
        """Set flux shadow price during balance analysis..
        
        Parameters
        ----------
        value: float
        """
        self.add_fba_info(f'fba_{rpCompound.__fba_str}', value)

    def set_fba_pfba_shadow_price(self, value: float) -> None:
        """Set flux shadow price during parcimonious balance analysis.
        
        Parameters
        ----------
        value: float
        """
        self.add_fba_info(f'pfba_{rpCompound.__fba_str}', value)

    def set_compartment(self, compartment: str) -> None:
        """Set compartment ID of the compound.
        
        Parameters
        ----------
        compartment: str
        """
        self.__compartment = compartment