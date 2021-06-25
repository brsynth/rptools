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

    thermo_str = 'thermo_standard_dg_formation'
    fba_biomass_str = rpObject.fba_biomass_str+'_shadow_price'
    fba_fraction_str = rpObject.fba_fraction_str+'_shadow_price'

    def __init__(
        self,
        id: str,
        smiles: str = '',
        inchi: str = '',
        inchikey: str = '',
        formula: str = '',
        name: str = '',
        logger: Logger = getLogger(__name__)
    ):
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
                **Compound._to_dict(self),
                **rpObject._to_dict(self),
                **self.__to_dict()
            }

    def __to_dict(self) -> Dict:
        return {
        }

    def get_thermo_standard_dg_formation(self) -> TypeVar:
        return self.get_thermo_info(rpCompound.thermo_str)

    def get_fba_biomass_shadow_price(self) -> TypeVar:
        return self.get_fba_info(rpCompound.fba_biomass_str)

    def get_fba_fraction_shadow_price(self) -> TypeVar:
        return self.get_fba_info(rpCompound.fba_fraction_str)

    def set_thermo_standard_dg_formation(self, value: float) -> None:
        self.set_thermo_info(rpCompound.thermo_str, value)

    def set_fba_biomass_shadow_price(self, value: float) -> None:
        self.set_fba_info(rpCompound.fba_biomass_shadow_price, value)

    def set_fba_fraction_shadow_price(self, value: float) -> None:
        self.set_fba_info(rpCompound.fba_fraction_shadow_price, value)