"""A class to represent a chemical reaction."""
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
    Union,
    TypeVar
)
from copy import deepcopy
from logging import (
    Logger,
    getLogger
)
from math import isnan
from chemlite import Reaction
from rptools.rplibs.rpObject import rpObject


class rpReaction(Reaction, rpObject):

    __default_fbc_lower = -10000
    __default_fbc_upper = 10000
    __default_fbc_units = 'mmol_per_gDW_per_hr'

    def __init__(
        self,
        id: str,
        ec_numbers: Union[List[str], str] = [],
        reactants: Dict[str, int] = {},
        products: Dict[str, int] = {},
        idx_in_path: int = -1,
        lower_flux_bound: float = __default_fbc_lower,
        upper_flux_bound: float = __default_fbc_upper,
        flux_bound_units: str = __default_fbc_units,
        reversible: bool = False,
        logger: Logger = getLogger(__name__)
    ):
        Reaction.__init__(
            self,
            id=id,
            ec_numbers=ec_numbers,
            reactants=reactants,
            products=products,
            logger=logger
        )
        rpObject.__init__(self)
        self.set_rp2_transfo_id(None)
        self.set_rule_id(None)
        self.set_tmpl_rxn_id(None)
        self.set_rule_score(float('nan'))
        self.set_idx_in_path(idx_in_path)
        self.set_fbc(
            l_value=lower_flux_bound,
            u_value=upper_flux_bound,
            units=flux_bound_units
        )
        self.set_reversible(reversible)

    ## OUT METHODS
    # def __repr__(self):
    #     return f'Reaction {self.get_name()}'

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
                **Reaction._to_dict(self),
                **rpObject._to_dict(self),
                **self.__to_dict()
            }

    def __to_dict(self) -> Dict:
        return {
            'rp2_transfo_id': self.get_rp2_transfo_id(),
            'rule_id': self.get_rule_id(),
            'tmpl_rxn_id': self.get_tmpl_rxn_id(),
            'rule_score': self.get_rule_score(),
            'idx_in_path': self.get_idx_in_path(),
            # 'fbc': deepcopy(self.get_fbc())
        }

    def __eq__(self, other) -> bool:
        if not isinstance(self, other.__class__):
            return False
        # Compare with some specific keys
        return all(
            self._to_dict().get(key) == other._to_dict().get(key)
            for key in [
                'ec_numbers',
                'reactants',
                'products'
            ]
        )

    ## READ METHODS
    def get_rp2_transfo_id(self) -> str:
        return self.__rp2_transfo_id

    def get_rule_id(self) -> str:
        return self.__rule_id

    def get_tmpl_rxn_id(self) -> str:
        return self.__tmpl_rxn_id

    def get_rule_score(self) -> float:
        return self.__rule_score

    def get_idx_in_path(self) -> int:
        return self.__idx_in_path
    
    def get_fbc(self) -> float:
        return self.__fbc

    def get_fbc_units(self) -> str:
        return self.__fbc_units

    def reversible(self) -> bool:
        return self.__reversible

    @staticmethod
    def get_default_fbc_units() -> str:
        return rpReaction.__default_fbc_units

    def get_fbc_lower(self) -> float:
        return self.__fbc_lower

    def get_fbc_upper(self) -> float:
        return self.__fbc_upper

    @staticmethod
    def get_default_fbc_lower() -> float:
        return rpReaction.__default_fbc_lower

    @staticmethod
    def get_default_fbc_upper() -> float:
        return rpReaction.__default_fbc_upper

    ## WRITE METHODS
    def set_rp2_transfo_id(self, transfo_id: str) -> None:
        self.__rp2_transfo_id = transfo_id

    def set_rule_id(self, rule_id: str) -> None:
        self.__rule_id = rule_id

    def set_tmpl_rxn_id(self, tmpl_rxn_id: str) -> None:
        self.__tmpl_rxn_id = tmpl_rxn_id

    def set_rule_score(self, rule_score: str) -> None:
        self.__rule_score = rule_score

    def set_idx_in_path(self, idx_in_path: str) -> None:
        self.__idx_in_path = idx_in_path

    def set_fbc_lower(self, value: float) -> None:
        self.__fbc_lower = value

    def set_fbc_upper(self, value: float) -> None:
        self.__fbc_upper = value

    def set_fbc_units(self, value: str) -> None:
        self.__fbc_units = value

    def set_fbc(
        self,
        l_value: float,
        u_value: float,
        units: str = __default_fbc_units
    ) -> None:
        self.__fbc = {}
        self.set_fbc_lower(l_value)
        self.set_fbc_upper(u_value)
        self.set_fbc_units(units)

    def set_reversible(self, value: bool) -> None:
        self.__reversible = value