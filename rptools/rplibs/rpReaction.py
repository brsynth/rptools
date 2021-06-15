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
    Union
)
from logging import (
    Logger,
    getLogger
)
from copy import deepcopy
from chemlite import Reaction


class rpReaction(Reaction):

    def __init__(
        self,
        id: str,
        ec_numbers: Union[List[str], str] = [],
        reactants: Dict[str, int] = {},
        products: Dict[str, int] = {},
        logger: Logger = getLogger(__name__)
    ):
        super().__init__(
            id=id,
            ec_numbers=ec_numbers,
            reactants=reactants,
            products=products,
            logger=logger
        )
        self.set_rp2_transfo_id(None)
        self.set_rule_id(None)
        self.set_tmpl_rxn_id(None)
        self.set_rule_score(float('nan'))

    ## OUT METHODS
    # def __repr__(self):
    #     return f'Reaction {self.get_name()}'

    def _to_dict(self) -> Dict:
        return {
            **super()._to_dict(),
            **self._infos_to_dict()
        }

    def _infos_to_dict(self) -> Dict:
        return {
            'rp2_transfo_id': self.get_rp2_transfo_id(),
            'rule_id': self.get_rule_id(),
            'tmpl_rxn_id': self.get_tmpl_rxn_id(),
            'rule_score': self.get_rule_score()
        }

    def __eq__(self, other) -> bool:
        if isinstance(self, other.__class__):
            return self._to_dict() == other._to_dict()
        return False

    ## READ METHODS
    def get_rp2_transfo_id(self) -> str:
        return self.__rp2_transfo_id

    def get_rule_id(self) -> str:
        return self.__rule_id

    def get_tmpl_rxn_id(self) -> str:
        return self.__tmpl_rxn_id

    def get_rule_score(self) -> float:
        return self.__rule_score

    ## WRITE METHODS
    def set_rp2_transfo_id(self, transfo_id: str) -> None:
        self.__rp2_transfo_id = transfo_id

    def set_rule_id(self, rule_id: str) -> None:
        self.__rule_id = rule_id

    def set_tmpl_rxn_id(self, tmpl_rxn_id: str) -> None:
        self.__tmpl_rxn_id = tmpl_rxn_id

    def set_rule_score(self, rule_score: str) -> None:
        self.__rule_score = rule_score
