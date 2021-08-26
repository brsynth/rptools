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
from .rpObject import rpObject


class rpReaction(Reaction, rpObject):
    """A class to implement a chemical reaction
    enriched with both FBA and thermodynamics informations,
    and with RP informations (from RetroPath2 and rp2paths).
    """

    __default_fbc_lower = -10000
    __default_fbc_upper = 10000
    __default_fbc_units = 'mmol_per_gDW_per_hr'

    @staticmethod
    def get_default_fbc_units() -> str:
        return rpReaction.__default_fbc_units
    @staticmethod
    def get_default_fbc_lower() -> float:
        return rpReaction.__default_fbc_lower
    @staticmethod
    def get_default_fbc_upper() -> float:
        return rpReaction.__default_fbc_upper

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
        miriam: Dict = {},
        logger: Logger = getLogger(__name__)
    ):
        """Create a rpReaction object with default settings.

        Parameters
        ----------
        id: str
            ID of the reaction
        ec_numbers: Union[List[str], str], optional
            (list of) String(s) to define Enzyme Commission Number
            of the reaction.
        reactants: Dict[str, int], optional
            Stoichiometric dictionary of reactants species
        products: Dict[str, int], optional
            Stoichiometric dictionary of products species
        idx_in_path: int, optional
            Index of the reaction within the metabolic pathway
        lower_flux_bound: float, optional
            Lower flux bound value
        upper_flux_bound: float, optional
            Upper flux bound value
        flux_bound_units: str, optional
            Flux bounds units
        reversible: bool, optional
            Tells if the reaction is reversible or not
        miriam: Dict, optional
            Extended informations as cross references
        logger : Logger, optional
        """
        if 'ec-code' in miriam:
            ec_numbers = miriam['ec-code']
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
        self.set_rule_score('NaN')
        self.set_idx_in_path(idx_in_path)
        self.set_fbc(
            l_bound=lower_flux_bound,
            u_bound=upper_flux_bound,
            units=flux_bound_units
        )
        self.set_selenzy({})
        self.set_reversible(reversible)
        self.set_miriam(miriam)

    ## OUT METHODS
    # def __repr__(self):
    #     return f'Reaction {self.get_name()}'

    def _to_dict(
        self,
        specific: bool = False
    ) -> Dict:
        """Get attributes as a dictionary.

        Parameters
        ----------
        specific: bool, optional
            If set to True, the returned dictionary will not
            contain attributes inherited from Reaction class.
        """
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
        """Returns a dictionary which contains attributes
        only from rpReaction class excluding inherited ones."""
        return {
            'rp2_transfo_id': self.get_rp2_transfo_id(),
            'rule_id': self.get_rule_id(),
            'tmpl_rxn_id': self.get_tmpl_rxn_id(),
            'rule_score': self.get_rule_score(),
            'idx_in_path': self.get_idx_in_path(),
            'selenzy': self.get_selenzy()
            # 'fbc': deepcopy(self.get_fbc())
        }

    def __eq__(self, other) -> bool:
        """Returns the equality between two rpReaction objects."""
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
        """Get the ID of chemical transformation that
        the current reaction has been generated from."""
        return self.__rp2_transfo_id

    def get_rule_id(self) -> str:
        """Get the ID of reaction rule that
        the chemical transformation of the
        current reaction has been generated from."""
        return self.__rule_id

    def get_tmpl_rxn_id(self) -> str:
        """Set the ID of the template (original) reaction
        that the reaction rule of the reaction rule of the
        current reaction has been generated from."""
        return self.__tmpl_rxn_id

    def get_rule_score(self) -> float:
        """Get the score of the reaction rule (from RetroPath2)."""
        return self.__rule_score

    def get_idx_in_path(self) -> int:
        """Get the index of the reaction within the metabolic pathway."""
        return self.__idx_in_path
    
    # def get_fbc(self) -> float:
    #     return self.__fbc

    def get_fbc_units(self) -> str:
        """Get flux bounds constraints units."""
        return self.__fbc_units

    def get_fbc_lower(self) -> float:
        """Get flux lower bound value."""
        return self.__fbc_lower

    def get_fbc_upper(self) -> float:
        """Get flux upper bound value."""
        return self.__fbc_upper

    def reversible(self) -> bool:
        """Tells if the reaction is reversible or not."""
        return self.__reversible

    def get_miriam(self) -> Dict:
        """Get extended informations as cross references."""
        return self.__miriam

    def get_selenzy(self) -> Dict[str, float]:
        """Same as get_selenzy_scores()."""
        return self.get_selenzy_scores()

    def get_selenzy_scores(self) -> Dict[str, float]:
        """Get selenzyme scores."""
        return self.__selenzy

    def get_selenzy_score(self, id: str) -> float:
        """Get selenzyme score for a specific UniProtID.
        
        Parameters
        ----------
        id: str,
            UniProtID to get the score for.
        """
        return self.get_selenzy_scores().get(id, None)

    ## WRITE METHODS
    def set_rp2_transfo_id(self, id: str) -> None:
        """Set the ID of chemical transformation that
        the current reaction has been generated from.
        
        Parameters
        ----------
        id: str
        """
        self.__rp2_transfo_id = id

    def set_rule_id(self, id: str) -> None:
        """Set the ID of reaction rule that
        the chemical transformation of the
        current reaction has been generated from.
        
        Parameters
        ----------
        id: str
        """
        self.__rule_id = id

    def set_tmpl_rxn_id(self, id: str) -> None:
        """Set the ID of the template (original) reaction
        that the reaction rule of the reaction rule of the
        current reaction has been generated from.
        
        Parameters
        ----------
        id: str
        """
        self.__tmpl_rxn_id = id

    def set_rule_score(self, score: float) -> None:
        """Set the score of the reaction rule of
        the chemical transformation of
        the current reaction has been generated from.
        
        Parameters
        ----------
        score: float
        """
        self.__rule_score = score

    def set_idx_in_path(self, index: int) -> None:
        """Set the position of the reaction within
        the metabolic pathway.
        
        Parameters
        ----------
        index: int
        """
        self.__idx_in_path = index

    def set_fbc_lower(self, bound: float) -> None:
        """Set the lower flux bound.
        
        Parameters
        ----------
        bound: float
        """
        self.__fbc_lower = bound

    def set_fbc_upper(self, bound: float) -> None:
        """Set the upper flux bound.
        
        Parameters
        ----------
        bound: float
        """
        self.__fbc_upper = bound

    def set_fbc_units(self, units: str) -> None:
        """Set the flux bounds units.
        
        Parameters
        ----------
        units: str
        """
        self.__fbc_units = units

    def set_fbc(
        self,
        l_bound: float,
        u_bound: float,
        units: str = __default_fbc_units
    ) -> None:
        """Set flux bound constraints.
        
        Parameters
        ----------
        l_bound: float
            Lower flux bound
        u_bound: float
            Upper flux bound
        units: str, optional
            Flux bounds constraints
        """
        # self.__fbc = {}
        self.set_fbc_lower(l_bound)
        self.set_fbc_upper(u_bound)
        self.set_fbc_units(units)

    def set_reversible(self, reversible: bool) -> None:
        """Set the reversibility of the reaction.
        
        Parameters
        ----------
        reversible: bool
        """
        self.__reversible = reversible

    def set_selenzy(self, scores: Dict[str, float]) -> None:
        """Same as set_selenzy_scores()."""
        self.set_selenzy_scores(scores)

    def set_selenzy_scores(self, scores: Dict[str, float]) -> None:
        """Set the selenzyme scores.
        
        Parameters
        ----------
        scores: Dict[str, float]
            A dictionary which links each UniProtID to a score
        """
        self.__selenzy = {}
        for uniprot_id, score in scores.items():
            self.add_selenzy_score(uniprot_id, score)

    def add_selenzy_score(self, id: str, score: float) -> None:
        """Add a selenzyme score.
        
        Parameters
        ----------
        id: str
            UniProtID
        score: float
            Score to link to the UniProtID
        """
        self.__selenzy[id] = score

    def set_miriam(self, miriam: Dict) -> None:
        """Set extended informations as cross references.
        
        Parameters
        ----------
        miriam: Dict
            A dictionary which contains extended informations
        """
        self.__miriam = deepcopy(miriam)

    def add_miriam(self, key: str, infos: TypeVar) -> None:
        """Add an extended information.
        
        Parameters
        ----------
        key: str
            Label of the information to add
        infos: TypeVar
            Extended informations to add
        """
        self.__miriam[key] = deepcopy(infos)
