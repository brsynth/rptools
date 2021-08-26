"""
Created on May 28 2021

@author: Joan Hérisson
"""

from unittest import TestCase
from copy import deepcopy
from rptools.rplibs import rpReaction


class Test_rpReaction(TestCase):

    def setUp(self):
        self.reactants = {
            "CMPD_0000000010": 1,
            "MNXM1": 1
        }
        self.products = {
            "CMPD_0000000003": 1,
            "MNXM13": 1
        }
        self.ec_numbers = [
            "4.1.1.63"
        ]
        self.id = "rxn"
        self.miriam = {'ec-code': ['4.1.1.63']}
        self.rxn = rpReaction(
            id=self.id,
            miriam=self.miriam,
            reactants=self.reactants,
            products=self.products
        )
        self.rp2_transfo_id = 'TRS_0_0_0'
        self.rule_id = 'RR-02-a0cc0be463ff412f-16-F'
        self.tmpl_rxn_id = 'MNXR96458'
        self.rule_score = 0.5982208769718989
        self.selenzy = {
            'UniProtID_1': 65.65,
            'UniProtID_2': 77.77,
        }
        self.inherited_dict = {
            'id': self.id,
            'reactants': self.reactants,
            'products': self.products,
            'ec_numbers': self.ec_numbers,
        }
        self.specific_dict = {
            'idx_in_path': -1,
            'rp2_transfo_id': self.rp2_transfo_id,
            'rule_id': self.rule_id,
            'tmpl_rxn_id': self.tmpl_rxn_id,
            'rule_score': self.rule_score,
            'selenzy': self.selenzy
        }

    ## READ METHODS
    def test_get_rp2_transfo_id(self):
        self.rxn.set_rp2_transfo_id(self.rp2_transfo_id)
        self.assertEqual(
            self.rxn.get_rp2_transfo_id(),
            self.rp2_transfo_id
        )

    def test_get_rule_id(self):
        self.rxn.set_rule_id(self.rule_id)
        self.assertEqual(
            self.rxn.get_rule_id(),
            self.rule_id
        )

    def test_get_tmpl_rxn_id(self):
        self.rxn.set_tmpl_rxn_id(self.tmpl_rxn_id)
        self.assertEqual(
            self.rxn.get_tmpl_rxn_id(),
            self.tmpl_rxn_id
        )

    def test_get_rule_score(self):
        self.rxn.set_rule_score(self.rule_score)
        self.assertEqual(
            self.rxn.get_rule_score(),
            self.rule_score
        )

    def test__to_dict(self):
        self.rxn.set_rp2_transfo_id(self.rp2_transfo_id)
        self.rxn.set_rule_id(self.rule_id)
        self.rxn.set_tmpl_rxn_id(self.tmpl_rxn_id)
        self.rxn.set_rule_score(self.rule_score)
        self.rxn.set_selenzy(self.selenzy)
        self.assertDictEqual(
            self.rxn._to_dict(),
            {
                **self.inherited_dict,
                **self.specific_dict
            }
        )
        self.assertDictEqual(
            self.rxn._to_dict(specific=True),
            self.specific_dict
        )
    
    def test___eq__(self):
        rxn = deepcopy(self.rxn)
        self.assertEqual(
            self.rxn,
            rxn
        )
        rxn.set_idx_in_path(2)
        self.assertEqual(
            self.rxn,
            rxn
        )
        rxn.add_reactant('reac', 2)
        # objects are not equal
        self.assertNotEqual(
            self.rxn,
            rxn
        )
        # objects are not the same type
        self.assertNotEqual(
            self.rxn,
            'rxn'
        )

    def test_get_default_fbc(self):
        fbc = {
            'units': rpReaction.get_default_fbc_units(),
            'lower': rpReaction.get_default_fbc_lower(),
            'upper': rpReaction.get_default_fbc_upper()
        }
        for id, val in fbc.items():
            with self.subTest(f'Testing get/set fbc_{id}()', id=id, val=val):
                getattr(self.rxn, f'set_fbc_{id}')(val)
                self.assertEqual(
                    getattr(self.rxn, f'get_fbc_{id}')(),
                    val
                )

    def test_reversible(self):
        self.rxn.set_reversible(False)
        self.assertFalse(self.rxn.reversible())
        self.rxn.set_reversible(True)
        self.assertTrue(self.rxn.reversible())

    def test_miriam(self):
        self.assertDictEqual(
            self.rxn.get_miriam(),
            self.miriam
        )

    def test_selenzy(self):
        self.rxn.set_selenzy(self.selenzy)
        self.assertDictEqual(
            self.rxn.get_selenzy(),
            self.selenzy
        )
        id = 'UniProtID_1'
        self.assertEqual(
            self.rxn.get_selenzy_score(id),
            self.selenzy[id]
        )

    def test_add_miriam(self):
        db = 'bigg'
        xref = 'bigg_ID'
        self.rxn.add_miriam(db, xref)
        self.assertDictEqual(
            self.rxn.get_miriam(),
            {
                **self.miriam,
                **{
                    db: xref
                }
            }
        )