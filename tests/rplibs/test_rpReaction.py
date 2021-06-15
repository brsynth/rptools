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
        self.rxn = rpReaction(
            id=self.id,
            ec_numbers=self.ec_numbers,
            reactants=self.reactants,
            products=self.products
        )
        self.rp2_transfo_id = 'TRS_0_0_0'
        self.rule_id = 'RR-02-a0cc0be463ff412f-16-F'
        self.tmpl_rxn_id = 'MNXR96458'
        self.idx_in_path = 1
        self.rule_score = 0.5982208769718989

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

    def test_get_idx_in_path(self):
        self.rxn.set_idx_in_path(self.idx_in_path)
        self.assertEqual(
            self.rxn.get_idx_in_path(),
            self.idx_in_path
        )

    def test_infos_to_dict(self):
        self.rxn.set_rp2_transfo_id(self.rp2_transfo_id)
        self.rxn.set_rule_id(self.rule_id)
        self.rxn.set_tmpl_rxn_id(self.tmpl_rxn_id)
        self.rxn.set_rule_score(self.rule_score)
        self.rxn.set_idx_in_path(self.idx_in_path)
        print(self.rxn._to_dict())
        self.assertDictEqual(
            self.rxn._to_dict(),
            {
                **{
                    'id': self.id,
                    'reactants': self.reactants,
                    'products': self.products,
                    'ec_numbers': self.ec_numbers,
                    'infos': {}
                },
                **{
                    'rp2_transfo_id': self.rp2_transfo_id,
                    'rule_id': self.rule_id,
                    'tmpl_rxn_id': self.tmpl_rxn_id,
                    'rule_score': self.rule_score,
                    'idx_in_path': self.idx_in_path
                }
            }
        )
    
    def test___eq__(self):
        rxn = deepcopy(self.rxn)
        self.assertEqual(
            self.rxn,
            rxn
        )
        rxn.add_info('test_info', 'test_data')
        self.assertNotEqual(
            self.rxn,
            rxn
        )
        self.assertNotEqual(
            self.rxn,
            'rxn'
        )
