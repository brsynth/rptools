"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

import logging
from unittest             import TestCase
from tempfile             import TemporaryDirectory
from rptools.rplibs       import rpCache
from rptools.rpcompletion import rp_completion
from rptools.rpcompletion.rpCompletion import update_rppaths, build_side_rxn, rxns_from_rules, rp2paths_to_dict
from os                   import path  as os_path
from os                   import stat  as os_stat
from io                   import open  as io_open
from json                 import load  as json_load


class Test_rpCompletion(TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, 'ERROR'))
        self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def test_rp_completion(self):
        with TemporaryDirectory() as temp_d:
            result = rp_completion(self.rpcache,
                                   self.rp2_pathways,
                                   self.rp2paths_compounds,
                                   self.rp2paths_pathways,
                                   temp_d,
                                   upper_flux_bound=999999,
                                   lower_flux_bound=0,
                                   max_subpaths_filter=10,
                                   pathway_id='rp_pathway',
                                   compartment_id='MNXC3',
                                   species_group_id='central_species',
                                   sink_species_group_id='rp_sink_species',
                                   pubchem_search=False,
                                   logger=self.logger)
            # Useless to sort files since smiles could be equivalent and not equal, then checksum will be different
            for file, size in self.files:
                self.assertTrue(os_path.isfile(os_path.join(temp_d, file)))
                self.assertEqual(os_stat(os_path.join(temp_d, file)).st_size, size)
         

    rpcache            = rpCache('file')
    data_path          = os_path.join(os_path.dirname(__file__), 'data' , 'lycopene')
    rp2_pathways       = os_path.join(data_path, '1-rp2_pathways.csv')
    rp2paths_compounds = os_path.join(data_path, '2-rp2paths_compounds.tsv')
    rp2paths_pathways  = os_path.join(data_path, '3-rp2paths_pathways.csv')

    files = [
    ('rp_1_11_sbml.xml',  32215),
    ('rp_1_1_sbml.xml',   32499),
    ('rp_1_6_sbml.xml',   32084),
    ('rp_2_12_sbml.xml',  32212),
    ('rp_2_22_sbml.xml',  32340),
    ('rp_2_2_sbml.xml',   32624),
    ('rp_3_10_sbml.xml',  33620),
    ('rp_3_131_sbml.xml', 34406),
    ('rp_3_132_sbml.xml', 34669),
    ('rp_3_140_sbml.xml', 33208),
    ('rp_3_1_sbml.xml',   34815),
    ('rp_3_261_sbml.xml', 34535),
    ('rp_3_262_sbml.xml', 34797),
    ('rp_3_270_sbml.xml', 33337),
    ('rp_3_2_sbml.xml',   35079),
    ]

    def test_update_rppaths(self):
        current_path_id = 2
        path_step = 3
        sub_path_step = 1
        rxn = {'rule_id'           : 'RR-02-ae41ec5771136ea5-14-F',
               'rule_ori_reac'     : 'MNXR113128',
               'rule_score'        : 0.7358363677022237,
               'left'              : {'CMPD_0000000001': 1},
               'right'             : {'TARGET_0000000001': 1},
               'path_id'           : current_path_id,
               'step'              : path_step,
               'transformation_id' : 'TRS_0_0_1'}
        rp_paths = update_rppaths({}, current_path_id, path_step, sub_path_step, rxn)
        self.assertEqual(rp_paths, {
                                    current_path_id: {
                                        path_step: {
                                            sub_path_step: rxn
                                        }
                                    }
                                    })


    def test_build_side_rxn(self):
        cid           = 'CMPD_0000000001'
        index         = 1
        deprecatedCID = {}
        self.assertDictEqual(build_side_rxn(str(index)+'.'+cid, deprecatedCID), {cid: index})


    def test_build_side_rxn_deprecatedCID_NoMatch(self):
        cid           = 'CMPD_0000000001'
        index         = 1
        deprecatedCID = {'CMPD_000000001': 'FOO'}
        self.assertDictEqual(build_side_rxn(str(index)+'.'+cid, deprecatedCID), {cid: index})


    def test_build_side_rxn_deprecatedCID_Match(self):
        cid           = 'CMPD_0000000001'
        index         = 1
        deprecatedCID = {'CMPD_0000000001': 'FOO'}
        self.assertDictEqual(build_side_rxn(str(index)+'.'+cid, deprecatedCID), {deprecatedCID[cid]: index})


    def test_rxns_from_rules(self):
        with open(os_path.join(self.data_path, 'refs', 'rxns_from_rules.json'), 'r') as read_file:
            data = json_load(read_file)
            self.assertDictEqual(rxns_from_rules('RR-02-ae41ec5771136ea5-14-F', self.rpcache.rr_reactions), data)


    def test_rp2paths_to_dict(self):
        with open(os_path.join(self.data_path, 'refs', 'rp2paths_to_dict.json'), 'r') as read_file:
            # object_ook is used to convert str keys into int keys as stored in rpCompletion functions
            data = json_load(read_file, object_hook=lambda d: {int(k) if k.lstrip('-').isdigit() else k: v for k, v in d.items()})
            self.assertDictEqual(rp2paths_to_dict(self.rp2paths_pathways,
                                                  self.rpcache.rr_reactions, self.rpcache.deprecatedCID_cid),
                                 data)