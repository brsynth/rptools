"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

import logging
from unittest             import TestCase
from tempfile             import TemporaryDirectory
from rptools.rplibs       import rpCache, rpSBML
from rptools.rpcompletion import rp_completion
from rptools.rpcompletion.rpCompletion import build_side_rxn, rp2paths_to_dict
from os                   import path  as os_path
from os                   import stat  as os_stat
from io                   import open  as io_open
from json                 import load  as json_load
from json                 import dumps as json_dumps


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
            rpsbml = rpSBML(os_path.join(temp_d, self.rp_1_11_xml))
            with open(self.rp_1_11_json, 'r') as f:
                self.assertDictEqual(rpsbml.toDict(), json_load(f))
                # self.assertEqual(os_stat(os_path.join(temp_d, file)).st_size, size)
         

    rpcache            = rpCache('file')
    data_path          = os_path.join(os_path.dirname(__file__), 'data' , 'lycopene')
    rp2_pathways       = os_path.join(data_path, '1-rp2_pathways.csv')
    rp2paths_compounds = os_path.join(data_path, '2-rp2paths_compounds.tsv')
    rp2paths_pathways  = os_path.join(data_path, '3-rp2paths_pathways.csv')
    test_file_pattern  = 'rp_1_11'
    rp_1_11_xml        = test_file_pattern+'_sbml.xml'
    rp_1_11_json       = os_path.join(data_path, 'refs', test_file_pattern+'.json')

    files = [
    ('rp_1_11_sbml.xml',  32217),
    ('rp_1_1_sbml.xml',   32501),
    ('rp_1_6_sbml.xml',   32086),
    ('rp_2_12_sbml.xml',  32214),
    ('rp_2_22_sbml.xml',  32342),
    ('rp_2_2_sbml.xml',   32626),
    ('rp_3_10_sbml.xml',  33622),
    ('rp_3_131_sbml.xml', 34408),
    ('rp_3_132_sbml.xml', 34671),
    ('rp_3_140_sbml.xml', 33210),
    ('rp_3_1_sbml.xml',   34817),
    ('rp_3_261_sbml.xml', 34537),
    ('rp_3_262_sbml.xml', 34799),
    ('rp_3_270_sbml.xml', 33339),
    ('rp_3_2_sbml.xml',   35081),
    ]

    # def test_update_rppaths(self):
    #     path_base_id = 2
    #     path_step = 3
    #     path_variant_idx = 1
    #     rxn = {'rule_id'           : 'RR-02-ae41ec5771136ea5-14-F',
    #            'rule_ori_reac'     : 'MNXR113128',
    #            'rule_score'        : 0.7358363677022237,
    #            'left'              : {'CMPD_0000000001': 1},
    #            'right'             : {'TARGET_0000000001': 1},
    #            'step'              : path_step,
    #            'transformation_id' : 'TRS_0_0_1'}
    #     rp_paths = update_rppaths({}, path_base_id, path_variant_idx, rxn)
    #     self.assertEqual(rp_paths, {
    #                                 path_base_id: {
    #                                     path_step: {
    #                                         path_variant_idx: rxn
    #                                     }
    #                                 }
    #                                 })


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


    # def test_rxns_from_rules(self):
    #     with open(os_path.join(self.data_path, 'refs', 'rxns_from_rules.json'), 'r') as read_file:
    #         data = json_load(read_file)
    #         self.assertDictEqual(rxns_from_rules('RR-02-ae41ec5771136ea5-14-F', self.rpcache.rr_reactions), data)


    def test_rp2paths_to_dict(self):
        with open(os_path.join(self.data_path, 'rp2paths_pathways.json'), 'r') as read_file:
            # object_hook is used to convert str keys into int keys as stored in rpCompletion functions
            data = json_load(read_file, object_hook=lambda d: {int(k) if k.lstrip('-').isdigit() else k: v for k, v in d.items()})
            self.assertDictEqual(rp2paths_to_dict(self.rp2paths_pathways,
                                                  self.rpcache.rr_reactions, self.rpcache.deprecatedCID_cid),
                                 data)