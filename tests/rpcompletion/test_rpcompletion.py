"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

import logging
from tempfile             import TemporaryDirectory
from rr_cache import rrCache
from rptools.rplibs       import rpSBML
from rptools.rpcompletion import rp_completion
from rptools.rpcompletion.rpCompletion import (
    build_side_rxn,
    rp2paths_to_dict
)
from os                   import path  as os_path
from os                   import (
    stat as os_stat,
    listdir
)
from pathlib import Path
from io                   import open  as io_open
from json                 import load  as json_load
from json                 import dumps as json_dumps
from unittest import TestCase
from brs_utils import (
    create_logger,
)


class Test_rpCompletion(TestCase):


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')


    def test_rp_completion(self):
        with TemporaryDirectory() as temp_d:
            temp_d = '/tmp/joan20'
            result = rp_completion(
                self.cache,
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
                logger=self.logger
            )
            # Useless to sort files since smiles could be equivalent and not equal, then checksum will be different
            for file in listdir(temp_d):
                self.assertEqual(
                    self.files[file],
                    Path(os_path.join(temp_d, file)).stat().st_size
                )
            rpsbml = rpSBML(os_path.join(temp_d, self.rpsbml_xml))
            # print(json_dumps(rpsbml.toDict(), indent=4))
            # self.assertTrue(False)
            # exit()
            with open(self.rpsbml_json, 'r') as f:
                self.assertDictEqual(rpsbml.toDict(), json_load(f))
                # self.assertEqual(os_stat(os_path.join(temp_d, file)).st_size, size)
         

    cache            = rrCache('file')
    data_path          = os_path.join(os_path.dirname(__file__), 'data' , 'lycopene')
    rp2_pathways       = os_path.join(data_path, '1-rp2_pathways.csv')
    rp2paths_compounds = os_path.join(data_path, '2-rp2paths_compounds.tsv')
    rp2paths_pathways  = os_path.join(data_path, '3-rp2paths_pathways.csv')
    test_file_pattern  = 'rp_002_0022'
    rpsbml_xml         = test_file_pattern+'_sbml.xml'
    rpsbml_json        = os_path.join(data_path, 'refs', test_file_pattern+'.json')

    files = {
        'rp_001_0011_sbml.xml': 32217,
        'rp_001_0001_sbml.xml': 32501,
        'rp_001_0006_sbml.xml': 32086,
        'rp_002_0012_sbml.xml': 32214,
        'rp_002_0022_sbml.xml': 32465,
        'rp_002_0002_sbml.xml': 32626,
        'rp_003_0001_sbml.xml': 34943,
        'rp_003_0002_sbml.xml': 35207,
        'rp_003_0010_sbml.xml': 33746,
        'rp_003_0131_sbml.xml': 34530,
        'rp_003_0132_sbml.xml': 34794,
        'rp_003_0140_sbml.xml': 33332,
        'rp_003_0261_sbml.xml': 34658,
        'rp_003_0262_sbml.xml': 34922,
        'rp_003_0270_sbml.xml': 33461,
    }

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


    def test_rp2paths_to_dict(self):
        with open(os_path.join(self.data_path, 'refs', 'rp2paths_pathways.json'), 'r') as read_file:
            # object_hook is used to convert str keys into int keys as stored in rpCompletion functions
            data = json_load(read_file, object_hook=lambda d: {int(k) if k.lstrip('-').isdigit() else k: v for k, v in d.items()})
            self.assertDictEqual(rp2paths_to_dict(self.rp2paths_pathways,
                                                  self.cache.rr_reactions, self.cache.deprecatedCID_cid),
                                 data)