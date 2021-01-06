"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

import logging
from unittest             import TestCase
from tempfile             import TemporaryDirectory
from rptools.rplibs       import rpCache
from rptools.rpcompletion import rp_completion
from os                   import path as os_path
from os                   import stat as os_stat


class Test_rpCompletion(TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, 'ERROR'))
        self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def test_rp2ToSBML(self):
        with TemporaryDirectory() as temp_d:
            result = rp_completion(self.cache,
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
         


    cache              = rpCache('file')
    data_path          = os_path.join(os_path.dirname(__file__), 'data' , 'lycopene')
    rp2_pathways       = os_path.join(data_path, '1-rp2_pathways.csv')
    rp2paths_compounds = os_path.join(data_path, '2-rp2paths_compounds.tsv')
    rp2paths_pathways  = os_path.join(data_path, '3-rp2paths_pathways.csv')

    files = [
    ('rp_1_11_sbml.xml',  32358),
    ('rp_1_1_sbml.xml',   32640),
    ('rp_1_6_sbml.xml',   32225),
    ('rp_2_12_sbml.xml',  32355),
    ('rp_2_22_sbml.xml',  32483),
    ('rp_2_2_sbml.xml',   32765),
    ('rp_3_10_sbml.xml',  33763),
    ('rp_3_131_sbml.xml', 34551),
    ('rp_3_132_sbml.xml', 34814),
    ('rp_3_140_sbml.xml', 33353),
    ('rp_3_1_sbml.xml',   34956),
    ('rp_3_261_sbml.xml', 34680),
    ('rp_3_262_sbml.xml', 34942),
    ('rp_3_270_sbml.xml', 33482),
    ('rp_3_2_sbml.xml',   35220),
    ]
