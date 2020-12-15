"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

import logging
from unittest       import TestCase
from tempfile       import TemporaryDirectory
from rptools.rplibs import inchikeyMIRIAM
from os             import path as os_path
from _main          import Main


class Test_inchikeyMIRIAM(TestCase):

    data_path = os_path.join(os_path.dirname(__file__), 'data')

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, 'ERROR'))
        self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def test_inchikeyMIRIAM(self):
        inchi = inchikeyMIRIAM()
        with TemporaryDirectory() as tempd:
            output_sbml = os_path.join(tempd, 'output.sbml')
            inchi.addInChiKey(os_path.join(self.data_path,'e_coli_model.sbml'), output_sbml)
            self.assertTrue(Main._check_file_hash(output_sbml, '0a26fa7dfc49480f87f79292af26fda39d10398a08f32ab4163de3908f814b61'))
