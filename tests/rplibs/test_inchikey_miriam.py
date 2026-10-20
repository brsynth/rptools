"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from os import path as os_path
from tempfile import NamedTemporaryFile

from brs_utils import extract_gz_to_string
from main_rplibs import Main_rplibs

from rptools.rplibs.inchikey_miriam import inchikeyMIRIAM


class Test_inchikeyMIRIAM(Main_rplibs):

    __test__ = False

    # def setUp(self):
    #     super().setUp()

    def test_inchikeyMIRIAM(self):

        inchi = inchikeyMIRIAM()

        with NamedTemporaryFile() as temp_f:

            inchi.addInChiKey(self.e_coli_model_path, temp_f.name)

            test = temp_f.read().decode("utf-8")

            ref = extract_gz_to_string(
                os_path.join(self.data_path, "output_inchikeyMIRIAM.sbml.gz")
            )

            self.assertListEqual(ref.split(), test.split())
