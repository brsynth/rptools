import logging
from unittest import TestCase
from os       import path as os_path
from rptools.rpfba.rpFBA import rp_fba, rp_fraction, rp_pfba
from rptools.rplibs      import rpSBML

class Test_rpFBA(TestCase):

    """
    @classmethod
    def setUpClass(self):
    """
    data_path = os_path.join(os_path.dirname(__file__), 'data')

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging, 'ERROR'))
        self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    def test_fba(self):
        rpsbml = rpSBML(os_path.join(self.data_path, 'merged.xml'))
        obj_value, rpsbml = rp_fba(rpsbml, 'Rxn_sink', logger=self.logger)
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 9.230769230769237)
        # make sure that the results are written to the file
        rpsbml_dict = rpsbml.toDict()
        self.assertAlmostEqual(rpsbml_dict['pathway']['brsynth']['fba_obj_Rxn_sink']['value'], 9.230769230769237)

    def test_fraction(self):
        rpsbml = rpSBML(os_path.join(self.data_path, 'merged.xml'))
        obj_value, rpsbml = rp_fraction(rpsbml, 'biomass', 1.0, 'Rxn_sink', 1.0, logger=self.logger)
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 2.3076923076923888)
        # make sure that the results are written to the file
        rpsbml_dict = rpsbml.toDict()
        self.assertAlmostEqual(rpsbml_dict['pathway']['brsynth']['fba_obj_Rxn_sink__restricted_biomass']['value'], 2.3076923076923888)
        self.assertAlmostEqual(rpsbml_dict['pathway']['brsynth']['fba_obj_biomass']['value'], 3.6794124272706443)

    def test_pfba(self):
        rpsbml = rpSBML(os_path.join(self.data_path, 'merged.xml'))
        obj_value, rpsbml = rp_pfba(rpsbml, 'Rxn_sink', logger=self.logger)
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 859.3846153846168)
        # make sure that the results are written to the file
        rpsbml_dict = rpsbml.toDict()
        self.assertAlmostEqual(rpsbml_dict['pathway']['brsynth']['fba_obj_Rxn_sink']['value'], 859.3846153846168)
