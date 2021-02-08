import logging
from unittest import TestCase
from os       import path as os_path
from rptools.rpfba.rpFBA import rp_fba, rp_fraction, rp_pfba
from rptools.rplibs      import rpSBML
from brs_utils import create_logger


class Test_rpFBA(TestCase):

    """
    @classmethod
    def setUpClass(self):
    """
    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    rpsbml = rpSBML(
        os_path.join(
            data_path,
            'merged.xml'
        )
    )


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')


    def test_fba(self):
        ref_score = 9.230769230769237
        obj_value, rpsbml = rp_fba(
                 rpsbml = self.rpsbml,
            reaction_id = 'rxn_sink',
                 logger = self.logger
        )
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(
            obj_value,
            ref_score
        )
        # make sure that the results are written to the file
        pathway = rpsbml.toDict()['pathway']['brsynth']
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_sink']['value'],
            ref_score
        )


    def test_fraction(self):
        ref_score = 2.3076923076923888
        obj_value, rpsbml = rp_fraction(
                rpsbml = self.rpsbml,
            src_rxn_id = 'biomass',
             src_coeff = 1.0,
            tgt_rxn_id = 'rxn_sink',
             tgt_coeff = 1.0,
                logger = self.logger
        )

        self.assertTrue(rpsbml)
        self.assertAlmostEqual(
            obj_value,
            ref_score
        )

        # make sure that the results are written to the file
        pathway = rpsbml.toDict()['pathway']['brsynth']
        print(pathway)
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_sink']['value'],
            ref_score
        )
        self.assertAlmostEqual(
            pathway['fba_obj_biomass']['value'],
            3.6794124272706443
        )


    def test_pfba(self):
        ref_score = 859.3846153846168
        obj_value, rpsbml = rp_pfba(
                 rpsbml = self.rpsbml,
            reaction_id = 'rxn_sink',
                 logger = self.logger
        )
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(
            obj_value,
            ref_score
        )
        # make sure that the results are written to the file
        pathway = rpsbml.toDict()['pathway']['brsynth']
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_sink']['value'],
            ref_score
        )
