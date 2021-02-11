from unittest import TestCase
from os import path as os_path
from rptools.rpfba.rpFBA import (
    rp_fba,
    rp_fraction,
    rp_pfba
)
from rptools.rplibs import rpSBML
from tempfile import (
    TemporaryDirectory,
    mkdtemp
)
from shutil    import rmtree
from brs_utils import (
    create_logger,
    extract_gz
)


class Test_rpFBA(TestCase):


    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    merged_path_gz = os_path.join(
        data_path,
        'merged.xml.gz'
    )
    rpsbml_path = os_path.join(
        data_path,
        'rp_1_1_sbml.xml'
    )
    e_coli_model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )

    rxn_tgt = 'rxn_target'


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        self.rpsbml = rpSBML(
            inFile = self.rpsbml_path,
            logger = self.logger
        )

        with TemporaryDirectory() as temp_d:
            self.e_coli_model_path = extract_gz(
                self.e_coli_model_path_gz,
                temp_d
            )
            self.merged_path = extract_gz(
                self.merged_path_gz,
                temp_d
            )
            self.rpsbml_model = rpSBML(
                inFile = self.e_coli_model_path,
                logger = self.logger
            )
            self.merged_rpsbml_1 = rpSBML(
                inFile = self.merged_path,
                logger = self.logger
            )
            self.merged_rpsbml_2, reactions_in_both = rpSBML.mergeModels(
                source_rpsbml = self.rpsbml,
                target_rpsbml = self.rpsbml_model,
                logger = self.logger
            )


    def test_fba_1(self):
        ref_score = 9.230769230769237
        obj_value, rpsbml = rp_fba(
                 rpsbml = self.merged_rpsbml_1,
            reaction_id = self.rxn_tgt,
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
            pathway['fba_obj_'+self.rxn_tgt]['value'],
            ref_score
        )


    def test_fba_2(self):
        ref_score = 5.144315545243618
        obj_value, rpsbml = rp_fba(
                 rpsbml = self.merged_rpsbml_2,
            reaction_id = self.rxn_tgt,
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
            pathway['fba_obj_'+self.rxn_tgt]['value'],
            ref_score
        )


    def test_fraction_1(self):
        ref_score = 2.3076923076923888
        obj_value, rpsbml = rp_fraction(
                rpsbml = self.merged_rpsbml_1,
            src_rxn_id = 'biomass',
             src_coeff = 1.0,
            tgt_rxn_id = self.rxn_tgt,
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
        self.assertAlmostEqual(
            pathway['fba_obj_'+self.rxn_tgt+'__restricted_biomass']['value'],
            ref_score
        )
        self.assertAlmostEqual(
            pathway['fba_obj_biomass']['value'],
            3.6794124272706443
        )


    def test_fraction_2(self):
        ref_score = 1.3296695186776557
        obj_value, rpsbml = rp_fraction(
                rpsbml = self.merged_rpsbml_2,
            src_rxn_id = 'biomass',
             src_coeff = 1.0,
            tgt_rxn_id = self.rxn_tgt,
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
        self.assertAlmostEqual(
            pathway['fba_obj_'+self.rxn_tgt+'__restricted_biomass']['value'],
            ref_score
        )
        self.assertAlmostEqual(
            pathway['fba_obj_biomass']['value'],
            0.7638744755010182
        )


    def test_pfba_1(self):
        ref_score = 859.3846153846168
        obj_value, rpsbml = rp_pfba(
                 rpsbml = self.merged_rpsbml_1,
            reaction_id = self.rxn_tgt,
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
            pathway['fba_obj_'+self.rxn_tgt]['value'],
            ref_score
        )

    def test_pfba_2(self):
        ref_score = 471.6390839533362
        obj_value, rpsbml = rp_pfba(
                 rpsbml = self.merged_rpsbml_2,
            reaction_id = self.rxn_tgt,
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
            pathway['fba_obj_'+self.rxn_tgt]['value'],
            ref_score
        )
