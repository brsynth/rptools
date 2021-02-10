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
    model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )

    def setUp(self):
        self.logger = create_logger(__name__, 'DEBUG')

        # Create persistent temp folder
        # to deflate compressed data file so that
        # it remains reachable outside of this method.
        # Has to remove manually it in tearDown() method 
        self.temp_d = mkdtemp()
        self.merged_path = extract_gz(
            self.merged_path_gz,
            self.temp_d
        )
        self.model_path = extract_gz(
            self.model_path_gz,
            self.temp_d
        )
        # objects below have to be created for each test instance
        # since some tests can modified them
        self.rpsbml = rpSBML(
            inFile = self.model_path,
            logger = self.logger
        )
        import json
        print(json.dumps(self.rpsbml.toDict(), indent=4))
        exit()


    def tearDown(self):
        rmtree(self.temp_d)


    def test_fba(self):
        ref_score = 9.230769230769237
        cobra_results = rp_fba(
                 rpsbml = self.rpsbml,
            reaction_id = 'rxn_target',
                 logger = self.logger
        )
        self.assertTrue(self.rpsbml)
        self.assertAlmostEqual(
            cobra_results.objective_value,
            ref_score
        )
        # make sure that the results are written to the file
        pathway = self.rpsbml.toDict()['pathway']['brsynth']
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_target']['value'],
            ref_score
        )


    def test_fraction(self):
        ref_score = 2.3076923076923888
        cobra_results = rp_fraction(
                rpsbml = self.rpsbml,
            src_rxn_id = 'biomass',
             src_coeff = 1.0,
            tgt_rxn_id = 'rxn_target',
             tgt_coeff = 1.0,
                logger = self.logger
        )

        self.assertTrue(self.rpsbml)
        self.assertAlmostEqual(
            cobra_results.objective_value,
            ref_score
        )

        # make sure that the results are written to the file
        pathway = rpsbml.toDict()['pathway']['brsynth']
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_target__restricted_biomass']['value'],
            ref_score
        )
        self.assertAlmostEqual(
            pathway['fba_obj_biomass']['value'],
            3.6794124272706443
        )


    def test_pfba(self):
        # ref_score = 859.3846153846168
        ref_score = 1761.107066440342
        cobra_results = rp_pfba(
                 rpsbml = self.rpsbml,
            reaction_id = 'rxn_target',
                 logger = self.logger
        )
        self.assertTrue(self.rpsbml)
        self.assertAlmostEqual(
            cobra_results.objective_value,
            ref_score
        )
        # make sure that the results are written to the file
        pathway = self.rpsbml.toDict()['pathway']['brsynth']
        import json
        print(json.dumps(self.rpsbml.toDict(), indent=4))
        exit()
        self.assertAlmostEqual(
            pathway['fba_obj_rxn_target']['value'],
            ref_score
        )
