from unittest import TestCase
from os import path as os_path
# from rptools.rpfba import rpFBA
from rptools.rpscore.rpScore import (
    minmax_score,
    score_from_reactions,
    score_from_pathway,
    compute_globalscore
)
from rptools.rplibs import rpSBML
# from tempfile import (
#     TemporaryDirectory,
#     mkdtemp
# )
# from shutil    import rmtree
from brs_utils import (
    create_logger,
    # extract_gz
)


class Test_rpScore(TestCase):

    __test__ = False

    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    rpsbml_path = os_path.join(
        data_path,
        'fba.sbml'
    )
    bounds = {
        'thermo': {
            'floor': 5000.0,
            'ceil' : -5000.0
        },
        'fba': {
            'floor': 0.0,
            'ceil' : 5.0
        },
        'max_rp_steps': 15
    }
    pathway_id = 'rp_pathway'


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        self.rpsbml = rpSBML(
            inFile = self.rpsbml_path,
            logger = self.logger
        )


    def test_minmax_score_in(self):
        value = 1.3
        floor = 0.5
        ceil = 2.0
        self.assertEqual(
            (value-floor) / (ceil-floor),
            minmax_score(
                value,
                floor,
                ceil
            )
        )


    def test_minmax_score_out_inf(self):
        value = 0.1
        floor = 0.5
        ceil = 2.0
        self.assertEqual(
            0.0,
            minmax_score(
                value,
                floor,
                ceil
            )
        )


    def test_minmax_score_out_sup(self):
        value = 2.1
        floor = 0.5
        ceil = 2.0
        self.assertEqual(
            1.0,
            minmax_score(
                value,
                floor,
                ceil
            )
        )


    def test_scores_from(self):
        scores = {}
        # score_names = ['dfG_prime_m', 'dfG_uncert', 'dfG_prime_o', 'rule_score', 'fba_obj_biomass', 'fba_obj_fraction']
        rpsbml_dict = self.rpsbml.toDict(self.pathway_id)
        rpsbml_dict, scores = score_from_reactions(
            rpsbml_dict,
            scores,
            self.bounds,
            self.pathway_id,
            self.logger
        )
        self.assertAlmostEqual(
            3.6794124272706443,
            rpsbml_dict['pathway']['brsynth']['fba_obj_biomass']['value']
        )
        self.assertAlmostEqual(
            2.3076923076923888,
            rpsbml_dict['pathway']['brsynth']['fba_obj_fraction']['value']
        )
        self.assertAlmostEqual(
            0.0,
            rpsbml_dict['reactions']['rxn_1']['brsynth']['fba_obj_biomass']['value']
        )
        self.assertAlmostEqual(
            0.0,
            rpsbml_dict['reactions']['rxn_1']['brsynth']['fba_obj_biomass']['value']
        )
        self.assertAlmostEqual(
            0.0,
            rpsbml_dict['reactions']['rxn_1']['brsynth']['norm_fba_obj_biomass']
        )
        self.assertAlmostEqual(
            0.4615384615384778,
            rpsbml_dict['reactions']['rxn_1']['brsynth']['norm_fba_obj_fraction']
        )
        rpsbml_dict = score_from_pathway(
            rpsbml_dict,
            scores,
            self.bounds,
            self.pathway_id,
            self.logger
        )
        self.assertAlmostEqual(
            0.7358824854541288,
            rpsbml_dict['pathway']['brsynth']['norm_fba_obj_biomass']
        )
        self.assertAlmostEqual(
            0.4615384615384778,
            rpsbml_dict['pathway']['brsynth']['norm_fba_obj_fraction']
        )
        self.assertAlmostEqual(
            1.0,
            rpsbml_dict['pathway']['brsynth']['norm_rule_score']
        )
        self.assertAlmostEqual(
            1.0,
            rpsbml_dict['pathway']['brsynth']['norm_steps']
        )

    def test_compute_globalscore(self):
        rpsbml_dict = compute_globalscore(self.rpsbml)
        self.assertAlmostEqual(
            0.7939839049864056,
            rpsbml_dict['pathway']['brsynth']['global_score']
        )
